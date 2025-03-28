#include "Simulation.h"
#include "ForcingData.h"
#include "InputManager.h"
#include "OutputManager.h"
#include "Timer.h"
#include "MeshFactory.h"
#include "NetworkMeshFactory.h"
#include "MaterialLibrary.h"
#include "NetworkMaterialLibrary.h"
#include "PhysicsState.h"
#include "DateTimeUtils.h"
#include <iostream>
#include <stdexcept>

/**
 * @brief Constructs a Simulation object with the specified JSON settings.
 * 
 * This constructor initializes the simulation by parsing the settings and
 * calling the internal initialize() function.
 * 
 * @param settings Parsed JSON settings.
 */
Simulation::Simulation(const JsonValue &settings)
    : settings_(settings)
{
    initialize();
}

/**
 * @brief Initializes the simulation.
 * 
 * This function reads the time control settings from the JSON configuration,
 * loads forcing data, creates and initializes the surface, subsurface, and network meshes,
 * maps plant parameters from the surface to subsurface cells, creates the mesh coupler,
 * instantiates the simulation models, and sets up the output manager.
 * 
 * @throws std::runtime_error If any required settings are missing or if data loading fails.
 */
void Simulation::initialize() {
	// Print the settings json data.
	settings_.print(settings_, 0);

    // Validate and read time control settings.
    if (!settings_.is_object() ||
        settings_.as_object().find("TimeControl") == settings_.as_object().end()) {
        throw std::runtime_error("Missing 'TimeControl' section in settings.");
    }
    
    const auto &timeControl = settings_.as_object().at("TimeControl");
    if (!timeControl.is_object() ||
        timeControl.as_object().find("StartDateTime") == timeControl.as_object().end() ||
        timeControl.as_object().find("EndDateTime") == timeControl.as_object().end() ||
        timeControl.as_object().find("TimeStep") == timeControl.as_object().end() ||
        timeControl.as_object().find("ResultsSaveInterval") == timeControl.as_object().end()) {
        throw std::runtime_error("Missing one or more required 'TimeControl' settings.");
    }
	
    // Parse simulation start and end times.
    std::string startDateTimeStr = timeControl.as_object().at("StartDateTime").as_string();
	startTime_ = parseDateTime(startDateTimeStr);
	std::string endDateTimeStr = timeControl.as_object().at("EndDateTime").as_string();
	std::time_t endTime = parseDateTime(endDateTimeStr);
    simulationLength_ = static_cast<double>(std::difftime(endTime, startTime_));
	timeStep_ = timeControl.as_object().at("TimeStep").as_double();
	resultsSaveInterval_ = timeControl.as_object().at("ResultsSaveInterval").as_double();
	std::cout << "timeStep: " << timeStep_ << std::endl;
    std::cout << "Simulation start: " << startDateTimeStr
              << "  End: " << endDateTimeStr
              << "  Duration (s): " << simulationLength_ << std::endl;

    // Convert start time to Excel serial date.
    startTimeExcel_ = timeToExcelDate(startTime_);

    const auto &files = settings_.as_object().at("Files").as_object();

    try {
        // Load environmental forcing data.
        auto forcingCSV = InputManager::loadCSV(files.at("AtmosphericForcing").as_string());
        forcingData_.load(forcingCSV, startTime_);

        // Create the surface mesh.
        surfaceMesh_ = MeshFactory::createSurfaceMesh(
            files.at("SurfaceGridVTK").as_string(),
            files.at("Materials2D").as_string(),
            files.at("InitialConditions2D").as_string(),
            files.at("BoundaryConditions2D").as_string());

        // Create the subsurface mesh.
        subsurfaceMesh_ = MeshFactory::createSubsurfaceMesh(
            files.at("SubsurfaceGridVTK").as_string(),
            files.at("Materials3D").as_string(),
            files.at("InitialConditions3D").as_string(),
            files.at("BoundaryConditions3D").as_string());

        // Map plant parameters from surface cells to corresponding subsurface cells.
        int numSurfaceCells = surfaceMesh_->cells.size();
        int numSubsurfaceCells = subsurfaceMesh_->cells.size();
        if (numSurfaceCells > 0) {
            int numLayers = numSubsurfaceCells / numSurfaceCells;
            std::cout << "Mapping plant parameters: " << numLayers
                      << " subsurface layers per surface cell." << std::endl;
            for (int i = 0; i < numSurfaceCells; i++) {
                // Get surface cell water state.
                SurfaceCell* surfCell = dynamic_cast<SurfaceCell*>(surfaceMesh_->cells[i].get());
                if (!surfCell || !surfCell->waterState)
                    continue;
                auto sws = std::dynamic_pointer_cast<SurfaceWaterState>(surfCell->waterState);
                if (!sws)
                    continue;
                // For each vertical layer under the current surface cell...
                for (int layer = 0; layer < numLayers; layer++) {
                    int idx = i + layer * numSurfaceCells;
                    if (idx < numSubsurfaceCells) {
                        VolumeCell* volCell = dynamic_cast<VolumeCell*>(subsurfaceMesh_->cells[idx].get());
                        if (!volCell || !volCell->waterState)
                            continue;
                        auto subWater = std::dynamic_pointer_cast<SubsurfaceWaterState>(volCell->waterState);
                        if (!subWater)
                            continue;
                        subWater->cropFactor = sws->cropFactor;
                        subWater->rootDepth = sws->rootDepth;
                        subWater->evaporationFraction = sws->evaporationFraction;
                    }
                }
            }
        }

        // Create the network mesh.
        networkMesh_ = NetworkMeshFactory::createNetworkMesh(
            files.at("NetworkJunctionVTK").as_string(),
            files.at("NetworkLinkVTK").as_string(),
            files.at("MaterialsNetJunc").as_string(),
            files.at("MaterialsNetLink").as_string(),
            files.at("InitialConditionsNetJunc").as_string(),
            files.at("BoundaryConditionsNetJunc").as_string(),
            files.at("InitialConditionsNetLink").as_string(),
            files.at("BoundaryConditionsNetLink").as_string());

        // Instantiate the mesh coupler.
        meshCoupler_ = std::make_unique<MeshCoupler>(surfaceMesh_.get(), subsurfaceMesh_.get(), networkMesh_.get());
		
        // Create simulation model objects and store them.
        models_.push_back(std::make_unique<SurfaceSnowModel>(surfaceMesh_.get()));
		models_.push_back(std::make_unique<NetworkFlowModel>(networkMesh_.get()));
        models_.push_back(std::make_unique<SurfaceFlowModel>(surfaceMesh_.get()));
        models_.push_back(std::make_unique<SubsurfaceFlowModel>(subsurfaceMesh_.get()));
		models_.push_back(std::make_unique<SubsurfaceHeatTransferModel>(subsurfaceMesh_.get()));
		
        // Initialize each model and set the mesh coupler.
        for (auto &model : models_) {
            model->initialize(settings_);
            if (auto netModel = dynamic_cast<NetworkFlowModel*>(model.get()))
                netModel->setMeshCoupler(meshCoupler_.get());
            else if (auto surfModel = dynamic_cast<SurfaceFlowModel*>(model.get()))
                surfModel->setMeshCoupler(meshCoupler_.get());
            else if (auto subModel = dynamic_cast<SubsurfaceFlowModel*>(model.get()))
                subModel->setMeshCoupler(meshCoupler_.get());
        }
		
		// Create the output manager with the specified results path.
		std::string resultsPath = files.at("ResultsCSV").as_string();
		outputManager_ = std::make_unique<OutputManager>(resultsPath);
		outputManager_->writeHeader();
    }
    catch (const std::exception &ex) {
        throw std::runtime_error("Error during input data loading: " + std::string(ex.what()));
    }
}

/**
 * @brief Processes a single simulation time step.
 *
 * This function queries the environmental forcing data for the current simulation time,
 * then runs each model for a step with the retrieved forcing. It prints status information to the console.
 */
void Simulation::processTimeStep() {
    // Query environmental forcing data.
    ForcingRecord fr = forcingData_.getForcingAtTime(currentTime_);
    double precipitation = fr.precipitation;  // in mm
    double pet = fr.pet;                      // in mm
    double temp = fr.airTemperature;          // in °C

    std::cout << "Time " << currentTime_ << " s => P=" << precipitation
              << " mm, PET=" << pet << " mm, T=" << temp << " C" << std::endl;

    // Run each model for one step.
    for (auto &model : models_) {
        model->runStep(timeStep_, fr);
    }
}

/**
 * @brief Outputs simulation results and writes VTK files if configured.
 *
 * This function computes water volumes from each model and then writes a CSV row of results.
 * It also writes VTK output files for the network, surface, and subsurface meshes if the respective
 * output folder paths are non-empty.
 *
 * @param currentTime The simulation time (in seconds) corresponding to the results.
 */
void Simulation::outputResults(double currentTime) {
    // Convert simulation time to Excel serial date (in days).
    double excelTime = startTimeExcel_ + (currentTime / 86400.0);

    // Compute water volumes from each model.
    double precipitationWaterVolume = 0.0;
    double petWaterVolume = 0.0;
    double upperStorageWaterVolume = 0.0;
    double evaporationWaterVolume = 0.0;
    double overlandWaterVolume = 0.0;
    double outfallWaterVolume = 0.0;
    double networkWaterVolume = 0.0;
    double infiltrationWaterVolume = 0.0;
    double transpirationWaterVolume = 0.0;
    double subsurfaceWaterVolume = 0.0;

    for (auto &model : models_) {
        if (auto sfm = dynamic_cast<SurfaceFlowModel*>(model.get())) {
            precipitationWaterVolume = sfm->getPrecipitationWaterVolume();
            petWaterVolume = sfm->getPetWaterVolume();
            upperStorageWaterVolume = sfm->getUpperStorageWaterVolume();
            evaporationWaterVolume = sfm->getEvaporationWaterVolume();
            overlandWaterVolume = sfm->getOverlandWaterVolume();
        } else if (auto ssm = dynamic_cast<SubsurfaceFlowModel*>(model.get())) {
            infiltrationWaterVolume = ssm->getInfiltrationWaterVolume();
            transpirationWaterVolume = ssm->getTranspirationWaterVolume();
            subsurfaceWaterVolume = ssm->getSubsurfaceWaterVolume();
        } else if (auto nfm = dynamic_cast<NetworkFlowModel*>(model.get())) {
            outfallWaterVolume = nfm->getOutfallWaterVolume();
            networkWaterVolume = nfm->getNetworkWaterVolume();
        }
    }
    
    // Write a CSV row of results.
    outputManager_->writeRow(
        static_cast<int>(currentTime),
        excelTime,
        precipitationWaterVolume,
        petWaterVolume,
        upperStorageWaterVolume, 
        evaporationWaterVolume, 
        overlandWaterVolume,
        outfallWaterVolume,
        networkWaterVolume,
        infiltrationWaterVolume, 
        transpirationWaterVolume,
        subsurfaceWaterVolume
    );
    
    // Print current storage values for inspection.
    std::cout << "Outputting results at t = " << currentTime << " s (Excel date = " << excelTime << ")" << std::endl
              << "Precipitation Water Volume = " << precipitationWaterVolume << " m^3" << std::endl
              << "Upper Storage Water Volume = " << upperStorageWaterVolume << " m^3" << std::endl
              << "Evaporation Water Volume = " << evaporationWaterVolume << " m^3" << std::endl
              << "Overland Water Volume = " << overlandWaterVolume << " m^3" << std::endl
              << "Outfall Water Volume = " << outfallWaterVolume << " m^3" << std::endl
              << "Network Water Volume = " << networkWaterVolume << " m^3" << std::endl
              << "Infiltration Water Volume = " << infiltrationWaterVolume << " m^3" << std::endl
              << "Transpiration Water Volume = " << transpirationWaterVolume << " m^3" << std::endl
              << "Subsurface Water Volume = " << subsurfaceWaterVolume << " m^3" << std::endl;

    // Write VTK output files if output folders are specified.
    const auto& filesObj = settings_.as_object().at("Files").as_object();
    std::string networkLinkFolder = filesObj.at("NetworkLinkVTKOutputFolder").as_string();
    std::string networkJuncFolder = filesObj.at("NetworkJunctionVTKOutputFolder").as_string();
    std::string surfaceFolder = filesObj.at("SurfaceVTKOutputFolder").as_string();
    std::string subsurfaceFolder = filesObj.at("SubsurfaceVTKOutputFolder").as_string();
    int fileNumber = static_cast<int>(currentTime / resultsSaveInterval_);
	
    if (!networkLinkFolder.empty() && !networkJuncFolder.empty()) {
        outputManager_->writeNetworkVTK(networkLinkFolder, networkJuncFolder, fileNumber, *networkMesh_);
    }
    if (!surfaceFolder.empty()) {
        outputManager_->writeSurfaceVTK(surfaceFolder, fileNumber, *surfaceMesh_);
    }
    if (!subsurfaceFolder.empty()) {
        outputManager_->writeSubsurfaceVTK(subsurfaceFolder, fileNumber, *subsurfaceMesh_);
    }
}

/**
 * @brief Runs the simulation.
 *
 * This function first outputs the initial state, then enters the main simulation loop.
 * At each time step, it processes the simulation step, outputs progress, and saves results
 * at specified intervals. Finally, it flushes and closes the output.
 */
void Simulation::run() {
    // Save initial state of the system.
    outputResults(currentTime_);
    
    std::cout << "Starting the simulation loop." << std::endl;
    Timer simulationTimer("Total simulation run time");
    int lastPercentage = 0;
    for (currentTime_ = 0.0; currentTime_ <= simulationLength_; currentTime_ += timeStep_) {
        // Process one simulation time step.
        processTimeStep();

        // Update and print simulation progress.
        int currentPercentage = static_cast<int>(((currentTime_ + timeStep_) * 100) / simulationLength_);
        if (currentPercentage > lastPercentage) {
            std::cout << "Simulation progress: " << currentPercentage << "%" << std::endl;
            lastPercentage = currentPercentage;
        }
        
        // Save output results if the current time matches the results save interval.
        if (resultsSaveInterval_ > 0.0 &&
            static_cast<int>(currentTime_ + timeStep_) % static_cast<int>(resultsSaveInterval_) == 0) {
            outputResults(currentTime_ + timeStep_);
        }
    }
    if (lastPercentage < 100) {
        std::cout << "Simulation progress: 100%" << std::endl;
    }
    // Flush and close the output file.
    outputManager_->close();
}
