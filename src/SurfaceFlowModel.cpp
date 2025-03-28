#include "SurfaceFlowModel.h"
#include "PhysicsState.h"
#include <iostream>
#include <cmath>
#include <stdexcept>

/**
 * @brief Initializes the SurfaceFlowModel using settings from the JSON configuration.
 *
 * This method reads the weather data interval from the "TimeControl" section and
 * the surface flow parameters from the "SurfaceFlowBasic" section. If any parameter
 * is missing, a default fallback value is used.
 *
 * @param settings A JSON object containing the simulation configuration.
 */
void SurfaceFlowModel::initialize(const JsonValue &settings) {
    // Load the weather data interval from the "TimeControl" section.
    if (settings.is_object() &&
        settings.as_object().find("TimeControl") != settings.as_object().end() &&
        settings.as_object().at("TimeControl").is_object() &&
        settings.as_object().at("TimeControl").as_object().find("WeatherDataInterval") != settings.as_object().at("TimeControl").as_object().end())
    {
        weatherDataInterval_ = settings.as_object().at("TimeControl").as_object().at("WeatherDataInterval").as_double();
    } else {
        weatherDataInterval_ = 120.0;
    }
	
    // Read surface flow model settings from the "SurfaceFlowBasic" section.
    if (settings.is_object() &&
        settings.as_object().find("WaterFlow") != settings.as_object().end() &&
        settings.as_object().at("WaterFlow").is_object() &&
        settings.as_object().at("WaterFlow").as_object().find("SurfaceFlowBasic") != settings.as_object().at("WaterFlow").as_object().end())
    {
		const JsonValue &surfaceFlowSettings = settings.as_object().at("WaterFlow").as_object().at("SurfaceFlowBasic");

        // Load Enabled flag.
        if (surfaceFlowSettings.is_object() &&
            surfaceFlowSettings.as_object().find("Enabled") != surfaceFlowSettings.as_object().end())
        {
            enabled_ = surfaceFlowSettings.as_object().at("Enabled").as_bool();
        } else {
            enabled_ = false;
        }
		
        if (surfaceFlowSettings.as_object().find("MaxIterationsSurface") != surfaceFlowSettings.as_object().end())
            maxIter_ = surfaceFlowSettings.as_object().at("MaxIterationsSurface").as_int();
        else
            maxIter_ = 100;
		
        if (surfaceFlowSettings.as_object().find("IterationCutThresholdSurface") != surfaceFlowSettings.as_object().end())
            iterThresh_ = surfaceFlowSettings.as_object().at("IterationCutThresholdSurface").as_double();
        else
            iterThresh_ = 1e-07;

        if (surfaceFlowSettings.as_object().find("BisectionIterationThresholdSurface") != surfaceFlowSettings.as_object().end())
            iterThreshBis_ = surfaceFlowSettings.as_object().at("BisectionIterationThresholdSurface").as_double();
        else
            iterThreshBis_ = 1e-07;

        if (surfaceFlowSettings.as_object().find("BisectionMaxIterationsSurface") != surfaceFlowSettings.as_object().end())
            maxIterBis_ = surfaceFlowSettings.as_object().at("BisectionMaxIterationsSurface").as_int();
        else
            maxIterBis_ = 100;

        if (surfaceFlowSettings.as_object().find("BisectionLeftDepthSurface") != surfaceFlowSettings.as_object().end())
            leftDepthBis_ = surfaceFlowSettings.as_object().at("BisectionLeftDepthSurface").as_double();
        else
            leftDepthBis_ = 0.0;
		
        if (surfaceFlowSettings.as_object().find("BisectionRightDepthSurface") != surfaceFlowSettings.as_object().end())
            rightDepthBis_ = surfaceFlowSettings.as_object().at("BisectionRightDepthSurface").as_double();
        else
            rightDepthBis_ = 10.0;
    } else {
        // Fallback defaults if the section is missing.
        enabled_ = false;
        maxIter_ = 100;
        iterThresh_ = 1e-07;
        iterThreshBis_ = 1e-07;
        maxIterBis_ = 100;
        leftDepthBis_ = 0.0;
        rightDepthBis_ = 10.0;
    }
}

/**
 * @brief Executes a simulation step for the surface flow model.
 *
 * This function applies environmental forcing (precipitation and potential PET)
 * to each surface cell, updating the upper storage and overland water depths. It
 * then calls the iterative solver to update the surface water state.
 *
 * @param timeStep The simulation time step (in seconds).
 * @param forcing A ForcingRecord containing the current forcing data.
 */
void SurfaceFlowModel::runStep(double timeStep, const ForcingRecord &forcing) {
    if (!enabled_)
        return;

    // Convert precipitation and potential PET values from mm to m over the forcing interval.
    double precipInput = (forcing.precipitation / weatherDataInterval_ * 0.001) * timeStep;
    double petPot      = (forcing.pet / weatherDataInterval_ * 0.001) * timeStep;

    // Process each cell in the surface mesh.
    for (auto &cellPtr : surfaceMesh_->cells) {
        SurfaceCell* cell = dynamic_cast<SurfaceCell*>(cellPtr.get());
        if (!cell || !cell->waterState || cell->material == 0)
            continue;
        auto state = std::dynamic_pointer_cast<SurfaceWaterState>(cell->waterState);
        if (!state)
            continue;

        double maxUpper   = state->upperStorage;
        double oldUpper   = state->upperStorageDepth;
        double oldOverland = state->waterDepth;

        // Add all precipitation to the upper storage first.
        double newUpper = oldUpper + precipInput;
        double excess   = 0.0;
        if (newUpper > maxUpper) {
            excess   = newUpper - maxUpper;
            newUpper = maxUpper;
        } else {
            excess = 0.0;
        }

        // Update the surface storage.
        state->upperStorageDepth = newUpper;
        double newOverland = oldOverland + excess;

        // Calculate potential evaporation from the cell.
        double totalEvapPot = state->evaporationFraction * petPot; // in meters

        // Evaporate from the upper storage.
        double evapFromUpper = 0.0;
        if (state->upperStorageDepth >= totalEvapPot) {
            evapFromUpper = totalEvapPot;
            state->upperStorageDepth -= totalEvapPot;
            totalEvapPot = 0.0;
        } else {
            evapFromUpper = state->upperStorageDepth;
            totalEvapPot -= state->upperStorageDepth;
            state->upperStorageDepth = 0.0;
        }

        // Evaporate the remaining potential from the overland water.
        double evapFromOverland = 0.0;
        if (newOverland >= totalEvapPot) {
            evapFromOverland = totalEvapPot;
            newOverland -= totalEvapPot;
            totalEvapPot = 0.0;
        } else {
            evapFromOverland = newOverland;
            newOverland = 0.0;
            totalEvapPot = 0.0;
        }

        // Save cumulative water balance results.
        double area = cell->computeArea(surfaceMesh_->points);
        state->cumPetVolume          += petPot * area;
        state->cumEvaporationVolume  += (evapFromUpper + evapFromOverland) * area;
        state->cumPrecipitationVolume += precipInput * area;

        // Save the new water depth state.
        state->waterDepth    = newOverland;
        state->waterDepthOld = state->waterDepth;
    }

    // Execute the iterative solver.
    iterate(timeStep);
}

/**
 * @brief Iteratively solves for the updated water depth in each surface cell.
 *
 * This function uses an OpenMP-parallelized loop to adjust the water depth
 * in each cell by applying a bisection method. The iteration continues until
 * the maximum deviation across all cells falls below a threshold or the maximum
 * number of iterations is reached.
 *
 * @param timeStep The simulation time step (in seconds).
 */
void SurfaceFlowModel::iterate(double timeStep) {
    double iterDev;
    int iterCount = 0;

    do {
        iterDev = 0.0;

        // Parallelize the loop over all cells.
        #pragma omp parallel
        {
            double threadMaxDev = 0.0;

            #pragma omp for nowait
            for (size_t i = 0; i < surfaceMesh_->cells.size(); i++) {
                auto& cellPtr = surfaceMesh_->cells[i];
                SurfaceCell* cell = dynamic_cast<SurfaceCell*>(cellPtr.get());
                if (!cell || cell->material == 0 || !cell->waterState)
                    continue;

                double newWaterDepth = bisection(timeStep, cell, leftDepthBis_, rightDepthBis_);
                auto surfWater = std::dynamic_pointer_cast<SurfaceWaterState>(cell->waterState);
                if (!surfWater)
                    continue;

                double dev = std::fabs(surfWater->waterDepth - newWaterDepth);
                threadMaxDev = std::max(threadMaxDev, dev);
                surfWater->waterDepth = newWaterDepth;
            }

            // In a critical section, update the global deviation.
            #pragma omp critical
            {
                iterDev = std::max(iterDev, threadMaxDev);
            }
        }
        iterCount++;
    } while (iterDev > iterThresh_ && iterCount < maxIter_);

    std::cout << "Surface flow iteration completed in " << iterCount 
              << " steps with final deviation: " << iterDev << std::endl;
}

/**
 * @brief Solves for water depth using a bisection method.
 *
 * This method adjusts the water depth between the given left and right bounds until
 * the residual (computed via calcResidual) is minimized. The final water depth is returned.
 *
 * @param timeStep The simulation time step (in seconds).
 * @param cell Pointer to the surface cell being updated.
 * @param left The lower bound for water depth.
 * @param right The upper bound for water depth.
 * @return double The computed water depth.
 */
double SurfaceFlowModel::bisection(double timeStep, SurfaceCell* cell, double left, double right) {
    int iterCountBis = 0;

    while (std::fabs(right - left) > 2.0 * iterThreshBis_ && iterCountBis < maxIterBis_) {
        double midpoint = (right + left) * 0.5;
        if (calcResidual(timeStep, cell, left) * calcResidual(timeStep, cell, midpoint) > 0.0) {
            left = midpoint;
        } else {
            right = midpoint;
        }
        iterCountBis++;
    }
    return (right + left) * 0.5;
}

/**
 * @brief Computes the residual of the water balance in a surface cell.
 *
 * The residual is based on the difference in water volume due to a change in water depth
 * (accounting for depression storage) and the fluxes exchanged with neighboring cells.
 *
 * @param timeStep The simulation time step (in seconds).
 * @param cell Pointer to the surface cell.
 * @param waterDepth The trial water depth.
 * @return double The water balance residual (m^3).
 */
double SurfaceFlowModel::calcResidual(double timeStep, SurfaceCell* cell, double waterDepth) {
    auto waterState = std::dynamic_pointer_cast<SurfaceWaterState>(cell->waterState);
    if (!waterState)
        return 0.0;

    double area = cell->computeArea(surfaceMesh_->points);
    std::vector<double> sideLengths = cell->computeSideLengths(surfaceMesh_->points);
    double waterVolume = (waterDepth - waterState->waterDepthOld) * area;
    int material = cell->material;

    if (material != 0) {
        double dischIn = 0.0;
        double dischOut = 0.0;

        // Adjust water depth by subtracting depression storage.
        double adjustedDepth = waterDepth - waterState->depressionStorage;
        if (adjustedDepth < 0.0)
            adjustedDepth = 0.0;

        // Compute hydraulic head in the cell.
        Point3D center = cell->computeCenter(surfaceMesh_->points);
        double elevation = center.z;
        double head = elevation + adjustedDepth;

        // Loop over neighbors (using neighborIndices stored in the cell).
        for (size_t i = 0; i < cell->neighborIndices.size(); i++) {
            int neighIndex = cell->neighborIndices[i];
            if (neighIndex < 0 || neighIndex >= static_cast<int>(surfaceMesh_->cells.size()))
                continue;
            SurfaceCell* neighCell = dynamic_cast<SurfaceCell*>(surfaceMesh_->cells[neighIndex].get());
            if (!neighCell || !neighCell->waterState)
                continue;
            auto neighWater = std::dynamic_pointer_cast<SurfaceWaterState>(neighCell->waterState);
            if (!neighWater)
                continue;

            // Do not allow water to move into undefined cells or between roofs and ground level.
            int neighMaterial = neighCell->material;
            if (neighMaterial == 0 || (material == 2 && neighMaterial != 2) ||
                (material != 2 && neighMaterial == 2))
                continue;

            // Get neighbor cell center and elevation.
            Point3D neighCenter = neighCell->computeCenter(surfaceMesh_->points);
            double elevationNeigh = neighCenter.z;

            // Adjust neighbor water depth by its depression storage.
            double neighDepthAdj = neighWater->waterDepth - neighWater->depressionStorage;
            if (neighDepthAdj < 0.0)
                neighDepthAdj = 0.0;
            double headNeigh = elevationNeigh + neighDepthAdj;

            // Compute slope between cells.
            double dx = neighCenter.x - center.x;
            double dy = neighCenter.y - center.y;
            double distance = std::sqrt(dx * dx + dy * dy);

            double slope = 0.0;
            if (distance > 0.0) {
                slope = (headNeigh - head) / distance;
            }

            // Determine effective water depth for flow.
            double waterDepthEff = 0.0;
            if (head > headNeigh) {
                if ((elevation - elevationNeigh) >= neighDepthAdj) {
                    waterDepthEff = adjustedDepth;
                } else {
                    waterDepthEff = head - headNeigh;
                }
            } else if (head < headNeigh) {
                if ((elevationNeigh - elevation) >= adjustedDepth) {
                    waterDepthEff = neighDepthAdj;
                } else {
                    waterDepthEff = headNeigh - head;
                }
            }

            // Calculate effective Manning's n as the arithmetic mean.
            double mannings = 0.5 * (waterState->mannings + neighWater->mannings);
            double velocity = 0.0;
            if (waterDepthEff > 0.0 && mannings > 0.0) {
                if (slope < 0.0)
                    velocity = 1.0 / mannings * std::sqrt(-slope) * std::pow(waterDepthEff, 2.0 / 3.0);
                else if (slope > 0.0)
                    velocity = -1.0 / mannings * std::sqrt(slope) * std::pow(waterDepthEff, 2.0 / 3.0);
            }

            // Compute flux between cells.
            double sideLength = sideLengths[i]; // assumed correct side length
            double flux = waterDepthEff * sideLength * velocity;

            // Update water volume based on flux and time step.
            waterVolume += flux * timeStep;
            if (flux < 0.0)
                dischOut -= flux;
            else
                dischIn += flux;
        }
    }

    return waterVolume;
}

/**
 * @brief Returns the cumulative precipitation water volume.
 *
 * @return double Total precipitation water volume (in m^3).
 */
double SurfaceFlowModel::getPrecipitationWaterVolume() const {
    double waterVolume = 0.0;
    for (auto &cellPtr : surfaceMesh_->cells) {
        SurfaceCell* cell = dynamic_cast<SurfaceCell*>(cellPtr.get());
        if (!cell || !cell->waterState || cell->material == 0)
            continue;
        auto sws = std::dynamic_pointer_cast<SurfaceWaterState>(cell->waterState);
        if (!sws)
            continue;
        waterVolume += sws->cumPrecipitationVolume;
    }
    return waterVolume;
}

/**
 * @brief Returns the cumulative potential PET water volume.
 *
 * @return double Total PET water volume (in m^3).
 */
double SurfaceFlowModel::getPetWaterVolume() const {
    double waterVolume = 0.0;
    for (auto &cellPtr : surfaceMesh_->cells) {
        SurfaceCell* cell = dynamic_cast<SurfaceCell*>(cellPtr.get());
        if (!cell || !cell->waterState || cell->material == 0)
            continue;
        auto sws = std::dynamic_pointer_cast<SurfaceWaterState>(cell->waterState);
        if (!sws)
            continue;
        waterVolume += sws->cumPetVolume;
    }
    return waterVolume;
}

/**
 * @brief Computes the total overland water volume based on cell water depths.
 *
 * The overland water volume is computed as the product of the water depth and the
 * cell area, summed over all cells.
 *
 * @return double The total overland water volume (in m^3).
 */
double SurfaceFlowModel::getOverlandWaterVolume() const {
    double waterVolume = 0.0;
    for (auto &cellPtr : surfaceMesh_->cells) {
        SurfaceCell* cell = dynamic_cast<SurfaceCell*>(cellPtr.get());
        if (!cell || !cell->waterState || cell->material == 0)
            continue;
        auto sws = std::dynamic_pointer_cast<SurfaceWaterState>(cell->waterState);
        if (!sws)
            continue;
        double area = cell->computeArea(surfaceMesh_->points);
        waterVolume += area * sws->waterDepth;
    }
    return waterVolume;
}

/**
 * @brief Returns the total upper storage water volume.
 *
 * The upper storage water volume is computed as the product of the upper storage depth
 * and the cell area, summed over all cells.
 *
 * @return double The total upper storage water volume (in m^3).
 */
double SurfaceFlowModel::getUpperStorageWaterVolume() const {
    double waterVolume = 0.0;
    for (auto &cellPtr : surfaceMesh_->cells) {
        SurfaceCell* cell = dynamic_cast<SurfaceCell*>(cellPtr.get());
        if (!cell || !cell->waterState || cell->material == 0)
            continue;
        auto sws = std::dynamic_pointer_cast<SurfaceWaterState>(cell->waterState);
        if (!sws)
            continue;
        double area = cell->computeArea(surfaceMesh_->points);
        waterVolume += area * sws->upperStorageDepth;
    }
    return waterVolume;
}

/**
 * @brief Returns the cumulative evaporation water volume.
 *
 * This value represents the total water volume lost to evaporation over time.
 *
 * @return double The cumulative evaporation water volume (in m^3).
 */
double SurfaceFlowModel::getEvaporationWaterVolume() const {
    double waterVolume = 0.0;
    for (auto &cellPtr : surfaceMesh_->cells) {
        SurfaceCell* cell = dynamic_cast<SurfaceCell*>(cellPtr.get());
        if (!cell || !cell->waterState || cell->material == 0)
            continue;
        auto sws = std::dynamic_pointer_cast<SurfaceWaterState>(cell->waterState);
        if (!sws)
            continue;
        waterVolume += sws->cumEvaporationVolume;
    }
    return waterVolume;
}
