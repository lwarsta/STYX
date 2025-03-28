/**
 * @file Simulation.h
 * @brief Declaration of the Simulation class.
 *
 * The Simulation class encapsulates the hydrological simulation by managing the surface,
 * subsurface, and network flow models, as well as the mesh coupler and output manager.
 */

#ifndef SIMULATION_H
#define SIMULATION_H

#include "json.h"
#include "SurfaceSnowModel.h"
#include "SurfaceFlowModel.h"
#include "SubsurfaceFlowModel.h"
#include "SubsurfaceHeatTransferModel.h"
#include "NetworkFlowModel.h"
#include "ForcingData.h"
#include "IModel.h"
#include "MeshCoupler.h"
#include "OutputManager.h"

/**
 * @brief The Simulation class encapsulates the hydrological simulation.
 *
 * This class initializes the meshes, models, forcing data, and output manager.
 * It then runs the simulation loop, processes each time step, and writes output results.
 */
class Simulation {
public:
    /**
     * @brief Constructs a Simulation instance with the given JSON settings.
     * 
     * @param settings The parsed JSON settings.
     */
    explicit Simulation(const JsonValue &settings);

    /**
     * @brief Runs the simulation.
     */
    void run();

private:
    /// Parsed JSON settings.
    JsonValue settings_;

    /// Unique pointer to the surface mesh.
    std::unique_ptr<SurfaceMesh> surfaceMesh_;
    /// Unique pointer to the subsurface mesh.
    std::unique_ptr<SubsurfaceMesh> subsurfaceMesh_;
    /// Unique pointer to the network mesh.
    std::unique_ptr<NetworkMesh> networkMesh_;
    /// Forcing data (e.g., precipitation, PET, temperature).
    ForcingData forcingData_;
    /// Unique pointer to the mesh coupler.
    std::unique_ptr<MeshCoupler> meshCoupler_;
    /// Vector of unique pointers to simulation models.
    std::vector<std::unique_ptr<IModel>> models_;
    /// Unique pointer to the output manager.
    std::unique_ptr<OutputManager> outputManager_;

    /// Total simulation length in seconds.
    double simulationLength_ = 0;
    /// Time step (seconds).
    double timeStep_ = 0;
    /// Interval (in seconds) at which simulation results are saved.
    double resultsSaveInterval_ = 0;
    /// Current simulation time (seconds).
    double currentTime_ = 0;
    /// Simulation start time (std::time_t).
    std::time_t startTime_ = 0;
    /// Excel serial date corresponding to the simulation start time.
    double startTimeExcel_ = 0.0;

    /**
     * @brief Initializes the simulation components.
     */
    void initialize();

    /**
     * @brief Processes one simulation time step.
     */
    void processTimeStep();

    /**
     * @brief Outputs simulation results.
     *
     * The simulation time (in seconds) is converted to an Excel serial date,
     * and various water volume balances are computed and written.
     * 
     * @param currentTime The current simulation time.
     */
    void outputResults(double currentTime);
	
	/**
	 * @brief Recursively prints a JsonValue.
	 *
	 * If the JsonValue is an object, it iterates over all key/value pairs;
	 * if it is an array, it prints each element. For basic types it prints
	 * the value.
	 *
	 * @param j The JsonValue to print.
	 * @param indent The current indentation level (used for pretty-printing).
	 */
	void printJson(const JsonValue &j, int indent = 0);
};

#endif // SIMULATION_H
