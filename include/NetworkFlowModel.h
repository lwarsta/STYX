/**
 * @file NetworkFlowModel.h
 * @brief Declaration of the NetworkFlowModel class.
 *
 * This class implements a network flow model for simulating water flow in network
 * components (junctions and links). It derives from the IModel interface and provides
 * methods for initialization, stepping through the simulation, and computing water volumes.
 */

#ifndef NETWORK_FLOW_MODEL_H
#define NETWORK_FLOW_MODEL_H

#include "IModel.h"
#include "NetworkMesh.h"
#include "MeshCoupler.h"

/**
 * @brief The NetworkFlowModel class.
 *
 * Implements a network flow solver that computes water flows in network elements.
 * It can be coupled with a MeshCoupler for water exchanges with other meshes.
 */
class NetworkFlowModel : public IModel {
public:
    /**
     * @brief Constructor.
     * @param mesh Pointer to the network mesh.
     */
    explicit NetworkFlowModel(NetworkMesh* mesh)
        : networkMesh_(mesh), 
          meshCoupler_(nullptr)
    {}

    /**
     * @brief Initializes the network flow model.
     * 
     * Reads model-specific settings from the provided JSON object.
     *
     * @param settings JSON settings for the model.
     */
    virtual void initialize(const JsonValue &settings) override;

    /**
     * @brief Executes a single simulation step.
     * 
     * Applies the network flow solver for the specified time step using the given forcing data.
     *
     * @param timeStep Simulation time step (seconds).
     * @param forcing Environmental forcing data for the current step.
     */
    virtual void runStep(double timeStep, const ForcingRecord &forcing) override;
	
    /**
     * @brief Returns the total water volume in the network.
     * @return Total network water volume in cubic meters.
     */
    double getNetworkWaterVolume();

    /**
     * @brief Returns the water volume discharged at outfall junctions.
     * @return Outfall water volume in cubic meters.
     */
    double getOutfallWaterVolume();

    /**
     * @brief Sets the MeshCoupler used for water exchange.
     *
     * @param coupler Pointer to a MeshCoupler object.
     */
    void setMeshCoupler(MeshCoupler* coupler) {
        meshCoupler_ = coupler;
    }

private:
    NetworkMesh* networkMesh_;    ///< Pointer to the network mesh.
    MeshCoupler* meshCoupler_;    ///< Pointer to the mesh coupler for water exchanges.

    // Model parameters and variables.
    bool enabled_;                ///< Flag indicating whether the network model is enabled.
    double weatherDataInterval_;  ///< Interval for weather data (seconds).
    double iterCutThresh_;        ///< Convergence threshold for the iterative solver.
    int maxIter_;                 ///< Maximum number of iterations.
    double iterCutThreshBis_;     ///< Convergence threshold for bisection iterations.
    int maxIterBis_;              ///< Maximum iterations for bisection.
    double leftDepthBis_;         ///< Left bound for bisection.
    double rightDepthBis_;        ///< Right bound for bisection.

    // Iterative solver functions.
    /**
     * @brief Iteratively solves the network flow equations.
     *
     * @param timeStep Simulation time step (seconds).
     */
    void iterate(double timeStep);

    /**
     * @brief Performs postprocessing on the network flow results.
     */
    void postprocess();

    // Helper functions.
    /**
     * @brief Reverts water levels to their previous values.
     */
    void revertHeads();

    /**
     * @brief Swaps current water levels with old water levels.
     */
    void swapHeads();

    /**
     * @brief Calculates the water depth in a link given a target area and the link diameter.
     *
     * @param target_area The target cross-sectional area (m²).
     * @param link_diameter The diameter of the link (m).
     * @return Calculated water depth (m).
     */
    double calcWaterDepthInLink(double target_area, double link_diameter);

    /**
     * @brief Computes the flow area and hydraulic radius for a link.
     *
     * @param waterDepth Water depth in the link (m).
     * @param linkDiameter Link diameter (m).
     * @param flowArea Output parameter for the computed flow area (m²).
     * @param hydrRad Output parameter for the computed hydraulic radius (m).
     */
    void getFlowAreaAndHydrRad(double waterDepth, double linkDiameter, double &flowArea, double &hydrRad);
};

#endif // NETWORK_FLOW_MODEL_H
