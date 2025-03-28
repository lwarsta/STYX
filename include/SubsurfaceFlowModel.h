/**
 * @file SubsurfaceFlowModel.h
 * @brief Declaration of the SubsurfaceFlowModel class for subsurface flow simulation.
 *
 * This file contains the declaration of the SubsurfaceFlowModel class that implements
 * the IModel interface for simulating water flow in the subsurface domain. It uses a
 * SubsurfaceMesh and a MeshCoupler to perform the computations.
 */

#ifndef SUBSURFACE_FLOW_MODEL_H
#define SUBSURFACE_FLOW_MODEL_H

#include "IModel.h"
#include "SubsurfaceMesh.h"
#include "ForcingData.h"
#include "MeshCoupler.h"
#include <nlohmann/json.hpp>
#include <memory>
#include <algorithm>

/**
 * @brief Class for simulating subsurface water flow.
 *
 * The SubsurfaceFlowModel class implements the IModel interface. It uses a SubsurfaceMesh
 * to represent the geometry and state of the subsurface and a MeshCoupler to manage exchanges
 * with other domains. The model is configured via JSON settings and processes forcing data
 * during each simulation time step.
 */
class SubsurfaceFlowModel : public IModel {
public:
    /**
     * @brief Constructs a SubsurfaceFlowModel.
     * @param mesh Pointer to a SubsurfaceMesh representing the subsurface.
     */
    explicit SubsurfaceFlowModel(SubsurfaceMesh* mesh)
        : subsurfaceMesh_(mesh),
          meshCoupler_(nullptr)
    {}

    /**
     * @brief Initializes the subsurface flow model with the given settings.
     * @param settings A JSON object containing simulation settings.
     */
    virtual void initialize(const JsonValue &settings) override;

    /**
     * @brief Executes a single simulation step.
     * @param timeStep The simulation time step (in seconds).
     * @param forcing Forcing data record to be applied during the time step.
     */
    virtual void runStep(double timeStep, const ForcingRecord &forcing) override;
    
    /**
     * @brief Sets the MeshCoupler used for exchange with other domains.
     * @param coupler Pointer to the MeshCoupler object.
     */
    void setMeshCoupler(MeshCoupler* coupler) { meshCoupler_ = coupler; }

    /**
     * @brief Returns the computed subsurface water volume.
     * @return Total subsurface water volume.
     */
    double getSubsurfaceWaterVolume() const;

    /**
     * @brief Returns the infiltration water volume in the subsurface.
     * @return Infiltration water volume.
     */
    double getInfiltrationWaterVolume() const;

    /**
     * @brief Returns the transpiration water volume computed from the subsurface.
     * @return Transpiration water volume.
     */
    double getTranspirationWaterVolume() const;

private:
    SubsurfaceMesh* subsurfaceMesh_;   ///< Pointer to the subsurface mesh.
    MeshCoupler* meshCoupler_;         ///< Pointer to the MeshCoupler for grid exchanges.
    
    // Model parameters and variables.
    bool enabled_;                     ///< Flag indicating if the model is enabled.
    double weatherDataInterval_;       ///< Weather data interval (seconds).
    int maxIter_;                      ///< Maximum iterations for the solver.
    double iterCutThresh_;             ///< Convergence threshold for the solver.

    // Iterative solver helper functions.
    /**
     * @brief Preprocesses the subsurface mesh prior to solving.
     * @param timeStep The simulation time step.
     * @param forcing Forcing data record.
     */
    void preprocess(double timeStep, const ForcingRecord &forcing);

    /**
     * @brief Iteratively solves the subsurface flow equations.
     * @param timeStep The simulation time step.
     */
    void iterate(double timeStep);
    
    // Helper functions implementing subsurface physics.
    /**
     * @brief Computes the moisture factor using van Genuchten parameters.
     * @param presHead The pressure head [m].
     * @param phMin The minimum pressure head [m].
     * @param phMax The maximum pressure head [m].
     * @param phWilt The wilting pressure head [m].
     * @return The moisture factor.
     */
    double computeMoistureFactor(double presHead, double phMin, double phMax, double phWilt);

    /**
     * @brief Applies transpiration effects to a vertical column of subsurface cells.
     * @param topCell Pointer to the top cell of the column.
     * @param timeStep The simulation time step.
     * @param pet_m Potential evapotranspiration (m).
     */
    void applyTranspirationToColumn(VolumeCell* topCell, double timeStep, double pet_m);

    /**
     * @brief Applies transpiration to the subsurface mesh.
     * @param timeStep The simulation time step.
     * @param pet_m Potential evapotranspiration (m).
     */
    void applyTranspiration(double timeStep, double pet_m);

    /**
     * @brief Calculates the unsaturated hydraulic conductivity for a cell.
     * @param cell Pointer to the VolumeCell.
     */
    void calcUnsatCond(VolumeCell* cell);

    /**
     * @brief Calculates intercell conductivity and flux properties.
     * @param cell Pointer to the VolumeCell.
     */
    void calcCondInter(VolumeCell* cell);

    /**
     * @brief Calculates the differential water capacity of a cell.
     * @param cell Pointer to the VolumeCell.
     */
    void calcDiffWatCap(VolumeCell* cell);

    /**
     * @brief Calculates the fluxes between a cell and its neighbors.
     * @param cell Pointer to the VolumeCell.
     */
    void calcFluxes(VolumeCell* cell);

    /**
     * @brief Calculates the flow velocities between a cell and its neighbors.
     * @param cell Pointer to the VolumeCell.
     */
    void calcFlowVel(VolumeCell* cell);
};

#endif // SUBSURFACE_FLOW_MODEL_H
