/**
 * @file SurfaceFlowModel.h
 * @brief Declaration of the SurfaceFlowModel class.
 *
 * This file declares the SurfaceFlowModel class which implements the IModel interface
 * for simulating surface water flow using a SurfaceMesh. It also provides functions to
 * query various water volume components.
 */

#ifndef SURFACE_FLOW_MODEL_H
#define SURFACE_FLOW_MODEL_H

#include "IModel.h"
#include "SurfaceMesh.h"
#include "ForcingData.h"
#include "MeshCoupler.h"
#include <iostream>
#include <typeinfo>

/**
 * @brief The SurfaceFlowModel class implements a simulation model for surface flow.
 *
 * The SurfaceFlowModel uses a SurfaceMesh to represent the domain and employs iterative
 * solvers to compute water depth changes over a simulation time step. It also supports
 * exchanging water with other meshes via a MeshCoupler.
 */
class SurfaceFlowModel : public IModel {
public:
    /**
     * @brief Constructs a SurfaceFlowModel with the specified surface mesh.
     * @param mesh Pointer to a SurfaceMesh that represents the 2D surface domain.
     */
    explicit SurfaceFlowModel(SurfaceMesh* mesh)
        : surfaceMesh_(mesh),
          meshCoupler_(nullptr)
    {}

    /**
     * @brief Initializes the model with settings from a JSON object.
     * @param settings A JSON object containing the configuration settings.
     */
    virtual void initialize(const JsonValue &settings) override;

    /**
     * @brief Executes a single simulation step for the surface flow model.
     * @param timeStep The time step duration in seconds.
     * @param forcing A ForcingRecord containing the environmental forcing data.
     */
    virtual void runStep(double timeStep, const ForcingRecord &forcing) override;

    /**
     * @brief Sets the MeshCoupler to be used for water exchange.
     * @param coupler Pointer to a MeshCoupler instance.
     */
    void setMeshCoupler(MeshCoupler* coupler) { meshCoupler_ = coupler; }

    /**
     * @brief Returns the precipitation water volume on the surface.
     * @return The precipitation water volume in cubic meters.
     */
    double getPrecipitationWaterVolume() const;

    /**
     * @brief Returns the potential PET water volume on the surface.
     * @return The potential PET water volume in cubic meters.
     */
    double getPetWaterVolume() const;

    /**
     * @brief Returns the overland water volume on the surface.
     * @return The overland water volume in cubic meters.
     */
    double getOverlandWaterVolume() const;

    /**
     * @brief Returns the upper storage water volume on the surface.
     * @return The upper storage water volume in cubic meters.
     */
    double getUpperStorageWaterVolume() const;

    /**
     * @brief Returns the evaporation water volume from the surface.
     * @return The evaporation water volume in cubic meters.
     */
    double getEvaporationWaterVolume() const;

private:
    SurfaceMesh* surfaceMesh_;      ///< Pointer to the surface mesh.
    MeshCoupler* meshCoupler_;      ///< Pointer to the mesh coupler for water exchange.

    // Model parameters and variables.
    bool   enabled_;                ///< Flag indicating if the model is enabled.
    double weatherDataInterval_;    ///< Weather data interval (in seconds).
    int    maxIter_;                ///< Maximum number of iterations.
    double iterThresh_;             ///< Convergence threshold for iterations.
    double iterThreshBis_;          ///< Convergence threshold for bisection.
    int    maxIterBis_;             ///< Maximum number of bisection iterations.
    double leftDepthBis_;           ///< Left depth bound for bisection.
    double rightDepthBis_;          ///< Right depth bound for bisection.

    // Iterative solver helper functions.
    /**
     * @brief Iterates the surface flow model for the specified time step.
     * @param timeStep The time step duration in seconds.
     */
    void iterate(double timeStep);

    /**
     * @brief Calculates the residual for the implicit surface flow equation.
     * @param timeStep The time step duration in seconds.
     * @param cell Pointer to the SurfaceCell for which the residual is calculated.
     * @param waterDepth The current water depth in the cell.
     * @return The computed residual value.
     */
    double calcResidual(double timeStep, SurfaceCell* cell, double waterDepth);

    /**
     * @brief Solves for the water depth in a cell using the bisection method.
     * @param timeStep The time step duration in seconds.
     * @param cell Pointer to the SurfaceCell for which the water depth is solved.
     * @param left The left bound for the water depth.
     * @param right The right bound for the water depth.
     * @return The computed water depth.
     */
    double bisection(double timeStep, SurfaceCell* cell, double left, double right);
};

#endif // SURFACE_FLOW_MODEL_H
