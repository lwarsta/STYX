#ifndef MESH_COUPLER_H
#define MESH_COUPLER_H

#include "SurfaceMesh.h"
#include "SubsurfaceMesh.h"
#include "NetworkMesh.h"

/**
 * @brief Handles water exchange computations between surface, subsurface, and network meshes.
 */
class MeshCoupler {
public:
    /**
     * @brief Constructs a MeshCoupler instance.
     * @param surface Pointer to the surface mesh.
     * @param subsurface Pointer to the subsurface mesh.
     * @param network Pointer to the network mesh.
     */
    MeshCoupler(SurfaceMesh* surface, SubsurfaceMesh* subsurface, NetworkMesh* network)
        : surfaceMesh(surface), 
          subsurfaceMesh(subsurface), 
          networkMesh(network)
    {}

    /**
     * @brief Pre-computation of water exchange fluxes between surface and network meshes.
     * @param timeStep Simulation timestep.
     */
    void preExchangeSurfaceNetwork(double timeStep);

    /**
     * @brief Post-computation adjustment after water exchange between surface and network meshes.
     * @param timeStep Simulation timestep.
     */
    void postExchangeSurfaceNetwork(double timeStep);

    /**
     * @brief Pre-computation of water exchange fluxes between surface and subsurface meshes.
     * @param timeStep Simulation timestep.
     */
    void preExchangeSurfaceSubsurface(double timeStep);

    /**
     * @brief Post-computation adjustment after water exchange between surface and subsurface meshes.
     * @param timeStep Simulation timestep.
     */
    void postExchangeSurfaceSubsurface(double timeStep);

    /**
     * @brief Static helper to compute water content from pressure head using van Genuchten parameters.
     *
     * Computes water content for a given subsurface cell based on pressure head, saturated water content,
     * residual water content, van Genuchten alpha, and van Genuchten N parameters.
     *
     * @param cell Volume cell pointer for which water content is calculated.
     */
    static void calcWatCont(VolumeCell* cell) {
        auto state = std::dynamic_pointer_cast<SubsurfaceWaterState>(cell->waterState);
        if (!state) return;

        if (state->presHead >= 0.0) {
            state->watCont = state->saturatedWaterContent + state->presHead * state->compressibility;
        } else {
            double tmp = 1.0 + std::pow(std::fabs(state->vanGenuchtenAlpha * state->presHead), state->vanGenuchtenN);
            state->watCont = state->residualWaterContent +
                             (state->saturatedWaterContent - state->residualWaterContent) / std::pow(tmp, state->vanGenuchtenM);
        }
    }

private:
    SurfaceMesh* surfaceMesh;         ///< Pointer to surface mesh
    SubsurfaceMesh* subsurfaceMesh;   ///< Pointer to subsurface mesh
    NetworkMesh* networkMesh;         ///< Pointer to network mesh

    /**
     * @brief Computes water flow from a surface cell to a network junction.
     * @param junction Target junction pointer.
     * @param cellSurf Source surface cell pointer.
     * @param timeStep Simulation timestep.
     */
    void compFlowFromCellToJunc(NetworkJunctionCell* junction, SurfaceCell* cellSurf, double timeStep);

    /**
     * @brief Updates the water content of a subsurface cell ensuring physical constraints.
     * @param cell Volume cell pointer.
     * @param watContNew Desired new water content value.
     * @return The updated water content after applying constraints.
     */
    double changeWatCont(VolumeCell* cell, double watContNew);
};

#endif // MESH_COUPLER_H
