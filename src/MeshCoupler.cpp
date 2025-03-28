#include "MeshCoupler.h"
#include <iostream>

/**
 * @brief Computes water exchange fluxes between surface and subsurface meshes before the main simulation step.
 *
 * This method calculates infiltration (surface to subsurface) and exfiltration (subsurface to surface)
 * fluxes based on hydraulic conditions. Fluxes are constrained by available water volume on each side,
 * maintaining physical realism.
 *
 * @param timeStep Simulation timestep in seconds.
 */
void MeshCoupler::preExchangeSurfaceSubsurface(double timeStep) {
    for (auto& cellSubPtr : subsurfaceMesh->cells) {
        auto cellSub = dynamic_cast<VolumeCell*>(cellSubPtr.get());

        // Validate subsurface cell.
        if (!cellSub || cellSub->material == 0 ||
            cellSub->gridConnection < 0 ||
            static_cast<size_t>(cellSub->gridConnection) >= surfaceMesh->cells.size() ||
            !cellSub->waterState) {
            continue;
        }

        auto stateSub = std::dynamic_pointer_cast<SubsurfaceWaterState>(cellSub->waterState);
        int idCellSurf = cellSub->gridConnection;

        auto cellSurf = dynamic_cast<SurfaceCell*>(surfaceMesh->cells[idCellSurf].get());

        // Validate surface cell.
        if (!cellSurf || cellSurf->material == 0 || !cellSurf->waterState) {
            continue;
        }

        auto stateSurf = std::dynamic_pointer_cast<SurfaceWaterState>(cellSurf->waterState);

        // Calculate hydraulic head of surface water.
        Point3D centerPointSurf = cellSurf->computeCenter(surfaceMesh->points);
        double hydrHeadSurf = stateSurf->waterDepth + centerPointSurf.z;

        // Geometric and hydraulic parameters.
        double distance = 0.5 * cellSub->computeThickness(subsurfaceMesh->points);  ///< Half-thickness, represents vertical distance between mesh centers.
        double areaSurf = cellSurf->computeArea(surfaceMesh->points);
        double volumeSub = cellSub->computeVolume(subsurfaceMesh->points);
        double fluxExchange = 0.0;

        // Handle infiltration from surface to subsurface.
        if (stateSub->presHead < 0.0 && hydrHeadSurf > stateSub->hydraulicHead) {
            // Compute infiltration flux based on Darcy's law.
            fluxExchange = stateSub->condUnsat * areaSurf * (hydrHeadSurf - stateSub->hydraulicHead) / distance;

            // Limit flux to available subsurface storage volume.
            double emptyVolume = std::max(0.0, (stateSub->saturatedWaterContent - stateSub->watCont) * volumeSub);
            fluxExchange = std::min(fluxExchange, emptyVolume / timeStep);

            // Limit flux to available surface water volume.
            double availableSurfaceWater = stateSurf->waterDepth * areaSurf;
            fluxExchange = std::min(fluxExchange, availableSurfaceWater / timeStep);

            // Update subsurface water content with constraints.
            double watContOld = stateSub->watCont;
            double watContNew = (watContOld * volumeSub + fluxExchange * timeStep) / volumeSub;
            double watContResult = changeWatCont(cellSub, watContNew);
            fluxExchange = (watContResult - watContOld) * volumeSub / timeStep;
        }
        // Handle exfiltration from subsurface to surface.
        else if (stateSub->presHead > 0.0 && hydrHeadSurf < stateSub->hydraulicHead) {
            // Compute exfiltration flux based on Darcy's law.
            fluxExchange = stateSub->condUnsat * areaSurf * (hydrHeadSurf - stateSub->hydraulicHead) / distance;

            // Temporarily lower hydraulic head by 1.0 m to estimate freely available water.
            Point3D centerPointSub = cellSub->computeCenter(subsurfaceMesh->points);
            double hydrHeadOriginal = stateSub->hydraulicHead;
            stateSub->hydraulicHead = centerPointSub.z - 1.0; ///< Temporarily lower by 1m to estimate drainable water.
            stateSub->presHead = stateSub->hydraulicHead - centerPointSub.z;

            calcWatCont(cellSub);
            double watContDry = stateSub->watCont;

            // Calculate drainable water volume.
            double soilWatVol = std::max(0.0, (stateSub->watContOld - watContDry) * volumeSub);

            // Restore original hydraulic conditions.
            stateSub->hydraulicHead = hydrHeadOriginal;
            stateSub->presHead = stateSub->hydraulicHead - centerPointSub.z;
            calcWatCont(cellSub);

            // Limit exfiltration flux to freely available water volume.
            fluxExchange = std::max(fluxExchange, -soilWatVol / timeStep);

            // Update subsurface water content after exfiltration.
            double watContOld = stateSub->watCont;
            double watContNew = (watContOld * volumeSub + fluxExchange * timeStep) / volumeSub;
            double watContResult = changeWatCont(cellSub, watContNew);
            fluxExchange = (watContResult - watContOld) * volumeSub / timeStep;
        }

        // Update state variables and flux tracking.
        stateSub->watContOld = stateSub->watCont;
        stateSub->hydraulicHeadOld = stateSub->hydraulicHead;
        stateSub->diffWatCap = stateSub->compressibility; ///< Retain for compatibility; evaluate necessity.
        stateSub->fluxExchangeSurfSub = fluxExchange;
        stateSub->cumExchangeSurfSub += fluxExchange * timeStep;
    }
}

/**
 * @brief Updates surface water depths after water exchange with the subsurface mesh.
 *
 * This function finalizes the surface water depths based on previously computed water fluxes
 * from the subsurface to surface mesh. It ensures that surface water depths do not become negative.
 *
 * @param timeStep Simulation timestep in seconds.
 */
void MeshCoupler::postExchangeSurfaceSubsurface(double timeStep) {
    for (auto& cellSubPtr : subsurfaceMesh->cells) {
        auto cellSub = dynamic_cast<VolumeCell*>(cellSubPtr.get());

        // Validate subsurface cell.
        if (!cellSub || cellSub->material == 0 ||
            cellSub->gridConnection < 0 ||
            static_cast<size_t>(cellSub->gridConnection) >= surfaceMesh->cells.size() ||
            !cellSub->waterState) {
            continue;
        }

        auto stateSub = std::dynamic_pointer_cast<SubsurfaceWaterState>(cellSub->waterState);
        int idCellSurf = cellSub->gridConnection;

        auto cellSurf = dynamic_cast<SurfaceCell*>(surfaceMesh->cells[idCellSurf].get());

        // Validate surface cell.
        if (!cellSurf || cellSurf->material == 0 || !cellSurf->waterState) {
            continue;
        }

        auto stateSurf = std::dynamic_pointer_cast<SurfaceWaterState>(cellSurf->waterState);

        // Update surface water depth based on subsurface-surface flux exchange.
        double areaSurf = cellSurf->computeArea(surfaceMesh->points);
        double fluxExchange = stateSub->fluxExchangeSurfSub;
        stateSurf->waterDepth -= fluxExchange * timeStep / areaSurf;

        // Prevent negative surface water depth due to numerical precision errors.
        if (stateSurf->waterDepth < 0.0) {
            stateSurf->waterDepth = 0.0;
        }
    }
}

/**
 * @brief Computes water exchange fluxes from surface cells to the network junctions before the main simulation step.
 *
 * This method transfers water from surface cells (including roof cells and cells directly overlying network junctions)
 * to the connected network junctions based on predefined connectivity and hydrological conditions.
 *
 * @param timeStep Simulation timestep in seconds.
 */
void MeshCoupler::preExchangeSurfaceNetwork(double timeStep) {
    std::vector<NetworkJunctionCell*>& junctions = networkMesh->junctions;

    // Transfer water from roof cells to designated junction outlets.
    for (auto& cellSurfPtr : surfaceMesh->cells) {
        auto cellSurf = dynamic_cast<SurfaceCell*>(cellSurfPtr.get());

        // Validate surface cell.
        if (!cellSurf || cellSurf->material == 0 || !cellSurf->waterState) {
            continue;
        }

        auto stateSurf = std::dynamic_pointer_cast<SurfaceWaterState>(cellSurf->waterState);

        // Verify valid outlet ID.
        if (stateSurf->outletId < 0 ||
            static_cast<size_t>(stateSurf->outletId) >= junctions.size()) {
            continue;
        }

        compFlowFromCellToJunc(junctions[stateSurf->outletId], cellSurf, timeStep);
    }

    // Transfer water from surface cells directly overlying junctions.
    for (size_t i = 0; i < junctions.size(); i++) {
        int gridConnection = junctions[i]->gridConnection;

        if (gridConnection >= 0 &&
            static_cast<size_t>(gridConnection) < surfaceMesh->cells.size()) {
            auto cellSurf = dynamic_cast<SurfaceCell*>(surfaceMesh->cells[gridConnection].get());

            compFlowFromCellToJunc(junctions[i], cellSurf, timeStep);
        }
    }
}

/**
 * @brief Updates surface water depths when network junction capacities are exceeded.
 *
 * If water levels in network junctions exceed their depth, excess water volume
 * is transferred back to the connected surface cells. The surface and junction
 * water depths are updated accordingly.
 *
 * @param timeStep Simulation timestep in seconds (currently unused but included for consistency).
 */
void MeshCoupler::postExchangeSurfaceNetwork(double timeStep) {
    (void)timeStep; // Parameter intentionally unused

    std::vector<NetworkJunctionCell*>& junctions = networkMesh->junctions;

    // Transfer excess water from junctions back to surface cells.
    for (size_t i = 0; i < junctions.size(); i++) {
        NetworkJunctionCell* junction = junctions[i];

        // Validate junction and its associated surface cell.
        if (!junction || junction->gridConnection < 0 ||
            static_cast<size_t>(junction->gridConnection) >= surfaceMesh->cells.size() ||
            !junction->waterState) {
            continue;
        }

        auto stateJunc = std::dynamic_pointer_cast<NetworkWaterState>(junction->waterState);
        double junctionWaterLevel = stateJunc->waterLevel;
        double junctionDepth = networkMesh->getJunctionDepth(junction);

        // Only proceed if junction water level exceeds the maximum depth.
        if (junctionWaterLevel > junctionDepth) {
            int surfaceCellId = junction->gridConnection;
            SurfaceCell* cellSurf = dynamic_cast<SurfaceCell*>(surfaceMesh->cells[surfaceCellId].get());

            if (!cellSurf || cellSurf->material == 0 || !cellSurf->waterState) {
                continue;
            }

            auto stateSurf = std::dynamic_pointer_cast<SurfaceWaterState>(cellSurf->waterState);
            double surfaceWaterDepth = stateSurf->waterDepth;
            double surfaceArea = cellSurf->computeArea(surfaceMesh->points);
            double surfaceWaterVolume = surfaceWaterDepth * surfaceArea;

            double junctionArea = networkMesh->getJunctionArea(junction->diameter);
            double excessWaterVolume = (junctionWaterLevel - junctionDepth) * junctionArea;

            // Update water volumes and depths.
            junctionWaterLevel = junctionDepth;
            stateJunc->waterLevel = junctionWaterLevel;
            stateJunc->waterLevelOld = junctionWaterLevel;

            double newSurfaceWaterVolume = surfaceWaterVolume + excessWaterVolume;
            stateSurf->waterDepth = newSurfaceWaterVolume / surfaceArea;
            stateSurf->waterDepthOld = stateSurf->waterDepth;
        }
    }
}

/**
 * @brief Computes water flux from a surface cell into a connected network junction.
 *
 * This method calculates water exchange based on surface cell water depth, slope,
 * and Manning's coefficient. The transfer is constrained by available surface water,
 * depression storage, and junction capacity.
 *
 * @param junction   Pointer to the network junction cell receiving water.
 * @param cellSurf   Pointer to the surface cell providing water.
 * @param timeStep   Simulation timestep in seconds.
 */
void MeshCoupler::compFlowFromCellToJunc(NetworkJunctionCell* junction, SurfaceCell* cellSurf, double timeStep) {
    if (!junction || !cellSurf || !junction->waterState || !cellSurf->waterState) {
        return;
    }

    auto stateJunc = std::dynamic_pointer_cast<NetworkWaterState>(junction->waterState);
    auto stateSurf = std::dynamic_pointer_cast<SurfaceWaterState>(cellSurf->waterState);

    if (!stateJunc || !stateSurf) {
        return;
    }

    // Calculate water available beyond depression storage.
    double availableWaterDepth = std::max(0.0, stateSurf->waterDepth - stateSurf->depressionStorage);

    // Compute geometry and slope parameters.
    double junctionArea = networkMesh->getJunctionArea(junction->diameter);
    double cellArea = cellSurf->computeArea(surfaceMesh->points);
    double cellSlope = cellSurf->computeAverageSlope(surfaceMesh->points, surfaceMesh->cells);
    double junctionWaterLevel = stateJunc->waterLevel;
    double junctionDepth = networkMesh->getJunctionDepth(junction);

    if (junctionDepth <= 0.0 || stateSurf->mannings <= 0.0 || junctionArea <= 0.0) {
        return;
    }

    // Calculate flow velocity and potential volume transferred.
    double cellVelocity = (1.0 / stateSurf->mannings) *
                          std::sqrt(cellSlope) *
                          std::pow(availableWaterDepth, 2.0 / 3.0);

    double potentialWaterVolume = sqrt(cellArea) * availableWaterDepth * cellVelocity * timeStep;

    // Limit actual water transfer based on constraints.
    double actualWaterVolume = std::min({
        potentialWaterVolume,
        (availableWaterDepth - stateSurf->depressionStorage) * cellArea,
        (junctionDepth - junctionWaterLevel) * junctionArea
    });

    // Update water depths based on the actual volume transferred.
    if (potentialWaterVolume > 0.0 && actualWaterVolume > 0.0) {
        stateSurf->waterDepth -= actualWaterVolume / cellArea;
        stateSurf->waterDepth = std::max(stateSurf->waterDepth, 0.0);
        stateSurf->waterDepthOld = stateSurf->waterDepth;

        stateJunc->waterLevel += actualWaterVolume / junctionArea;
        stateJunc->waterLevelOld = stateJunc->waterLevel;
    }
}

/**
 * @brief Adjusts soil water content within realistic physical constraints.
 *
 * This method ensures that soil water content remains between residual and saturated
 * water contents. It recalculates pressure and hydraulic head accordingly.
 *
 * @param cell       Pointer to the subsurface volume cell to update.
 * @param watContNew Proposed new water content value.
 * @return Adjusted water content value.
 */
double MeshCoupler::changeWatCont(VolumeCell* cell, double watContNew) {
    auto state = std::dynamic_pointer_cast<SubsurfaceWaterState>(cell->waterState);

    if (!state || state->residualWaterContent <= 0.0 || state->saturatedWaterContent <= state->residualWaterContent ||
        state->vanGenuchtenN <= 1.0 || state->vanGenuchtenAlpha <= 0.0) {
        std::cerr << "-> Error in MeshCoupler::changeWatCont (invalid state parameters)." << std::endl;
        return state ? state->watCont : 0.0;
    }

    const double presHeadMax = 10.0;  ///< Max allowed absolute value for pressure head [m]
    const double epsilon = 1e-6;      ///< Small margin to maintain water content above residual

    // Enforce upper saturation limit.
    if (watContNew >= state->saturatedWaterContent) {
        // Keep current state, no additional water added.
    }
    // Enforce lower limit just above residual water content.
    else if (watContNew <= state->residualWaterContent + epsilon) {
        state->watCont = state->residualWaterContent + epsilon;
        state->presHead = -presHeadMax;
    } else {
        // Normal condition: compute updated water content and pressure head.
        state->watCont = watContNew;
        state->presHead = -std::pow(
            std::pow((state->saturatedWaterContent - state->residualWaterContent) /
                     (state->watCont - state->residualWaterContent), 1.0 / state->vanGenuchtenM) - 1.0,
            1.0 / state->vanGenuchtenN) / state->vanGenuchtenAlpha;

        // Limit pressure head to realistic bounds.
        if (state->presHead < -presHeadMax) {
            state->presHead = -presHeadMax;
            state->watCont = state->residualWaterContent + epsilon;
        }
    }

    // Compute and update hydraulic head.
    Point3D centerPoint = cell->computeCenter(subsurfaceMesh->points);
    state->hydraulicHead = state->presHead + centerPoint.z;

    return state->watCont;
}
