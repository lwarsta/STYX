/**
 * @file SubsurfaceFlowModel.cpp
 * @brief Implements the subsurface flow model for simulating water flow in the subsurface domain.
 *
 * This file defines the member functions of the SubsurfaceFlowModel class. The model is
 * initialized using a JSON configuration and processes forcing data during each simulation step.
 */

#include "SubsurfaceFlowModel.h"
#include "PhysicsState.h"
#include <iostream>
#include <cmath>
#include <stdexcept>

/**
 * @brief Initializes the subsurface flow model using the provided JSON settings.
 *
 * This function reads configuration parameters from the JSON settings. It checks for the
 * "SubsurfaceFlowBasic" section under "WaterFlow" to determine if the model is enabled and
 * to load parameters such as maximum iterations and convergence thresholds. It also loads
 * the weather data interval from the "TimeControl" section. Finally, it initializes basic cell
 * properties needed for the initial water balance computation.
 *
 * @param settings A JSON object containing simulation settings.
 */
void SubsurfaceFlowModel::initialize(const JsonValue &settings) {
	// Check if the settings object contains the "WaterFlow" section and its "SubsurfaceFlowBasic" sub‐section.
    if (settings.is_object() &&
        settings.as_object().find("WaterFlow") != settings.as_object().end() &&
        settings.as_object().at("WaterFlow").is_object() &&
        settings.as_object().at("WaterFlow").as_object().find("SubsurfaceFlowBasic") != settings.as_object().at("WaterFlow").as_object().end())
    {
		const JsonValue& subsurfSettings = settings.as_object().at("WaterFlow").as_object().at("SubsurfaceFlowBasic");

        // Load the Enabled flag; default to true if not present.
        if (subsurfSettings.is_object() &&
            subsurfSettings.as_object().find("Enabled") != subsurfSettings.as_object().end())
        {
            enabled_ = subsurfSettings.as_object().at("Enabled").as_bool();
        } else {
            enabled_ = true;
        }

        // If the model is disabled, print a message and return.
        if (!enabled_) {
            std::cout << "SubsurfaceFlowModel is disabled." << std::endl;
            return;
        }

        // Load other model parameters with fallback defaults.
        if (subsurfSettings.as_object().find("MaxIterationsSubsurface") != subsurfSettings.as_object().end())
            maxIter_ = subsurfSettings.as_object().at("MaxIterationsSubsurface").is_int() ?
                       subsurfSettings.as_object().at("MaxIterationsSubsurface").as_int() :
                       static_cast<int>(subsurfSettings.as_object().at("MaxIterationsSubsurface").as_double());
        else
            maxIter_ = 50;

        if (subsurfSettings.as_object().find("IterationCutThresholdSubsurface") != subsurfSettings.as_object().end())
            iterCutThresh_ = subsurfSettings.as_object().at("IterationCutThresholdSubsurface").as_double();
        else
            iterCutThresh_ = 1e-08;
    } else {
        // Default values if the section is missing.
        enabled_ = false;
        maxIter_ = 50;
        iterCutThresh_ = 1e-08;
    }

    // Load weather data interval from the TimeControl section.
    if (settings.is_object() &&
        settings.as_object().find("TimeControl") != settings.as_object().end() &&
        settings.as_object().at("TimeControl").is_object() &&
        settings.as_object().at("TimeControl").as_object().find("WeatherDataInterval") != settings.as_object().at("TimeControl").as_object().end())
    {
        weatherDataInterval_ = settings.as_object().at("TimeControl").as_object().at("WeatherDataInterval").as_double();
    } else {
        weatherDataInterval_ = 120.0;
    }

    // Initiate cell properties needed for the initial water balance computation.
    for (auto &cellPtr : subsurfaceMesh_->cells) {
        auto cell = dynamic_cast<VolumeCell*>(cellPtr.get());
        if (!cell)
            continue;
        auto state = std::dynamic_pointer_cast<SubsurfaceWaterState>(cell->waterState);
        if (!state)
            continue;
        Point3D center = cell->computeCenter(subsurfaceMesh_->points);
        state->presHead = state->hydraulicHead - center.z;
        MeshCoupler::calcWatCont(cell);
        state->watContOld = state->watCont;
        calcDiffWatCap(cell);
    }
}

/**
 * @brief Executes a single simulation step for subsurface flow.
 *
 * This function applies the forcing data to the subsurface mesh by preprocessing the cells,
 * exchanging water with the surface (if a MeshCoupler is set), performing the iterative solver,
 * and then exchanging water again after the iteration.
 *
 * @param timeStep The simulation time step (in seconds).
 * @param forcing Forcing data record containing PET, precipitation, and temperature.
 */
void SubsurfaceFlowModel::runStep(double timeStep, const ForcingRecord &forcing) {
    if (!enabled_)
        return;

    preprocess(timeStep, forcing);

    if (meshCoupler_) {
        meshCoupler_->preExchangeSurfaceSubsurface(timeStep);
    }

    iterate(timeStep);

    if (meshCoupler_) {
        meshCoupler_->postExchangeSurfaceSubsurface(timeStep);
    }
}

/**
 * @brief Preprocesses the subsurface mesh before solving.
 *
 * This function updates the basic properties of each cell, such as unsaturated conductivity
 * and water content, then computes intercell properties (conductivity, fluxes, and flow velocities)
 * and applies transpiration effects using the potential evapotranspiration (PET) value.
 *
 * @param timeStep The simulation time step (in seconds).
 * @param forcing Forcing data record containing PET and other environmental data.
 */
void SubsurfaceFlowModel::preprocess(double timeStep, const ForcingRecord &forcing) {
    // Update basic cell properties.
    for (auto &cellPtr : subsurfaceMesh_->cells) {
        auto cell = dynamic_cast<VolumeCell*>(cellPtr.get());
        if (!cell)
            continue;
        auto state = std::dynamic_pointer_cast<SubsurfaceWaterState>(cell->waterState);
        if (!state)
            continue;
        calcUnsatCond(cell);
        state->hydraulicHeadOld = state->hydraulicHead;
        state->watContOld = state->watCont;
    }

    // Compute intercell properties.
    for (auto &cellPtr : subsurfaceMesh_->cells) {
        auto cell = dynamic_cast<VolumeCell*>(cellPtr.get());
        if (!cell)
            continue;
        calcCondInter(cell);
        calcFluxes(cell);
        calcFlowVel(cell);
    }

    applyTranspiration(timeStep, forcing.pet);
}

/**
 * @brief Iteratively solves for the hydraulic head in each subsurface cell.
 *
 * This function uses an implicit iterative method to update the hydraulic head in each
 * cell of the subsurface mesh until the maximum change (deviation) is below a specified
 * convergence threshold or the maximum number of iterations is reached. The computation
 * is parallelized with OpenMP.
 *
 * @param timeStep The simulation time step in seconds.
 */
void SubsurfaceFlowModel::iterate(double timeStep) {
    double iterDev = 0.0;
    int iterCount = 0;
    do {
        iterDev = 0.0;
        
        // Parallelize the loop over all cells.
        #pragma omp parallel
        {
            double threadMaxDev = 0.0;
            
            #pragma omp for nowait
            for (size_t i = 0; i < subsurfaceMesh_->cells.size(); i++) {
                VolumeCell* cell = dynamic_cast<VolumeCell*>(subsurfaceMesh_->cells[i].get());
                if (!cell || cell->material == 0)
                    continue;
                auto state = std::dynamic_pointer_cast<SubsurfaceWaterState>(cell->waterState);
                if (!state)
                    continue;
                
                // Base of the implicit equation.
                double volume = cell->computeVolume(subsurfaceMesh_->points);
                double numerator = state->diffWatCap * volume / timeStep * state->hydraulicHeadOld;
                double denominator = state->diffWatCap * volume / timeStep;
                Point3D center = cell->computeCenter(subsurfaceMesh_->points);
                
                // Sum contributions from neighbors.
                for (size_t j = 0; j < cell->neighborIndices.size(); j++) {
                    int neighIndex = cell->neighborIndices[j];
                    if (neighIndex < 0 || static_cast<size_t>(neighIndex) >= subsurfaceMesh_->cells.size())
                        continue;
                    VolumeCell* neigh = dynamic_cast<VolumeCell*>(subsurfaceMesh_->cells[neighIndex].get());
                    if (!neigh || !neigh->waterState || neigh->material == 0)
                        continue;
                    auto neighborState = std::dynamic_pointer_cast<SubsurfaceWaterState>(neigh->waterState);
                    if (neighborState) {
                        numerator -= state->fluxes[j] * neighborState->hydraulicHead;
                        denominator -= state->fluxes[j];
                    }
                }
                
                // Include transpiration.
                numerator -= state->transpFlux;
                
                // Compute new hydraulic head.
                double newHead = (denominator != 0.0) ? numerator / denominator : state->hydraulicHeadOld;
                double cellDev = std::fabs(state->hydraulicHead - newHead);
                threadMaxDev = std::max(threadMaxDev, cellDev);
                
                // Update cell state.
                state->hydraulicHead = newHead;
                state->presHead = state->hydraulicHead - center.z;
                MeshCoupler::calcWatCont(cell);
                calcDiffWatCap(cell);
            }
            
            // Update global deviation.
            #pragma omp critical
            {
                iterDev = std::max(iterDev, threadMaxDev);
            }
        }
        iterCount++;
    } while (iterDev > iterCutThresh_ && iterCount < maxIter_);
    
    std::cout << "Subsurface flow iteration completed in " << iterCount 
              << " steps with final deviation: " << iterDev << std::endl;
}

/**
 * @brief Computes the moisture factor for transpiration using pressure head.
 *
 * This function calculates a moisture factor based on the effective pressure head
 * relative to specified minimum, maximum, and wilting thresholds.
 *
 * @param presHead The pressure head (m).
 * @param phMin The minimum pressure head (m).
 * @param phMax The maximum pressure head (m).
 * @param phWilt The wilting pressure head (m).
 * @return double The computed moisture factor (0 to 1).
 */
double SubsurfaceFlowModel::computeMoistureFactor(double presHead, double phMin, double phMax, double phWilt) {
    if (presHead > 0.0)
        presHead = 0.0;
    if (presHead < phWilt)
        presHead = phWilt;
    if (presHead >= phMin && presHead <= phMax)
        return 1.0;
    else if (presHead > phMax)
        return 1.0 - (presHead - phMax) / (0.0 - phMax);
    else
        return 1.0 - (presHead - phMin) / (phWilt - phMin);
}

/**
 * @brief Applies transpiration flux down a column of subsurface cells.
 *
 * Starting from a top cell, this function computes the transpiration flux for each cell
 * in the column until the cumulative depth reaches the root depth. The transpiration flux
 * is computed based on a moisture factor, root factor, and plant parameters.
 *
 * @param topCell The top cell of the subsurface column.
 * @param timeStep The simulation time step (s).
 * @param pet_m The potential evapotranspiration (m).
 */
void SubsurfaceFlowModel::applyTranspirationToColumn(VolumeCell* topCell, double timeStep, double pet_m) {
    double cumulativeDepth = 0.0;
    VolumeCell* currentCell = topCell;
    
    while (currentCell != nullptr) {
        auto state = std::dynamic_pointer_cast<SubsurfaceWaterState>(currentCell->waterState);
        if (!state)
            break;
        // Determine the fraction of the cell within the root zone.
        double thickness = currentCell->computeThickness(subsurfaceMesh_->points);
        double rootDepth = state->rootDepth;
        double rootFactor = 0.0;
        if (cumulativeDepth + thickness <= rootDepth)
            rootFactor = thickness / rootDepth;
        else if (cumulativeDepth < rootDepth && cumulativeDepth + thickness > rootDepth)
            rootFactor = (rootDepth - cumulativeDepth) / rootDepth;
        // Compute the moisture factor.
        double phMin = state->pressHeadMin;
        double phMax = state->pressHeadMax;
        double phWilt = state->pressHeadWilt;
        double moistureFactor = computeMoistureFactor(state->presHead, phMin, phMax, phWilt);
        // Compute transpiration volume.
        std::vector<double> faceAreas = currentCell->computeFaceAreas(subsurfaceMesh_->points);
        double cropFactor = state->cropFactor;
        double evapFraction = state->evaporationFraction;
        double cellTranspVol = (1.0 - evapFraction) * moistureFactor * rootFactor * cropFactor * 0.001 * pet_m / weatherDataInterval_ * faceAreas[0] * timeStep;
        state->transpFlux = cellTranspVol / timeStep;
        state->cumTranspVolume += cellTranspVol;
        // Update cumulative depth.
        cumulativeDepth += thickness;
        if (cumulativeDepth >= rootDepth)
            break;
        // Move to the bottom neighbor.
        int bottomNeighborIdx = currentCell->neighborIndices.back();
        if (bottomNeighborIdx == -1)
            break;
        currentCell = dynamic_cast<VolumeCell*>(subsurfaceMesh_->cells[bottomNeighborIdx].get());
    }
}

/**
 * @brief Applies transpiration flux across the subsurface mesh.
 *
 * This function iterates over all subsurface cells and applies transpiration to each
 * column starting from top cells (identified by a missing neighbor at the top edge).
 *
 * @param timeStep The simulation time step (s).
 * @param pet_m The potential evapotranspiration (m).
 */
void SubsurfaceFlowModel::applyTranspiration(double timeStep, double pet_m) {
    for (size_t i = 0; i < subsurfaceMesh_->cells.size(); i++) {
        VolumeCell* cell = dynamic_cast<VolumeCell*>(subsurfaceMesh_->cells[i].get());
        if (!cell || cell->material == 0)
            continue;
        // A top cell is identified by neighborIndices[0] == -1.
        if (!cell->neighborIndices.empty() && cell->neighborIndices[0] == -1) {
            applyTranspirationToColumn(cell, timeStep, pet_m);
        }
    }
}

/**
 * @brief Calculates unsaturated hydraulic conductivity using the van Genuchten model.
 *
 * This function computes the relative saturation and applies the van Genuchten formula to
 * update the unsaturated conductivity in the cell's water state.
 *
 * @param cell A pointer to a volume cell.
 */
void SubsurfaceFlowModel::calcUnsatCond(VolumeCell* cell) {
    auto state = std::dynamic_pointer_cast<SubsurfaceWaterState>(cell->waterState);
    if (!state)
        return;

    // Compute relative saturation.
    double sat = 1.0;
    double delta = state->saturatedWaterContent - state->residualWaterContent;
    if (delta > 0.0)
        sat = (state->watCont - state->residualWaterContent) / delta;
    sat = std::clamp(sat, 0.0, 1.0);

    // Compute m = 1 - 1/n.
    if (state->vanGenuchtenN > 1.0)
        state->vanGenuchtenM = 1.0 - 1.0 / state->vanGenuchtenN;
    else
        state->vanGenuchtenM = 0.0;

    double tmp = 1.0;
    if (state->vanGenuchtenM > 0.0)
        tmp = 1.0 - std::pow(sat, 1.0 / state->vanGenuchtenM);
    tmp = (tmp > 0.0) ? 1.0 - std::pow(tmp, state->vanGenuchtenM) : 1.0;

    state->condUnsat = state->saturatedHydraulicConductivity * std::sqrt(sat) * tmp * tmp;
}

/**
 * @brief Calculates intercell conductivity between a cell and its neighbors.
 *
 * For each face, the intercell conductivity is computed as the geometric mean of the unsaturated
 * conductivities of the cell and the neighboring cell.
 *
 * @param cell A pointer to a volume cell.
 */
void SubsurfaceFlowModel::calcCondInter(VolumeCell* cell) {
    auto state = std::dynamic_pointer_cast<SubsurfaceWaterState>(cell->waterState);
    if (!state)
        return;

    size_t numNeighbors = cell->neighborIndices.size();
    if (state->condInter.size() < numNeighbors)
        state->condInter.resize(numNeighbors, 0.0);

    for (size_t i = 0; i < cell->neighborIndices.size(); i++) {
        int neighIndex = cell->neighborIndices[i];
        if (neighIndex < 0 || static_cast<size_t>(neighIndex) >= subsurfaceMesh_->cells.size())
            continue;
        VolumeCell* neighborCell = dynamic_cast<VolumeCell*>(subsurfaceMesh_->cells[neighIndex].get());
        if (neighborCell && neighborCell->material > 0 && neighborCell->waterState) {
            auto neighborState = std::dynamic_pointer_cast<SubsurfaceWaterState>(neighborCell->waterState);
            if (neighborState) {
                state->condInter[i] = std::sqrt(state->condUnsat * neighborState->condUnsat);
            }
        }
    }
}

/**
 * @brief Calculates the differential water capacity of a cell.
 *
 * This function computes the differential water capacity as the change in water content divided
 * by the change in hydraulic head. If changes are too small, it defaults to the cell's compressibility.
 *
 * @param cell A pointer to a volume cell.
 */
void SubsurfaceFlowModel::calcDiffWatCap(VolumeCell* cell) {
    auto state = std::dynamic_pointer_cast<SubsurfaceWaterState>(cell->waterState);
    if (!state)
        return;
    double hydraulicHeadDelta = state->hydraulicHead - state->hydraulicHeadOld;
    double watContDelta = state->watCont - state->watContOld;

    if (std::fabs(hydraulicHeadDelta) > 1e-12 && std::fabs(watContDelta) > 1e-12)
        state->diffWatCap = watContDelta / hydraulicHeadDelta;
    else
        state->diffWatCap = state->compressibility;
}

/**
 * @brief Calculates fluxes between a cell and its neighbors.
 *
 * For each face, the flux is computed as the negative product of intercell conductivity and face area,
 * divided by the distance to the neighboring cell.
 *
 * @param cell A pointer to a volume cell.
 */
void SubsurfaceFlowModel::calcFluxes(VolumeCell* cell) {
    auto state = std::dynamic_pointer_cast<SubsurfaceWaterState>(cell->waterState);
    if (!state)
        return;

    std::vector<double> faceAreas = cell->computeFaceAreas(subsurfaceMesh_->points);
    size_t numNeighbors = cell->neighborIndices.size();
    if (state->fluxes.size() < numNeighbors)
        state->fluxes.resize(numNeighbors, 0.0);

    for (size_t i = 0; i < numNeighbors; i++) {
        double distance = cell->getDistanceToNeighbor(i, subsurfaceMesh_->points, subsurfaceMesh_->cells);
        if (distance > 0.0 && i < faceAreas.size()) {
            state->fluxes[i] = -state->condInter[i] * faceAreas[i] / distance;
        }
    }
}

/**
 * @brief Calculates flow velocities between a cell and its neighbors.
 *
 * This function computes the flow velocity for each face based on the intercell conductivity,
 * average water content, and the difference in hydraulic head between the cell and its neighbor.
 *
 * @param cell A pointer to a volume cell.
 */
void SubsurfaceFlowModel::calcFlowVel(VolumeCell* cell) {
    auto state = std::dynamic_pointer_cast<SubsurfaceWaterState>(cell->waterState);
    if (!state)
        return;

    size_t numNeighbors = cell->neighborIndices.size();
    if (state->velInter.size() < numNeighbors)
        state->velInter.resize(numNeighbors, 0.0);

    for (size_t i = 0; i < numNeighbors; i++) {
        size_t neighIndex = cell->neighborIndices[i];
        if (neighIndex < subsurfaceMesh_->cells.size()) {
            VolumeCell* neighborCell = dynamic_cast<VolumeCell*>(subsurfaceMesh_->cells[neighIndex].get());
            if (neighborCell && neighborCell->waterState) {
                auto neighborState = std::dynamic_pointer_cast<SubsurfaceWaterState>(neighborCell->waterState);
                if (!neighborState)
                    continue;
                double watContInter = 0.5 * (state->watCont + neighborState->watCont);
                double distance = cell->getDistanceToNeighbor(i, subsurfaceMesh_->points, subsurfaceMesh_->cells);
                if (watContInter > 0.0 && distance > 0.0) {
                    state->velInter[i] = state->condInter[i] / watContInter *
                        (neighborState->hydraulicHead - state->hydraulicHead) / distance;
                } else {
                    state->velInter[i] = 0.0;
                }
            }
        }
    }
}

/**
 * @brief Computes the total subsurface water volume.
 *
 * The total water volume is the sum over all cells of the product of cell volume and water content.
 *
 * @return double Total subsurface water volume in cubic meters.
 */
double SubsurfaceFlowModel::getSubsurfaceWaterVolume() const {
    double waterVolume = 0.0;
    for (auto &cellPtr : subsurfaceMesh_->cells) {
        VolumeCell* cell = dynamic_cast<VolumeCell*>(cellPtr.get());
        if (!cell || !cell->waterState || cell->material == 0)
            continue;
        auto ssw = std::dynamic_pointer_cast<SubsurfaceWaterState>(cell->waterState);
        if (!ssw)
            continue;
        double volume = cell->computeVolume(subsurfaceMesh_->points);
        waterVolume += volume * ssw->watCont;
    }
    return waterVolume;
}

/**
 * @brief Computes the cumulative transpiration water volume.
 *
 * @return double Total transpiration water volume in cubic meters.
 */
double SubsurfaceFlowModel::getTranspirationWaterVolume() const {
    double waterVolume = 0.0;
    for (auto &cellPtr : subsurfaceMesh_->cells) {
        VolumeCell* cell = dynamic_cast<VolumeCell*>(cellPtr.get());
        if (!cell || !cell->waterState || cell->material == 0)
            continue;
        auto ssw = std::dynamic_pointer_cast<SubsurfaceWaterState>(cell->waterState);
        if (!ssw)
            continue;
        waterVolume += ssw->cumTranspVolume;
    }
    return waterVolume;
}

/**
 * @brief Computes the cumulative infiltration water volume.
 *
 * @return double Total infiltration water volume in cubic meters.
 */
double SubsurfaceFlowModel::getInfiltrationWaterVolume() const {
    double waterVolume = 0.0;
    for (auto &cellPtr : subsurfaceMesh_->cells) {
        VolumeCell* cell = dynamic_cast<VolumeCell*>(cellPtr.get());
        if (!cell || !cell->waterState || cell->material == 0)
            continue;
        auto ssw = std::dynamic_pointer_cast<SubsurfaceWaterState>(cell->waterState);
        if (!ssw)
            continue;
        waterVolume += ssw->cumExchangeSurfSub;
    }
    return waterVolume;
}
