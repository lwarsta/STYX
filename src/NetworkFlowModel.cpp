/**
 * @file NetworkFlowModel.cpp
 * @brief Implements the NetworkFlowModel class for network flow simulation.
 */

#include "NetworkFlowModel.h"
#include "PhysicsState.h"
#include <iostream>
#include <cmath>
#include <algorithm>

/**
 * @brief Initializes the network flow model with settings from the JSON configuration.
 *
 * The function reads parameters from the "WaterFlow"->"NetworkFlowBasic" section.
 * If the section is missing, fallback default values are used.
 *
 * @param settings The JSON object containing configuration settings.
 */
void NetworkFlowModel::initialize(const JsonValue &settings) {
    // Check if the settings object contains the "WaterFlow" section and its "NetworkFlowBasic" subsection.
    if (settings.is_object() &&
        settings.as_object().find("WaterFlow") != settings.as_object().end() &&
        settings.as_object().at("WaterFlow").is_object() &&
        settings.as_object().at("WaterFlow").as_object().find("NetworkFlowBasic") != settings.as_object().at("WaterFlow").as_object().end())
    {
		const JsonValue& networkFlowSettings = settings.as_object().at("WaterFlow").as_object().at("NetworkFlowBasic");

        // Load the Enabled flag; default to true if not present.
        if (networkFlowSettings.as_object().find("Enabled") != networkFlowSettings.as_object().end())
            enabled_ = networkFlowSettings.as_object().at("Enabled").as_bool();
        else
            enabled_ = true;

        // Load the iteration cut threshold.
        if (networkFlowSettings.as_object().find("IterationCutThresholdNetwork") != networkFlowSettings.as_object().end())
            iterCutThresh_ = networkFlowSettings.as_object().at("IterationCutThresholdNetwork").as_double();
        else
            iterCutThresh_ = 1e-07;

        // Load the maximum iterations.
        if (networkFlowSettings.as_object().find("MaxIterationsNetwork") != networkFlowSettings.as_object().end())
            maxIter_ = networkFlowSettings.as_object().at("MaxIterationsNetwork").is_int() ?
                networkFlowSettings.as_object().at("MaxIterationsNetwork").as_int() :
                static_cast<int>(networkFlowSettings.as_object().at("MaxIterationsNetwork").as_double());
        else
            maxIter_ = 100;
		
        // Load the bisection iteration threshold.
        if (networkFlowSettings.as_object().find("BisectionIterationThresholdNetwork") != networkFlowSettings.as_object().end())
            iterCutThreshBis_ = networkFlowSettings.as_object().at("BisectionIterationThresholdNetwork").as_double();
        else
            iterCutThreshBis_ = 1e-08;

        // Load the bisection maximum iterations.
        if (networkFlowSettings.as_object().find("BisectionMaxIterationsNetwork") != networkFlowSettings.as_object().end())
            maxIterBis_ = networkFlowSettings.as_object().at("BisectionMaxIterationsNetwork").is_int() ?
                networkFlowSettings.as_object().at("BisectionMaxIterationsNetwork").as_int() :
                static_cast<int>(networkFlowSettings.as_object().at("BisectionMaxIterationsNetwork").as_double());
        else
            maxIterBis_ = 100;

        // Load bisection left depth.
        if (networkFlowSettings.as_object().find("BisectionLeftDepthNetwork") != networkFlowSettings.as_object().end())
            leftDepthBis_ = networkFlowSettings.as_object().at("BisectionLeftDepthNetwork").as_double();
        else
            leftDepthBis_ = 0.0;

        // Load bisection right depth.
        if (networkFlowSettings.as_object().find("BisectionRightDepthNetwork") != networkFlowSettings.as_object().end())
            rightDepthBis_ = networkFlowSettings.as_object().at("BisectionRightDepthNetwork").as_double();
        else
            rightDepthBis_ = 10.0;
    } else {
        // Fallback defaults if the section is missing.
        enabled_ = false;
        iterCutThresh_ = 1e-07;
        maxIter_ = 100;
        iterCutThreshBis_ = 1e-08;
        maxIterBis_ = 100;
        leftDepthBis_ = 0.0;
        rightDepthBis_ = 10.0;
    }
}

/**
 * @brief Executes a simulation step for the network flow model.
 *
 * The function first updates the previous water levels for junctions and links.
 * Then, using a sub time stepping scheme, it exchanges water with the surface (if a mesh coupler is set),
 * calls the iterative solver (iterate), and then performs postprocessing.
 *
 * @param timeStep Simulation time step (seconds).
 * @param forcing Environmental forcing data (currently unused).
 */
void NetworkFlowModel::runStep(double timeStep, const ForcingRecord &forcing) {
    if (!enabled_)
        return;
    
    (void) forcing; // Forcing is not used in the current implementation.
    
    // Store previous state for network junctions.
    for (auto &cellPtr : networkMesh_->cells) {
        auto junction = dynamic_cast<NetworkJunctionCell*>(cellPtr.get());
        if (!junction || !junction->waterState || junction->material == 0)
            continue;
        auto state = std::dynamic_pointer_cast<NetworkWaterState>(junction->waterState);
        if (!state)
            continue;

        state->waterLevelOld = state->waterLevel;
    }

    // Store previous state for network links.
    for (auto &cellPtr : networkMesh_->cells) {
        auto link = dynamic_cast<NetworkLinkCell*>(cellPtr.get());
        if (!link || !link->waterState)
            continue;

        auto state = std::dynamic_pointer_cast<NetworkWaterState>(link->waterState);
        if (!state)
            continue;

        state->waterLevelOld = state->waterLevel;
    }
    
    // Run the model using a sub time stepping scheme.
    double subTimeAccum = 0.0;
    double subTimeStep = 1.0; // magic number for sub time step

    while (subTimeAccum < timeStep) {
        if (meshCoupler_) {
            meshCoupler_->preExchangeSurfaceNetwork(subTimeStep);
        }
        
        iterate(subTimeStep);
        postprocess();
        
        if (meshCoupler_) {
            meshCoupler_->postExchangeSurfaceNetwork(subTimeStep);
        }

        subTimeAccum += subTimeStep;
    }
}

/**
 * @brief Performs the iterative solver for the network flow model over a sub time step.
 *
 * This function subdivides the given time step into smaller sub‐steps if needed, iteratively
 * solving for the new water levels at junctions and links until the solution converges
 * (or a maximum number of iterations is reached). It also updates the state of the mesh
 * by reverting and swapping water levels if necessary.
 *
 * @param subTimeStep The length of the time step (in seconds) to be processed.
 */
void NetworkFlowModel::iterate(double subTimeStep) {
    // Obtain a reference to the junction vector from the network mesh.
    std::vector<NetworkJunctionCell*>& junctions = networkMesh_->junctions;

    double timeLoc = 0.0;
    double subTimeStepFrac;
    const double errorThresh = 0.0001;        // Convergence threshold (magic number)
    int iterations;
    const int iterationsMax = 1000;           // Maximum iterations allowed (magic number)
    const double waterDepthMinThresh = -0.01;   // Minimum allowed water depth (magic number)
    double waterDepthMin;
    double volumeError = 0.0;

    // Break the sub time step into smaller parts if required.
    do {
        iterations = 0;
        double sysWaterVolumeOld = getNetworkWaterVolume();
        subTimeStepFrac = 1.0;

        // Ensure the sub time step does not exceed the remaining time.
        if (subTimeStepFrac > subTimeStep) {
            subTimeStepFrac = subTimeStep;
        }
        if (subTimeStepFrac > subTimeStep - timeLoc && (subTimeStep - timeLoc) > 0.0) {
            subTimeStepFrac = subTimeStep - timeLoc;
        }

        // Begin error iteration loop.
        do {
            // Loop over each junction.
            for (size_t i = 0; i < junctions.size(); i++) {
                NetworkJunctionCell* junc = junctions[i];
                auto state = std::dynamic_pointer_cast<NetworkWaterState>(junc->waterState);
                if (!state)
                    continue;

                // Retrieve junction properties and compute the hydraulic head.
                double areaJunc = networkMesh_->getJunctionArea(junc->diameter);
                double elevBottJunc = junc->bottom.z;
                double depthWaterOldJunc = state->waterLevelOld;
                double headOldJunc = elevBottJunc + depthWaterOldJunc;

                // Build temporary vectors of neighbor junctions and corresponding links.
                std::vector<NetworkJunctionCell*> neighJuncs;
                std::vector<NetworkLinkCell*> linksNeigh;
                for (auto link : junc->connectedLinks) {
                    NetworkJunctionCell* other = (link->junction1 == junc) ? link->junction2 : link->junction1;
                    if (other) {
                        neighJuncs.push_back(other);
                        linksNeigh.push_back(link);
                    }
                }

                // Loop over each neighbor junction.
                for (size_t j = 0; j < neighJuncs.size(); j++) {
                    NetworkJunctionCell* neigh = neighJuncs[j];
                    auto neighState = std::dynamic_pointer_cast<NetworkWaterState>(neigh->waterState);
                    if (!neighState)
                        continue;

                    // Retrieve neighbor junction properties.
                    double areaJuncNeigh = networkMesh_->getJunctionArea(neigh->diameter);
                    double elevBottJuncNeigh = neigh->bottom.z;
                    double depthWaterOldJuncNeigh = neighState->waterLevelOld;
                    double headOldJuncNeigh = elevBottJuncNeigh + depthWaterOldJuncNeigh;

                    // Retrieve link properties connecting the junctions.
                    NetworkLinkCell* link = linksNeigh[j];
                    double diamLink = link->diameter;
                    double slopeLink = networkMesh_->getLinkSlope(link);
                    double lengthFlatLink = networkMesh_->getLinkLateralDistance(link);
                    double lengthLink = networkMesh_->getLinkLength(link);
                    double areaLink = networkMesh_->getLinkCrossectionArea(link);
                    Point3D centrePointLink = networkMesh_->getLinkCentrePoint(link);
                    double elevCentrePointLink = centrePointLink.z;
                    auto linkState = std::dynamic_pointer_cast<NetworkWaterState>(link->waterState);
                    double waterDepthOldLink = (linkState ? linkState->waterLevelOld : 0.0);
                    double headOldLink = elevCentrePointLink + waterDepthOldLink;
                    double filledAreaOldLink = 0.0;
                    double hydraulicRadOldLink = 0.0;
                    getFlowAreaAndHydrRad(waterDepthOldLink, diamLink, filledAreaOldLink, hydraulicRadOldLink);
                    double filledVolumeOldLink = lengthLink * filledAreaOldLink;
                    double fullVolumeLink = lengthLink * areaLink;
                    double freeVolumeOldLink = fullVolumeLink - filledVolumeOldLink;
                    double mannN = link->roughness;

                    // Determine which end of the link is associated with the current junction.
                    int idLnkEnd = (link->junction1 == junc) ? 0 : 1;
                    double elevBottLink, elevBottLinkNeigh;
                    if (idLnkEnd == 0) {
                        elevBottLink = link->end1.z;
                        elevBottLinkNeigh = link->end2.z;
                    } else {
                        elevBottLink = link->end2.z;
                        elevBottLinkNeigh = link->end1.z;
                    }

                    // Unpressurized flow computation.
                    if (freeVolumeOldLink > 0.0 || headOldJunc < elevBottLink + diamLink ||
                        headOldJuncNeigh < elevBottLinkNeigh + diamLink) {
                        if (headOldJunc > headOldLink) {  // Flow from junction to link.
                            double depthWaterOldJuncEff = headOldJunc - elevBottLink;
                            if (depthWaterOldJuncEff < 0.0)
                                depthWaterOldJuncEff = 0.0;
                            double filledArea = 0.0;
                            double hydraulicRad = 0.0;
                            getFlowAreaAndHydrRad(depthWaterOldJuncEff, diamLink, filledArea, hydraulicRad);
                            double velocity = 1.0 / mannN * sqrt(slopeLink) * pow(hydraulicRad, 2.0 / 3.0);
                            double discharge = velocity * filledArea;
                            double volume = discharge * subTimeStepFrac;
                            double deltaWaterDepth = volume / areaJunc;
                            state->waterLevel -= deltaWaterDepth;
                            double depthWaterLinkLoc = linkState->waterLevel;
                            double filledAreaLinkLoc = 0.0;
                            double hydraulicRadLinkLoc = 0.0;
                            getFlowAreaAndHydrRad(depthWaterLinkLoc, diamLink, filledAreaLinkLoc, hydraulicRadLinkLoc);
                            double volWaterLinkLoc = filledAreaLinkLoc * lengthLink;
                            double freeVolumeLinkLoc = fullVolumeLink - volWaterLinkLoc;
                            double volumeExcessLoc = volume - freeVolumeLinkLoc;
                            if (volumeExcessLoc > 0.0) {
                                volume = freeVolumeLinkLoc;
                                double deltaWaterDepthNeigh = volumeExcessLoc / areaJuncNeigh;
                                neighState->waterLevel += deltaWaterDepthNeigh;
                            }
                            filledAreaLinkLoc = (volWaterLinkLoc + volume) / lengthLink;
                            double depthLinkLoc = calcWaterDepthInLink(filledAreaLinkLoc, 0.5 * diamLink);
                            linkState->waterLevel = depthLinkLoc;
                        } else {  // Flow from link to junction.
                            double velocity = 1.0 / mannN * sqrt(slopeLink) * pow(hydraulicRadOldLink, 2.0 / 3.0);
                            double discharge = velocity * filledAreaOldLink;
                            double volume = discharge * subTimeStepFrac;
                            double deltaWaterDepth = volume / areaJunc;
                            state->waterLevel += deltaWaterDepth;
                            double depthWaterLinkLoc = linkState->waterLevel;
                            double filledAreaLinkLoc = 0.0;
                            double hydraulicRadLinkLoc = 0.0;
                            getFlowAreaAndHydrRad(depthWaterLinkLoc, diamLink, filledAreaLinkLoc, hydraulicRadLinkLoc);
                            double volWaterLinkLoc = filledAreaLinkLoc * lengthLink - volume;
                            filledAreaLinkLoc = volWaterLinkLoc / lengthLink;
                            double depthLinkLoc = calcWaterDepthInLink(filledAreaLinkLoc, 0.5 * diamLink);
                            linkState->waterLevel = depthLinkLoc;
                        }
                    }
                    // Pressurized flow: compute flow directly between wells.
                    else {
                        if (headOldJunc > headOldJuncNeigh && headOldJunc > elevBottLink + diamLink) {
                            double slope = (headOldJunc - headOldJuncNeigh) / lengthFlatLink;
                            double hydraulicRad = M_PI * diamLink / areaLink;
                            double velocity = 1.0 / mannN * sqrt(slope) * pow(hydraulicRad, 2.0 / 3.0);
                            double discharge = velocity * areaLink;
                            double volume = discharge * subTimeStepFrac;
                            double deltaWaterDepth = volume / areaJuncNeigh;
                            state->waterLevel -= deltaWaterDepth;
                        } else if (headOldJunc <= headOldJuncNeigh && headOldJuncNeigh > elevBottLinkNeigh + diamLink) {
                            double slope = (headOldJuncNeigh - headOldJunc) / lengthFlatLink;
                            double hydraulicRad = M_PI * diamLink / areaLink;
                            double velocity = 1.0 / mannN * sqrt(slope) * pow(hydraulicRad, 2.0 / 3.0);
                            double discharge = velocity * areaLink;
                            double volume = discharge * subTimeStepFrac;
                            double deltaWaterDepth = volume / areaJunc;
                            state->waterLevel += deltaWaterDepth;
                        }
                    }
                } // end neighbor loop
            } // end loop over junctions

            // Compute the minimum water level among all junctions and links.
            waterDepthMin = std::numeric_limits<double>::max();
            for (const auto &junction : networkMesh_->junctions) {
                auto state = std::dynamic_pointer_cast<NetworkWaterState>(junction->waterState);
                if (state && state->waterLevel < waterDepthMin) {
                    waterDepthMin = state->waterLevel;
                }
            }
            for (auto link : networkMesh_->links) {
                auto state = std::dynamic_pointer_cast<NetworkWaterState>(link->waterState);
                if (state && state->waterLevel < waterDepthMin) {
                    waterDepthMin = state->waterLevel;
                }
            }

            // Compute the volume error during the sub time step.
            double sysWaterVolume = getNetworkWaterVolume();
            volumeError = fabs(sysWaterVolume - sysWaterVolumeOld);

            // Check if the solution is accurate enough; if not, revert state and reduce sub time step.
            if (volumeError > errorThresh || waterDepthMin < waterDepthMinThresh) {
                subTimeStepFrac *= 0.5;
                revertHeads();
            }
            iterations++;
        } while ((volumeError > errorThresh || waterDepthMin < waterDepthMinThresh) &&
                 iterations < iterationsMax);

        // Update the mesh state: swap new water levels to old values.
        swapHeads();

        timeLoc += subTimeStepFrac;
    } while (timeLoc < subTimeStep);
}

/**
 * @brief Post-processes the network flow model by removing water from outfall cells.
 *
 * This function iterates over all network junction cells. For those with a type value of 1 (indicating an outfall),
 * the water volume (based on the junction’s effective area) is accumulated into cumOutfallVolume and the water levels 
 * (both current and previous) are reset to zero.
 */
void NetworkFlowModel::postprocess() {
    for (auto &cellPtr : networkMesh_->cells) {
        // Remove water from outfall junctions.
        if (NetworkJunctionCell* junction = dynamic_cast<NetworkJunctionCell*>(cellPtr.get())) {
            if (junction->type == 1) {
                auto state = std::dynamic_pointer_cast<NetworkWaterState>(junction->waterState);
                if (!state) continue;
                double area = networkMesh_->getJunctionArea(junction->diameter);
                state->cumOutfallVolume += state->waterLevel * area;
                state->waterLevel = 0.0;
                state->waterLevelOld = 0.0;
            }
        }
    }
}

/**
 * @brief Reverts the water levels for network junctions and links to their previous values.
 *
 * This function resets the current water level (waterLevel) to the stored previous water level (waterLevelOld)
 * for both junctions and links.
 */
void NetworkFlowModel::revertHeads() {
    // Revert water levels for network junctions.
    for (auto junction : networkMesh_->junctions) {
        auto state = std::dynamic_pointer_cast<NetworkWaterState>(junction->waterState);
        if (state) {
            state->waterLevel = state->waterLevelOld;
        }
    }
    // Revert water levels for network links.
    for (auto link : networkMesh_->links) {
        auto state = std::dynamic_pointer_cast<NetworkWaterState>(link->waterState);
        if (state) {
            state->waterLevel = state->waterLevelOld;
        }
    }
}

/**
 * @brief Swaps the current water levels into the previous water level variables.
 *
 * After a time step, this function updates each network cell so that the current water level is saved
 * as the previous water level (waterLevelOld) for use in the next time step.
 */
void NetworkFlowModel::swapHeads() {
    // Swap water levels for network junctions.
    for (auto junction : networkMesh_->junctions) {
        auto state = std::dynamic_pointer_cast<NetworkWaterState>(junction->waterState);
        if (state) {
            state->waterLevelOld = state->waterLevel;
        }
    }
    // Swap water levels for network links.
    for (auto link : networkMesh_->links) {
        auto state = std::dynamic_pointer_cast<NetworkWaterState>(link->waterState);
        if (state) {
            state->waterLevelOld = state->waterLevel;
        }
    }
}

/**
 * @brief Computes the water depth in a link using a bisection method.
 *
 * Given a target cross-sectional flow area and the link radius, this function computes the water depth
 * that produces a cross-sectional area equal to the target area. The method uses a bisection algorithm
 * with a fixed tolerance.
 *
 * @param targetArea The desired flow area.
 * @param radius The radius of the link (half the diameter).
 * @return double The computed water depth in the link.
 */
double NetworkFlowModel::calcWaterDepthInLink(double targetArea, double radius) {
    const double tolerance = 1e-7;
    double heightLow = 0.0;
    double heightHigh = 2.0 * radius; // Maximum water depth is the full diameter.
    double heightMid;

    while (heightHigh - heightLow > tolerance) {
        heightMid = 0.5 * (heightLow + heightHigh);
        double theta = 2 * acos((radius - heightMid) / radius);
        double flowArea = 0.5 * radius * radius * (theta - sin(theta));

        if (flowArea < targetArea) {
            heightLow = heightMid;
        } else {
            heightHigh = heightMid;
        }
    }
    return 0.5 * (heightHigh + heightLow);
}

/**
 * @brief Computes the total water volume in the network mesh.
 *
 * The water volume in junctions is computed as the product of their effective area (approximated as a circle)
 * and their current water level. The volume in links is computed as the product of the filled cross-sectional
 * area and the link length.
 *
 * @return double The total network water volume.
 */
double NetworkFlowModel::getNetworkWaterVolume() {
    double waterVolume = 0.0;

    // Compute water volume stored in junctions.
    for (auto junction : networkMesh_->junctions) {
        if (!junction->waterState)
            continue;
        auto state = std::dynamic_pointer_cast<NetworkWaterState>(junction->waterState);
        if (!state)
            continue;
        double area = networkMesh_->getJunctionArea(junction->diameter);
        waterVolume += state->waterLevel * area;
    }

    // Compute water volume stored in links.
    for (auto link : networkMesh_->links) {
        if (!link->waterState)
            continue;
        auto state = std::dynamic_pointer_cast<NetworkWaterState>(link->waterState);
        if (!state)
            continue;
        double diameter = link->diameter;
        double waterDepth = state->waterLevel;
        double filledArea = 0.0;
        double hydraulicRad = 0.0;
        getFlowAreaAndHydrRad(waterDepth, diameter, filledArea, hydraulicRad);
        double length = networkMesh_->getLinkLength(link);
        waterVolume += length * filledArea;
    }

    return waterVolume;
}

/**
 * @brief Computes the total outfall water volume accumulated in the network.
 *
 * Only junctions with a type value of 1 (indicating outfall) contribute to the outfall volume.
 *
 * @return double The total outfall water volume.
 */
double NetworkFlowModel::getOutfallWaterVolume() {
    double waterVolume = 0.0;
    for (auto junction : networkMesh_->junctions) {
        if (!junction->waterState)
            continue;
        if (junction->type == 1) {
            auto state = std::dynamic_pointer_cast<NetworkWaterState>(junction->waterState);
            if (!state) continue;
            waterVolume += state->cumOutfallVolume;
        }
    }
    return waterVolume;
}

/**
 * @brief Computes the flow area and hydraulic radius for a given water depth and link diameter.
 *
 * Uses a circular cross-sectional model. For water depths less than the full diameter, the area is computed 
 * using a circular segment.
 *
 * @param waterDepth The water depth in the link.
 * @param linkDiameter The diameter of the link.
 * @param flowArea Output parameter to store the computed flow area.
 * @param hydrRad Output parameter to store the computed hydraulic radius.
 */
void NetworkFlowModel::getFlowAreaAndHydrRad(double waterDepth,
                                  double linkDiameter,
                                  double &flowArea,
                                  double &hydrRad) {
    double radius = 0.5 * linkDiameter;
    if (radius <= 0.0 || waterDepth <= 0.0) {
        flowArea = 0.0;
        hydrRad = 0.0;
    } else if (waterDepth >= 2.0 * radius) {
        // Full pipe: area and hydraulic radius for a full circle.
        flowArea = M_PI * radius * radius;
        hydrRad = radius / 2.0;
    } else {
        double theta = 2.0 * acos((radius - waterDepth) / radius);
        flowArea = 0.5 * radius * radius * (theta - sin(theta));
        double perimeter = radius * theta;
        hydrRad = flowArea / perimeter;
    }
}
