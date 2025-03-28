/**
 * @file PhysicsState.h
 * @brief Definitions of physics state classes used in the simulation.
 *
 * This file defines the abstract base class for physics state objects,
 * and concrete classes for surface water, subsurface water, and network water.
 */

#ifndef PHYSICS_STATE_H
#define PHYSICS_STATE_H

/**
 * @brief Abstract base class for physics state objects.
 */
class PhysicsState {
public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~PhysicsState() = default;
};

/**
 * @brief Physics state for surface water cells (2D).
 *
 * This class holds the water state variables for surface water cells,
 * including water depths, boundary conditions, material properties, and cumulative mass balance values.
 */
class SurfaceWaterState : public PhysicsState {
public:
    double waterDepth = 0.0;              ///< Initial water depth [m].
    double waterDepthOld = 0.0;           ///< Water depth from the previous time step [m].
    double upperStorageDepth = 0.0;       ///< Upper storage water depth [m].

    // Boundary condition fields:
    int outletId = -1;                  ///< Outlet identifier.
    double distanceToOutlet = 0.0;      ///< Distance to the outlet [m].
    int sinkId = -1;                    ///< Sink identifier.

    // Material properties:
    double mannings = 0.0;              ///< Manning's roughness coefficient.
    double depressionStorage = 0.0;     ///< Depression storage [m].
    double upperStorage = 0.0;          ///< Maximum upper storage capacity [m].
    double evaporationFraction = 0.0;   ///< Fraction of water evaporated from the surface.
    double cropFactor = 0.0;            ///< Crop factor affecting water uptake.
    double rootDepth = 0.0;             ///< Root depth [m].

    // Mass balance related variables:
    double cumEvaporationVolume = 0.0;  ///< Cumulative evaporation volume [m^3].
    double cumPrecipitationVolume = 0.0;  ///< Cumulative precipitation volume [m^3].
    double cumPetVolume = 0.0;          ///< Cumulative potential evapotranspiration volume [m^3].
};

/**
 * @brief Physics state for subsurface water cells (3D).
 *
 * This class stores the hydraulic state of subsurface cells,
 * along with soil hydraulic properties, derived parameters (via the van Genuchten model),
 * intercell flux variables, transpiration variables, and exchange variables between grids.
 */
class SubsurfaceWaterState : public PhysicsState {
public:
    // Primary state.
    double hydraulicHead = 0.0;         ///< Hydraulic head [m].
    double hydraulicHeadOld = 0.0;      ///< Hydraulic head from the previous time step [m].

    // Soil hydraulic properties (from the material library).
    double saturatedHydraulicConductivity = 0.0; ///< Saturated hydraulic conductivity [m/s].
    double compressibility = 0.0;                ///< Soil compressibility [1/m].
    double saturatedWaterContent = 0.0;          ///< Saturated water content [m³/m³].
    double residualWaterContent = 0.0;           ///< Residual water content [m³/m³].
    double vanGenuchtenAlpha = 0.0;              ///< van Genuchten α parameter [1/m].
    double vanGenuchtenN = 0.0;                  ///< van Genuchten n parameter (must be > 1).
    double dryWeight = 0.0;                      ///< Dry weight.
    double heatCapacity = 0.0;                   ///< Heat capacity.
    double heatConductivity = 0.0;               ///< Heat conductivity.
    double heatConductivityMultiplier = 0.0;     ///< Heat conductivity multiplier.
    double dispersivityTrans = 0.0;              ///< Transverse dispersivity.
    double dispersivityLong = 0.0;               ///< Longitudinal dispersivity.

    // Derived parameters computed during simulation.
    double vanGenuchtenM = 0.0;                  ///< Computed as 1 - 1/vanGenuchtenN.
    double presHead = 0.0;                       ///< Pressure head [m].
    double watCont = 0.0;                        ///< Water content [m³/m³].
    double watContOld = 0.0;                     ///< Previous water content [m³/m³].
    double diffWatCap = 0.0;                     ///< Differential water capacity [m³/m³/m].
    double condUnsat = 0.0;                      ///< Unsaturated hydraulic conductivity [m/s].

    // Intercell variables.
    std::vector<double> condInter;               ///< Intercell conductivity for each face.
    std::vector<double> fluxes;                  ///< Fluxes across each face.
    std::vector<double> velInter;                ///< Flow velocities across each face.

    /**
     * @brief Initializes the intercell variables for a given number of neighbors.
     * @param numNeighbors Number of neighboring cells (faces).
     */
    void initializeIntercellVariables(size_t numNeighbors) {
        condInter.resize(numNeighbors, 0.0);
        fluxes.resize(numNeighbors, 0.0);
        velInter.resize(numNeighbors, 0.0);
    }

    // Transpiration related variables.
    double cropFactor = 0.0;                     ///< Crop factor (from surface material properties).
    double rootDepth = 0.0;                      ///< Root depth (from surface material properties) [m].
    double evaporationFraction = 0.0;            ///< Evaporation fraction (from surface material properties).
    double pressHeadMin = -5.0;                  ///< Minimum pressure head [m].
    double pressHeadMax = -0.1;                  ///< Maximum pressure head [m].
    double pressHeadWilt = -150.0;               ///< Pressure head at wilting point [m].
    double transpFlux = 0.0;                     ///< Transpiration flux [m/s].
    double cumTranspVolume = 0.0;                ///< Cumulative transpiration volume [m^3].

    // Exchange between grids.
    double fluxExchangeSurfSub = 0.0;            ///< Flux exchanged between surface and subsurface [m/s].
    double cumExchangeSurfSub = 0.0;             ///< Cumulative exchange volume between surface and subsurface [m^3].
};

/**
 * @brief Physics state for network components (e.g., wells, pipes).
 *
 * This class holds the water level and cumulative outfall volume
 * for network components such as junctions and pipes.
 */
class NetworkWaterState : public PhysicsState {
public:
    double waterLevel = 0.0;      ///< Water level in a junction or pipe [m].
    double waterLevelOld = 0.0;   ///< Water level from the previous time step [m].
    double cumOutfallVolume = 0.0;///< Cumulative outfall volume (valid for junctions) [m^3].
};

#endif // PHYSICS_STATE_H
