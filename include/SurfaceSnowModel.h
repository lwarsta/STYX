#ifndef SURFACE_SNOW_MODEL_H
#define SURFACE_SNOW_MODEL_H

#include "IModel.h"
#include "SurfaceMesh.h"
#include "ForcingData.h"
#include "MeshCoupler.h"
#include <iostream>
#include <typeinfo>

class SurfaceSnowModel : public IModel {
public:
    /**
     * @brief Constructs a SurfaceSnowModel with the specified surface mesh.
     * @param mesh Pointer to a SurfaceMesh that represents the 2D surface domain.
     */
    explicit SurfaceSnowModel(SurfaceMesh* mesh)
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
     * @brief Sets the MeshCoupler to be used for heat exchange.
     * @param coupler Pointer to a MeshCoupler instance.
     */
    void setMeshCoupler(MeshCoupler* coupler) { meshCoupler_ = coupler; }

private:
    SurfaceMesh* surfaceMesh_;      ///< Pointer to the surface mesh.
    MeshCoupler* meshCoupler_;      ///< Pointer to the mesh coupler for heat exchange.

    // Common model parameters and variables.
    bool   enabled_;                ///< Flag indicating if the model is enabled.
    double weatherDataInterval_;    ///< Weather data interval (in seconds).
	
	// Snow processes.
	void preprocess();
	double computeWetBulbTemperature(
		double heatCapAirLoc,
		double atmPresLoc,
		double latHeatSubLoc,
		double tempAirLoc,
		double humRelLoc); 
	double SVPW(double tempAirLoc);
	double SVPI(double tempAirLoc);
	//double getSaturationVaporPressure(double temperature);
	void iterate();
	void postprocess();
	void compCanopAerRes(double &windSpeed);
	void compSnowAlbedo(double timeStep, double tempAirLoc, double precSnowLoc);
	void computeSnowAlbedo(double timeStep, double tempAir, double precSnow);
	void compSnowState(
		double timeStep,
		double precWat,
		double precSnow,
		double tempAir,
		double tempWetBulb,
		double windSpeed,
		double humRel,
		double pressAtm,
		double radShortIn,
		double radLongIn);
	void compSurfBal(
		double tempAir, 
		double humRel,
		double pressAtm,
		double windSpeed,
		double precWat, 
		double precSnow, 
		double radShortIn, 
		double &radShortOut, 
		double radLongIn, 
		double &radLongOut, 
		double &fluxEnergyPrec, 
		double &fluxEnergySens, 
		double &fluxEnergyLat,
		double *fluxEnergyCond,
		double &fluxWaterEvap);
	double compSurfEnBal(
		double tempAir,
		double humRel,
		double pressAtm,
		double windSpeed,
		double precWat,
		double precSnow,
		double tempSurfLoc, 
		double radShortIn, 
		double &radShortOut, 
		double radLongIn, 
		double &radLongOut, 
		double &fluxEnergyPrec, 
		double &fluxEnergySens, 
		double &fluxEnergyLat, 
		double *fluxEnergyCond,
		double &fluxWaterEvap);
	void compPrecEnerFlux(double tempAir, double precSnow, double precWat, double &fluxEnergyPrec);
	void compTurbFluxes(double tempSurfLoc, double tempAir, double humRel, double pressAtm, double windSpeed, double &fluxEnergyLat, double &fluxEnergySens, double &fluxWaterEvap);
	void compTurbCoeff(double windSpeed, double tempAir,  double tempSurfLoc, double &traCofSensHeat, double &traCofLatHeat);
	void compCondFluxes(double tempSurfLoc, double *fluxEnergyCond);
	void compSnowTemperature();
	void compMeltRates(double timeStep, double precWat, double *fluxWaterMelt/*, double *watLiqDiff*/);
	void compSnowDensity(int layer, double timeStep, double precSnow, double precWat, double tempWetBulb, double &fluxWaterEvap, double *watLiqDiff);
	
	// Forcing data.
	double timeStep = 0.0;				// Time step length of the simulation [s] - may need to be converted ?
	double atmPres = 101325.0;			// Athmospheric pressure [Pa]
	double tempAir = 10.0;				// Air temperature [oC]
	double humRel = 98.7;           	// Relative humidity [%]
	double windSpeed = 1.0;				// Wind speed [m s-1]
	double radShortIn = 1.0;        	// Short wave radiation [MJ m-2 h-1]
	double radLongIn = 1.0;         	// Longwave radiation [MJ m-2 h-1]
	
	// Computed variables.
	double sumTime = 0.0;				// Cumulative simulation time [s]
	double tempWetBulb = 0.0;			// Wet-bulb temperature [oC]
	double depthSnow = 0.0;				// Total snow pack depth [m]
	double condSnowAvr = 0.0;			// Average snow pack thermal conductivity []
	double canAeroRes = 0.0;			// Canopy aerodynamic resistance [s m-1]
	double snowAge = 0.0;				// Snow age in albedo computing [?] FIX THIS!
	double albedo = 0.0;				// Snow surface albedo [-]
	double watEqSnow[2];				// Snow water equivalent [mm] CHECK THIS!
	double watLiq[2];					// The amount of liquid water in the snow layers [mm] CHECK THIS!
	double densSnow[2];					// Density of the snow layers [kg m-3]
	double energySnow[2];				// Energy of the snow layers [?] FIX THIS
	double precWat = 0.0;				// Liquid water amount of precipitation [mm] CHECK THIS!
	double precSnow = 0.0;				// Snow water equivalent of precipitation [mm] CHECK THIS!
	double precCum[2];					// Cumulative precipitation [mm]  CHECK THIS!
	double fluxEvapOut = 0.0;			// Evaporation water flux [mm/s] CHECK THIS!
	double fluxWaterOut = 0.0;			// Liquid water out flux of the snow pack [mm/s] CHECK THIS!
	double snowFallCum = 0.0;			// Cumulative snow fall on the snow pack [mm] CHECK THIS!
	double snowFallPer = 0.0;			// Snow fall period length [s] CHECK THIS!
	double tempSnow[2];					// Temperature of the snow layers [oC]
	double fluxEnergySoil = 0.0;		// Energy flux between top soil layers and snow pack [?] FIX THIS!
	double tempSurf = 0.0;				// snow surface Temperature [oC]
	
	// Parameters.
	double iterConvThr = 0.0001;		// Interation convergence threshold [oC]
	double iterCountMax = 100;			// Maximum number of interations [-]
	double albedoMax = 0.85;			// Max albedo [-]
	double albedoMin = 0.2;				// Min albedo [-]
	double albPrecLim = 0.002;			// Snow albedo parameter, precipitation limit [mm] CHECK UNIT!
	double albTimeLim = 2.0;			// Snow albedo parameter, time limit [h]
	double atmStabCof = 1.0;			// Atmospheric stability weighing factor [-]
	double condSatSnow = 20.0;			// Saturated hydraulic conductivity of wet snow [m h-1]
	double condSnowPara[2];				// Snow thermal conductivity parameter [W m-1 oC-1] 
										// Snow thermal conductivity parameter [W m-1 oC-1]
	double convCofLatHeat = 0.0;		// Windless convection coefficient for latent heat flux [kJ m-2 h-1 Pa-1]
	double convCofSensHeat = 7.2;		// Windless convection coefficient for sensible heat flux [kJ m-2 h-1 oC-1]
	double densSnowLim = 0.15;			// Density limit multiplier of snow (0.15 [g cm-3]) FIX THIS!
	double densSnowMax = 500.0;			// Maximum density of snow [kg m-3]
	double depthSnowRef = 0.005;		// Snow surface roughness height [m]
	double emisSnow = 0.99;				// Emissivity of snow [-]
	double refHeightTemp = 2.0;			// Reference height of measured temperature and relative humidity [m]
	double refHeightWind = 2.0;			// Reference height of measured wind [m]
	double richNumLimMax = 0.16;		// Maximum limit of the estimated Richardson number [-]
	double richNumLimMin = -10.0;		// Minimum limit of the estimated Richardson number [-]
	double snoWatEqTopMax = 0.02;		// Maximum snow water equivalent of the upper snow layer [m]
	double snowLiqCap = 0.1;			// Liquid water holding capacity of snow [m3 m-3]
	double snowMetaPara[6];				// Snow metamorphism parameters C1 [cm-1 h-1]
										// Snow metamorphism parameters C2 [cm3 g-1]
										// Snow metamorphism parameters C3 [h-1]
										// Snow metamorphism parameters C4 [K-1]
										// Snow metamorphism parameters C5 [?] FIX THIS!
										// Snow metamorphism parameters Tc [oC]
	double tempSnowMaxChange = 5.0;		// Maximum temperature change of snow in an iteration round [oC]
	double tmpSnowThresMin = 0.0;		// Snowfall/mixed-rain-snow temperature lower limit [oC]
	double tmpSnowThresMax = 2.0;		// Snowfall/mixed-rain-snow temperature upper limit [oC]
	
	// Constants.
	double bolzConst = 5.67E-8;			// Stefan-Boltzmann constant [W m-2 oK-4]
	double densIce = 917.0;				// Density of ice [kg m-3]
	double densWater = 1000.0;			// Density of water [kg m-3]
	double gravAcc = 9.81;				// Acceleration of gravity [m s-2]
	double heatCapAir = 1.005;			// Heat capacity of air (Specific Heat of Moist Air) 1.013 [kJ kg-1 oC-1]
	double heatCapIce = 2.09;			// Heat capacity of ice [kJ kg-1 oC-1]
	double heatCapWater = 4.18;			// Heat capacity of water [kJ kg-1 oC-1]
	double idealGasConst = 287.0;		// Ideal gas constant for dry air [J kg- 1 K-1]
	double latHeatFus = 333.7;			// Latent heat of fusion [kJ kg-1]
	double latHeatSub = 2834.0;			// Latent heat of sublimation [kJ kg-1]
	double tempZeroKel = 273.15;		// Degrees Kelvin at freezing point [K]
	double vonKarmConst = 0.41;			// Von Karman constant [-]
};

#endif // SURFACE_SNOW_MODEL_H
