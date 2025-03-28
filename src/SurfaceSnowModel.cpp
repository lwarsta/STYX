#include "SurfaceSnowModel.h"
#include "PhysicsState.h"
#include <iostream>
#include <cmath>
#include <stdexcept>

/**
 * @brief Initializes the SurfaceSnowModel using settings from the JSON configuration.
 *
 * This method reads the weather data interval from the "TimeControl" section and
 * the surface snow parameters from the "SurfaceSnow" section. If any parameter
 * is missing, a default fallback value is used.
 *
 * @param settings A JSON object containing the simulation configuration.
 */
void SurfaceSnowModel::initialize(const JsonValue &settings) {
	(void)settings;
	
	watEqSnow[0] = 0.0;
	watEqSnow[1] = 0.0;
	watLiq[0] = 0.0;
	watLiq[1] = 0.0;
	densSnow[0] = 0.0;
	densSnow[1] = 0.0;
	energySnow[0] = 0.0;
	energySnow[1] = 0.0;
	precCum[0] = 0.0;
	precCum[1] = 0.0;
	condSnowPara[0] = 0.021;
	condSnowPara[1] = 2.51;
	snowMetaPara[0] = 0.01;
	snowMetaPara[1] = 21.0;
	snowMetaPara[2] = 0.01;
	snowMetaPara[3] = 0.04;
	snowMetaPara[4] = 2.0;
	snowMetaPara[5] = 0.0;
	tempSnow[0] = tempSnow[1] = 0.0;
}

void SurfaceSnowModel::runStep(double timeStep, const ForcingRecord &forcing) {
	(void)timeStep;
	(void)forcing;
	
	preprocess();
	iterate();
	postprocess();
}

void SurfaceSnowModel::preprocess() {
	// Compute wet bulb temperature.
	tempWetBulb = computeWetBulbTemperature(heatCapAir, atmPres, latHeatSub, tempAir, humRel);
}

double SurfaceSnowModel::computeWetBulbTemperature(double heatCapAirLoc, double atmPresLoc, double latHeatSubLoc, double tempAirLoc, double humRelLoc) {
	double e0 = 0.0;
	double ed = 0.0;
	double heatCapAirConv =  heatCapAirLoc / 4.1868;								// 1.005/4.1868 cal/g/oC, see table 3.5, ASCE, this could be done only once?
	double pa = atmPresLoc / 100.0;													// e.g. 1008 mb, see table 2.4 in ASCE Manual
	double lambda = latHeatSubLoc / 4.1868;											// 677.0  cal/g, this could be done only once? - check that parameter is correct

	if (tempAirLoc < -90.0 || humRelLoc < 0.0) {
		std::cout << "\nWet bulb temperature estimation fails";
	}
	else {
		e0 = SVPW(tempAirLoc) / 100.0; // mb
		ed = e0 * humRelLoc  / 100.0; // mb
	}

	double td = -6141.9 / (log(ed / (3.5558E10)));
	double gamma = heatCapAirConv * pa / 0.622 / lambda;							//  mb/oC eq. 7.13 ASCE
	double tmp = (td + tempAirLoc + 273.15) / 2.0 - 273.15 + 237.3;
	double delta = 4098.0 * e0 / (tmp * tmp);										// eq. 7.15 ASCE
	return (gamma * (tempAirLoc + 273.15) + delta * td) / (delta + gamma) - 273.15;	//  eq. 7.19 ASCE
}

double SurfaceSnowModel::SVPW(double tempAirLoc) {
	// Calculates the vapour pressure at a specified temperature over water
	//  using polynomial from Lowe (1977).Taken from UEB (Tarboton and Luce, 1996)	
	const double a = 6.107799961;
	const double b = 0.4436518521;
	const double c = 0.01428945805;
	const double d = 0.0002650648471;
	const double e = 3.031240936e-06;
	const double f = 2.034080948e-08;
	const double h = 6.136820929e-11;

	return 100.0 * (a + tempAirLoc * (b + tempAirLoc * (c + tempAirLoc * (d + tempAirLoc * (e + tempAirLoc * (f + tempAirLoc * h)))))); // Pa
}

double SurfaceSnowModel::SVPI(double tempAirLoc) {
	// Calculates the vapour pressure at a specified temperature over water
	//  using polynomial from Lowe (1977).Taken from UEB (Tarboton and Luce, 1996)	
	const double a = 6.109177956;
	const double b = 0.503469897;
	const double c = 0.01886013408;
	const double d = 0.0004176223716;
	const double e = 5.82472028e-06;
	const double f = 4.838803174e-08;
	const double h = 1.838826904e-10;
	
	return 100.0 * (a + tempAirLoc * (b + tempAirLoc * (c + tempAirLoc * (d + tempAirLoc * (e + tempAirLoc * (f + tempAirLoc * h)))))); // Pa
}

void SurfaceSnowModel::iterate() {
    for (auto &cellPtr : surfaceMesh_->cells) {
        auto cell = dynamic_cast<SurfaceCell*>(cellPtr.get());
        if (!cell)
            continue;
		
		compCanopAerRes(windSpeed);
		compSnowAlbedo(timeStep, tempAir, precSnow);
		compSnowState(timeStep, precWat, precSnow, tempAir, tempWetBulb, windSpeed, humRel, atmPres, radShortIn, radLongIn);
	}
}

void SurfaceSnowModel::compCanopAerRes(double &windSpeed)
{
	if (windSpeed > 0.0) {
		// ! Case 3 or 4: no overstory or overstory is snow-covered
		//double UUnder = UOver; // tuulen nopeus 2 m
		// wind at a height of 2.0 m above ground (snow input reference height)
		//double USnow = UUnder;
		//double ZrefSnow = 2.0; // mittauskorkeus tuulelle 

		//ground is covered with snow 

		//If wind is measured higher than T/RH, estimate wind speed at the height of Zrefsnow
		//tämä silloin jos tuuli on mitattu korkeammalla kuin T & RH

		// Velocity is fixed to the proper elevation
		if (refHeightWind > refHeightTemp) {
			// friction velocity at the top of the canopy
			double uFriction = windSpeed * vonKarmConst / log((refHeightWind - depthSnow) / depthSnowRef);
			// wind speed at the height of Zref 
			windSpeed = uFriction / vonKarmConst * log(refHeightTemp / depthSnowRef);
			//USnow = UOver;
		}
		// Else in the function is a bit odd???????????????????????????
		//else
		//{
		//	// friction velocity 
		//	// VonKarman = 0.4
		//	// Zref = mittauskorkeus maan yläpuolella 2 m
		//	// Z0Snow = lumen karkeuskorkeus , likimain 0.005 m
		//	// Sdepth = lumen syvyys
		//	
		//	Ufriction = UOver * VonKarman / dDenom;
		//}
		// logarithmic part from VegPara%Zref to depthSnow + VegPara%Z0Snow
		// Zref used here instead of ZrefU, is it ok?
		double logPart = log((refHeightTemp - depthSnow) / depthSnowRef);
		canAeroRes = 1.0 / (vonKarmConst * vonKarmConst * windSpeed) * logPart * logPart;
	}
	else {
		canAeroRes = 1.0E20;
	}
}

void SurfaceSnowModel::compSnowAlbedo(double timeStep, double tempAirLoc, double precSnowLoc)
{
	// Snow on ground
	if (watEqSnow[0] > 0.0) {
		double oldAge = snowAge;

		// Update snow surface age in hours
		snowAge += timeStep;

		// Update accumulated event snowfall
		if (precSnowLoc > 0.0 && tempAirLoc < tmpSnowThresMin) {
			snowFallCum = snowFallCum + precSnowLoc * timeStep;
			snowFallPer = 0.0;
		}

		// Update duration of a dry period between snowfall events
		if (precSnowLoc <= 0.0 /*snowDepthMin*/) {
			snowFallPer += timeStep;
		}

		// Check if snowfall event is over.
		if (snowFallPer > albTimeLim) {
			snowFallCum = 0.0;
		}

		// Check if new event have started and a sufficient amount of snowfall has occurred.
		if (snowFallCum > albPrecLim) {
			snowAge = 0.0;
			albedo = albedoMax;
		}

		// YHDISTÄ YLLÄ OLEVAAN!
		if (snowAge <= 0.0 /*dSnowFallLim*/) {
			oldAge = 0.0;
		}

		// Decrease of snow albebo in cold conditions
		// Unit is days in age variable in power function
		if (tempAirLoc < 0.0) {
			double power = pow(snowAge / 24.0, 0.58);
			albedo = albedo + albedoMax * pow(0.94, power);
			power = pow(oldAge / 24.0, 0.58);
			albedo = albedo - albedoMax * pow(0.94, power);
		}
		// Decrease of snow albedo in warm (> 0 oC) conditions
		else {
			double power = pow(snowAge / 24.0, 0.46);
			albedo = albedo + albedoMax * pow(0.82, power);
			power = pow(oldAge / 24.0, 0.46);
			albedo = albedo - albedoMax * pow(0.82, power);
		}

		// Fix minimum albedo
		if (albedo < albedoMin) {
			albedo = albedoMin;
		}
	}
	// Snow cover is not present on the ground
	else {
		albedo = albedoMin;
	} 
}

void SurfaceSnowModel::compSnowState(
	double timeStep,
	double precWat,
	double precSnow,
	double tempAir,
	double tempWetBulb,
	double windSpeed,
	double humRel,
	double pressAtm,
	double radShortIn,
	double radLongIn) {
	// Energy and water fluxes.
	double fluxEnergyMelt[2];
	fluxEnergyMelt[0] = fluxEnergyMelt[1] = 0.0;
	//double energySoil = 0.0; // not used
	double fluxWaterEvap = 0.0;
	double fluxEnergyCond[2];
	fluxEnergyCond[0] = fluxEnergyCond[1] = 0.0;
	double fluxEnergyPrec = 0.0;
	double fluxEnergySens = 0.0;
	double fluxEnergyLat = 0.0;
	double fluxWaterMelt[2];
	fluxWaterMelt[0] = fluxWaterMelt[1] = 0.0;
	double waterLiqDiff[2];
	waterLiqDiff[0] = waterLiqDiff[1] = 0.0;
	
	// Compute snow temperatures in the two layers.
	compSnowTemperature();

	// Derive snow depth.
	depthSnow = 0.0;
	if (densSnow[0] > 0.0 && densSnow[1] > 0.0) {
		depthSnow = (watEqSnow[0] - watLiq[0]) * densWater / densSnow[0] + (watEqSnow[1] - watLiq[1]) * densWater / densSnow[1];
	}

	// Estimate snow surface energy fluxes and snow surface temperature. 
	double radShortOut = 0.0;
	double radLongOut = 0.0;
	compSurfBal(tempAir, humRel, pressAtm, windSpeed, precWat, precSnow, radShortIn, radShortOut, radLongIn, radLongOut, fluxEnergyPrec, fluxEnergySens, fluxEnergyLat, fluxEnergyCond, fluxWaterEvap);

	// Estimate melt water discharged out of the two snow layers.
	compMeltRates(timeStep, precWat, fluxWaterMelt/*, waterLiqDiff*/);

	// Compute energy advected with melt outflow.
	for (int i = 0; i < 2; i++) {
		fluxEnergyMelt[i] = densWater * latHeatFus * fluxWaterMelt[i];
	}
	
	// Estimate snow layer densities.
	compSnowDensity(0, timeStep, precSnow, precWat, tempWetBulb, fluxWaterEvap, waterLiqDiff);

	if (watEqSnow[1] > 0.0) {
		compSnowDensity(1, timeStep, precSnow, precWat, tempWetBulb, fluxWaterEvap, waterLiqDiff);
	}
	else  {
		densSnow[1] = densSnow[0];
	}
	
	// Derive energy balances.


	// Two snow layers present.
	double minQinThinSnow = 0.0;

	if (watEqSnow[1] > 0.0) {
		// Energy is negative in cold snow and positive in wet snow.
		energySnow[0] += timeStep * (radShortIn - radShortOut + radLongIn - radLongOut + fluxEnergyPrec + fluxEnergySens + fluxEnergyLat - fluxEnergyMelt[0] - fluxEnergyCond[1]);
		energySnow[1] += timeStep * (fluxEnergyMelt[0] - fluxEnergyMelt[1] + fluxEnergyCond[1] - fluxEnergySoil);
		
		// Check small snow layer problem (small layer heat capacity).
		if (watEqSnow[1] < 0.001) { // snow layer less than 1 mm, Magic number! 
			if (energySnow[0] < 0 && watEqSnow[0] > 0.0) {
				minQinThinSnow = energySnow[0] / watEqSnow[0] * watEqSnow[1];
			}
			// Thin bottom layer cannot reside cold content per unit depth
			// that is less than cold content in top layer per unit depth
			if (energySnow[1] < minQinThinSnow) {
				energySnow[0] += energySnow[1] - minQinThinSnow;
				energySnow[1] = minQinThinSnow;
			}
		}
	}
	else {
		// One snow layer present
		if (watEqSnow[1] <= 0.0 /*dSnoDepthMin*/ && watEqSnow[0] > 0.0) {
			// A TEST BY LASSI. WHEN SNOW LAYER IS THIN, SOIL HEAT FLUX IS NOT CONSIDERED ANYMORE (INSTABILITIES OCCURED).
			//if (watEqSnow[0] < 0.005)
			//{
			//	fluxEnergySoil = 0.0;
			//}

			energySnow[0] += timeStep * (radShortIn - radShortOut + radLongIn - radLongOut + fluxEnergyPrec + fluxEnergySens + fluxEnergyLat - fluxEnergyMelt[0] - fluxEnergySoil);
			energySnow[1] = 0.0;

			// Check small snow layer problem
			if (watEqSnow[0] < 0.001) { // Magic number!   
				//if (energySoil < 0.0 && dDepthSoil > 0.0)
				//{
				//	minQinThinSnow = energySoil / dDepthSoil * watEqSnow[0];
				//}

				minQinThinSnow = 0.0;

				if (energySnow[0] < minQinThinSnow/* && timeStep > 0.0*/) {
					//dFlxEnSoil += (energySnow[0] - minQinThinSnow) / timeStep;
					energySnow[0] = minQinThinSnow;
				}
			}
		}
		// No snow cover
		else { 
			// Balance over soil surface is disregarded
			energySnow[0] = 0.0;
			energySnow[1] = 0.0;
		}
	}
	
	// Derive mass balance
	// Meltwater is discharged out of the system at the snow/soil interface
	// Upper snow layer
	// There is no snow prior to this time
	if (watEqSnow[0] <= 0.0) {
		// Is there enough mass input to start snow accumulation?
		if (watEqSnow[0] + timeStep * (precWat + precSnow - fluxWaterMelt[0] - fluxWaterEvap) > 0.0 && precSnow > 0.0) {
			watEqSnow[0] += timeStep * (precWat + precSnow - fluxWaterMelt[0] - fluxWaterEvap);

			// Insignificant snow accumulaton is disregarded
			if (watEqSnow[0] < 1.0E-10 && timeStep > 0.0) { // magic numbers! 
				fluxWaterMelt[0] += watEqSnow[0] / timeStep;
				fluxEnergyMelt[0] += watEqSnow[0] * latHeatFus * densWater / timeStep;
				watEqSnow[0] = 0.0;
				watLiq[0] = 0.0;
			}

		}
	}
	// Snow cover is present, check if snow in the top layer disappears all at once
	else  {
		if (watEqSnow[0] + timeStep * (precWat + precSnow - fluxWaterMelt[0] - fluxWaterEvap) > 0.0) {
			watEqSnow[0] += timeStep * (precWat + precSnow - fluxWaterMelt[0] - fluxWaterEvap); 
			// Check if the energy in top snow layer is large enough to melt it all at once
			if (energySnow[0] > watEqSnow[0] * latHeatFus * densWater && timeStep > 0.0) {
				fluxWaterMelt[0] += watEqSnow[0] / timeStep;
				fluxEnergyMelt[0] += watEqSnow[0] * latHeatFus * densWater / timeStep;
				watEqSnow[0] = 0.0;
				watLiq[0] = 0.0;
			}
			// Ignore insignificant amount of snow
			if (watEqSnow[0] < 1.0E-10 && timeStep > 0.0) { // Magic numbers! 
				fluxWaterMelt[0] += watEqSnow[0] / timeStep;
				fluxEnergyMelt[0] += watEqSnow[0] * latHeatFus * densWater / timeStep;
				watEqSnow[0] = 0.0;
				watLiq[0] = 0.0;
			}
		}
		else {
			// All snow in top layer disappeared
			if (timeStep * fluxWaterEvap < watEqSnow[0] + timeStep * (precWat + precSnow) && timeStep > 0.0) {
				// Top layer is lost through both evaporation (if positive) and melt
				// values below are replaced totally, is it ok?
				fluxWaterMelt[0] = watEqSnow[0] / timeStep + precWat + precSnow - fluxWaterEvap;
				fluxEnergyMelt[0] = fluxWaterMelt[0] * latHeatFus * densWater;
			}
			// top layer evaporates all at once
			else if (timeStep > 0.0) {
				fluxWaterEvap = watEqSnow[0] / timeStep + precWat + precSnow;
			}
			watEqSnow[0] = 0.0;
			watLiq[0] = 0.0;
		}
	} 
	// if no snow prior to this time
	// Lower snow layer:
	// Check if the lower snow layer is lost entirely
	if (watEqSnow[1] + timeStep * (-fluxWaterMelt[1]) > 0.0) {
		watEqSnow[1] += timeStep * (fluxWaterMelt[0] - fluxWaterMelt[1]);
		if (watEqSnow[1] < 1.0E-10 && timeStep > 0.0) { // Magic number? 
			fluxWaterMelt[1] += fluxWaterMelt[0] + watEqSnow[1] / timeStep;
			fluxEnergyMelt[1] += watEqSnow[1] * latHeatFus * densWater / timeStep;
			watEqSnow[1] = 0.0;
			watLiq[1] = 0.0;
		}
	}
	//! lower snow layer is lost
	else {
		if (watEqSnow[1] > 0.0 && timeStep > 0.0) {
			fluxWaterMelt[1] = fluxWaterMelt[0] + watEqSnow[1] / timeStep;
			fluxEnergyMelt[1] = fluxWaterMelt[1] * latHeatFus * densWater;
			watEqSnow[1] = 0.0;
		}
		// the lower snow layer does not exist
		else {
			fluxWaterMelt[1] = fluxWaterMelt[0];
			watEqSnow[1] = 0.0;
			watLiq[1] = 0.0;
		}
	}
	
	// Snow layer mass and energy redistribution
	// Two layers present
	if (watEqSnow[0] + watEqSnow[1] > snoWatEqTopMax) {
		// Top snow layer has gained mass
		if (watEqSnow[0] > snoWatEqTopMax && snoWatEqTopMax > 0.0) {
			energySnow[1] += (watEqSnow[0] - snoWatEqTopMax) / watEqSnow[0] * energySnow[0];
			energySnow[0] -= (watEqSnow[0] - snoWatEqTopMax) / watEqSnow[0] * energySnow[0];
			watEqSnow[1] += watEqSnow[0] - snoWatEqTopMax;
			watEqSnow[0] = snoWatEqTopMax;
		}
		// Top snow layer has lost mass
		else if (watEqSnow[1] > 0.0) {
			// This includes the case when top layer has melted and has to be reformed
			energySnow[0] += (watEqSnow[0] - snoWatEqTopMax) / watEqSnow[1] * energySnow[1];
			energySnow[1] -= (watEqSnow[0] - snoWatEqTopMax) / watEqSnow[1] * energySnow[1];
			watEqSnow[1] += watEqSnow[0] - snoWatEqTopMax;
			watEqSnow[0] = snoWatEqTopMax;
		}
	}
	// Single layer present
	else {
		watEqSnow[0] += watEqSnow[1];
		energySnow[0] += energySnow[1];
		watEqSnow[1] = 0.0;
		energySnow[1] = 0.0;
	}

	// Save precipitation
	precCum[0] += precWat * timeStep;
	precCum[1] += precSnow * timeStep;

	// Save outflux to overland domain
	fluxEvapOut += fluxWaterEvap * timeStep;
	fluxWaterOut += fluxWaterMelt[1] * timeStep;

	
	// Test print variables - use only a single core
	sumTime += timeStep;
}

void SurfaceSnowModel::compSurfBal(
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
	double &fluxWaterEvap
	) {
	// Evaluate snow surface Temperature, when ground is covered with snow.
	double tempSurfLoc;

	if (depthSnow > 0.0) {
		// Initial value for snow surface temperature 
		tempSurfLoc = tempAir;
		double tempSurfLocOld;
		int iterCount = 0;
		double tempDelta = -1.0E-1; // Magic number! -1.0E-1

		// Iteration loop for estimating snow surface temperature and energy fluxes
		// The root is solved with the Newton-Raphson method
		do {
			// Kysy tästä?
			//if (iterCount == iIterCountMax - 50)
			//{
			//	tempDelta = -tempDelta;    
			//}

			// Compute the new temperature
			tempSurfLocOld = tempSurfLoc;
			double res = compSurfEnBal(tempAir, humRel, pressAtm, windSpeed, precWat, precSnow, tempSurfLoc, radShortIn, radShortOut, radLongIn, radLongOut, fluxEnergyPrec, fluxEnergySens, fluxEnergyLat, fluxEnergyCond, fluxWaterEvap);
			double func = res;
			res = compSurfEnBal(tempAir, humRel, pressAtm, windSpeed, precWat, precSnow, tempSurfLoc + tempDelta, radShortIn, radShortOut, radLongIn, radLongOut, fluxEnergyPrec, fluxEnergySens, fluxEnergyLat, fluxEnergyCond, fluxWaterEvap);
			if (tempDelta != 0.0) {
				double deriv = (res - func) / tempDelta;

				// Calculate a new estimate of the snow surface temperature
				if (deriv != 0) {
					tempSurfLoc = tempSurfLoc - 0.5 * func / deriv; // Magic number 0.5?
				}
			}
			
			// Next iteration round
			iterCount = iterCount + 1;
		} while (fabs(tempSurfLoc - tempSurfLocOld) > iterConvThr && iterCount < iterCountMax);

		// Temperature is forced to 0 oC during snowmelt
		if (tempSurfLoc > 0.0 && watEqSnow[0] > 0.0) {
			tempSurfLoc = 0.0;
		}

		// Make a correction in a case of an unrealistic temperature fluctuation
		if (fabs(tempSurfLoc - tempAir) > 30.0) { // Magic number? 
			tempSurfLoc = tempAir;
		}
	}
	// No snow, soil surface temperature is taken as air temperature.
	else {
		tempSurfLoc = tempAir;
	}
	
	// Update surface energy fluxes after estimating surface temperature dTempSurf
	// NOTE: In the case of a bare ground, energy fluxes are computed using dTempSurf = tempAir,
	// but energy fluxes across (soil) surface are not balanced in this case.
	double fun = compSurfEnBal(tempAir, humRel, pressAtm, windSpeed, precWat, precSnow, tempSurfLoc, radShortIn, radShortOut, radLongIn, radLongOut, fluxEnergyPrec, fluxEnergySens, fluxEnergyLat, fluxEnergyCond, fluxWaterEvap);
	(void)fun; // check this

	// Save the solved surface temperature to the class variable
	tempSurf = tempSurfLoc;
}

double SurfaceSnowModel::compSurfEnBal(
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
	double &fluxWaterEvap
) {
	// Reflected outgoing short-wave radiation
	radShortOut = radShortIn * albedo;	
	
	// Estimate outgoing long-wave radiation
	double temp = std::max(tempSurfLoc, -tempZeroKel) + tempZeroKel;
	radLongOut = bolzConst * emisSnow * temp * temp * temp * temp + (1.0 - emisSnow) * radLongIn;
	
	// Estimate energy advected with precipitation
	compPrecEnerFlux(tempAir, precSnow, precWat, fluxEnergyPrec);

	// Estimate turbulent fluxes of latent and sensible heat
	compTurbFluxes(tempSurfLoc, tempAir, humRel, pressAtm, windSpeed, fluxEnergyLat, fluxEnergySens, fluxWaterEvap);

	// Estimate conducted heat, siirrä tämä ennen turbulent fluxes laskentaa
	compCondFluxes(tempSurfLoc, fluxEnergyCond);
	
	return radShortIn - radShortOut + radLongIn - radLongOut + fluxEnergyPrec + fluxEnergySens + fluxEnergyLat - fluxEnergyCond[0];
}

void SurfaceSnowModel::compPrecEnerFlux(double tempAir, double precSnow, double precWat, double &fluxEnergyPrec) {
	double tempAirMin = 0.0;
	double tempAirMax = 0.0;
	
	if (tempAir < 0.0) {
		tempAirMin = tempAir;
	}
	if (tempAir > 0.0) {
		tempAirMax = tempAir;
	}
	
	fluxEnergyPrec = precSnow * densWater * heatCapIce * tempAirMin + precWat * (latHeatFus * densWater + heatCapWater * densWater * tempAirMax);
}

void SurfaceSnowModel::compTurbFluxes(double tempSurfLoc, double tempAir, double humRel, double pressAtm, double windSpeed, double &fluxEnergyLat, double &fluxEnergySens, double &fluxWaterEvap) {
	if (depthSnow > 0.0) { 
		// Estimate turbulent exchange coefficients
		double traCofSensHeat = 0.0;
		double traCofLatHeat = 0.0;
		
		compTurbCoeff(windSpeed, tempAir, tempSurfLoc, traCofSensHeat, traCofLatHeat);
		
		// VAIHDA dPresAirBase KORJATTU VERSIO PASCALEINA - LASKETTU MALLISSA
		double dDensAir = pressAtm / (idealGasConst * (tempAir + tempZeroKel)); //! kg/m3 UEB, why is this computed here?
		//! Sensible heat flux:
		//! The suggested windless coefficient of exchange is from Jordan 1991: 2.0 W/m2/oK

		//traCofSensHeat ~600 m/d
		fluxEnergySens = (traCofSensHeat * dDensAir * heatCapAir + convCofSensHeat) * (tempAir - tempSurfLoc); //! kJ/m2/h
		
		//! UEB model code claims that measurement is usually taken with respect to water
		//! Leena Tanskanen from Vaisala OYj assures that relative humidity is always
		//!  measured with respect to water. What does this mean? Does it mean that there are
		//!  two Esat(T)-functions in cold temperatures? Here, in any case, sat vapor pressure 
		//!  over water is used both in cold temperatures and in temp above 0 oC }

		double humAirSpec = SVPW(tempAir) * humRel / 100.0; //! Pa  RH = relative humidity????????????? eair
		double humSurfSpec; // esurf
		
		if (tempSurfLoc < 0.0) {
			humSurfSpec = SVPI(tempSurfLoc);
		}
		else {
			humSurfSpec = SVPW(tempSurfLoc);
		}

		double latHeatSubLoc;
		if (tempSurfLoc > 0.1) { // Magic number?
			latHeatSubLoc = latHeatSub - latHeatFus;
		}
		else {
			if (tempSurfLoc < 0.0) {
				latHeatSubLoc = latHeatSub;
			}
			else {
				latHeatSubLoc = latHeatSub - tempSurfLoc * 10.0 * latHeatFus; // Magic number?
			}
		}

		//! Latent heat flux:
		//!  The suggested windless coefficient of exchange is in Jordan 1991:
		//!  2.0 W/m2/mb, but this correction is suggested not to be used, see Jordan 1992 
		fluxEnergyLat = ((traCofLatHeat * latHeatSubLoc) * 0.622 / (idealGasConst * (tempAir + tempZeroKel)) + convCofLatHeat) * (humAirSpec - humSurfSpec);

		// Sublimation / evaporation water flux
		fluxWaterEvap = -fluxEnergyLat / (densWater * latHeatSubLoc);
	}
	else { //! no snow
		//! Fluxes are not computed during summer, because surface temperature computation 
		//! during the summer is complex: We should separate canopy temperatures and 
		//! soil surface T, as well as the efect of humus layer

		//!		Kh = 1.0/(Rau) * 3600.0 * 24.0 ! m/d
		//!		AirDns = P0 / (Rd * (TAir + tempZeroKel)) ! kg/m3 UEB
		//! Sensible heat flux:
		//!		QH = (Kh * AirDns * heatCapAir) * (Ta - Tsrf) ! kJ/m2/d
		//!		QE = -Eau * (WaterDns * (Hv - Hf))
		//! Fluxes from surface are set to zero
		fluxEnergyLat = 0.0;
		fluxEnergySens = 0.0;
		fluxWaterEvap = 0.0;
	}
}

void SurfaceSnowModel::compTurbCoeff(double windSpeed, double tempAir,  double tempSurfLoc, double &traCofSensHeat, double &traCofLatHeat) {
	//! Assume logarithmic wind profile and neutral atmospheric stability
	double kneutral = 0.0;
	if (canAeroRes != 0.0) {
		kneutral = 1.0 / canAeroRes; //! m/s
		kneutral = kneutral * 3600.0; //! m/h
	}

	double kadj;
	if (windSpeed > 0.0) {
		//! Correction for atmospheric stability is carried out for the total rss+ru+ra

		//! Estimate Richardson number 
		//! NOTE: unlimited iteration of surface temperature may lead to iteration
		//!   values lower than -273 K
		double richNum = gravAcc * (tempAir - tempSurfLoc) * (refHeightTemp - depthSnow) / (windSpeed * windSpeed * (std::max(0.5 * (tempAir + tempSurfLoc), -tempZeroKel) + tempZeroKel));

		if (richNum > richNumLimMax) {
			richNum = richNumLimMax;
		}
		if (richNum < richNumLimMin) {
			richNum = richNumLimMin;
		}

		//! This is the procedure based on Storck (2000)
		double riu = 1.0 / (log(refHeightTemp / depthSnowRef) + 5.0);
		if (riu > richNumLimMax) {
			riu = richNumLimMax;
		}
		if (richNum > riu) {
			richNum = riu;
		}

		double raAdj; 
		if (richNum < 0.0) { //! unstable 
			raAdj = canAeroRes / (pow(1.0 - 16.0 * richNum, 0.5));
		}
		else {
			if (richNum > 0.0 && richNumLimMax != 0.0 && richNumLimMax != richNum) { //! stable
				double tmp = 1.0 - richNum / richNumLimMax;
				raAdj = canAeroRes / (tmp * tmp);
			}
			else {
				raAdj = canAeroRes;
			}
		}

		kadj = 3600.0 / raAdj; //! m/h
	}
	else {
		kadj = kneutral;
	}

	// These are the same????????????????????????????????????
	traCofSensHeat = kneutral + atmStabCof * (kadj - kneutral); //! Turbulent transfer coefficient for sensible heat flux
	traCofLatHeat = kneutral + atmStabCof * (kadj - kneutral); //! Turbulent transfer coefficient for latent heat flux
}

void SurfaceSnowModel::compCondFluxes(double tempSurfLoc, double *fluxEnergyCond) {
	// Estimate snow thermal conductivity
	double tmp0 = 0.0;
	double tmp1 = 0.0;
	double condSnow[2];
	condSnow[0] = condSnow[1] = 0.0;
	if (densWater > 0.0) {
		tmp0 = densSnow[0] / densWater;
		tmp1 = densSnow[1] / densWater;
	}

	// Anderson (1976) equations for snow thermal conductivity [kJ h-1 m-1 oC-1].
	condSnow[0] = condSnowPara[0] + condSnowPara[1] * tmp0 * tmp0;
	condSnow[1] = condSnowPara[0] + condSnowPara[1] * tmp1 * tmp1;
	
	// compute conductive fluxes
	// Two snow layers present
	if (watEqSnow[1] > 0.0) {
		double depth[2];
		depth[0] = depth[1] = 0.0;
		
		if (densSnow[0] > 0.0) {
			depth[0] = 0.5 * (watEqSnow[0] - watLiq[0]) * densWater / densSnow[0];
		}

		if (densSnow[1] > 0.0) {
			depth[1] = 0.5 * (watEqSnow[1] - watLiq[1]) * densWater / densSnow[1];
		}
		
		double dCondAver = 0.0;
		if (condSnow[0] > 0.0 && condSnow[1] > 0.0 && depth[0] + depth[1] > 0.0) {
			dCondAver = (depth[0] + depth[1]) / (depth[0] / condSnow[0] + depth[1] / condSnow[1]);
		}

		if (depth[0] > 0.0) {
			fluxEnergyCond[0] = condSnow[0] * (tempSurfLoc - tempSnow[0]) / depth[0];
		}

		if (depth[0] + depth[1] > 0.0) {
			fluxEnergyCond[1] = dCondAver * (tempSnow[0] - tempSnow[1]) / (depth[0] + depth[1]);
		}

		// Save snow variables
		depthSnow = 2.0 * (depth[0] + depth[1]); // where does the multiplier (2) come from?
		condSnowAvr = dCondAver;
		//condSnowLower = condSnow[1];
		//depthSnowLower = depth[1];
	}
	// One snow layer present
	else {
		if (watEqSnow[1] <= 0.0 /*snowDepthMin*/ && watEqSnow[0] > 0.0) {		
			double depth = 0.0;
			
			if (densSnow[0] > 0.0) {
				depth = 0.5 * (watEqSnow[0] - watLiq[0]) * densWater / densSnow[0];
			}
			
			if (depth > 0.0) {
				fluxEnergyCond[0] = condSnow[0] * (tempSurfLoc - tempSnow[0]) / depth;
			}

			fluxEnergyCond[1] = 0.0;

			// Save snow variables
			depthSnow = 2.0 * depth; // where does the multiplier (2) come from?
			condSnowAvr = condSnow[0];
			//condSnowLower = condSnow[0];
			//depthSnowLower = dDepth;
		}
		// No snow 
		else {
			fluxEnergyCond[0] = 0.0;
			fluxEnergyCond[1] = 0.0;
			
			// Save snow variables
			depthSnow = 0.0;
			condSnowAvr = 0.0;
			//condSnowLower = 0.0;
			//depthSnowLower = 0.0;
		}
	}
}

void SurfaceSnowModel::compSnowTemperature() {
	// Loop over the two snow layers
	for (int i = 0; i < 2; i++) {
		if (watEqSnow[i] > 0.0) {
			// energySnow is equivalent to water that is needed to freeze and heat snowpack
			double tempSnowOld = tempSnow[i];

			if (energySnow[i] > 0.0) {
				tempSnow[i] = 0.0;
			}
			else if (watEqSnow[i] > 0.0 && heatCapIce > 0.0 && densWater > 0.0) {
				tempSnow[i] = energySnow[i] / watEqSnow[i] / heatCapIce / densWater;
			}

			// It is assured that max temperature change in snow is tempSnowMaxChange oC
			if (tempSnow[i] > tempSnowOld + tempSnowMaxChange) {
				tempSnow[i] = tempSnowOld + tempSnowMaxChange;
			}
			if (tempSnow[i] < tempSnowOld - tempSnowMaxChange) {
				tempSnow[i] = tempSnowOld - tempSnowMaxChange;
			}
		}
		// No snow cover on the ground
		else {
			tempSnow[i] = 0.0;
		}
	}
}

void SurfaceSnowModel::compMeltRates(double timeStep, double precWat, double *fluxWatMelt/*, double *waterLiqDiff*/) {
	// Snow cover on ground
	if (watEqSnow[0] > 0.0) {
		fluxWatMelt[0] = fluxWatMelt[1] = 0.0;
		
		//for (int i = 0; i < 2; i++)
		//{
		//	fluxWatMelt[i] = 0.0;
		//	//waterLiqDiff[i] = watLiq[i];
		//}

		for (int i = 0; i < 2; i++) {
			// Compute the energy advected with melt outflow from the top snow layer to the layer below
			double energyMeltTop = 0.0;

			if (watEqSnow[i] > 0.0 && i == 1) {
				energyMeltTop = timeStep * fluxWatMelt[0] * densWater * latHeatFus;
			}

			// There is liquid water held in snow layer.
			if (watLiq[i] > 0.0 || energyMeltTop + energySnow[i] > 0.0 ) { // snowE
				// Estimate liquid water mass fraction
				double watLiqFrac = 0.0;

				if (watEqSnow[i] < 1.0E-10) { // Magic number! 
					watLiqFrac = 0.99; // Magic number!
				}
				else {
					if (energyMeltTop + energySnow[i] > 0.0) {
						if ((energyMeltTop + energySnow[i]) / (densWater * latHeatFus) < watEqSnow[i] && densWater * latHeatFus > 0.0) {
							watLiq[i] = (energyMeltTop + energySnow[i]) / (densWater * latHeatFus);
						}
						else {
							watLiq[i] = watEqSnow[i];
						}
					}
					else {
						watLiq[i] = 0.0;
					}
					
					if (watLiq[i] < 0.0) {
						watLiq[i] = 0.0;
					}
					
					if (watEqSnow[i] > 0.0 && watLiq[i] / watEqSnow[i] > 1.0) {
						watLiq[i] = watEqSnow[i]; //! All is liquid
					}
					
					// Change in liquid water content due to energy content change
					//waterLiqDiff[i] = watLiq[i] - waterLiqDiff[i];

					// Estimate liquid water mass fraction
					if (watEqSnow[i] > 0.0) {
						watLiqFrac = watLiq[i] / watEqSnow[i];
					}
					
					if (watLiqFrac > 0.99) {
						watLiqFrac = 0.99;
					}
				}
				
				// Estimate snow layer saturation (volumetric)
				double dSnowSat = 0.0;
				if (1.0 - watLiqFrac > 0.0 && watLiqFrac / (1.0 - watLiqFrac) < snowLiqCap) {
					dSnowSat = 0.0;
				}
				else if (1.0 - watLiqFrac > 0.0 && densSnow[i] > 0.0 && densIce > 0.0 && densWater / densSnow[i] - densWater / densIce - snowLiqCap > 0.0) {
					dSnowSat = (watLiqFrac / (1.0 - watLiqFrac) - snowLiqCap) / (densWater / densSnow[i] - densWater / densIce - snowLiqCap);
				}
				if (dSnowSat > 1.0) {
					dSnowSat = 1.0;
				}
				
				// check this
				// Estimate melt outflow (unit mm / 1 h)
				fluxWatMelt[i] = condSatSnow * dSnowSat * dSnowSat * dSnowSat; //! m/h
				if (fluxWatMelt[i] > condSatSnow) {
					fluxWatMelt[i] = condSatSnow;
				}
				// Melt water volume exceeds available liquid water
				if (timeStep > 0.0 && fluxWatMelt[i] > watLiq[i] / timeStep/* - snowDepthMin*/) {
					fluxWatMelt[i] = watLiq[i] / timeStep;
					watLiq[i] = 0.0;
				}
				else {
					watLiq[i] = watLiq[i] - fluxWatMelt[i] * timeStep;
				}
			}
			// Snowpack temperature is negative, no liquid water
			else {
				watLiq[i] = 0.0;
				//waterLiqDiff[i] = watLiq[i] - waterLiqDiff[i];
			}
		}
		// If only one snow layer present, set melt outflow variables equal
		if (watEqSnow[1] <= 0.0 /* snowDepthMin*/) {
			fluxWatMelt[1] = fluxWatMelt[0];
		}
	}
	 // If there is no snow on the ground
	else {
		fluxWatMelt[0] = precWat;
		fluxWatMelt[1] = fluxWatMelt[0];
	}
}

void SurfaceSnowModel::compSnowDensity(int layer, double timeStep, double precSnow, double precWat, double tempWetBulb, double &fluxWatEvap, double *watLiqDiff) {
	//densSnowLim = 0.15; //dnslimit = parameters[63];
	//densSnowMax = 500.0; //MaxSnowDns = parameters[32];	
	
	//! Snow methamorphism routine is from Anderson (1976).
	//! Snow compaction estimate is based on the weight of the snow layer above 
	//! plus 50% of the current snow layer weight 
	//! Set constants
	double snowMetaParaLoc;
	if (watLiq[layer] > 0.0) {
		snowMetaParaLoc = snowMetaPara[4]; // local variable version!!!!!!!!!!!!!!!!
	}
	else {
		snowMetaParaLoc = 1.0; //! 1.0 -- 2.0, Multiplied by this if liquid water is present
	}

	//! Wet bulb temperature used as a temperature estimate in Anderson (1976) equation
	double newSnowDensity;
	if (tempWetBulb + tempZeroKel > 258.16) {
		newSnowDensity = 1000.0 * (0.05 + 0.0017 * (tempWetBulb + tempZeroKel - 258.16) * sqrt(tempWetBulb + tempZeroKel - 258.16));
	}
	else {
		newSnowDensity = 50.0; //! kg/m3
	}

	if (watEqSnow[layer] > 0.0) {
		//! Only snowfall and evaporation from snow (during snowfall) results in density change, not rain
		//! Estimate evaporation from snow
		double evapSnow; // = watEqSnow[layer] * 0.5;

		if (layer == 0) {
			evapSnow = watEqSnow[layer] * 0.5;
		}
		else {
			evapSnow = watEqSnow[layer] * 0.5 + watEqSnow[layer - 1];
		}

		fluxWatEvap = fluxWatEvap - precWat;
		if (fluxWatEvap < 0.0) {
			fluxWatEvap = 0.0; //! It is assumed that condensation does not change snow density
		}
		if (precSnow > 0.0 && precSnow > fluxWatEvap && layer == 0) { // Warningm the same variable used here!!!!!!!!!!!!!
			depthSnow = (watEqSnow[layer] - watLiq[layer]) * densWater / densSnow[layer] + timeStep * (precSnow - fluxWatEvap) * densWater / newSnowDensity;
			densSnow[layer] = (watEqSnow[layer] - watLiq[layer] + timeStep * (precSnow - fluxWatEvap)) / depthSnow * densWater;
		}

		// Density change due to liquid water freezing in snow
		if (watLiqDiff[layer] < 0.0 && watEqSnow[layer] > 0.0) {
			// Refreezing of liquid water increases snow density.
			densSnow[layer] = (watEqSnow[layer] - watLiq[layer]) / ((watEqSnow[layer] - watLiq[layer] + watLiqDiff[layer]) / densSnow[layer] - watLiqDiff[layer] / densIce);
		}
		
		// Snow density increases with time.
		// Density change due to compaction (following Anderson 1976).
		double tempSnowLoc = tempSnow[layer];
		if (tempSnowLoc > 0.0) {
			tempSnowLoc = 0.0;
		}

		// Density change due to destructive metamorphism (following Anderson 1976).
		double Change1 = densSnow[layer] * timeStep * (snowMetaPara[0] * exp(-0.08 * (snowMetaPara[5] - tempSnowLoc)) * (100.0 * evapSnow) * exp(-snowMetaPara[1] * densSnow[layer] / 1000.0));
		double Change2;
		if (densSnow[layer] > 1000.0 * densSnowLim) {
			Change2 = densSnow[layer] * timeStep * (snowMetaPara[2] * exp(-snowMetaPara[3] * (snowMetaPara[5] - tempSnowLoc)) * exp(-46.0 * (densSnow[layer] / 1000.0 - densSnowLim))) * snowMetaParaLoc;
		}
		else {
			Change2 = densSnow[layer] * timeStep * (snowMetaPara[2] * exp(-snowMetaPara[3] * (snowMetaPara[5] - tempSnowLoc))) * snowMetaParaLoc;
		}
		densSnow[layer] = densSnow[layer] + Change1 + Change2;
	}
	// Density changes due to precipitation input
	else {
		densSnow[layer] = newSnowDensity;
	}

	if (densSnow[layer] > densSnowMax) {
		densSnow[layer] = densSnowMax;
	}
}

void SurfaceSnowModel::postprocess() {
	
}
