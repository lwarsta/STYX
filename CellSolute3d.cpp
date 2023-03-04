#include "CellSolute3d.h"

CellSolute3d::CellSolute3d()
{
    geom = 0;
    water = 0;
	conc = 0.0;
	concOld = 0.0;
	mass = 0.0;
	massOld = 0.0;
	//disp[0] = disp[1] = disp[2] = 0.0;
	//for (int i = 0; i < 6; i++)
	//{
	//	dispFlux[i] = 0.0;
	//	dispInter[i] = 0.0;
	//	advFlux[i][0] = advFlux[i][1] = 0.0;
	//}
	alphaCofLong = 0.0;
	alphaCofTrans = 0.0;
	molDiffRate = 0.0;
	decayRate[0] = decayRate[1] = 0.0;
	decayRateInit = 0.0;
	decayTarget[0] = decayTarget[1] = -1;
	adsParam[0] = adsParam[1] = adsParam[2] = 0.0;
	ads = 0.0;
	adsMult = 0.0;
	infConc = 0.0;
	infMass = 0.0;
	drainMass = 0.0;
	decayMassOut[0] = decayMassOut[1] = 0.0;
	decayMassIn[0] = decayMassIn[1] = 0.0;
}

void CellSolute3d::setNeighPointers(
	std::vector<CellSolute3d*> neighPointersNew)
{
    neighPointers = neighPointersNew;
    disp.assign(neighPointers.size(), 0.0); // THIS SHOULD ONLY SIZE OF THREE CURRENTLY!
    dispFlux.assign(neighPointers.size(), 0.0);
    dispInter.assign(neighPointers.size(), 0.0);
	for (size_t i = 0; i < neighPointers.size(); i++)
	{ 
		std::vector<double> innerVec;
		innerVec.assign(2, 0.0);
		advFlux.assign(neighPointers.size(), innerVec);
	}
}

void CellSolute3d::calcAdsMult()
{
	if (water != 0 && water->getWatCont() > 0.0 && geom != 0 && conc > 0.0)
	{
		// Freundlich isotherm.
		if ((int)adsParam[0] == 0)
		{
			adsMult = 1.0 + (water->getDryWeight() / water->getWatCont()) * 
				(adsParam[1] * adsParam[2] * pow(conc, adsParam[2] - 1.0));
		}
		// Langmuir isotherm.
		else if ((int)adsParam[0] == 1)
		{
			if (adsParam[1] >= 0 && adsParam[2] >= 0)
			{
				adsMult = 1.0 + (water->getDryWeight() / water->getWatCont()) * 
					((adsParam[1] * adsParam[2]) / pow(1 + adsParam[1] * conc, 2));
			}
		}
		else
		{
			adsMult = 1.0;
		}
	}
	else
	{
		adsMult = 1.0;
	}
}
/*
void CellSolute3d::calcAds()
{

	if (water != 0 && water->getWatCont() > 0.0 && geom != 0)
	{
		// Freundlich isotherm.
		if ((int)adsParam[0] == 0)
		{
			ads = adsParam[1] * pow(conc, adsParam[2]) * conc * geom->getVolume() * 
			    water->getDryWeight();
		}
		// Langmuir isotherm.
		else if ((int)adsParam[0] == 1)
		{

			ads = (adsParam[1] * adsParam[2] * conc) / (1 + adsParam[1] * conc);
		}
		else
		{
			ads = 0.0;
		}
	}
}
*/
void CellSolute3d::calcDisp()
{
	if (water != 0)
	{
		// Calculate mean velocity.
		double velMean = 0.0;

		for (size_t i = 0; i < neighPointers.size(); i++)
		{
			velMean += water->getVelInt(i) * water->getVelInt(i);
		}

		velMean = sqrt(velMean);

		if (velMean > 0.0)
		{
			// Calculate average directional velocity components.
			// Arithmetic mean. SHOULD THE SUM BE SUBSTRACTION INSTEAD?
			double velMeanX = 0.5 * (fabs(water->getVelInt(2)) + 
				fabs(water->getVelInt(4))); // ADAPT FOR GENERIC GEOMETRY - DIRECTIONS CAN VARY!
			double velMeanY = 0.5 * (fabs(water->getVelInt(1)) + 
				fabs(water->getVelInt(3))); // ADAPT FOR GENERIC GEOMETRY - DIRECTIONS CAN VARY!
			double velMeanZ = 0.5 * (fabs(water->getVelInt(0)) + 
				fabs(water->getVelInt(5))); // ADAPT FOR GENERIC GEOMETRY - DIRECTIONS CAN VARY!
			// Calculate directional dispersion components.
			disp[0] = alphaCofLong * velMeanX * velMeanX / velMean + alphaCofTrans * 
				velMeanY * velMeanY / velMean  + alphaCofTrans * velMeanZ * velMeanZ / 
				velMean + molDiffRate;
			disp[1] = alphaCofLong * velMeanY * velMeanY / velMean + alphaCofTrans * 
				velMeanX * velMeanX / velMean  + alphaCofTrans * velMeanZ * velMeanZ / 
				velMean + molDiffRate;
			disp[2] = alphaCofLong * velMeanZ * velMeanZ / velMean + alphaCofTrans * 
				velMeanX * velMeanX / velMean  + alphaCofTrans * velMeanY * velMeanY / 
				velMean + molDiffRate;
		}
		else
		{
			disp[0] = disp[1] = disp[2] = molDiffRate;
		}
	}
	/*
    if (water != 0)
	{
		// Calculate mean velocity.
		double velMean = 0.0;
		for (size_t i = 0; i < neighPointers.size(); i++)
		{
			velMean += water->getVelInt(i) * water->getVelInt(i);
		}
		velMean = sqrt(velMean);
		// Calculate directional dispersion components.
        if (velMean > 0.0)
		{
            for (size_t i = 0; i < neighPointers.size(); i++)
            {
                disp.at(i) = alphaCofLong * water->getVelInt(i) * water->getVelInt(i) / velMean;
                for (size_t j = 0; j < neighPointers.size(); j++)
                {
                    if (j != i)
                    {
                        disp.at(i) += alphaCofTrans * water->getVelInt(j) * 
						    water->getVelInt(j) / (neighPointers.size() - 1) / velMean;
                    }
                }
                disp.at(i) += molDiffRate;
            }
		}
		else
		{
            for (size_t i = 0; i < neighPointers.size(); i++)
            {
                disp.at(i) = molDiffRate;
            }
		}
	}
	*/
}

void CellSolute3d::calcDispInter()
{
	int dir = -1;

	for (size_t i = 0; i < neighPointers.size(); i++)
	{
		if (i == 0 || i == 5) // ADAPT FOR GENERIC GEOMETRY - DIRECTIONS CAN VARY!
		{
			dir = 2;
		}
		else if (i == 1 || i == 3) // ADAPT FOR GENERIC GEOMETRY - DIRECTIONS CAN VARY!
		{
			dir = 1;
		}
		else if (i == 2 || i == 4) // ADAPT FOR GENERIC GEOMETRY - DIRECTIONS CAN VARY!
		{
			dir = 0;
		}
		if (neighPointers.at(i) != 0)
		{
			// Arithmetic mean.
			dispInter[i] = 0.5 * (neighPointers.at(i)->getDisp(dir) + disp[dir]);
			// Geometric mean.
			//dispInter[i] = sqrt(cellNeigh[i]->getD(dir) * dD[dir]);
		}
	}
}

void CellSolute3d::calcDispFluxes()
{
	if (water != 0 && geom != 0)
	{
		for (size_t i = 0; i < neighPointers.size(); i++)
		{
			if (neighPointers.at(i) != 0)
			{
				CellWater3d * neighWater = neighPointers.at(i)->getWater();
				if (neighWater != 0)
				{
					// Arithmetic mean.
					double watContAvg = 0.5 * (water->getWatCont() + 
						neighWater->getWatCont());
					if (geom->getDistance(i) > 0.0)
					{
						dispFlux[i] = dispInter[i] * geom->getFaceArea(i) * 
							watContAvg / geom->getDistance(i);
					}
				}
			}
		}
	}
}

void CellSolute3d::calcAdvFluxes()
{
	if (geom != 0 && water != 0)
	{
		for (size_t i = 0; i < neighPointers.size(); i++)
		{
			double velInt = water->getVelInt(i);
			// Advection out of the cell.
			if (velInt < 0.0)
			{
				advFlux[i][0] = 0.0;
				advFlux[i][1] = geom->getFaceArea(i) * velInt * water->getWatCont();
			}
			// Advection into the cell.
			else
			{
				if (neighPointers.at(i) != 0)
				{
					CellWater3d * neighWater = neighPointers.at(i)->getWater();
					if (neighWater != 0)
					{
						advFlux[i][0] = geom->getFaceArea(i) * velInt * neighWater->getWatCont();
						advFlux[i][1] = 0.0;
					}
				}
			}
		}
	}
}

double CellSolute3d::getConcNeigh(int index)
{
	if (neighPointers.at(index) != 0)
	{
		return neighPointers.at(index)->getConc();
	}
	else
	{
		return 0.0;
	}
}

double CellSolute3d::getConcNeighOld(int index)
{
	if (neighPointers.at(index) != 0)
	{
		return neighPointers.at(index)->getConcOld();
	}
	else
	{
		return 0.0;
	}
}

void CellSolute3d::convertMassToConc()
{
	if (water != 0 && geom != 0)
	{
		if (water->getWatCont() * geom->getVolume() > 0.0)
		{
			conc = mass / (water->getWatCont() * geom->getVolume());
			mass = 0.0;
		}
		else
		{
			conc = 0.0;
		}
	}
}

void CellSolute3d::convertConcToMass()
{
	if (water != 0 && geom != 0)
	{
		if (conc * water->getWatCont() * geom->getVolume() > 0.0) // is this needed?
		{
			// CHECK IF MASS TERM ON THE RIGHT SIDE COULD BE REMOVED?
			mass = mass + conc * water->getWatCont() * geom->getVolume();
		}
	}
}
/*
// NOT READY!
void CellSolute3d::calcInfConc()
{
	if (water != 0 && surfSolute != 0)
	{
		double infVol = water->getFluxInf();
		CellSurfWaterAgri * surfWater = surfSolute->getWater();
		double dSurfVol = surfWater->getWatVol();
		concInf = 0.0;
		if (infVol >= 0.0 && dSurfVol >= 0.0 && infVol + dSurfVol > 0.0 && surfSolute->getMass() > 0.0)
		{
			// Fix the influx concentration if the influx is taking more mass from the surface than there is available
			concInf = surfSolute->getMass() / (infVol + dSurfVol);

			if (concInf * infVol * water->getPoreFrac() > surfSolute->getMass())
			{
				concInf = surfSolute->getMass() / (infVol * water->getPoreFrac());
			}
			// If the surface solute concentration drops below a certain level (kg m-3) set it as zero.
			// If this is not done, adsorption can cause the concentration to start to approach zero but not reach it.
			if (concInf < 0.000000000001)
			{
				concInf = 0.0;
			}
		}
	}
}
*/
