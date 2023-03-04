#include "CellHeat3d.h"

CellHeat3d::CellHeat3d()
{
    geom = 0;
    water = 0;
	C = 0.0;
	K = 0.0;
	oldT = 0.0;
	T = 0.0;
	U = 0.0;
	fracSoil = 0.0;
	fracWat = 0.0;
	fracIce = 0.0;
	fracWatTot = 0.0;
	fracAir = 0.0;
	densWat = 0.0;
	densIce = 0.0;
	densAir = 0.0;
	capSoil = 0.0;
	capWat = 0.0;
	capIce = 0.0;
	capAir = 0.0;
	condSoil = 0.0;
	condWat = 0.0;
	condIce = 0.0;
	condAir = 0.0;
	condMultSoil = 0.0;
	condMultWat = 0.0;
	condMultIce = 0.0;
	condMultAir = 0.0;
	latHeatFus = 0.0;
	multFreez = 0.0;
	tempFreez = 0.0;
	gradMax = 15.0;
	surfDiffFlux = 0.0;
}

void CellHeat3d::setNeighPointers(std::vector<CellHeat3d*> neighPointersNew)
{
    neighPointers = neighPointersNew;
	fluxDiff.assign(neighPointers.size(), 0.0);
	interK.assign(neighPointers.size(), 0.0);
	gradLimiters.assign(neighPointers.size(), 1.0);
	for (size_t i = 0; i < neighPointers.size(); i++)
	{
		std::vector<double> innerVec;
		innerVec.assign(2, 0.0);
		fluxConv.assign(neighPointers.size(), innerVec);
	}
}

void CellHeat3d::calcCond()
{
	double condSum = condMultSoil * fracSoil +
		             condMultWat * fracWat +
		             condMultIce * fracIce +
		             condMultAir * fracAir;

	if (condSum > 0.0)
	{
		K = (condSoil * condMultSoil * fracSoil + 
			 condWat * condMultWat * fracWat + 
			 condIce * condMultIce * fracIce + 
			 condAir * condMultAir * fracAir) / condSum;
	}
	else
	{
		K = 0.0;
	}
}

void CellHeat3d::calcInterCond()
{
	for (size_t i = 0; i < neighPointers.size(); i++)
	{
		// IS THIS NEEDED: && neighPointers.at(i)->getK() >= 0.0 && K >= 0.0
		if (neighPointers.at(i) != 0 && 
			neighPointers.at(i)->getK() >= 0.0 && 
			K >= 0.0)
		{
			interK.at(i) = sqrt(neighPointers.at(i)->getK() * K);
		}
	}
}

void CellHeat3d::calcDiffFluxes()
{
	if (geom != 0)
	{
		for (size_t i = 0; i < neighPointers.size(); i++)
		{
			if (geom->getDistance(i) > 0.0)
			{
				// A test with a gradient limiter!
				// THIS WAS CHANGED FROM - TO + !?
				fluxDiff.at(i) = interK.at(i) * geom->getFaceArea(i) / 
					geom->getDistance(i); // * gradLimiters.at(i);
			}
		}
	}
}

void CellHeat3d::calcConvFluxes()
{
	if (water != 0 && geom != 0)
	{
		for (size_t i = 0; i < neighPointers.size(); i++)
		{
			double flux = water->getVelInt(i) * geom->getFaceArea(i) * densWat * capWat;
			if (flux < 0.0)
			{
				fluxConv.at(i).at(0) = 0.0;
				fluxConv.at(i).at(1) = flux;
			}
			else
			{
				fluxConv.at(i).at(0) = flux;
				fluxConv.at(i).at(1) = 0.0;
			}
		}
	}
}

void CellHeat3d::calcC()
{
	if (water != 0)
	{
		// CHANGED: water->getPorosity() -> water->getWatContSat()
		fracSoil = 1.0 - water->getWatContSat();
		if (fracSoil < 0.0)
		{
			fracSoil = 0.0;
		}
		fracAir = 1.0 - fracSoil - fracWat - fracIce;
		if (fracAir < 0.0)
		{
			fracAir = 0.0;
		}
		if (T < 0.0 && multFreez > 0.0)
		{
			// CHANGED: geom->getSoilDensity() -> water->getDryWeight()
			C = fracSoil * water->getDryWeight() * capSoil +
				fracWat * densWat * capWat + 
				fracIce * densIce * capIce + 
				fracAir * densAir * capAir + 
				fracWat * densWat * latHeatFus / multFreez;
		}
		else
		{
			// CHANGED: geom->getSoilDensity() -> water->getDryWeight()
			C = fracSoil * water->getDryWeight() * capSoil + 
				fracWat * densWat * capWat + 
				fracIce * densIce * capIce + 
				fracAir * densAir * capAir;
		}
	}
}

void CellHeat3d::calcU()
{
	// Calculate melt and frozen water fractions.
	calcWatAndIceFracs();
	if (geom != 0)
	{
		// CHANGED: geom->getSoilDensity() -> water->getDryWeight()
		U = (fracSoil * water->getDryWeight() * capSoil * T + 
			 fracWat * densWat * capWat * T + 
			 fracIce * densIce * capIce * T + 
			 fracAir * densAir * capAir * T - 
			 fracIce * densIce * latHeatFus) * geom->getVolume();
	}
}

void CellHeat3d::calcDiffEnergy(double timeStep)
{
	// Exchange of energy between neighbour cells in frozen conditions.
	double deltaU = 0.0;
	for (size_t i = 0; i < neighPointers.size(); i++)
	{
		if (neighPointers.at(i) != 0)
		{
			// SHOULD OLD TEMPERATURE VALUE USED HERE?
			deltaU -= -fluxDiff[i] * timeStep * (T - neighPointers.at(i)->getT());
		}
	}
	// Add change of energy to the cell energy.
	U += deltaU;
}

void CellHeat3d::calcConvEnergy(double timeStep)
{
	// Exchange of energy between neighbour cells in frozen conditions.
	double deltaU = 0.0;
	for (size_t i = 0; i < neighPointers.size(); i++)
	{
		// One of the components is zero.
		if (neighPointers.at(i) != 0)
		{
			// SHOULD OLD TEMPERATURE VALUE USED HERE?
			deltaU += fluxConv.at(i).at(0) * timeStep * neighPointers.at(i)->getT();
			deltaU += fluxConv.at(i).at(1) * timeStep * T;
		}
	}
	// Add change of energy to the cell energy.
	U += deltaU;
}

void CellHeat3d::calcWatAndIceFracs()
{
	// Difference between water and ice is neglected
	if (T < tempFreez && multFreez > 0.0)
	{
		fracIce = fracWatTot * (1.0 - exp(-(tempFreez - T) / multFreez));
	}
	else
	{
		fracIce = 0.0;
	}
	fracWat = fracWatTot - fracIce;
}

double CellHeat3d::calcResidual(double newT)
{
	T = newT;
	calcWatAndIceFracs();
	if (water != 0 && geom != 0)
	{
		// CHANGED: geom->getSoilDensity() -> water->getDryWeight()
		return (fracSoil * water->getDryWeight() * capSoil * newT + 
			    fracWat * densWat * capWat * newT + 
			    fracIce * densIce * capIce * newT + 
			    fracAir * densAir * capAir * newT - 
			    fracIce * densIce * latHeatFus) * geom->getVolume() - U;
	}
}

void CellHeat3d::compGradLimiters()
{
	for (size_t i = 0; i < neighPointers.size(); i++)
	{
		if (neighPointers.at(i) != 0 && geom != 0 && geom->getDistance(i) > 0.0)
		{
			double grad = fabs((getT() - neighPointers.at(i)->getT()) / 
				geom->getDistance(i));
			if (grad > gradMax && grad > 0.0)
			{
				gradLimiters.at(i) = gradMax / grad;
			}
			else
			{
				gradLimiters.at(i) = 1.0;
			}
		}
		else
		{
			gradLimiters.at(i) = 1.0;
		}
	}
}

void CellHeat3d::calcSurfDiffFlux()
{
	if (geom != 0)
	{
		/*
		// UNDER DEVELOPMENT.
		double depthSnow = 0.5 * atmos->getDepthSnow();
		if (depthSnow < 0.0)
		{
			depthSnow = 0.0;
		}
		double condSnow = atmos->getCondSnow();
		double depthSoil = 0.5 * geom->getDim(2);
		double condAver = 0.0;
		// Compute heat flux between snow cover/air and top subsurface soil layer. 
		surfDiffFlux = 0.0;

		if (depthSoil + depthSnow > 0.0)
		{
			surfDiffFlux = (condSnow * depthSnow / (depthSoil + depthSnow) + K * 
				depthSoil / (depthSoil + depthSnow)) * geom->getFaceArea(0) * 
				(atmos->getTemp() - T) / (depthSoil + depthSnow);
		}
		*/
	}
}
