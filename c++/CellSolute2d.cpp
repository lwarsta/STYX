#include "CellSolute2d.h"

CellSolute2d::CellSolute2d()
{
    geom = 0;
    water = 0;
	conc = 0.0;
	concOld = 0.0;
	mass = 0.0;
	massOld = 0.0;
	decayRate = 0.0;
	decayTarget = -1;
	adsParam[0] = adsParam[1] = adsParam[2] = 0.0;
	depos = 0.0;
	deposWetRate = 0.0;
	deposDryRate = 0.0;
}

void CellSolute2d::convertMassToConc()
{
	if (water != 0 && geom != 0 && water->getWaterDepth() > 0.0)
	{
		conc = mass / (water->getWaterDepth() * geom->getArea());
		mass = 0.0;
	}
	else
	{
		conc = 0.0;
	}
}

void CellSolute2d::convertConcToMass()
{
	if (water != 0 && geom != 0 && water->getWaterDepth() > 0.0)
	{
		mass = mass + conc * water->getWaterDepth() * geom->getArea();
		massOld = concOld * water->getWaterDepth() * geom->getArea(); // Should new or old water depth be used?
	}
	else
	{
		mass = 0.0;
		massOld = 0.0;
	}
}
