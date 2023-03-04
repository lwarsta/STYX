#include "CellWater3d.h"

CellWater3d::CellWater3d()
{
    geom = 0;
    hydrHead = 0.0;
    hydrHeadOld = 0.0;
    presHead = 0.0;
	watContRes = 0.0;
	watContSat = 0.0;
	watCont = 0.0;
	watContOld = 0.0;
	diffWatCap = 0.0;
	alpha = 0.0;
	m = 0.0;
	n = 0.0;
	dryWeight = 0.0;
	compress = 0.0;
	condSat = 0.0;
	condUnsat = 0.0;
	fluxDrain = 0.0;
	drainFactor = 0.0;
	dischDrain = 0.0;
	fluxInf = 0.0;
	infVolume = 0.0;
	petCellFact = 0.0;
	fluxEvap = 0.0;
	evapVolume = 0.0;
	fracIce = 0.0;
}

void CellWater3d::calcPresHead()
{
    if (geom != 0)
	{
        Vertex centerPoint = geom->getCentrePoint();
        presHead = hydrHead - centerPoint.z;
    }
}

void CellWater3d::assignGeom(CellGeom3d *geomNew)
{
    geom = geomNew;

    // Resize vectors depicting intercell properties.
    if (geom != 0)
	{
        fluxes.assign(geom->getNumOfNeigh(), 0.0);
        condInter.assign(geom->getNumOfNeigh(), 0.0);
        velInter.assign(geom->getNumOfNeigh(), 0.0);
    }
}

void CellWater3d::setN(double nNew)
{
	if (nNew > 1.0)
	{
		n = nNew;
		m = 1.0 - 1.0 / n;
	}
	else {
		std::cout << "\n-> Error, n value must be higher than 1.0 (CellWater3d).";
	}
}

void CellWater3d::calcWatCont()
{
	if (presHead >= 0.0)
	{
		watCont = watContSat + presHead * compress;
	}
	else
	{
		double tmp = 1.0 + pow(fabs(alpha * presHead), n);
		watCont = watContRes + (watContSat - watContRes) / pow(tmp, m);
	}
}

void CellWater3d::calcUnsatCond()
{
	// Calculate relative saturation ( 0 <= Sr <= 1 )
	double saturation = 1.0;

	if (watContSat - watContRes > 0.0)
	{
		saturation = (watCont - watContRes) / (watContSat - watContRes);
	}

	// Fix saturation in case of under or over saturation.
	if (saturation > 1.0 || saturation < 0.0)
	{
		saturation = 1.0; // Check the latter conditional?
	}

	// Incomplete beta function?
	double tmp = 1.0;

	if (m > 0.0)
	{
		tmp = 1.0 - pow(saturation, 1.0 / m);
	}

	if (tmp > 0)
	{
		tmp = 1.0 - pow(tmp, m);
	}
	else
	{
		tmp = 1.0;
	}

	// Calculate relative conductivity - note exponent l = 0.5 (pore connectivity)
	// Calculate unsaturated conductivity.
	condUnsat = condSat * sqrt(saturation) * tmp * tmp;
}

void CellWater3d::calcCondInter()
{
    for (size_t i = 0; i < condInter.size(); i++)
	{
        // Compute using geometric mean.
        if (neighPointers.at(i) != 0)
		{
            condInter.at(i) = sqrt(neighPointers.at(i)->getCondUnsat() * condUnsat);
        }
        //// Compute using arithmetic mean.
        //if (neighPointers.at(i) != 0) {
        //    condInter.at(i) = 0.5 * (neighPointers.at(i)->getCondUnsat() + 
		//		condUnsat);
        //}
    }
}

void CellWater3d::calcDiffWatCap()
{
	double deltaHydrHead = hydrHead - hydrHeadOld;
	double deltaWatCont = watCont - watContOld;

	if (deltaHydrHead != 0.0 && deltaWatCont != 0.0)
	{
		diffWatCap = deltaWatCont / deltaHydrHead;
	}
	else
	{
		// Use compressibility in steady-state conditions.
		diffWatCap = compress;
	}
}

void CellWater3d::calcFluxes()
{
	for (size_t i = 0; i < fluxes.size(); i++)
	{
		if (geom != 0 && neighPointers.at(i) != 0)
		{ // is this the best way to check this?
			if (geom->getDistance(i) > 0.0)
			{
				// The signs could be reversed to make the equations more clear?
				fluxes.at(i) = -condInter.at(i) * geom->getFaceArea(i) / 
					geom->getDistance(i);
			}
		}
	}
}

void CellWater3d::calcFlowVel()
{
	for (size_t i = 0; i < fluxes.size(); i++)
	{
		velInter.at(i) = 0.0;

		if (geom != 0 && neighPointers.at(i) != 0)
		{
			// Should geometric mean used instead of arithmetic mean?
			double watContInter = 0.5 * (watCont +  neighPointers.at(i)->getWatCont());

			if (watContInter > 0.0 && geom != 0 && geom->getDistance(i) > 0.0)
			{
				velInter.at(i) = condInter.at(i) / watContInter * 
					(neighPointers.at(i)->getHydrHead() - hydrHead) / (geom->getDistance(i));
			}
		}
	}
}

double CellWater3d::getRelSat()
{
	if (watContSat - watContRes > 0.0)
	{
		double watContRel = (watCont - watContRes) / (watContSat - watContRes);

		if (watContRel > 1.0)
		{
			watContRel = 1.0;
		}
		else if (watContRel < 0.0)
		{
			watContRel = 0.0;
		}
		return watContRel;
	}
	else
	{
		return 0.0;
	}
}

void CellWater3d::calcFluxDrain()
{
    //std::cout << "\n" << id;
    //std::cout << "\n" << drainFactor;
    //std::cout << "\n" << condUnsat;
    //std::cout << "\n";
    fluxDrain = drainFactor * condUnsat;
}

void CellWater3d::calcFluxEvap(double pet)
{
	if (pet > 0.0 && petCellFact > 0.0)
	{
        double evapHeadMax = -0.1;
        double evapHeadMin = -5.0;
        double evapHeadWilt = -150.0;
		double petMult = 1.0;
		double presEff = presHead;

		// fix pressure
		if (presEff > 0.0)
		{
			presEff = 0.0;
		}
		if (presEff < evapHeadWilt)
		{
			presEff = evapHeadWilt;
		}

		// Optimal case.
		if (presEff >= evapHeadMin && presEff <= evapHeadMax )
		{
			petMult = 1.0;
		}
		// Soil is too wet - wet end seems to hinder et too much.
		else if (presEff > evapHeadMax)
		{
			petMult = 1.0;
			//petMult = 1.0 - (presEff - evapHeadMax) / (0.0 - evapHeadMax);
		}
		// Soil is too dry.
		else
		{
			petMult = 1.0 - (presEff - evapHeadMin) / (evapHeadWilt - evapHeadMin);
		}

		fluxEvap = petMult * pet * petCellFact;
	}
	else
	{
		fluxEvap = 0.0;
	}
}
