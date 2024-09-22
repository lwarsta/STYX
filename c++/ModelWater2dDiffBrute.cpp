#include "ModelWater2dDiffBrute.h"

void ModelWater2dDiffBrute::configure(int iterStopNew, double iterCutThreshNew, double implicityNew, double timeStepNew, 
                                      double iterThreshBisNew, int iterStopBisNew, double threshLeftBisNew, double threshRightBisNew)
{
	iterStop = iterStopNew;
	iterCutThresh = iterCutThreshNew;
	implicity = implicityNew;
	timeStep = timeStepNew;
    iterThreshBis = iterThreshBisNew;
    iterStopBis = iterStopBisNew;
    threshLeftBis = threshLeftBisNew;
    threshRightBis = threshRightBisNew;
}

void ModelWater2dDiffBrute::run(Grid2d &grid2d, Grid3d &grid3d)
{
	preprocess(grid2d, grid3d);
    iterate(grid2d, grid3d);
    postprocess(grid2d, grid3d);
}

void ModelWater2dDiffBrute::preprocess(Grid2d &grid2d, Grid3d &grid3d)
{
	std::vector<CellWater2d> * waterCells2d = grid2d.get_water_cells();
    // Calculate cell properties.
	//#pragma omp parallel for
	for (size_t i = 0; i < waterCells2d->size(); i++)
	{
        if (waterCells2d->at(i).getWaterDepth() < 0.0)
        {
            waterCells2d->at(i).setWaterDepth(0.0);
        }

		waterCells2d->at(i).swap();
	}
}

double ModelWater2dDiffBrute::calcResidual(CellWater2d& water, double waterDepth)
{	
	// Calculate water volume in the cell.
	CellGeom2d* geom = water.getGeom();
	double waterVolume = (waterDepth - water.getWaterDepthOld()) * geom->getArea();
	int material = geom->getMaterial();

	// Bypass building and inactive cells.
	//if (material != 2 && material != 0)
	if (material != 0)
	{
		// Initialise local inter cell variables.
		std::vector<double> fluxInter;
		std::vector<double> velocityInter;
		fluxInter.assign(geom->getNumOfNeigh(), 0.0);
		velocityInter.assign(geom->getNumOfNeigh(), 0.0);

		// Decrease water depth by flow threshold depth.
		waterDepth -= water.getDeprStor();

		// Water depth cannot decrease below zero.
		if (waterDepth < 0.0)
		{
			waterDepth = 0.0;
		}

		// Hydraulic head in the cell.
		Vertex centerPoint = geom->getCentrePoint();
		double elevation = centerPoint.z;
		double head = elevation + waterDepth;

		// Compute overland flow between cells.
		double dischIn = 0.0;
		double dischOut = 0.0;
		for (int i = 0; i < geom->getNumOfNeigh(); i++)
		{
			CellWater2d * waterNeigh = water.getNeigh(i);

			// Check for a border cell.
			if (waterNeigh != 0)
			{
				// Compute neighbour cell elevation.
				CellGeom2d * geomNeigh = waterNeigh->getGeom();
				int materialNeigh = geomNeigh->getMaterial();

				// Bypass building.
				//if (materialNeigh != 2 && materialNeigh != 0)
				if (((material == 2 && materialNeigh == 2) or (material != 2 && materialNeigh != 2)) && materialNeigh != 0)
				//if (materialNeigh != 0)
				{
					Vertex centerPointNeigh = geomNeigh->getCentrePoint();
					double elevationNeigh = centerPointNeigh.z;

					// Flow threshold might be different in the neighbour cell?
					double waterDepthNeigh = waterNeigh->getWaterDepth() - waterNeigh->getDeprStor();
					if (waterDepthNeigh < 0.0)
					{
						waterDepthNeigh = 0.0;
					}
					double headNeigh = elevationNeigh + waterDepthNeigh;

					// Calculate gradient. Should the distance exclude elevation?
					double slope = 0.0;
					if (geom->getDistance(i) > 0.0)
					{
						slope = (headNeigh - head) / geom->getDistance(i);
					}

					// Calculate effective flow depth.
					double waterDepthEff = 0.0;
					if (head > headNeigh)
					{
						if (elevation - elevationNeigh >= waterDepthNeigh)
						{
							waterDepthEff = waterDepth;
						}
						else
						{
							waterDepthEff = head - headNeigh;
						}
					}
					else if (head < headNeigh) // equality is discarded
					{
						if (elevationNeigh - elevation >= waterDepth)
						{
							waterDepthEff = waterDepthNeigh;
						}
						else
						{
							waterDepthEff = headNeigh - head;
						}
					}

					// Calculate Mannings n between cells (with arithmetic mean).
					double mannN = 0.5 * (water.getMannN() + waterNeigh->getMannN());
					// Calculate velocity.
					double velocity = 0.0;
					if (waterDepthEff > 0.0 && mannN > 0.0)
					{
						// Cube root hack could be tested here. 
						if (slope < 0.0)
						{
							velocity = 1.0 / mannN * sqrt(-slope) * pow(waterDepthEff, 2.0 / 3.0);
						}
						else if (slope > 0.0) // equality is discarded
						{
							velocity = -1.0 / mannN * sqrt(slope) * pow(waterDepthEff, 2.0 / 3.0);
						}
					}

					// Save velocity and flux and update water volume.
					velocityInter.at(i) = -velocity;
					double flux = waterDepthEff * geom->getSideLength(i) * velocity;
					fluxInter.at(i) = -flux;
					waterVolume += flux * timeStep;
					if (flux < 0.0)
					{
						dischOut -= flux;
					}
					else
					{
						dischIn += flux;
					}
				}
			}
		}
		// Save average discharge in the cell.
		water.setAvgDisch(0.5 * (dischIn + dischOut));
	}

	return waterVolume;
}

double ModelWater2dDiffBrute::bisection(CellWater2d& water, double left, double right)
{
	int iterCountBis = 0;

	while (fabs(right - left) > 2.0 * iterThreshBis && iterCountBis < iterStopBis) // iterCutThresh -> iterThreshBis
	{
		// Calculate midpoint of domain.
		double midpoint = (right + left) * 0.5;

		// Find f(midpoint).
		if (calcResidual(water, left) * calcResidual(water, midpoint) > 0.0)
		{
			// Throw away left half.
			left = midpoint;
		}
		else
		{
			// Throw away right half.
			right = midpoint;
		}

		iterCountBis++;
	}
	return (right + left) * 0.5;
}

void ModelWater2dDiffBrute::iterate(Grid2d &grid2d, Grid3d &grid3d)
{
	// The main iteration loop.
    std::vector<CellWater2d>* waterCells2d = grid2d.get_water_cells();
	double iterDev;
	int iterCount = 0;
	do
    {
		iterDev = 0.0;
		//omp_set_dynamic(0);
		//omp_set_num_threads(4);
		//#pragma omp parallel for reduction(max:iterDev)
		#pragma omp parallel for
		for (int i = 0; i < waterCells2d->size(); i++)
        {
			// Bypass inactive cells.
			CellGeom2d * geom2d = waterCells2d->at(i).getGeom();
			if (geom2d->getMaterial() == 0)
			{
				continue;
			}
			// Compute new water depth.
            double waterDepthNew = bisection(waterCells2d->at(i), threshLeftBis, threshRightBis);
            // Save iteration deviation.
			double iterDevNew = fabs(waterCells2d->at(i).getWaterDepth() - waterDepthNew);
			// OpenMP double-checked locking.
			if (iterDevNew > iterDev)
            {
				#pragma omp critical
				if (iterDevNew > iterDev) iterDev = iterDevNew;
            }
            // Save the new water depth.
            waterCells2d->at(i).setWaterDepth(waterDepthNew);
		}
		iterCount++;
	} while (iterDev > iterCutThresh && iterCount < iterStop);
}

void ModelWater2dDiffBrute::postprocess(Grid2d &grid2d, Grid3d &grid3d)
{
	//std::vector<CellWater2d> * waterCells2d = grid2d.getWaterCells();
    //for (size_t i = 0; i < waterCells2d->size(); i++)
    //{
    //
    //}
}
