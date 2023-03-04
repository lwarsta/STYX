#include "ModelWater3dBrute.h"

void ModelWater3dBrute::configure(int iterStopNew, double iterCutThreshNew, double implicityNew, double timeStepNew)
{
	iterStop = iterStopNew;
	iterCutThresh = iterCutThreshNew;
	implicity = implicityNew;
	timeStep = timeStepNew;
}

void ModelWater3dBrute::run(Grid2d &grid2d, Grid3d &grid3d)
{
    preprocess(grid2d, grid3d);
    iterate(grid2d, grid3d);
    postprocess(grid2d, grid3d);
}

void ModelWater3dBrute::preprocess(Grid2d &grid2d, Grid3d &grid3d)
{
    std::vector<CellWater3d> * waterCells3d = grid3d.get_water_cells();
    std::vector<CellWater2d> * waterCells2d = grid2d.get_water_cells();

	// Calculate cell properties.
    //#pragma omp parallel for
	for (size_t i = 0; i < waterCells3d->size(); i++)
	{
        waterCells3d->at(i).swap();
        waterCells3d->at(i).calcUnsatCond();
        waterCells3d->at(i).calcFluxDrain();

        // Calculate infiltration / exfiltration between 2d and 3d cells.
        CellGeom3d * geom3d = waterCells3d->at(i).getGeom();
        int cell2dId = geom3d->getGridConnection();
        if (cell2dId >= 0 && cell2dId < waterCells2d->size())
        {
            Vertex centre3d = geom3d->getCentrePoint();
            CellGeom2d * geom2d =  waterCells2d->at(cell2dId).getGeom();
            Vertex centre2d = geom2d->getCentrePoint();
            double waterDepth = waterCells2d->at(cell2dId).getWaterDepth();
            double hydrHeadSurf = waterDepth + centre2d.z;
            Vertex vec = geom3d->createVector(centre2d, centre3d);
            double distance = geom3d->computeVectorLength(vec);

            // Infiltration.
            double fluxInfExf = 0.0;
            if (waterCells3d->at(i).getPresHead() < 0.0 && hydrHeadSurf > waterCells3d->at(i).getHydrHead())
            {
                // Compute flux from surface domain to subsurface domain.
                fluxInfExf = waterCells3d->at(i).getCondSat() * geom2d->getArea() * (hydrHeadSurf - waterCells3d->at(i).getHydrHead()) / distance;

                // Check that infiltrating water fits into the available pore space.
                double emptySpace = (waterCells3d->at(i).getWatContSat() - waterCells3d->at(i).getWatCont()) * geom3d->getVolume();
                if (emptySpace < 0.0)
                {
                    emptySpace = 0.0;
                }
                if (fluxInfExf * timeStep > emptySpace)
                {
                    fluxInfExf = emptySpace / timeStep;
                }

                // Check that there is enough water on the surface.
                if (hydrHeadSurf > waterCells3d->at(i).getHydrHead() && fluxInfExf * timeStep > waterDepth * geom2d->getArea())
                {
                    fluxInfExf = waterDepth * geom2d->getArea() / timeStep;
                }
            }
            // Exfiltration. THE DRY END WATER CONTENT COULD BE PRECOMPUTED IN INITIALIZATION PHASE?
            else if (waterCells3d->at(i).getPresHead() > 0.0 && hydrHeadSurf < waterCells3d->at(i).getHydrHead())
            {
                // Compute flux from subsurface domain to surface domain.
                fluxInfExf = waterCells3d->at(i).getCondUnsat() * geom2d->getArea() * (hydrHeadSurf - waterCells3d->at(i).getHydrHead()) / distance;

                // Cap the exfiltration flux to the available water volume.
                double watCont = waterCells3d->at(i).getWatCont();

                // Set the freely available water to empty at a suction of 1 m.
                double hydrHeadNew = centre3d.z - 1.0;
                waterCells3d->at(i).setHydrHead(hydrHeadNew);
                waterCells3d->at(i).calcPresHead();
                waterCells3d->at(i).calcWatCont();
                double watContDry = waterCells3d->at(i).getWatCont();
                double soilWatVol = (watCont - watContDry) * geom3d->getVolume();
                if (soilWatVol < 0.0)
                {
                    soilWatVol = 0.0;
                }
                if (-fluxInfExf * timeStep > soilWatVol)
                {
                    fluxInfExf = -soilWatVol / timeStep;
                }

                // Restore old hydraulic and pressure head values and recompute water content values.
                waterCells3d->at(i).setHydrHead(waterCells3d->at(i).getHydrHeadOld());
                waterCells3d->at(i).calcPresHead();
                waterCells3d->at(i).calcWatCont(); // Probably not needed
            }
            waterCells3d->at(i).setFluxInf(fluxInfExf);
        }
	}

	// Calculate inter cell properties.
	for (size_t i = 0; i < waterCells3d->size(); i++)
	{
        waterCells3d->at(i).calcCondInter();
        waterCells3d->at(i).calcFluxes();
        waterCells3d->at(i).calcFlowVel(); // is this needed here?
	}
}

void ModelWater3dBrute::iterate(Grid2d &grid2d, Grid3d &grid3d)
{
	// The main iteration loop.
	std::vector<CellWater3d> * waterCells3d = grid3d.get_water_cells();
	double iterDev;
	int iterCount = 0;
	do
    {
		iterDev = 0.0;
		//omp_set_dynamic(0);
		//omp_set_num_threads(4);
		//#pragma omp parallel for reduction(max:iterDev)
        #pragma omp parallel for
		for (int i = 0; i < waterCells3d->size(); i++)
        {
            // Bypass inactive cells.
            CellGeom3d* geom = waterCells3d->at(i).getGeom();
            if (geom->getMaterial() == 0)
            {
                continue;
            }

            // Compute new water depth.
            double numerator = waterCells3d->at(i).getDiffWatCap() * geom->getVolume() / timeStep * waterCells3d->at(i).getHydrHeadOld();
            double denominator = waterCells3d->at(i).getDiffWatCap() * geom->getVolume() / timeStep;

            // Include fluxes between neighbor cells into the solution.
            for (size_t j = 0; j < geom->getNumOfNeigh(); j++)
            {
                CellWater3d * neigh = waterCells3d->at(i).getNeigh(j);
                if (neigh != 0)
                {
                    numerator -= waterCells3d->at(i).getFlux(j) * neigh->getHydrHead();
                    denominator -= waterCells3d->at(i).getFlux(j);
                }
            }

            // Include fluxes to drains. Pressure is currently set to atmospheric pressure.
            Vertex centerPoint = geom->getCentrePoint();
            double hydrHeadDrain = centerPoint.z + 0.0;
            //if (waterCells3d->at(i).getHydrHead() > hydrHeadDrain)
            if (waterCells3d->at(i).getHydrHeadOld() > hydrHeadDrain)
            {
                numerator -= waterCells3d->at(i).getFluxDrain() * (waterCells3d->at(i).getHydrHeadOld() - hydrHeadDrain);
                //numerator += waterCells3d->at(i).getFluxDrain() * hydrHeadDrain;
                //denominator += waterCells3d->at(i).getFluxDrain();
                //waterCells3d->at(i).setDrainVol(waterCells3d->at(i).getFluxDrain() * (waterCells3d->at(i).getHydrHead() - hydrHeadDrain) * timeStep);
                waterCells3d->at(i).setDrainVol(waterCells3d->at(i).getFluxDrain() * (waterCells3d->at(i).getHydrHeadOld() - hydrHeadDrain) * timeStep);
            }
            else
            {
                waterCells3d->at(i).setDrainVol(0.0);
            }

            // Include infiltration/exfiltration flux.
            numerator += waterCells3d->at(i).getFluxInf();
            waterCells3d->at(i).setInfVol(waterCells3d->at(i).getFluxInf() * timeStep);

            // Include evapotranspiration flux.
            numerator -= waterCells3d->at(i).getFluxEvap();
            waterCells3d->at(i).setEvapVol(waterCells3d->at(i).getFluxEvap() * timeStep);

            // Compute the new head. Check that denominator is not zero!
            double hydrHeadNew = numerator / denominator;

            // Save iteration deviation.
            double iterDevNew = fabs(waterCells3d->at(i).getHydrHead() - hydrHeadNew );

            // OpenMP double-checked locking.
            if (iterDevNew > iterDev)
            {
                #pragma omp critical
                if (iterDevNew > iterDev) iterDev = iterDevNew;
            }

            // Set the new head.
            waterCells3d->at(i).setHydrHead((1.0 - implicity) * waterCells3d->at(i).getHydrHead() + implicity * hydrHeadNew);

            // Update differential water capacity.
            // Previously computed in a separate loop outside this loop.
            waterCells3d->at(i).calcPresHead();
            waterCells3d->at(i).calcWatCont();
            waterCells3d->at(i).calcDiffWatCap();
		}
		iterCount++;
	} while (iterDev > iterCutThresh && iterCount < iterStop);
}

void ModelWater3dBrute::postprocess(Grid2d &grid2d, Grid3d &grid3d)
{
    std::vector<CellWater3d> * waterCells3d = grid3d.get_water_cells();
    std::vector<CellWater2d> * waterCells2d = grid2d.get_water_cells();
    //#pragma omp parallel for
    for (size_t i = 0; i < waterCells3d->size(); i++)
    {
        // Remove infiltrating water from the surface cell.
        CellGeom3d * geom3d = waterCells3d->at(i).getGeom();
        size_t cell2dId = geom3d->getGridConnection();
        if (cell2dId >= 0 && cell2dId < waterCells2d->size())
        {
            double waterDepth = waterCells2d->at(cell2dId).getWaterDepth();
            double infVol = waterCells3d->at(i).getInfVol();
            CellGeom2d * geom2d =  waterCells2d->at(cell2dId).getGeom();
            waterCells2d->at(cell2dId).setWaterDepth(waterDepth - infVol / geom2d->getArea());
            // TEMPORARY.
            if (waterCells2d->at(cell2dId).getWaterDepth() < 0.0)
            {
                waterCells2d->at(cell2dId).setWaterDepth(0.0);
            }
        }
    }
}
