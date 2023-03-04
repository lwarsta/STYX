#include "ModelWater3dTri.h"

void ModelWater3dTri::configure(int iterStopNew, double iterCutThreshNew, double implicityNew, double timeStepNew, std::vector < std::vector<int>> gridMapNew)
{
	iterStop = iterStopNew;
	iterCutThresh = iterCutThreshNew;
	implicity = implicityNew;
	timeStep = timeStepNew;
    gridMap = gridMapNew;
}

void ModelWater3dTri::run(Grid2d &grid2d, Grid3d &grid3d)
{
    preprocess(grid2d, grid3d);
    iterate(grid2d, grid3d);
    postprocess(grid2d, grid3d);
}

void ModelWater3dTri::preprocess(Grid2d &grid2d, Grid3d &grid3d)
{
    std::vector<CellWater3d> * waterCells3d = grid3d.get_water_cells();
    std::vector<CellWater2d> * waterCells2d = grid2d.get_water_cells();

	// Calculate cell properties.
	for (size_t i = 0; i < waterCells3d->size(); i++)
	{
        waterCells3d->at(i).swap();
        waterCells3d->at(i).calcUnsatCond();
        waterCells3d->at(i).calcFluxDrain();

        // Calculate infiltration/exfiltration between 2d and 3d cells.
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
            double fluxInf = waterCells3d->at(i).getCondSat() * geom2d->getArea() * (hydrHeadSurf - waterCells3d->at(i).getHydrHead()) / distance;

            // Check that infiltrating water fits into the available pore space.
            double emptySpace = (waterCells3d->at(i).getWatContSat() - waterCells3d->at(i).getWatCont()) * geom3d->getVolume();
            if (emptySpace < 0.0)
            {
                emptySpace = 0.0;
            }
            if (fluxInf * timeStep > emptySpace)
            {
                fluxInf = emptySpace / timeStep;
            }

            // Check that there is water on the surface.
            if (hydrHeadSurf > waterCells3d->at(i).getHydrHead() && fluxInf * timeStep > waterDepth * geom2d->getArea())
            {
                fluxInf = waterDepth * geom2d->getArea() / timeStep;
            }
            waterCells3d->at(i).setFluxInf(fluxInf);
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

// MOVE THIS TO A MODELTRIBASE CLASS!
void ModelWater3dTri::triDiag(
    std::vector<double>& triA,
    std::vector<double>& triB,
    std::vector<double>& triC,
    std::vector<double>& triD,
    std::vector<double>& triG,
    std::vector<double>& triU,
    std::vector<double>& triV,
    std::vector<double>& triX)
{
    int iEquations = (int)triA.size(); // this is not required

    // Check that all the vectors are the same length?
    triV[0] = triB[0];
    triG[0] = triC[0] / triV[0];
    triU[0] = triD[0] / triV[0];
    for (int i = 1; i < iEquations; i++)
    {
        triV[i] = triB[i] - triA[i] * triG[i - 1];
        triG[i] = triC[i] / triV[i];
        triU[i] = (triD[i] - triA[i] * triU[i - 1]) / triV[i];
    }
    triX[iEquations - 1] = triU[iEquations - 1];
    for (int i = iEquations - 2; i >= 0; i--)
    {
        triX[i] = triU[i] - triG[i] * triX[i + 1];
    }
}

void ModelWater3dTri::iterate(Grid2d &grid2d, Grid3d &grid3d)
{
	std::vector<CellWater3d> * waterCells3d = grid3d.get_water_cells();
    double iterDev;
    int iterCount = 0;
    
    do
    {
        // Solve hydraulic heads in a column of cells.
        iterDev = 0.0;
        //omp_set_dynamic(0);
        //omp_set_num_threads(4);
        #pragma omp parallel for
        for (int m = 0; m < gridMap.size(); m++)
        {
            // Solve hydraulic heads in a column of cells.
            std::vector<double>triA;
            std::vector<double>triB;
            std::vector<double>triC;
            std::vector<double>triD;
            std::vector<double>triG;
            std::vector<double>triU;
            std::vector<double>triV;
            std::vector<double>triX;
            triA.assign(gridMap.at(m).size(), 0.0);
            triB.assign(gridMap.at(m).size(), 0.0);
            triC.assign(gridMap.at(m).size(), 0.0);
            triD.assign(gridMap.at(m).size(), 0.0);
            triG.assign(gridMap.at(m).size(), 0.0);
            triU.assign(gridMap.at(m).size(), 0.0);
            triV.assign(gridMap.at(m).size(), 0.0);
            triX.assign(gridMap.at(m).size(), 0.0);

            for (size_t n = 0; n < gridMap.at(m).size(); n++)
            {
                // Get index from the map.
                int ind = gridMap.at(m).at(n);
                CellGeom3d * geom = waterCells3d->at(ind).getGeom();
                triB.at(n) = waterCells3d->at(ind).getDiffWatCap() * geom->getVolume() / timeStep;
                triD.at(n) = waterCells3d->at(ind).getDiffWatCap() * geom->getVolume() / timeStep * waterCells3d->at(ind).getHydrHeadOld();

                // Vertical flux at the top
                triA.at(n) += waterCells3d->at(ind).getFlux(0);
                triB.at(n) -= waterCells3d->at(ind).getFlux(0);

                // Vertical flux at the bottom
                triC.at(n) += waterCells3d->at(ind).getFlux(5);
                triB.at(n) -= waterCells3d->at(ind).getFlux(5);

                // Horizontal fluxes between cell columns.
                for (size_t k = 1; k < geom->getNumOfNeigh() - 1; k++)
                {
                    CellWater3d * neigh = waterCells3d->at(ind).getNeigh(k);
                    if (neigh != 0)
                    {
                        triD.at(n) -= waterCells3d->at(ind).getFlux(k) * neigh->getHydrHead();
                        triB.at(n) -= waterCells3d->at(ind).getFlux(k);
                    }
                }

                // Include infiltration flux.
                triD.at(n) += waterCells3d->at(ind).getFluxInf();
                waterCells3d->at(ind).setInfVol(waterCells3d->at(ind).getFluxInf() * timeStep);

                // Include fluxes to drains. Pressure is currently set to atmospheric pressure.
                Vertex centerPoint = geom->getCentrePoint();
                double hydrHeadDrain = centerPoint.z + 0.0;
                //if (waterCells3d->at(i).getHydrHead() > hydrHeadDrain)
                if (waterCells3d->at(ind).getHydrHeadOld() > hydrHeadDrain)
                {
                    triD.at(n) -= waterCells3d->at(ind).getFluxDrain() * (waterCells3d->at(ind).getHydrHeadOld() - hydrHeadDrain);
                    //triD.at(i) += waterCells3d->at(i).getFluxDrain() * hydrHeadDrain;
                    //triC.at(i) += waterCells3d->at(i).getFluxDrain();
                    //waterCells3d->at(i).setDrainVol(waterCells3d->at(i).getFluxDrain() * (waterCells3d->at(i).getHydrHead() - hydrHeadDrain) * timeStep);
                    waterCells3d->at(ind).setDrainVol(waterCells3d->at(ind).getFluxDrain() * (waterCells3d->at(ind).getHydrHeadOld() - hydrHeadDrain) * timeStep);
                }
                else
                {
                    waterCells3d->at(ind).setDrainVol(0.0);
                }

                // Include evapotranspiration flux.
                triD.at(n) -= waterCells3d->at(ind).getFluxEvap();
                waterCells3d->at(ind).setEvapVol(waterCells3d->at(ind).getFluxEvap() * timeStep);
            }

            // Solve values in the column.
            triDiag(triA, triB, triC, triD, triG, triU, triV, triX);

            // Compute deviation.
            for (size_t n = 0; n < gridMap.at(m).size(); n++)
            {
                // Save deviation
                int ind = gridMap.at(m).at(n);
                double iterDevNew = fabs(waterCells3d->at(ind).getHydrHead() - triX.at(n));

                // OpenMP double-checked locking.
                if (iterDevNew > iterDev)
                {
                    #pragma omp critical
                    if (iterDevNew > iterDev) iterDev = iterDevNew;
                }
            }

            // Save new hydraulic head values to the cell column.
            for (size_t n = 0; n < gridMap.at(m).size(); n++)
            {
                int ind = gridMap.at(m).at(n);
                waterCells3d->at(ind).setHydrHead(triX.at(n));
                waterCells3d->at(ind).calcPresHead();
                waterCells3d->at(ind).calcWatCont();
                waterCells3d->at(ind).calcDiffWatCap();
            }
        }

        iterCount++;
    } while (iterDev > iterCutThresh && iterCount < iterStop);
}

void ModelWater3dTri::postprocess(Grid2d &grid2d, Grid3d &grid3d)
{
    std::vector<CellWater3d> * waterCells3d = grid3d.get_water_cells();
    std::vector<CellWater2d> * waterCells2d = grid2d.get_water_cells();
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
        }
    }
}
