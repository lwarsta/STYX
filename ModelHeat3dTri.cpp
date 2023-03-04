#include "ModelHeat3dTri.h"

void ModelHeat3dTri::configure(int iterStopNew, double iterCutThreshNew, 
    double implicityNew, double timeStepNew, double boundTempNew,
    double iterThreshBisNew, int iterStopBisNew, double threshLeftBisNew, 
    double threshRightBisNew, std::vector < std::vector<int>> gridMapNew)
{
    iterStop = iterStopNew;
    iterCutThresh = iterCutThreshNew;
    implicity = implicityNew;
    timeStep = timeStepNew;
    boundTemp = boundTempNew;
    iterThreshBis = iterThreshBisNew;
    iterStopBis = iterStopBisNew;
    threshLeftBis = threshLeftBisNew;
    threshRightBis = threshRightBisNew;
    gridMap = gridMapNew;
}

void ModelHeat3dTri::run(Grid2d &grid2d, Grid3d &grid3d)
{
    preprocess(grid2d, grid3d);
    iterate(grid2d, grid3d);
    postprocess(grid2d, grid3d);
}

void ModelHeat3dTri::preprocess(Grid2d &grid2d, Grid3d &grid3d)
{
    std::vector<CellHeat2d>* heatCells2d = grid2d.get_heat_cells();
    std::vector<CellHeat3d>* cells_heat_3d = grid3d.get_heat_cells();

    // Calculate cell properties.
    for (size_t i = 0; i < cells_heat_3d->size(); i++)
    {
        cells_heat_3d->at(i).calcFracWatTot();
        cells_heat_3d->at(i).swap();
        cells_heat_3d->at(i).calcC();
        cells_heat_3d->at(i).calcCond();
        cells_heat_3d->at(i).calcSurfDiffFlux();
    }

    // Calculate inter cell properties.
    for (int i = 0; i < cells_heat_3d->size(); i++)
    {
        cells_heat_3d->at(i).calcInterCond();
        cells_heat_3d->at(i).calcDiffFluxes();
    }
}

void ModelHeat3dTri::triDiag(
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

void ModelHeat3dTri::iterate(Grid2d& grid2d, Grid3d& grid3d)
{
    std::vector<CellHeat3d>* cells_heat_3d = grid3d.get_heat_cells();

    // Compute explicit solution.
    std::vector<double> tempExp;
    tempExp.assign(cells_heat_3d->size(), 0.0);

    for (size_t i = 0; i < cells_heat_3d->size(); i++)
    {
        CellGeom3d* geom = cells_heat_3d->at(i).getGeom();

        // Define numerator and denominator.
        double numerator = cells_heat_3d->at(i).getC() * geom->getVolume() /
            timeStep * cells_heat_3d->at(i).getOldT();
        double denominator = cells_heat_3d->at(i).getC() * geom->getVolume() /
            timeStep;

        // Include fluxes between neighbor cells into the solution.
        for (size_t k = 0; k < geom->getNumOfNeigh(); k++)
        {
            CellHeat3d* neigh = cells_heat_3d->at(i).getNeigh(k);

            if (neigh != 0)
            {
                // Compute dispersive fluxes between cells.
                numerator += cells_heat_3d->at(i).getDiffFlux(k) * 
                    (neigh->getOldT() - cells_heat_3d->at(i).getOldT());

                // Compute advective fluxes between cells.
                numerator += cells_heat_3d->at(i).getConvFlux(k, 0) * 
                    neigh->getOldT();
                numerator += cells_heat_3d->at(i).getConvFlux(k, 1) * 
                    cells_heat_3d->at(i).getOldT();
            }
        }

        // Heat flux between soil surface and top soil layer.
        numerator += cells_heat_3d->at(i).getSurfDiffFlux();

        // Top boundary condition - TEMPORARY
        if (i == 0)
        {
            denominator = 1.0;
            numerator = boundTemp; // boundary condition here (e.g. 10.0).
        }

        // Bottom boundary condition - TEMPORARY
        //if (i == cells_heat_3d->size() - 1)
        //{
        //    denominator = 1.0;
        //    numerator = 5.0; // boundary condition here.
        //}
        
        // Compute the new concentration.
        double tempExpNew = 0.0;
        if (denominator != 0.0)
        {
            tempExpNew = numerator / denominator;
        }
        if (tempExpNew < 0.0)
        {
            tempExpNew = 0.0;
        }

        // Set the new concentration.
        tempExp.at(i) = tempExpNew;

        if (implicity <= 0.0)
        {
            cells_heat_3d->at(i).setT(tempExpNew);
        }
    }

    // Compute implicit solution.
    double iterDev;
    int iterCount = 0;

    do
    {
        iterDev = 0.0;
        //omp_set_dynamic(0);
        //omp_set_num_threads(4);
        #pragma omp parallel for
        for (int m = 0; m < gridMap.size(); m++)
        {
            // Solve heat transport in a column of cells.
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
                CellGeom3d* geom = cells_heat_3d->at(ind).getGeom();
                triB.at(n) = cells_heat_3d->at(ind).getC() * 
                    geom->getVolume() / timeStep;
                triD.at(n) = cells_heat_3d->at(ind).getC() * 
                    geom->getVolume() / timeStep * 
                    cells_heat_3d->at(ind).getOldT();

                // Vertical fluxes at the top
                triA.at(n) += cells_heat_3d->at(ind).getDiffFlux(0);
                triB.at(n) -= cells_heat_3d->at(ind).getDiffFlux(0);
                triA.at(n) -= cells_heat_3d->at(ind).getConvFlux(0, 0);
                triB.at(n) -= cells_heat_3d->at(ind).getConvFlux(0, 1);

                // Vertical fluxes at the bottom
                triC.at(n) += cells_heat_3d->at(ind).getDiffFlux(5);
                triB.at(n) -= cells_heat_3d->at(ind).getDiffFlux(5);
                triC.at(n) -= cells_heat_3d->at(ind).getConvFlux(5, 0);
                triB.at(n) -= cells_heat_3d->at(ind).getConvFlux(5, 1);

                // Horizontal fluxes between cell columns
                for (size_t j = 1; j < geom->getNumOfNeigh() - 1; j++)
                {
                    CellHeat3d* neigh = cells_heat_3d->at(ind).getNeigh(j);

                    if (neigh != 0)
                    {
                        triD.at(n) += cells_heat_3d->at(ind).getDiffFlux(j) *
                            neigh->getT();
                        triB.at(n) += cells_heat_3d->at(ind).getDiffFlux(j);
                        triD.at(n) += cells_heat_3d->at(ind).getConvFlux(j, 0) * 
                            neigh->getT();
                        triB.at(n) -= cells_heat_3d->at(ind).getConvFlux(j, 1);
                    }
                }

                // Heat flux between soil surface and top soil layer.
                triD.at(n) += cells_heat_3d->at(ind).getSurfDiffFlux();

                // Top boundary condition.
                //if (i == 0)
                //{
                //    dTriA.(i) = dTriC.(i) = 0.0;
                //    dTriB.(i) = 1.0;
                //    dTriD.(i) = boundTemp;
                //}
            }

            // Solve values in the column.
            triDiag(triA, triB, triC, triD, triG, triU, triV, triX);

            // Compute deviation.
            for (size_t n = 0; n < gridMap.at(m).size(); n++)
            {
                // Save deviation
                int ind = gridMap.at(m).at(n);
                double temp_new = (1.0 - implicity) * tempExp.at(ind) + 
                    implicity * triX.at(n);
                double iterDevNew = fabs(cells_heat_3d->at(ind).getT() - 
                    temp_new);

                // OpenMP double-checked locking.
                if (iterDevNew > iterDev)
                {
                    #pragma omp critical
                    if (iterDevNew > iterDev) iterDev = iterDevNew;
                }

                // Save new temperature value.
                cells_heat_3d->at(ind).setT(temp_new);
            }
        }

        iterCount++;
    } while (iterDev > iterCutThresh && iterCount < iterStop);
}

void ModelHeat3dTri::postprocess(Grid2d& grid2d, Grid3d& grid3d)
{
    std::vector<CellHeat3d>* cells_heat_3d = grid3d.get_heat_cells();
    for (size_t i = 0; i < cells_heat_3d->size(); i++)
    {
        // Calc. U [kJ].
        // test this without the threshold (1.0) with a short time step
        double tempThreshold = 1000.0;

        if (cells_heat_3d->at(i).getT() < tempThreshold)
        {
            cells_heat_3d->at(i).calcU();
        }

        // Calc. diff. energy fluxes.
        if (cells_heat_3d->at(i).getT() < tempThreshold && i != 0 && 
            i != cells_heat_3d->size() - 1)
        {
            cells_heat_3d->at(i).calcDiffEnergy(timeStep);
        }

        // Calc. conv. energy fluxes.
        if (cells_heat_3d->at(i).getT() < tempThreshold && i != 0 && 
            i != cells_heat_3d->size() - 1)
        {
            cells_heat_3d->at(i).calcConvEnergy(timeStep);
        }

        // Solve the new temperature.
        if (cells_heat_3d->at(i).getT() < tempThreshold && i != 0 && 
            i != cells_heat_3d->size() - 1)
        {
            cells_heat_3d->at(i).setT(bisection(cells_heat_3d->at(i), 
                threshLeftBis, threshRightBis));
        }

        // Move water between frozen and liquid phases in the water cell.
        CellWater3d* water = cells_heat_3d->at(i).getWater();
        water->setFracIce(cells_heat_3d->at(i).getFracIce());

        // Change energy state of the surface cell.
        //CellSurfAtmosAgri* atmos = heat->at(iInd).getAtmos();
        //if (atmos != 0)
        //{
        //    CellSubGeomAgri* geom = heat->at(iInd).getGeom();
        //    atmos->setSoilDiffFlux(heat->at(iInd).getSurfDiffFlux() / 
        //        geom->getSideArea(0));
        //}
    }
}

double ModelHeat3dTri::bisection(CellHeat3d& heat, double left, double right)
{
    while (fabs(right - left) > 2.0 * iterThreshBis)
    {
        //Calculate midpoint of domain.
        double midpoint = 0.5 * (right + left);

        //Find f(midpoint).
        if (heat.calcResidual(left) * 
            heat.calcResidual(midpoint) > 0.0)
        {
            //Throw away left half.
            left = midpoint;
        }
        else
        {
            //Throw away right half.
            right = midpoint;
        }
    }
    return 0.5 * (right + left);
}
