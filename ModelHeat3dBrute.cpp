#include "ModelHeat3dBrute.h"

void ModelHeat3dBrute::configure(int iterStopNew, double iterCutThreshNew, 
    double implicityNew, double timeStepNew, double boundTempNew, 
    double iterThreshBisNew, int iterStopBisNew, double threshLeftBisNew,
    double threshRightBisNew)
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
}

void ModelHeat3dBrute::run(Grid2d &grid2d, Grid3d &grid3d)
{
    preprocess(grid2d, grid3d);
    iterate(grid2d, grid3d);
    postprocess(grid2d, grid3d);
}

void ModelHeat3dBrute::preprocess(Grid2d &grid2d, Grid3d &grid3d)
{
    std::vector<CellHeat2d> * cells_heat_2d = grid2d.get_heat_cells();
	std::vector<CellHeat3d> * cells_heat_3d = grid3d.get_heat_cells();

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

void ModelHeat3dBrute::iterate(Grid2d &grid2d, Grid3d &grid3d)
{
    std::vector<CellHeat3d> * cells_heat_3d = grid3d.get_heat_cells();

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

        // Compute the new temperature.
        double tempExpNew = 0.0;

        if (denominator != 0.0)
        {
            tempExpNew = numerator / denominator;
        }
        if (tempExpNew < 0.0)
        {
            tempExpNew = 0.0;
        }

        // Set the new temperature.
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
        //omp_set_num_threads(2);
        #pragma omp parallel for
		for (int i = 0; i < cells_heat_3d->size(); i++)
		{
            CellGeom3d * geom = cells_heat_3d->at(i).getGeom();

            // Define numerator and denominator.
            double numerator = cells_heat_3d->at(i).getC() * 
                geom->getVolume() / timeStep * cells_heat_3d->at(i).getOldT();
            double denominator = cells_heat_3d->at(i).getC() * 
                geom->getVolume() / timeStep;

            // Include fluxes between neighbor cells into the solution.
            for (size_t j = 0; j < geom->getNumOfNeigh(); j++)
            {
                CellHeat3d * neigh = cells_heat_3d->at(i).getNeigh(j);

                if (neigh != 0)
                {
                    numerator += cells_heat_3d->at(i).getDiffFlux(j) * 
                        neigh->getT();
                    denominator += cells_heat_3d->at(i).getDiffFlux(j);
                    numerator += cells_heat_3d->at(i).getConvFlux(j, 0) * 
                        neigh->getT();
                    denominator -= cells_heat_3d->at(i).getConvFlux(j, 1);
                }
            }

            // Heat flux between soil surface and top soil layer.
            numerator += cells_heat_3d->at(i).getSurfDiffFlux();

            // Tom boundary condition - TEMPORARY
            if (i == 0)
            {
                denominator = 1.0;
                numerator = boundTemp; // boundary condition here (e.g. 10.0).
            }

            /*
            // Bottom boundary condition - TEMPORARY
            if (i == cells_heat_3d->size() - 1)
            {
                denominator = 1.0;
                numerator = 5.0; // boundary condition here.
            }
            */

            // Compute the new temperature. Check that denominator is not zero!
            double TNew = numerator / denominator;

            // Save iteration deviation.
            double iterDevNew = fabs(cells_heat_3d->at(i).getT() - TNew);

            if (iterDevNew > iterDev)
            {
                iterDev = iterDevNew;
            }

            // Save the new temperature.
            cells_heat_3d->at(i).setT((1.0 - implicity) * tempExp.at(i) + 
                implicity * TNew);
		}

		iterCount++;
	} while (iterDev > iterCutThresh && iterCount < iterStop);
}

void ModelHeat3dBrute::postprocess(Grid2d &grid2d, Grid3d &grid3d)
{
    std::vector<CellHeat3d> * cells_heat_3d = grid3d.get_heat_cells();
    //std::cout << "\n8<--------------------------------------";
    for (size_t i = 0; i < cells_heat_3d->size(); i++)
    {
        // Calc. U [kJ].
        // CHECK WHY THE THRESHOLD VALUE IS USED BELOW!
        // test this without the threshold (1.0) with a short time step
        double tempThreshold = 1000.0;

        if (cells_heat_3d->at(i).getT() < tempThreshold)
        {
            cells_heat_3d->at(i).calcU();
        }

        /*
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
        */
        // Move water between frozen and liquid phases in the water cell.
        //CellWater3d * water = cells_heat_3d->at(i).getWater();
        //water->setFracIce(cells_heat_3d->at(i).getFracIce());
        //std::cout << "\n" << i << ": " << cells_heat_3d->at(i).getT();
        
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

double ModelHeat3dBrute::bisection(CellHeat3d & heat, double left, 
    double right)
{
    while (fabs(right - left) > 2.0 * iterThreshBis)
    {
        //Calculate midpoint of domain.
        double midpoint = 0.5 * (right + left);

        //Find f(midpoint).
        if (heat.calcResidual(left) * heat.calcResidual(midpoint) > 0.0)
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
