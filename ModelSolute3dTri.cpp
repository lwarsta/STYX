#include "ModelSolute3dTri.h"

void ModelSolute3dTri::configure(int iterStopNew, double iterCutThreshNew, 
    double implicityNew, double timeStepNew, int numOfSolutesNew, 
    std::vector < std::vector<int>> gridMapNew)
{
	iterStop = iterStopNew;
	iterCutThresh = iterCutThreshNew;
	implicity = implicityNew;
	timeStep = timeStepNew;
	numOfSolutes = numOfSolutesNew;
    gridMap = gridMapNew;
}

void ModelSolute3dTri::run(Grid2d &grid2d, Grid3d &grid3d)
{
    preprocess(grid2d, grid3d);
    iterate(grid2d, grid3d);
    postprocess(grid2d, grid3d);
}

void ModelSolute3dTri::preprocess(Grid2d &grid2d, Grid3d &grid3d)
{
    std::vector<CellSolute2d> * soluteCells2d = grid2d.get_solute_cells();
	std::vector<CellSolute3d> * soluteCells3d = grid3d.get_solute_cells();
	size_t cellsPerSolute2d = 0;

	if (numOfSolutes > 0)
    {
        cellsPerSolute2d = soluteCells2d->size() / numOfSolutes;
    }

	size_t cellsPerSolute3d = 0;

	if (numOfSolutes > 0)
    {
        cellsPerSolute3d = soluteCells3d->size() / numOfSolutes;
    }

	for (size_t i = 0; i < numOfSolutes; i++)
	{
        // Calculate cell properties.
        for (size_t j = 0; j < cellsPerSolute3d; j++)
        {
            int ind_3d = i * cellsPerSolute3d + j;
            soluteCells3d->at(ind_3d).convertMassToConc();
            soluteCells3d->at(ind_3d).swap();
            soluteCells3d->at(ind_3d).calcDisp();
            soluteCells3d->at(ind_3d).calcAdvFluxes();
            soluteCells3d->at(ind_3d).calcAdsMult();
            // Compute solute infiltration flux.
            CellGeom3d * geom3d = soluteCells3d->at(ind_3d).getGeom();
            int cellGeom2dInd = geom3d->getGridConnection();

            if (cellGeom2dInd >= 0)
            {
                int ind_2d = i * cellsPerSolute2d + cellGeom2dInd;
                CellGeom2d * geom2d =  soluteCells2d->at(ind_2d).getGeom();
                CellWater2d * water2d =  soluteCells2d->at(ind_2d).getWater();
                CellWater3d * water3d = soluteCells3d->at(ind_3d).getWater();
                double infWatVol = water3d->getFluxInf() * timeStep;
                double surfWatVol = water2d->getWaterDepth() * 
                    geom2d->getArea();
                double concInf = 0.0;

                if (infWatVol >= 0.0 && surfWatVol >= 0.0 && 
                    soluteCells2d->at(ind_2d).getMass() > 0.0)
                {
                    // Fix the influx concentration if the influx 
                    // is taking more mass from the surface than is available.
                    concInf = soluteCells2d->at(ind_2d).getMass() / 
                        (infWatVol + surfWatVol);

                    if (concInf * infWatVol > 
                        soluteCells2d->at(ind_2d).getMass())
                    {
                        concInf = soluteCells2d->at(ind_2d).getMass() / 
                            infWatVol;
                    }

                    // If the surface solute concentration drops below 
                    // a certain level (kg m-3) set it as zero.
                    // If this is not done, adsorption can cause concentration
                    // to start to approach zero but not reach it.
                    //if (concInf < 0.000000000001)
                    //{
                    //    concInf = 0.0;
                    //}
                    soluteCells3d->at(ind_3d).setInfConc(concInf);
                }
                else
                {
                    soluteCells3d->at(ind_3d).setInfConc(0.0);
                }
            }
        }

        // Calculate inter cell properties.
        for (int j = 0; j < cellsPerSolute3d; j++)
        {
            int indexSol3d = i * cellsPerSolute3d + j;
            soluteCells3d->at(indexSol3d).calcDispInter();
            soluteCells3d->at(indexSol3d).calcDispFluxes();
        }
	}
}

void ModelSolute3dTri::triDiag(
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

void ModelSolute3dTri::iterate(Grid2d& grid2d, Grid3d& grid3d)
{
    //std::vector<CellSolute2d> * soluteCells2d = grid2d.getSoluteCells();
    std::vector<CellSolute3d>* soluteCells3d = grid3d.get_solute_cells();
    size_t cellsPerSolute = 0;

    if (numOfSolutes > 0)
    {
        cellsPerSolute = soluteCells3d->size() / numOfSolutes;
    }

    int numOfDecPaths = 2;

    for (size_t i = 0; i < numOfSolutes; i++)
    {
        // Compute explicit solution.
        std::vector<double> concExp;
        concExp.assign(cellsPerSolute, 0.0);

        #pragma omp parallel for
        for (int j = 0; j < cellsPerSolute; j++)
        {
            int ind = i * cellsPerSolute + j;
            CellGeom3d* geom = soluteCells3d->at(ind).getGeom();
            CellWater3d* water = soluteCells3d->at(ind).getWater();
            double watVol = water->getWatCont() * geom->getVolume();

            // Define numerator and denominator.
            double numerator = soluteCells3d->at(ind).getAdsMult() * watVol /
                timeStep * soluteCells3d->at(ind).getConcOld();
            double denominator = soluteCells3d->at(ind).getAdsMult() * 
                watVol / timeStep;

            // Include fluxes between neighbor cells into the solution.
            for (size_t k = 0; k < geom->getNumOfNeigh(); k++)
            {
                CellSolute3d* neigh = soluteCells3d->at(ind).getNeigh(k);

                if (neigh != 0)
                {
                    // Compute dispersive fluxes between cells.
                    numerator += soluteCells3d->at(ind).getDispFlux(k) *
                        (neigh->getConcOld() - 
                            soluteCells3d->at(ind).getConcOld()) /
                        neigh->getAdsMult() * 
                        soluteCells3d->at(ind).getAdsMult();
                    
                    // Compute advective fluxes between cells.
                    numerator += soluteCells3d->at(ind).getAdvFlux(k, 0) * 
                        neigh->getConcOld() / neigh->getAdsMult() * 
                        soluteCells3d->at(ind).getAdsMult();
                    numerator += soluteCells3d->at(ind).getAdvFlux(k, 1) * 
                        soluteCells3d->at(ind).getConcOld() / 
                        neigh->getAdsMult() * 
                        soluteCells3d->at(ind).getAdsMult();
                }
            }
            // Include infiltration and exfiltration fluxes into the solution.
            if (water->getFluxInf() > 0.0 && 
                soluteCells3d->at(ind).getInfConc() > 0.0)
            {
                numerator += water->getFluxInf() * 
                    soluteCells3d->at(ind).getInfConc();
            }

            // Include drain flux into the solution. 
            // Prevent solute flux from the drain into the soil
            if (water->getDrainVol() > 0.0 && 
                soluteCells3d->at(ind).getConcOld() > 0.0)
            {
                numerator -= water->getDrainVol() / timeStep * 
                    soluteCells3d->at(ind).getConcOld();
                soluteCells3d->at(ind).setDrainMass(water->getDrainVol() * 
                    soluteCells3d->at(ind).getConcOld() / 
                    soluteCells3d->at(ind).getAdsMult());
            }
            else
            {
                soluteCells3d->at(ind).setDrainMass(0.0);
            }

            // Include decay input into the solution.
            for (size_t k = 0; k < numOfSolutes; k++)
            {
                for (int l = 0; l < numOfDecPaths; l++)
                {
                    int indexInput = k * cellsPerSolute + j;

                    if (ind == indexInput)
                    {
                        continue;
                    }

                    int decTarg = 
                        soluteCells3d->at(indexInput).getDecayTarget(l);
                    int cellId = soluteCells3d->at(indexInput).getId();

                    if (decTarg == i && 
                        cellId != soluteCells3d->at(ind).getId())
                    {
                        double fluxDecayIn = 
                            soluteCells3d->at(indexInput).getFluxDecay(l) / 
                            soluteCells3d->at(indexInput).getAdsMult() * 
                            soluteCells3d->at(ind).getAdsMult();
                        numerator += fluxDecayIn;
                        soluteCells3d->at(ind).setDecayMassIn(l, 
                            fluxDecayIn * timeStep);
                    }
                }
            }

            // Include solute decay into the solution.
            for (int k = 0; k < numOfDecPaths; k++)
            {
                numerator -= soluteCells3d->at(ind).getFluxDecay(k);
                double dec_mass_out = soluteCells3d->at(ind).getFluxDecay(k) * 
                    timeStep / soluteCells3d->at(ind).getAdsMult();
                soluteCells3d->at(ind).setDecayMassOut(k, dec_mass_out);
            }

            // Compute the new concentration.
            double concExpNew = 0.0;

            if (denominator != 0.0)
            {
                concExpNew = numerator / denominator;
            }

            if (concExpNew < 0.0)
            {
                concExpNew = 0.0;
            }

            // Set the new concentration.
            concExp.at(j) = concExpNew;

            if (implicity <= 0.0)
            {
                soluteCells3d->at(ind).setConc(concExpNew);
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
                    // Empty vector position - is this needed?
                    int ind = i * cellsPerSolute + gridMap.at(m).at(n);
                    CellGeom3d* geom = soluteCells3d->at(ind).getGeom();
                    CellWater3d* water = soluteCells3d->at(ind).getWater();
                    double watVol = water->getWatCont() * geom->getVolume();
                    triB.at(n) = soluteCells3d->at(ind).getAdsMult() * watVol /
                        timeStep;
                    triD.at(n) = soluteCells3d->at(ind).getAdsMult() * watVol /
                        timeStep * soluteCells3d->at(ind).getConcOld();

                    // Vertical dispersion and advection at the top
                    triA.at(n) -= soluteCells3d->at(ind).getDispFlux(0);
                    triB.at(n) += soluteCells3d->at(ind).getDispFlux(0);
                    triA.at(n) -= soluteCells3d->at(ind).getAdvFlux(0, 0);
                    triB.at(n) -= soluteCells3d->at(ind).getAdvFlux(0, 1);

                    // Vertical dispersion and advection at the bottom
                    triC.at(n) -= soluteCells3d->at(ind).getDispFlux(5);
                    triB.at(n) += soluteCells3d->at(ind).getDispFlux(5);
                    triC.at(n) -= soluteCells3d->at(ind).getAdvFlux(5, 0);
                    triB.at(n) -= soluteCells3d->at(ind).getAdvFlux(5, 1);

                    // Horizontal fluxes between cell columns
                    for (size_t k = 1; k < geom->getNumOfNeigh() - 1; k++)
                    {
                        CellSolute3d* neigh = soluteCells3d->at(ind).getNeigh(k);

                        if (neigh != 0)
                        {
                            triD.at(n) += soluteCells3d->at(ind).getDispFlux(k) *
                                neigh->getConc() / neigh->getAdsMult() * 
                                soluteCells3d->at(ind).getAdsMult();
                            triB.at(n) += soluteCells3d->at(ind).getDispFlux(k);
                            triD.at(n) += soluteCells3d->at(ind).getAdvFlux(k, 0) *
                                neigh->getConc() / neigh->getAdsMult() * 
                                soluteCells3d->at(ind).getAdsMult();
                            triB.at(n) -= soluteCells3d->at(ind).getAdvFlux(k, 1);
                        }
                    }

                    // Infiltration and exfiltration.
                    if (water->getFluxInf() > 0.0 && 
                        soluteCells3d->at(ind).getInfConc() > 0.0)
                    {
                        triD.at(n) += water->getFluxInf() * 
                            soluteCells3d->at(ind).getInfConc();
                    }

                    // Drain flux.
                    // Prevent solute flux from the drain into the soil.
                    if (water->getDrainVol() > 0.0 && 
                        soluteCells3d->at(ind).getConcOld() > 0.0)
                    {
                        triD.at(n) -= water->getDrainVol() / timeStep * 
                            soluteCells3d->at(ind).getConcOld();
                        soluteCells3d->at(ind).setDrainMass(
                            water->getDrainVol() * 
                            soluteCells3d->at(ind).getConcOld() / 
                            soluteCells3d->at(ind).getAdsMult());
                    }
                    else
                    {
                        soluteCells3d->at(ind).setDrainMass(0.0);
                    }

                    // Include decay input into the solution.
                    for (size_t k = 0; k < numOfSolutes; k++)
                    {
                        for (int l = 0; l < numOfDecPaths; l++)
                        {
                            int ind_in = k * cellsPerSolute + 
                                gridMap.at(m).at(n); //j -> gridMap.at(m).at(n)

                            if (ind == ind_in)
                            {
                                continue;
                            }

                            int decTarg = soluteCells3d->at(ind_in).getDecayTarget(l);
                            int cellId = soluteCells3d->at(ind_in).getId();

                            if (decTarg == i && 
                                cellId != soluteCells3d->at(ind).getId())
                            {
                                double fluxDecayIn = soluteCells3d->at(ind_in).getFluxDecay(l) /
                                    soluteCells3d->at(ind_in).getAdsMult() * 
                                    soluteCells3d->at(ind).getAdsMult();
                                triD.at(n) += fluxDecayIn;
                                soluteCells3d->at(ind).setDecayMassIn(l, 
                                    fluxDecayIn * timeStep);
                            }
                        }
                    }

                    // Include solute decay into the solution.
                    for (int k = 0; k < numOfDecPaths; k++)
                    {
                        triD.at(n) -= soluteCells3d->at(ind).getFluxDecay(k);
                        soluteCells3d->at(ind).setDecayMassOut(k, 
                            soluteCells3d->at(ind).getFluxDecay(k) * timeStep /
                            soluteCells3d->at(ind).getAdsMult());
                    }
                }

                // Solve values in the column.
                triDiag(triA, triB, triC, triD, triG, triU, triV, triX);

                // Compute deviation and new concentration values.
                for (size_t n = 0; n < gridMap.at(m).size(); n++)
                {
                    // Save deviation.
                    int ind = i * cellsPerSolute + gridMap.at(m).at(n);
                    int ind_exp = gridMap.at(m).at(n);
                    double conc_new = (1.0 - implicity) * concExp.at(ind_exp) +
                        implicity * triX.at(n);
                    double iterDevNew = fabs(soluteCells3d->at(ind).getConc() -
                        conc_new);

                    // OpenMP double-checked locking.
                    if (iterDevNew > iterDev)
                    {
                        #pragma omp critical
                        if (iterDevNew > iterDev) iterDev = iterDevNew;
                    }

                    // Save new concentration value.
                    soluteCells3d->at(ind).setConc(conc_new);
                }
            }
            
            iterCount++;
        } while (iterDev > iterCutThresh && iterCount < iterStop);
    }
}

void ModelSolute3dTri::postprocess(Grid2d &grid2d, Grid3d &grid3d)
{
    std::vector<CellSolute2d> * soluteCells2d = grid2d.get_solute_cells();
    std::vector<CellSolute3d> * soluteCells3d = grid3d.get_solute_cells();
	size_t cellsPerSolute2d = 0;

	if (numOfSolutes > 0)
    {
        cellsPerSolute2d = soluteCells2d->size() / numOfSolutes;
    }

	size_t cellsPerSolute3d = 0;

	if (numOfSolutes > 0)
    {
        cellsPerSolute3d = soluteCells3d->size() / numOfSolutes;
    }

	for (size_t i = 0; i < numOfSolutes; i++)
	{
        for (size_t j = 0; j < cellsPerSolute3d; j++)
        {
            int indexSol3d = i * cellsPerSolute3d + j;
            CellGeom3d * geom =  soluteCells3d->at(indexSol3d).getGeom();
            CellWater3d * water = soluteCells3d->at(indexSol3d).getWater();

            // Infiltration.
            if (water->getFluxInf() > 0.0 && 
                soluteCells3d->at(indexSol3d).getInfConc() > 0.0)
            {
                double infMass = 0.0;
                infMass = timeStep * water->getFluxInf() * 
                    soluteCells3d->at(indexSol3d).getInfConc() / 
                    soluteCells3d->at(indexSol3d).getAdsMult();
                soluteCells3d->at(indexSol3d).setInfMass(infMass);
                CellGeom3d * geom3d = soluteCells3d->at(indexSol3d).getGeom();
                int cellGeom2dInd = geom3d->getGridConnection();

                if (cellGeom2dInd >= 0)
                {
                    int indexSol2d = i * cellsPerSolute2d + cellGeom2dInd;
                    soluteCells2d->at(indexSol2d).setMass(
                        soluteCells2d->at(indexSol2d).getMass() - infMass);
                }
            }
            else
            {
                soluteCells3d->at(indexSol3d).setInfMass(0.0);
            }

            // Convert solute concentration back to mass.
            soluteCells3d->at(indexSol3d).convertConcToMass();
        }
	}
}

