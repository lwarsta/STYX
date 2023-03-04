#include "AnalyticalModels.h"

void AnalyticalModels::run_comparisons()
{
    run_1d_flow_comparison();
    run_1d_transport_comparison();
}

void AnalyticalModels::run_1d_flow_comparison()
{
    // Initialize parameters.
    double wat_cont_res = 0.15;
    double wat_cont_sat = 0.45;
    double compress = 0.0001;
    double vgAlpha = 7.0;
    double vgN = 2.0;
    double dryWeight = 1680.0;
    double timeStep = 0.05;
    double dx = 5.0;
    double h_r = -100.0;
    double h_s = 0.0;
    double condSat = 10.0;
    double L = 205.0;
    int N = 41;
    double simTimeMax = -5.0 * L * L * ((wat_cont_sat - wat_cont_res) / 
        (6.0 * h_r * condSat));
    
    // Load data, tokenize data and build the 2d mesh.
    std::string mesh2dfilePathIn = "mesh2d_analytical_tests_1d_vert_41.vtk";
    std::vector < std::vector<std::string> > tokens2d;
    tokens2d = fileIO.load_and_tokenize_file(mesh2dfilePathIn, ' ');
    Grid2d grid2d;
    grid2d.build_grid(tokens2d);
    grid2d.create_water_cells();

    // Load data, tokenize data and build the 3d mesh.
    std::string mesh3dfilePathIn = "mesh3d_analytical_tests_1d_vert_41.vtk";
    std::vector < std::vector<std::string> > tokens3d;
    tokens3d = fileIO.load_and_tokenize_file(mesh3dfilePathIn, ' ');
    Grid3d grid3d;
    grid3d.build_grid(tokens3d);
    grid3d.create_water_cells();

    // Parametrise the mesh.
    std::vector<CellWater3d>* waterCells3d = grid3d.get_water_cells();

    for (size_t i = 0; i < waterCells3d->size(); i++)
    {
        waterCells3d->at(i).setId(i);
        waterCells3d->at(i).setWatContRes(wat_cont_res);
        waterCells3d->at(i).setWatContSat(wat_cont_sat);
        waterCells3d->at(i).setCondSat(condSat);
        waterCells3d->at(i).setCompress(compress);
        waterCells3d->at(i).setAlpha(vgAlpha);
        waterCells3d->at(i).setN(vgN);
        waterCells3d->at(i).setDryWeight(dryWeight);
        waterCells3d->at(i).setHydrHead(0.0);
        waterCells3d->at(i).calcPresHead();
        waterCells3d->at(i).calcWatCont();
        waterCells3d->at(i).calcDiffWatCap();
    }

    for (size_t i = 0; i < waterCells3d->size(); i++)
    {
        waterCells3d->at(i).swap();
        waterCells3d->at(i).calcUnsatCond();
    }

    for (size_t i = 0; i < waterCells3d->size(); i++)
    {
        waterCells3d->at(i).calcCondInter();
        waterCells3d->at(i).calcFluxes();
        waterCells3d->at(i).calcFlowVel();
    }

    // Compute initial pressure heads in the grid.
    for (int i = 0; i < N; i++)
    {
        double x = i * dx + 0.5 * dx;
        double pressHead = h_r * (1.0 - 1.0 / 6.0 * (x - L) / L * (x - L) / L);
        waterCells3d->at(i).setHydrHead(pressHead); // elevation head is zero

    }

    // Run the simulation.
    double simTime = 0.0;

    while (simTime < simTimeMax)
    {
        // Run the model.
        //modelWater3dBrute.configure(1000, 0.000000000000000001, 1.0, timeStep);
        //modelWater3dBrute.run(grid2d, grid3d);
        flow_model(grid3d, timeStep, 0.000000000000000001, 1000);

        // Set bounadry conditions.
        double pressBound = h_r * (1.0 - 1.0 / 
            (6.0 - 5.0 * simTime / simTimeMax));

        //std::cout << "\nPressure boundary: " << pressBound;
        waterCells3d->at(0).setHydrHead(pressBound);
        waterCells3d->at(0).calcPresHead();

        double press_head_sat = 0.0;
        double press_head_res = -100.0;
        double watContNew = wat_cont_res + (wat_cont_sat - wat_cont_res) / 
            (press_head_sat - press_head_res) * (pressBound - press_head_res);
        waterCells3d->at(0).setWatCont(watContNew);
        waterCells3d->at(0).swap();
        waterCells3d->at(0).calcDiffWatCap();
        simTime += timeStep;
    }

    // Compute the results.
    double mse = 0.0;

    for (int i = 0; i < N; i++)
    {
        double x = 0.5 * dx + i * dx;
        double press_analytic = h_r * (1.0 - ((x - L) / L) * ((x - L) / L)) /
            (6.0 - 5.0 * simTime / simTimeMax);
        double press_model = waterCells3d->at(i).getPresHead();
        mse += (press_analytic - press_model) * (press_analytic - press_model);
        //std::cout << "\n" << x << "," << press_analytic << "," << press_model;
    }

    std::cout << "\n-> 1d analytical water flow test MSE: " << mse / N;
}

void AnalyticalModels::flow_model(Grid3d& grid3d, double timeStep, 
    double iterCutThresh, int iterStop)
{
    std::vector<CellWater3d>* waterCells3d = grid3d.get_water_cells();

    // Calculate cell properties.
    for (size_t i = 0; i < waterCells3d->size(); i++)
    {
        waterCells3d->at(i).swap();
        // Compute unsaturated hydraulic conductivity.
        double press_head = waterCells3d->at(i).getHydrHead();
        double press_head_sat = 0.0;
        double press_head_res = -100.0;
        double relCond = (press_head - press_head_res) / 
            (press_head_sat - press_head_res);
        double condSat = waterCells3d->at(i).getCondSat();
        waterCells3d->at(i).setCondUnsat(relCond * condSat);
    }

    // Calculate inter cell properties.
    for (size_t i = 0; i < waterCells3d->size(); i++)
    {
        waterCells3d->at(i).calcCondInter();
        waterCells3d->at(i).calcFluxes();
    }

    // Iterate new solution.
    double iterDev;
    int iterCount = 0;

    do
    {
        iterDev = 0.0;
        for (size_t i = 0; i < waterCells3d->size(); i++)
        {
            CellGeom3d* geom = waterCells3d->at(i).getGeom();
            double numerator = waterCells3d->at(i).getDiffWatCap() * 
                geom->getVolume() / 
                timeStep * waterCells3d->at(i).getHydrHeadOld();
            double denominator = waterCells3d->at(i).getDiffWatCap() * 
                geom->getVolume() / timeStep;
            
            // Include fluxes between neighbor cells into the solution.
            for (size_t j = 0; j < geom->getNumOfNeigh(); j++)
            {
                CellWater3d* neigh = waterCells3d->at(i).getNeigh(j);
                if (neigh != 0)
                {
                    numerator -= waterCells3d->at(i).getFlux(j) * 
                        neigh->getHydrHead();
                    denominator -= waterCells3d->at(i).getFlux(j);
                }
            }

            // Compute the new head. Check that denominator is not zero!
            double hydrHeadNew = numerator / denominator;

            // Save iteration deviation.
            double iterDevNew = fabs(waterCells3d->at(i).getHydrHead() - 
                hydrHeadNew);
            if (iterDevNew > iterDev)
            {
                iterDev = iterDevNew;
            }

            // Set the new head.
            waterCells3d->at(i).setHydrHead(hydrHeadNew);

            // Update differential water capacity.
            waterCells3d->at(i).calcPresHead();
            double wat_cont_sat = waterCells3d->at(i).getWatContSat();
            double wat_cont_res = waterCells3d->at(i).getWatContRes();
            double press_head = hydrHeadNew;
            double press_head_sat = 0.0;
            double press_head_res = -100.0;
            double watContNew = wat_cont_res + (wat_cont_sat - wat_cont_res) / 
                (press_head_sat - press_head_res) * 
                (press_head - press_head_res);
            //std::cout << "\nWat. cont.: " << watContNew;
            waterCells3d->at(i).setWatCont(watContNew);
            waterCells3d->at(i).calcDiffWatCap();
        }

        iterCount++;
        //std::cout << "\nIter. count: " << iterCount;
        //break;
    } while (iterDev > iterCutThresh && iterCount < iterStop);
}

double AnalyticalModels::calcWatCont(double wat_cont_sat, double wat_cont_res,
    double press_head, double press_head_sat, double press_head_res)
{
    return wat_cont_res + (wat_cont_sat - wat_cont_res) / 
        (press_head_sat - press_head_res) * (press_head - press_head_res);
}

double AnalyticalModels::calcRelCond(double press_head, double press_head_sat,
    double press_head_res)
{
    return (press_head - press_head_res) / (press_head_sat - press_head_res);
}

void AnalyticalModels::run_1d_transport_comparison()
{
    std::cout << "\nRunnning a 1d analytical transport test:";

    // Set parameters.
    double porosity = 0.25; //  0.25; 1.0
    double c_in = 100.0;
    double v = 0.1;
    double d = 0.06; // 0.06 - 0.5 * 0.1 * 1 = 0.01 when numerical dispersion is removed
    double kb = 0.001;
    double timeStep = 0.05;
    double simTimeMax = 200.0;
    double dx = 1.0;
    double l = 40.5;
    int n = 41;
    double dispLong = 0.1414211125; // this leads to dispersion of 0.01
    double dispTrans = 0.01414211125;
    double condSat = 0.025;
    double vgAlpha = 7.0;
    double vgN = 2.0;
    double dryWeight = 1680.0;
    double compress = 0.0001;
    double watContRes = 0.0;
    int numOfSolutes = 1;
    int modelType = 2;

    // Load data, tokenize data and build the 2d mesh.
    std::string mesh2dfilePathIn = "mesh2d_analytical_tests_1d_vert_41.vtk";
    std::vector < std::vector<std::string> > tokens2d;
    tokens2d = fileIO.load_and_tokenize_file(mesh2dfilePathIn, ' ');
    Grid2d grid2d;
    grid2d.build_grid(tokens2d);
    grid2d.create_water_cells();
    grid2d.create_solute_cells(1);

    // Load data, tokenize data and build the 3d mesh.
    std::string mesh3dfilePathIn = "mesh3d_analytical_tests_1d_vert_41.vtk";
    std::vector < std::vector<std::string> > tokens3d;
    tokens3d = fileIO.load_and_tokenize_file(mesh3dfilePathIn, ' ');
    Grid3d grid3d;
    grid3d.build_grid(tokens3d);
    grid3d.create_water_cells();
    grid3d.create_solute_cells(numOfSolutes);

    // Parametrise the mesh.
    std::vector<CellWater3d>* waterCells3d = grid3d.get_water_cells();

    for (size_t i = 0; i < waterCells3d->size(); i++)
    {
        waterCells3d->at(i).setId(i);
        waterCells3d->at(i).setWatContRes(watContRes);
        waterCells3d->at(i).setWatContSat(porosity);
        waterCells3d->at(i).setCondSat(condSat);
        waterCells3d->at(i).setCompress(compress);
        waterCells3d->at(i).setAlpha(vgAlpha);
        waterCells3d->at(i).setN(vgN);
        waterCells3d->at(i).setDryWeight(dryWeight);
        waterCells3d->at(i).setHydrHead(-0.5 - (double)i);
        waterCells3d->at(i).calcPresHead(); 
        waterCells3d->at(i).calcWatCont();
        waterCells3d->at(i).calcDiffWatCap();
    }

    for (size_t i = 0; i < waterCells3d->size(); i++)
    {
        waterCells3d->at(i).swap();
        waterCells3d->at(i).calcUnsatCond();
    }

    for (size_t i = 0; i < waterCells3d->size(); i++)
    {
        waterCells3d->at(i).calcCondInter();
        waterCells3d->at(i).calcFluxes();
        waterCells3d->at(i).calcFlowVel();
    }

    std::vector<CellSolute3d>* soluteCells3d = grid3d.get_solute_cells();

    for (size_t i = 0; i < soluteCells3d->size(); i++)
    {
        soluteCells3d->at(i).setId(i);
        soluteCells3d->at(i).setDispCofs(dispLong, dispTrans);
        soluteCells3d->at(i).setMolDiffRate(0.0);
        soluteCells3d->at(i).setAdsParam(0, -1);
        soluteCells3d->at(i).setAdsParam(1, 0);
        soluteCells3d->at(i).setAdsParam(2, 1);
        soluteCells3d->at(i).setDecayRate(0, kb);
        soluteCells3d->at(i).setDecayTarget(0, -1);
        soluteCells3d->at(i).setDecayRate(1, 0);
        soluteCells3d->at(i).setDecayTarget(1, -1);
    }

    // Calculate cell properties.
    for (size_t i = 0; i < n; i++)
    {
        soluteCells3d->at(i).calcDisp();
        soluteCells3d->at(i).calcAdvFluxes();
        soluteCells3d->at(i).calcAdsMult();
    }

    // Calculate inter cell properties.
    for (int i = 0; i < n; i++)
    {
        soluteCells3d->at(i).calcDispInter();
        soluteCells3d->at(i).calcDispFluxes();
    }

    /*
    // Print values.
    for (size_t i = 0; i < n; i++)
    {
        CellGeom3d* geom = waterCells3d->at(i).getGeom();
        std::cout << "\n" << geom->getId();
        Vertex cp = geom->getCentrePoint();
        std::cout << ",(" << cp.x << ";" << cp.y << ";" << cp.z << ")";
        std::cout << "," << waterCells3d->at(i).getHydrHead();
        std::cout << "," << waterCells3d->at(i).getPresHead();
        std::cout << "," << waterCells3d->at(i).getWatCont();
        std::cout << "," << waterCells3d->at(i).getCondUnsat();
        std::cout << "," << waterCells3d->at(i).getVelInt(0);
        std::cout << "," << waterCells3d->at(i).getVelInt(5);
        std::cout << "," << soluteCells3d->at(i).getAdvFlux(0, 0);
        std::cout << "," << soluteCells3d->at(i).getAdvFlux(0, 1);
        std::cout << "," << soluteCells3d->at(i).getAdvFlux(1, 0);
        std::cout << "," << soluteCells3d->at(i).getAdvFlux(1, 1);
        std::cout << "," << soluteCells3d->at(i).getAdvFlux(2, 0);
        std::cout << "," << soluteCells3d->at(i).getAdvFlux(2, 1);
        std::cout << "," << soluteCells3d->at(i).getAdvFlux(3, 0);
        std::cout << "," << soluteCells3d->at(i).getAdvFlux(3, 1);
        std::cout << "," << soluteCells3d->at(i).getAdvFlux(4, 0);
        std::cout << "," << soluteCells3d->at(i).getAdvFlux(4, 1);
        std::cout << "," << soluteCells3d->at(i).getAdvFlux(5, 0);
        std::cout << "," << soluteCells3d->at(i).getAdvFlux(5, 1);
        std::cout << "," << soluteCells3d->at(i).getDispFlux(0);
        std::cout << "," << soluteCells3d->at(i).getDispFlux(1);
        std::cout << "," << soluteCells3d->at(i).getDispFlux(2);
        std::cout << "," << soluteCells3d->at(i).getDispFlux(3);
        std::cout << "," << soluteCells3d->at(i).getDispFlux(4);
        std::cout << "," << soluteCells3d->at(i).getDispFlux(5);
    }
    */
    // Set initial conditions.
    soluteCells3d->at(0).setConc(c_in);
    soluteCells3d->at(0).setMass(0.0);
    soluteCells3d->at(0).convertConcToMass();
    soluteCells3d->at(soluteCells3d->size() - 1).setConc(0.0);
    soluteCells3d->at(soluteCells3d->size() - 1).setMass(0.0);
    soluteCells3d->at(soluteCells3d->size() - 1).convertConcToMass();

    // Run the simulation.
    double simTime = 0.0;

    while (simTime < simTimeMax)
    {
        // Run the model.
        if (modelType == 1)
        {
            modelSolute3dBrute.configure(1000, 0.000000000000000001, 1.0, 
                timeStep, 1);
            modelSolute3dBrute.run(grid2d, grid3d);
        }
        //else if (modelType == 2)
        //{
        //    modelSolute3dTri.configure(1000, 0.000000000000000001, 1.0, 
        //        timeStep, 1);
        //    modelSolute3dTri.run(grid2d, grid3d);
        //}

        // Set boundary conditions.
        soluteCells3d->at(0).setConc(c_in);
        soluteCells3d->at(0).setMass(0.0);
        soluteCells3d->at(0).convertConcToMass();
        soluteCells3d->at(soluteCells3d->size() - 1).setConc(0.0);
        soluteCells3d->at(soluteCells3d->size() - 1).setMass(0.0);
        soluteCells3d->at(soluteCells3d->size() - 1).convertConcToMass();
        simTime += timeStep;
        //std::cout << "\nCompleted " << simTime / simTimeMax * 100.0 << " %";
    }
    // Compute the results.
    double mse = 0.0;

    for (int i = 0; i < n; i++)
    {
        double z = 0.5 * dx + i * dx;
        double conc_analytic = transport1d(c_in, v, d, kb, simTimeMax, z);
        double conc_model = soluteCells3d->at(i).getConc();
        mse += (conc_analytic - conc_model) * (conc_analytic - conc_model);
        //std::cout << "\n" << x << "," << conc_analytic << "," << conc_model;
    }

    std::cout << "\n-> 1d analytical solute transport test MSE: " << mse / n;
}

double AnalyticalModels::transport1d(double c_in, double v, double d, 
    double kb, double t, double x)
{
    // Problem 1 / Sun / 3.2.14
    if (v * x / d > 10000.0)
    {
        return 0.0;
    }
    else
    {
        return c_in / 2.0 * exp(0.5 * v * x / d) * (exp(-0.5 * x / d * 
            sqrt(v * v + 4.0 * kb * d)) * 
            erfc((x - sqrt(v * v + 4.0 * kb * d) * t) / sqrt(4.0 * d * t)) + 
            exp(0.5 * x / d * sqrt(v * v + 4.0 * kb * d)) * 
            erfc((x + sqrt(v * v + 4.0 * kb * d) * t) / sqrt(4.0 * d * t)));
    }
}
