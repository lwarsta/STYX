#include "Framework.h"

Framework::Framework()
{

}

Framework::~Framework()
{

}

int Framework::initialize(std::string pathToSettings)
{
    // Load settings csv file and extract parameter values.
    std::cout << "Loading settings file and extracting settings: " << 
        pathToSettings.c_str() << "\n";
    std::vector < std::vector<std::string> > settings_data;
    settings_data = fileIO.load_and_tokenize_file(pathToSettings, ',');
    
    if (settings_data.size() == 0)
    {
        std::cout << "-> Settings file was not found or empty\n";
        return 1;
    }

    settings.initialize(settings_data);

    // Run analytical model tests.
    if (settings.get_int("run_tests"))
    {
        std::cout << "Running analytical model tests:\n";
        AnalyticalModels analytical_models;
        analytical_models.run_comparisons();
    }
    else
    {
        std::cout << "Bypassing analytical model tests:\n";
    }

    // Set OpenMP parallel threads.
    omp_set_dynamic(0);
    omp_set_num_threads(settings.get_int("num_of_threads"));
    
    // Load and tokenize network junction geometry.
    std::cout << "Loading network junction geometry:\n";
    std::vector < std::vector<std::string> > tokens_junc;
    tokens_junc = fileIO.load_and_tokenize_file(
        settings.get_str("network_junc_file_path_in"), ' ');

    if (tokens_junc.size() == 0)
    {
        std::cout << "-> Junction geometry file was not found or empty\n";
        std::string path = settings.get_str("network_junc_file_path_in");
        std::cout << path.c_str() << "\n";
        return 1;
    }

    // Load and tokenize network link geometry.
    std::cout << "Loading network link geometry:\n";
    std::vector < std::vector<std::string> > tokens_link;
    tokens_link = fileIO.load_and_tokenize_file(
        settings.get_str("network_link_file_path_in"), ' ');

    if (tokens_link.size() == 0)
    {
        std::cout << "-> Link geometry file was not found or empty\n";
        std::string path = settings.get_str("network_link_file_path_in");
        std::cout << path.c_str() << "\n";
        return 1;
    }

    // Build the network.
    std::cout << "Initializing network geometry:\n";
    network.build_network(tokens_junc, tokens_link);

    // Load and tokenize data and build the 2d mesh.
    std::cout << "Initializing 2d geometry cells:\n";
    std::vector < std::vector<std::string> > tokens2d;
    tokens2d = fileIO.load_and_tokenize_file(
        settings.get_str("mesh_2d_file_path_in"), ' ');
    
    if (tokens2d.size() == 0)
    {
        std::cout << "-> 2d geometry file was not found or empty\n";
        std::string path = settings.get_str("mesh_2d_file_path_in");
        std::cout << path.c_str() << "\n";
        return 1;
    }

    grid2d.build_grid(tokens2d);
    
    // Load and tokenize data and build the 3d mesh.
    std::cout << "Initializing 3d geometry cells:\n";
    std::vector < std::vector<std::string> > tokens3d;
    tokens3d = fileIO.load_and_tokenize_file(
        settings.get_str("mesh_3d_file_path_in"), ' ');
    
    if (tokens3d.size() == 0)
    {
        std::cout << "-> 3d geometry file was not found or empty\n";
        return 1;
    }

    grid3d.build_grid(tokens3d);

    // Load network materials, initial conditions and boundary conditions.
    std::vector<std::vector<std::string>> materials_net_junc;
    materials_net_junc = fileIO.load_and_tokenize_file(
        settings.get_str("materials_net_junc_path"), ',');

    if (materials_net_junc.size() == 0)
    {
        std::cout << "-> net junc material file was not found or empty\n";
        return 1;
    }

    std::vector<std::vector<std::string>> init_cond_net_junc;
    init_cond_net_junc = fileIO.load_and_tokenize_file(
        settings.get_str("init_cond_net_junc_path"), ',');

    if (init_cond_net_junc.size() == 0)
    {
        std::cout << "-> net junc initial conditions file was not found or empty\n";
        return 1;
    }

    std::vector < std::vector<std::string> > bound_cond_net_junc;
    bound_cond_net_junc = fileIO.load_and_tokenize_file(
        settings.get_str("bound_cond_net_junc_path"), ',');

    if (bound_cond_net_junc.size() == 0)
    {
        std::cout << "-> net junc boundary condition file was not found or empty\n";
        return 1;
    }

    std::vector<std::vector<std::string>> materials_net_link;
    materials_net_link = fileIO.load_and_tokenize_file(
        settings.get_str("materials_net_link_path"), ',');

    if (materials_net_link.size() == 0)
    {
        std::cout << "-> net link material file was not found or empty\n";
        return 1;
    }

    std::vector<std::vector<std::string>> init_cond_net_link;
    init_cond_net_link = fileIO.load_and_tokenize_file(
        settings.get_str("init_cond_net_link_path"), ',');

    if (init_cond_net_link.size() == 0)
    {
        std::cout << "-> net link initial conditions file was not found or empty\n";
        return 1;
    }

    std::vector < std::vector<std::string> > bound_cond_net_link;
    bound_cond_net_link = fileIO.load_and_tokenize_file(
        settings.get_str("bound_cond_net_link_path"), ',');

    if (bound_cond_net_link.size() == 0)
    {
        std::cout << "-> net link boundary condition file was not found or empty\n";
        return 1;
    }

    // Load 2d materials, initial conditions and boundary conditions.
    std::cout << "Loading 2d materials, initial conditions and" 
                 "boundary conditions:\n";
    std::vector < std::vector<std::string> > materials_2d;
    materials_2d = fileIO.load_and_tokenize_file(
        settings.get_str("materials_2d_path"), ',');
    
    if (materials_2d.size() == 0)
    {
        std::cout << "-> 2d material file was not found or empty\n";
        return 1;
    }
    
    std::vector < std::vector<std::string> > init_cond_2d;
    init_cond_2d = fileIO.load_and_tokenize_file(
        settings.get_str("init_cond_2d_path"), ',');
    
    if (init_cond_2d.size() == 0)
    {
        std::cout << "-> 2d initial conditions file was not found or empty\n";
        return 1;
    }
    
    std::vector < std::vector<std::string> > bound_cond_2d;
    bound_cond_2d = fileIO.load_and_tokenize_file(
        settings.get_str("bound_cond_2d_path"), ',');
    
    if (bound_cond_2d.size() == 0)
    {
        std::cout << "-> 2d boundary condition file was not found or empty\n";
        return 1;
    }

    // Load 3d materials, initial conditions and boundary conditions.
    std::cout << "Loading 3d materials, initial conditions and" 
                 "boundary conditions:\n";
    std::vector < std::vector<std::string> > materials_3d;
    materials_3d = fileIO.load_and_tokenize_file(
        settings.get_str("materials_3d_path"), ',');
    
    if (materials_3d.size() == 0)
    {
        std::cout << "-> 3d material file was not found or empty\n";
        return 1;
    }

    std::vector < std::vector<std::string> > init_cond_3d;
    init_cond_3d = fileIO.load_and_tokenize_file(
        settings.get_str("init_cond_3d_path"), ',');
    
    if (init_cond_3d.size() == 0)
    {
        std::cout << "-> 3d initial conditions file was not found or empty\n";
        return 1;
    }

    std::vector < std::vector<std::string> > bound_cond_3d;
    bound_cond_3d = fileIO.load_and_tokenize_file(
        settings.get_str("bound_cond_3d_path"), ',');
    
    if (bound_cond_3d.size() == 0)
    {
        std::cout << "-> 3d boundary condition file was not found or empty\n";
        return 1;
    }
    
    // Create and initialize 2d water cells.
    std::cout << "Creating and initializing 2d water cells:\n";
    grid2d.create_water_cells();
    grid2d.init_water_cells(settings, materials_2d, bound_cond_2d, 
        init_cond_2d);

    // Create and initialize the water network.
    std::cout << "Creating and initializing the water network:\n";
    network.create_water_network_items();
    network.init_water_network(settings, materials_net_junc, bound_cond_net_junc,
                               init_cond_net_junc, materials_net_link, 
                               bound_cond_net_link, init_cond_net_link);
    
    // Create and initialize 3d water cells.
    std::cout << "Creating and initializing 3d water cells:\n";
    grid3d.create_water_cells();
    grid3d.init_water_cells(settings, materials_3d, bound_cond_3d, 
        init_cond_3d);

    // Create and initialize 2d heat cells.
    std::cout << "Creating and initializing 2d heat cells:\n";
    grid2d.create_heat_cells();
    grid2d.init_heat_cells(settings, materials_2d, bound_cond_2d, 
        init_cond_2d);
    
    // Create and initialize 3d heat cells.
    std::cout << "Creating and initializing 3d heat cells:\n";
    grid3d.create_heat_cells();
    grid3d.init_heat_cells(settings, materials_3d, bound_cond_3d, 
        init_cond_3d);
    
    // Load solute library from file.
    std::cout << "Loading solute library file:\n";
    std::vector < std::vector<std::string> > soluteLib;
    soluteLib = fileIO.load_and_tokenize_file(
        settings.get_str("solute_properties_path"), ',');

    if (soluteLib.size() == 0)
    {
        std::cout << "-> Solute library file was not found or empty\n";
        return 1;
    }

    // Define simulated species either via PHREEQ or via solute library file.
    std::cout << "Defining simulated species:\n";
    int num_of_species = 0;

    if (settings.get_int("use_phreeqcrm") == 1)
    {
        phreeqcWrapper.initialize(settings, bound_cond_2d, bound_cond_3d, 
            grid2d, grid3d);
        num_of_species = phreeqcWrapper.get_num_of_species();
        species_names = phreeqcWrapper.get_species_names();
    }
    else
    {
        num_of_species = (int)(soluteLib.size() - 1);
        for (size_t i = 0; i < soluteLib.size() - 1; i++)
        {
            //std::cout << "Solute: " << soluteLib.at(i + 1).at(0) << "\n";
            species_names.push_back(soluteLib.at(i + 1).at(0));
        }
    }

    settings.set_value("num_of_species", num_of_species);
    
    // Create and initialize 2d solute cells.
    std::cout << "Creating and initializing 2d solute cells:\n";
    grid2d.create_solute_cells(num_of_species);
    grid2d.init_solute_cells(settings, materials_2d, bound_cond_2d, 
        init_cond_2d, soluteLib, species_names);
    
    // Create and initialize 3d solute cells.
    std::cout << "Creating and initializing 3d solute cells:\n";
    grid3d.create_solute_cells(num_of_species);
    std::vector<double> * phreeqc_conc = phreeqcWrapper.get_phreeqc_conc();
    std::map<std::string, double> * mol_weight_map;
    mol_weight_map = phreeqcWrapper.get_mol_weight_map();
    // Currently initial concentrations are derived from boundary conditions 
    // in a bit convoluted way.
    std::vector<double> conc_init;
    std::vector<CellSolute2d>* cells_solute_2d = grid2d.get_solute_cells();
    std::vector<CellGeom2d>* cells_geom_2d = grid2d.get_geom_cells();
    size_t num_of_cells_2d = cells_geom_2d->size();
    std::vector<CellGeom3d>* cells_geom_3d = grid3d.get_geom_cells();
    size_t num_of_cells_3d = cells_geom_3d->size();
    int use_phreeqcrm = settings.get_int("use_phreeqcrm");

    for (size_t i = 0; i < num_of_species; i++)
    {
        size_t index_2d = i * num_of_cells_2d + 0;

        for (size_t j = 0; j < num_of_cells_3d; j++)
        {
            double conc_new = 0.0;

            if (use_phreeqcrm == 1)
            {
                // Unit conversion [g l-1] to [g m-3] done below (1000.0 *). 
                conc_new = 1000.0 * phreeqcWrapper.get_bound_cond(i);
            }
            else
            {
                conc_new = cells_solute_2d->at(index_2d).getDeposWetRate();
            }

            conc_init.push_back(conc_new);
        }
    }

    grid3d.init_solute_cells(settings, materials_3d, bound_cond_3d, 
        init_cond_3d, soluteLib, species_names, conc_init);

	// Save atmospheric settings and load atmospheric forcing data.
    std::cout << "Saving settings and loading atmospheric forcing data:\n";
    atmosControl.set_weather_data_interval(settings.get_double("weather_data_interval"));
    atmosControl.set_precip_fact(settings.get_double("area_fact"));
    atmosControl.set_pet_fact(settings.get_double("pet_fact"));
    std::vector < std::vector<std::string> > atmos_forcing_data;
    atmos_forcing_data = fileIO.load_and_tokenize_file(
        settings.get_str("atmos_forcing_data_path"), ',');
    
    if (atmos_forcing_data.size() == 0)
    {
        std::cout << "-> Atmospheric forcing data file not found or empty\n";
        return 1;
    }

    atmosControl.set_atmos_forcing_data(atmos_forcing_data);
    atmosControl.print_atmos_forcing_data();

    // Load 3d grid map for tridiagonal solvers.
    std::cout << "Loading 3d grid map for tridiagonal solvers:\n";
    std::vector < std::vector<std::string> > gridMapStr;
    gridMapStr = fileIO.load_and_tokenize_file(
        settings.get_str("grid_map_file_path_in"), ',');
    
    if (gridMapStr.size() == 0)
    {
        std::cout << "-> Grid map file was not found or empty\n";
        return 1;
    }

    for (size_t j = 0; j < gridMapStr.size(); j++)
    {
        std::vector<int> columnIndices;

        for (size_t i = 0; i < gridMapStr.at(j).size(); i++)
        {
            columnIndices.push_back(atoi(gridMapStr.at(j).at(i).c_str()));
        }

        gridMap.push_back(columnIndices);
    }

    return 0;
}

int Framework::run()
{
    // Initialize local variables.
    int num_of_species = settings.get_int("num_of_species");
    int sim_time_max = settings.get_int("sim_time_max"); // should this be double?
    int use_phreeqcrm = settings.get_int("use_phreeqcrm");
    double area_fact = settings.get_double("area_fact");
    std::vector<CellWater2d> * cells_water_2d = grid2d.get_water_cells();
    std::vector<CellWater3d> * cells_water_3d = grid3d.get_water_cells();
    std::vector<CellHeat3d> * cells_heat_3d = grid3d.get_heat_cells();
    std::vector<CellSolute2d> * cells_solute_2d = grid2d.get_solute_cells();
    std::vector<CellSolute3d> * cells_solute_3d = grid3d.get_solute_cells();
	size_t cellsPerSolute2d = 0;

	if (num_of_species > 0)
    {
        cellsPerSolute2d = cells_solute_2d->size() / num_of_species;
    }

	size_t cellsPerSolute3d = 0;

	if (num_of_species > 0)
    {
        cellsPerSolute3d = cells_solute_3d->size() / num_of_species;
    }

    double sim_time = 0.0; // this was changed from int to double
	double precipVolCum = 0.0;
	double evap_vol_cum = 0.0;
    double transp_vol_cum = 0.0;
	double infWatCum = 0.0;
    double drainVolCum = 0.0;
    double drainVolCum5min = 0.0;
    double sinkVolCum = 0.0;
    double outfall_vol_cum = 0.0;
    double outfall_vol = 0.0;
    std::vector<double> infMassCum;
    infMassCum.assign(num_of_species, 0.0);
    std::vector<double> drainMassCum;
    drainMassCum.assign(num_of_species, 0.0);
    std::vector<double> drainMassCum5Min;
    drainMassCum5Min.assign(num_of_species, 0.0);
    std::vector<double> decayMassOutCum;
    decayMassOutCum.assign(num_of_species, 0.0);
    std::vector<double> decayMassInCum;
    decayMassInCum.assign(num_of_species, 0.0);
    double timePrint = 0;
    int progPerc = 0;
    int gridOutCount = 1;
    std::string fileName2d = "mesh2d_out_";
    std::string fileName3d = "mesh3d_out_";
    std::string ext = ".vtk";
    std::vector<std::string> header;
    header.push_back("time [h]");
    header.push_back("time step [s]");
    header.push_back("cell id [-]");
    header.push_back("material [-]");
    header.push_back("hydraulic head [m]");
    header.push_back("pressure head [m]");
    header.push_back("water content [m3/m3]");
    header.push_back("network water volume [m3]");
    header.push_back("upper storage water volume [m3]");
    header.push_back("surface water volume [m3]");
    header.push_back("subsurface water volume [m3]");
    header.push_back("precipitation cum. [m3]");
    header.push_back("Outfall cum. [m3]");
    header.push_back("Outfall [l/s]");
    header.push_back("evaporation cum. [m3]");
    header.push_back("transpiration cum. [m3]");
    header.push_back("water infiltration cum. [m3]");
    header.push_back("drain discharge cum. [m3]");
    header.push_back("drain discharge cum. [l/5 min]");
    header.push_back("discharge to surface sinks [m3]");

    for (size_t i = 0; i < num_of_species; i++)
	{
        std::string solName = species_names.at(i);
        header.push_back("surface solute mass [g] (" + solName + ")");
        header.push_back("subsurface solute mass [g] (" + solName + ")");
        header.push_back("solute infiltration cum. [g] (" + solName + ")");
        header.push_back("drain mass cum. [g] (" + solName + ")");
        header.push_back("drain mass cum. [g/5 min] (" + solName + ")");
        header.push_back("decay mass in cum. [g] (" + solName + ")");
        header.push_back("decay mass out cum. [g] (" + solName + ")");
	}

    results.push_back(header);

    // Run the simulation.
    while (sim_time < sim_time_max)
    {
        // Initialize local variables.
        double time_step = settings.get_double("time_step");
        
        // Get atmospheric forcing data.
        double precip = atmosControl.get_precip();
        double pet = atmosControl.get_pet();
        double airTemp = atmosControl.get_air_temp();

        // Add precipitation to the overland domain.
        for (size_t i = 0; i < cells_water_2d->size(); i++)
        {
            CellGeom2d* geom2d = cells_water_2d->at(i).getGeom();

            if (geom2d->getMaterial() == 0)
            {
                continue;
            }

            double precip_depth = precip * time_step;
            precipVolCum += precip_depth * geom2d->getArea();
            double up_stor_thresh = cells_water_2d->at(i).get_up_stor_thresh();
            double up_stor_depth = cells_water_2d->at(i).get_up_stor_depth();
            up_stor_depth += precip_depth;

            if (up_stor_depth >= up_stor_thresh) {
                precip_depth = up_stor_depth - up_stor_thresh;
                up_stor_depth = up_stor_thresh;
            }
            else if (up_stor_depth < up_stor_thresh) {
                precip_depth = 0.0;
            }

            cells_water_2d->at(i).set_up_stor_depth(up_stor_depth);
            double surf_wat_depth = cells_water_2d->at(i).getWaterDepth();
            cells_water_2d->at(i).setWaterDepth(surf_wat_depth + precip_depth);
            cells_water_2d->at(i).swap();
        }

        // Remove water from overland domain that discharges into sinks.
        for (size_t i = 0; i < cells_water_2d->size(); i++)
        {
            cells_water_2d->at(i).removeWatBySinks(time_step);
        }
        
        // Remove water due to evaporation from surfaces and transpiration from soil due to plants.
        for (size_t i = 0; i < cells_water_3d->size(); i++)
        {
            CellGeom3d* geom3d = cells_water_3d->at(i).getGeom();
            int cell_id = geom3d->getGridConnection();
            
            if (geom3d->getMaterial() == 0)
            {
                continue;
            }

            if (cell_id >= 0 && cell_id < cells_water_2d->size())
            {
                // Extract evaporation from the surface domain.
                CellGeom2d* geom2d = cells_water_2d->at(cell_id).getGeom();
                double cell_area = geom2d->getArea();
                double up_stor_depth = cells_water_2d->at(cell_id).get_up_stor_depth();
                double evap_frac = cells_water_2d->at(cell_id).get_evap_frac();
                double evap_pot = evap_frac * pet * time_step;
                
                if (up_stor_depth >= evap_pot) {
                    up_stor_depth -= evap_pot;
                    cells_water_2d->at(cell_id).set_up_stor_depth(up_stor_depth);
                    evap_vol_cum += evap_pot * cell_area;
                    evap_pot = 0.0;
                }
                else if (up_stor_depth < evap_pot) {
                    evap_pot -= up_stor_depth;
                    evap_vol_cum += up_stor_depth * cell_area;
                    up_stor_depth = 0.0;
                    cells_water_2d->at(cell_id).set_up_stor_depth(up_stor_depth);
                }
                if (cells_water_2d->at(cell_id).getWaterDepth() >= evap_pot) {
                    cells_water_2d->at(cell_id).setWaterDepth(cells_water_2d->at(cell_id).getWaterDepth() - evap_pot);
                    evap_vol_cum += evap_pot * cell_area;
                    evap_pot = 0.0;
                }
                else if (cells_water_2d->at(cell_id).getWaterDepth() < evap_pot) {
                    evap_pot -= cells_water_2d->at(cell_id).getWaterDepth();
                    evap_vol_cum += cells_water_2d->at(cell_id).getWaterDepth() * cell_area;
                    cells_water_2d->at(cell_id).setWaterDepth(0.0);
                }
                
                // Extract transpiration from the soil.
                double crop_factor = cells_water_2d->at(cell_id).get_crop_factor();
                double transp_pot = ((1.0 - evap_frac) * pet * time_step) * crop_factor + evap_pot;
                double root_depth = cells_water_2d->at(cell_id).get_root_depth();
                
                if (crop_factor > 0.0 && root_depth > 0.0 && geom2d != 0) {
                    std::vector<Vertex> faceCentrePoints = geom3d->getFaceCentrePoints();
                    Vertex cell_top_vert = faceCentrePoints.at(0);
                    CellWater3d * cell_water_3d = & cells_water_3d->at(i);
                    double depth_sum = 0.0;
                
                    while (cell_water_3d != 0) {
                        geom3d = cell_water_3d->getGeom();
                        size_t num_of_neigh = geom3d->getNumOfNeigh();
                        faceCentrePoints = geom3d->getFaceCentrePoints();
                        Vertex cell_bott_vert = faceCentrePoints.at(num_of_neigh-1);
                        double cell_root_length = 0.0;
                    
                        if (depth_sum + cell_top_vert.z - cell_bott_vert.z > root_depth) {
                            cell_root_length = root_depth - depth_sum;
                        }
                        else {
                            cell_root_length = cell_top_vert.z - cell_bott_vert.z;
                        }

                        double press_head_eff = cell_water_3d->getPresHead();
                        double press_head_max = -0.1;
                        double press_head_min = -5.0;
                        double press_head_wilt = -150.0;
                        double moist_factor = 0.0;

                        // Fix out of bounds pressure head.
                        if (press_head_eff > 0.0)
                        {
                            press_head_eff = 0.0;
                        }
                        if (press_head_eff < press_head_wilt)
                        {
                            press_head_eff = press_head_wilt;
                        }

                        // Soil wetness is in optimal range.
                        if (press_head_eff >= press_head_min && press_head_eff <= press_head_max)
                        {
                            moist_factor = 1.0;
                        }
                        // Soil is too wet.
                        else if (press_head_eff > press_head_max)
                        {
                            moist_factor = 1.0 - (press_head_eff - press_head_max) / (0.0 - press_head_max);
                        }
                        // Soil is too dry.
                        else
                        {
                            moist_factor = 1.0 - (press_head_eff - press_head_min) / (press_head_wilt - press_head_min);
                        }

                        double cell_transp_vol = transp_pot * cell_root_length / root_depth * moist_factor * cell_area;
                        double wat_vol_old = geom3d->getVolume() * cell_water_3d->getWatCont();
                        double wat_cont_new = (wat_vol_old - cell_transp_vol) / geom3d->getVolume();
                        cell_water_3d->changeWatCont(wat_cont_new);
                        cell_water_3d->swap();
                        transp_vol_cum += cell_transp_vol;

                        // Continue to the next cell.
                        depth_sum += cell_top_vert.z - cell_bott_vert.z;
                        cell_top_vert = cell_bott_vert;
                        cell_water_3d = cell_water_3d->getNeigh(num_of_neigh - 1);
                    }
                }
            }
        }
        
        // Add solute deposition on the surface.
        for (size_t i = 0; i < num_of_species; i++)
        {
            for (size_t j = 0; j < cellsPerSolute2d; j++)
            {
                size_t ind = i * cellsPerSolute2d + j;
                double massExisting = cells_solute_2d->at(ind).getMass();
                CellGeom2d * geom2d =  cells_solute_2d->at(ind).getGeom();
                double massInPrecip = 0.0;

                if (use_phreeqcrm == 1)
                {
                    // Unit conversion [g l-1] to [g m-3] (1000.0 * ).
                    massInPrecip = 1000.0 * precip * geom2d->getArea() * 
                        phreeqcWrapper.get_bound_cond(i);
                }
                else
                {
                    massInPrecip = precip * geom2d->getArea() * 
                        cells_solute_2d->at(ind).getDeposWetRate();
                }

                cells_solute_2d->at(ind).setMass(massExisting + time_step *
                    cells_solute_2d->at(ind).getDeposDryRate() * 
                    area_fact * geom2d->getArea() + massInPrecip);
            }
        }
        
        // Run brute force network flow model with 
        // iterative diffusion wave simplification.
        if (settings.get_int("wat_flow_net_solver") == 1)
        {
            modelWaterNetDiffBrute.configure(
                settings.get_int("iter_stop_wat_net"),
                settings.get_double("iter_cut_thresh_wat_net"),
                settings.get_double("implicity_wat_net"),
                time_step,
                settings.get_double("bis_iter_cut_thresh_wat_net"),
                settings.get_int("bis_iter_stop_wat_net"),
                settings.get_double("bis_iter_cut_left_wat_net"),
                settings.get_double("bis_iter_cut_right_wat_net"));
            modelWaterNetDiffBrute.run(grid2d, network);
        }

        // Run explicit network flow model.
        else if (settings.get_int("wat_flow_net_solver") == 2)
        {
            modelWaterNetExplicit.configure(
                settings.get_int("iter_stop_wat_net"),
                settings.get_double("iter_cut_thresh_wat_net"),
                settings.get_double("implicity_wat_net"),
                time_step,
                settings.get_double("bis_iter_cut_thresh_wat_net"),
                settings.get_int("bis_iter_stop_wat_net"),
                settings.get_double("bis_iter_cut_left_wat_net"),
                settings.get_double("bis_iter_cut_right_wat_net"));
            modelWaterNetExplicit.run(grid2d, network);
        }
        
        // Run brute force 2d overland flow model with 
        // diffusion wave simplification.
        if (settings.get_int("wat_flow_2d_solver") == 1)
        {
            modelWater2dDiffBrute.configure(
                settings.get_int("iter_stop_wat_2d"),
                settings.get_double("iter_cut_thresh_wat_2d"),
                settings.get_double("implicity_wat_2d"),
                time_step, 
                settings.get_double("bis_iter_cut_thresh_wat_2d"),
                settings.get_int("bis_iter_stop_wat_2d"),
                settings.get_double("bis_iter_cut_left_wat_2d"),
                settings.get_double("bis_iter_cut_right_wat_2d"));
            modelWater2dDiffBrute.run(grid2d, grid3d);
        }

        // Run brute force solver subsurface flow model.
        if (settings.get_int("wat_flow_solver") == 1)
        {
            modelWater3dBrute.configure(
                settings.get_int("iter_stop_wat"),
                settings.get_double("iter_cut_thresh_wat"),
                settings.get_double("implicity_wat"),
                time_step);
            modelWater3dBrute.run(grid2d, grid3d);
        }
        // Run tridiagonal matrix algorithm based subsurface flow model.
        else if (settings.get_int("wat_flow_solver") == 2)
        {
            modelWater3dTri.configure(
                settings.get_int("iter_stop_wat"),
                settings.get_double("iter_cut_thresh_wat"),
                settings.get_double("implicity_wat"),
                time_step, 
                gridMap);
            modelWater3dTri.run(grid2d, grid3d);
        }

        // Run brute force based subsurface heat transport model.
        if (settings.get_int("heat_trans_solver") == 1)
        {
            modelHeat3dBrute.configure(
                settings.get_int("iter_stop_heat"), 
                settings.get_double("iter_cut_thresh_heat"), 
                settings.get_double("implicity_heat"), 
                time_step, 
                airTemp,
                settings.get_double("bis_iter_cut_thresh_heat"), 
                settings.get_int("bis_iter_stop_heat"), 
                settings.get_double("bis_iter_cut_left_heat"), 
                settings.get_double("bis_iter_cut_right_heat"));
            modelHeat3dBrute.run(grid2d, grid3d);
        }
        // Run tridiagonal matrix algorithm based heat transport model.
        else if (settings.get_int("heat_trans_solver") == 2)
        {
            modelHeat3dTri.configure(
                settings.get_int("iter_stop_heat"), 
                settings.get_double("iter_cut_thresh_heat"), 
                settings.get_double("implicity_heat"), 
                time_step, 
                airTemp,
                settings.get_double("bis_iter_cut_thresh_heat"), 
                settings.get_int("bis_iter_stop_heat"), 
                settings.get_double("bis_iter_cut_left_heat"), 
                settings.get_double("bis_iter_cut_right_heat"), 
                gridMap);
            modelHeat3dTri.run(grid2d, grid3d);
        }

        if (num_of_species > 0)
        {
            // Run brute force solver subsurface solute transport model.
            if (settings.get_int("sol_trans_solver") == 1)
            {
                modelSolute3dBrute.configure(
                    settings.get_int("iter_stop_sol"),
                    settings.get_double("iter_cut_thresh_sol"),
                    settings.get_double("implicity_sol"),
                    time_step,
                    num_of_species);
                modelSolute3dBrute.run(grid2d, grid3d);
            }
            // Run tridiagonal matrix algorithm based solute transport model.
            else if (settings.get_int("sol_trans_solver") == 2)
            {
                modelSolute3dTri.configure(
                    settings.get_int("iter_stop_sol"),
                    settings.get_double("iter_cut_thresh_sol"),
                    settings.get_double("implicity_sol"),
                    time_step, 
                    num_of_species, 
                    gridMap);
                modelSolute3dTri.run(grid2d, grid3d);
            }
        }

        // Advance time.
        atmosControl.advance_time(time_step);
        timePrint += time_step;
        sim_time += time_step;
        settings.set_value("sim_time", sim_time);

        // Compute new progress value.
        int progPercNew = (int)((double)sim_time / 
            settings.get_double("sim_time_max") * 100.0 + 0.5);

        // Run PHREEQCRM model.
        if (use_phreeqcrm == 1)
        {
            phreeqcWrapper.run(settings, grid2d, grid3d);
        }

        // Update progress.
        if (progPercNew > progPerc)
        {
            std::cout << "Completed: " << progPercNew << "%\n";
            progPerc = progPercNew;
        }

        // Process 2d grid water balance.
        for (size_t i = 0; i < cells_water_2d->size(); i++)
        {
            sinkVolCum += cells_water_2d->at(i).getWatVolRemovedBySinks();
        }

        // Process 3d grid water balance.
        for (size_t i = 0; i < cells_water_3d->size(); i++)
        {
            drainVolCum += cells_water_3d->at(i).getDrainVol();
            drainVolCum5min += 1000.0 * cells_water_3d->at(i).getDrainVol();
            infWatCum += cells_water_3d->at(i).getInfVol();
            //et_vol_cum += cells_water_3d->at(i).getEvapVol();
        }

        // Process 3d solute mass balance.
        for (size_t i = 0; i < num_of_species; i++)
        {
            for (size_t j = 0; j < cellsPerSolute3d; j++)
            {
                size_t ind = i * cellsPerSolute3d + j;
                infMassCum.at(i) += cells_solute_3d->at(ind).
                    getInfMass();
                drainMassCum5Min.at(i) += cells_solute_3d->at(ind).
                    getDrainMass();
                drainMassCum.at(i) += cells_solute_3d->at(ind).
                    getDrainMass();
                decayMassOutCum.at(i) += cells_solute_3d->at(ind).
                    getDecayMassOut(0);
                decayMassOutCum.at(i) += cells_solute_3d->at(ind).
                    getDecayMassOut(1);
                decayMassInCum.at(i) += cells_solute_3d->at(ind).
                    getDecayMassIn(0);
                decayMassInCum.at(i) += cells_solute_3d->at(ind).
                    getDecayMassIn(1);
            }
        }

        // Save results.
        if (timePrint >= settings.get_double("time_print_thresh"))
        {
            // Compute water volume in the network.
            double wat_vol_net = 0.0;
            std::vector<JuncWater>* water_juncs = network.get_water_juncs();

            for (size_t i = 0; i < water_juncs->size(); i++)
            {
                JuncGeom* geom_junc = water_juncs->at(i).get_geom();
                wat_vol_net += water_juncs->at(i).get_water_depth() * geom_junc->get_area();
            }

            // outfall_vol_cum
            for (size_t i = 0; i < water_juncs->size(); i++)
            {
                outfall_vol_cum += water_juncs->at(i).get_outfall_volume();
                outfall_vol += 1000.0 * water_juncs->at(i).get_outfall_volume() / settings.get_double("time_print_thresh");
                water_juncs->at(i).set_outfall_volume(0.0);
            }
 
            // Compute water volume in the 2d grid upper storage.
            double waterVolume2dUpper = 0.0;

            for (size_t i = 0; i < cells_water_2d->size(); i++)
            {
                CellGeom2d* geom = cells_water_2d->at(i).getGeom();

                if (geom->getMaterial() != 0) {
                    waterVolume2dUpper += cells_water_2d->at(i).get_up_stor_depth() *
                        geom->getArea();
                }
            }

            // Compute water volume in the 2d grid.
            double waterVolume2d = 0.0;

            for (size_t i = 0; i < cells_water_2d->size(); i++)
            {
                CellGeom2d * geom =  cells_water_2d->at(i).getGeom();

                if (geom->getMaterial() != 0) {
                    waterVolume2d += cells_water_2d->at(i).getWaterDepth() *
                        geom->getArea();
                }
            }

            // Compute water volume in the 3d grid.
            double waterVolume3d = 0.0;

            for (size_t i = 0; i < cells_water_3d->size(); i++)
            {
                CellGeom3d * geom =  cells_water_3d->at(i).getGeom();

                if (geom->getMaterial() != 0) {
                    waterVolume3d += cells_water_3d->at(i).getWatCont() *
                        geom->getVolume();
                }
            }

            // Compute solute mass in the 2d grid.
            std::vector<double> soluteMass2d;
            soluteMass2d.assign(num_of_species, 0.0);

            for (size_t i = 0; i < num_of_species; i++)
            {
                for (size_t j = 0; j < cellsPerSolute2d; j++)
                {
                    size_t index = i * cellsPerSolute2d + j;
                    soluteMass2d.at(i) += cells_solute_2d->at(index).getMass();
                }
            }

            // Compute solute mass in the 3d grid.
            std::vector<double> soluteMass3d;
            soluteMass3d.assign(num_of_species, 0.0);

            for (size_t i = 0; i < num_of_species; i++)
            {
                for (size_t j = 0; j < cellsPerSolute3d; j++)
                {
                    size_t index = i * cellsPerSolute3d + j;
                    soluteMass3d.at(i) += cells_solute_3d->at(index).getMass();
                }
            }

            // Save a results row.
            std::vector<std::string> resultRow;
            resultRow.push_back(std::to_string(sim_time / 3600.0));
            resultRow.push_back(std::to_string(time_step));
            resultRow.push_back(std::to_string(cells_water_3d->at(
                settings.get_int("cell_index")).getId()));
            resultRow.push_back(std::to_string(cells_water_3d->at(
                settings.get_int("cell_index")).getMaterial()));
            resultRow.push_back(std::to_string(cells_water_3d->at(
                settings.get_int("cell_index")).getHydrHead()));
            resultRow.push_back(std::to_string(cells_water_3d->at(
                settings.get_int("cell_index")).getPresHead()));
            resultRow.push_back(std::to_string(cells_water_3d->at(
                settings.get_int("cell_index")).getWatCont()));
            resultRow.push_back(std::to_string(wat_vol_net));
            resultRow.push_back(std::to_string(waterVolume2dUpper));
            resultRow.push_back(std::to_string(waterVolume2d));
            resultRow.push_back(std::to_string(waterVolume3d));
            resultRow.push_back(std::to_string(precipVolCum));
            resultRow.push_back(std::to_string(outfall_vol_cum));
            resultRow.push_back(std::to_string(outfall_vol));
            outfall_vol = 0.0;
            resultRow.push_back(std::to_string(evap_vol_cum));
            resultRow.push_back(std::to_string(transp_vol_cum));
            resultRow.push_back(std::to_string(infWatCum));
            resultRow.push_back(std::to_string(drainVolCum));
            resultRow.push_back(std::to_string(drainVolCum5min));
            resultRow.push_back(std::to_string(sinkVolCum));

            // TEMPORARILY ADD WATER DEPTHS INTO THE RESULTS.
            for (size_t i = 0; i < water_juncs->size(); i++)
            {
                resultRow.push_back(std::to_string(water_juncs->at(i).get_water_depth()));
            }

            for (size_t i = 0; i < num_of_species; i++)
            {
                resultRow.push_back(std::to_string(soluteMass2d.at(i)));
                resultRow.push_back(std::to_string(soluteMass3d.at(i)));
                resultRow.push_back(std::to_string(infMassCum.at(i)));
                resultRow.push_back(std::to_string(drainMassCum.at(i)));
                resultRow.push_back(std::to_string(drainMassCum5Min.at(i)));
                resultRow.push_back(std::to_string(decayMassInCum.at(i)));
                resultRow.push_back(std::to_string(decayMassOutCum.at(i)));
            }

            results.push_back(resultRow);
            drainVolCum5min = 0.0;

            for (size_t i = 0; i < num_of_species; i++)
            {
                drainMassCum5Min.at(i) = 0.0;
            }

            // Save vtk output grids.
            if (sim_time >= settings.get_int("time_print_vtk_start") && 
                sim_time <= settings.get_int("time_print_vtk_end") && 
                gridOutCount < 1000) // MAGIC number here
            {
                std::stringstream ss;
                ss << std::setw(3) << std::setfill('0') << gridOutCount;
                std::string outNum = ss.str();

                // Write the 2d mesh into a vtk file.
                if (settings.get_str("output_surf_vtk_path") != "-")
                {
                    std::string mesh2dFilePathOut = settings.get_str(
                        "output_surf_vtk_path") + outNum + ext;
                    std::string strMesh2d = grid2d.parse_unstruct_vtk_mesh();
                    fileIO.save_grid(mesh2dFilePathOut, strMesh2d);
                }

                // Write the 3d mesh into a vtk file.
                if (settings.get_str("output_subs_vtk_path") != "-")
                {
                    std::string mesh3dFilePathOut = settings.get_str(
                        "output_subs_vtk_path") + outNum + ext;
                    std::string strMesh3d = grid3d.parse_unstruct_vtk_mesh();
                    fileIO.save_grid(mesh3dFilePathOut, strMesh3d);
                }
                
                gridOutCount += 1;
            }

            // Update counters.
            timePrint = 0;
        }
    }

    return 0;
}
void Framework::saveOutputVTK(std::string path2d, std::string path3d)
{
    // Write the 2d mesh into a vtk file.
    std::cout << "Writing the 2d mesh into a vtk file:\n";
    std::string strMesh2d = grid2d.parse_unstruct_vtk_mesh();
    fileIO.save_grid(path2d, strMesh2d);

    // Write the 3d mesh into a vtk file.
    std::cout << "Writing the 3d mesh into a vtk file:\n";
    std::string strMesh3d = grid3d.parse_unstruct_vtk_mesh();
    fileIO.save_grid(path3d, strMesh3d);
}

int Framework::finalize()
{
    // Save results to disk.
    std::cout << "Saving results to disk:\n";
    fileIO.save_results(settings.get_str("output_csv_path"), results);

    return 0;
}
