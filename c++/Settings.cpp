#include "Settings.h"

Settings::Settings()
{
    /*
    runTests = true;
    numOfThreads = 0;
    simTimeMax = 0;
    timeStep = 0;
    timePrintThresh = 0;
    timePrintVtkStart = 0;
    timePrintVtkEnd = 0;
    watFlow2dSolver = 0;
    iterStopWat2d = 0;
    iterCutThreshWat2d = 0.0;
    implicityWat2d = 0.0;
    bisIterCutThreshWat2d = 0.0;
    bisIterStopWat2d = 0;
    bisIterCutLeftWat2d = 0.0;
    bisIterCutRightWat2d = 0.0;
    watFlowSolver = 0;
    cellIndex = 0;
    iterStopWat = 0;
    iterCutThreshWat = 0.0;
    implicityWat = 0.0;
    areaFact = 0.0;
    petFact = 0.0;
    heatTransSolver = 0;
    iterStopHeat = 0;
    iterCutThreshHeat = 0.0;
    implicityHeat = 0.0;
    bisIterCutThreshHeat = 0.0;
    bisIterStopHeat = 0;
    bisIterCutLeftHeat = 0.0;
    bisIterCutRightHeat = 0.0;
    latHeatFus = 0.0;
    tempFreez = 0.0;
    multFreez = 0.0;
    densAir = 0.0;
    densIce = 0.0;
    densWat = 0.0;
    capAir = 0.0;
    capIce = 0.0;
    capWat = 0.0;
    condAir = 0.0;
    condIce = 0.0;
    condWat = 0.0;
    condMultAir = 0.0;
    condMultIce = 0.0;
    condMultWat = 0.0;
    tempSoilInit = 0.0;
    usePhreeqcrm = false;
    solTransSolver = 0;
    numOfSolutes = 0;
    iterStopSol = 0;
    iterCutThreshSol = 0.0;
    implicitySol = 0.0;
    molDiffRate = 0.0;
    mesh2dfilePathIn = "";
    mesh3dfilePathIn = "";
    gridMapfilePathIn = "";
    materials2dPath = "";
    materials3dPath = "";
    initCond2dPath = "";
    initCond3dPath = "";
    boundCond2dPath = "";
    boundCond3dPath = "";
    atmosForcingdataPath = "";
    solutePropertiesPath = "";
    phreeqcInputPath = "";
    phreeqcChemDbPath = "";
    phreeqcSpecMapPath = "";
    phreeqcMolWeightMapPath = "";
    outputCsvPath = "";
    outputSurfVtkPath = "";
    outputSubsVtkPath = "";
    */
}

Settings::~Settings()
{
}

void Settings::initialize(std::vector < std::vector<std::string> > settings)
{
    // Improve this by getting map keys directly from the input file.
    settings_str["run_tests"] =                   settings.at(1).at(1);
    settings_str["num_of_threads"] =              settings.at(2).at(1);
    settings_str["sim_time_max"] =                settings.at(3).at(1);
    settings_str["time_step"] =                   settings.at(4).at(1);
    settings_str["cell_index"] =                  settings.at(5).at(1);
    settings_str["time_print_thresh"] =           settings.at(6).at(1);
    settings_str["time_print_vtk_start"] =        settings.at(7).at(1);
    settings_str["time_print_vtk_end"] =          settings.at(8).at(1);
    settings_str["wat_flow_2d_solver"] =          settings.at(9).at(1);
    settings_str["iter_stop_wat_2d"] =            settings.at(10).at(1);
    settings_str["iter_cut_thresh_wat_2d"] =      settings.at(11).at(1);
    settings_str["implicity_wat_2d"] =            settings.at(12).at(1);
    settings_str["bis_iter_cut_thresh_wat_2d"] =  settings.at(13).at(1);
    settings_str["bis_iter_stop_wat_2d"] =        settings.at(14).at(1);
    settings_str["bis_iter_cut_left_wat_2d"] =    settings.at(15).at(1);
    settings_str["bis_iter_cut_right_wat_2d"] =   settings.at(16).at(1);
    settings_str["wat_flow_solver"] =             settings.at(17).at(1);
    settings_str["iter_stop_wat"] =               settings.at(18).at(1);
    settings_str["iter_cut_thresh_wat"] =         settings.at(19).at(1);
    settings_str["implicity_wat"] =               settings.at(20).at(1);
    settings_str["area_fact"] =                   settings.at(21).at(1);
    settings_str["pet_fact"] =                    settings.at(22).at(1);
    settings_str["heat_trans_solver"] =           settings.at(23).at(1);
    settings_str["iter_stop_heat"] =              settings.at(24).at(1);
    settings_str["iter_cut_thresh_heat"] =        settings.at(25).at(1);
    settings_str["implicity_heat"] =              settings.at(26).at(1);
    settings_str["bis_iter_cut_thresh_heat"] =    settings.at(27).at(1);
    settings_str["bis_iter_stop_heat"] =          settings.at(28).at(1);
    settings_str["bis_iter_cut_left_heat"] =      settings.at(29).at(1);
    settings_str["bis_iter_cut_right_heat"] =     settings.at(30).at(1);
    settings_str["lat_heat_fus"] =                settings.at(31).at(1);
    settings_str["temp_freez"] =                  settings.at(32).at(1);
    settings_str["mult_freez"] =                  settings.at(33).at(1);
    settings_str["dens_air"] =                    settings.at(34).at(1);
    settings_str["dens_ice"] =                    settings.at(35).at(1);
    settings_str["dens_wat"] =                    settings.at(36).at(1);
    settings_str["cap_air"] =                     settings.at(37).at(1);
    settings_str["cap_ice"] =                     settings.at(38).at(1);
    settings_str["cap_wat"] =                     settings.at(39).at(1);
    settings_str["cond_air"] =                    settings.at(40).at(1);
    settings_str["cond_ice"] =                    settings.at(41).at(1);
    settings_str["cond_wat"] =                    settings.at(42).at(1);
    settings_str["cond_mult_air"] =               settings.at(43).at(1);
    settings_str["cond_mult_ice"] =               settings.at(44).at(1);
    settings_str["cond_mult_wat"] =               settings.at(45).at(1);
    settings_str["temp_soil_init"] =              settings.at(46).at(1);
    settings_str["use_phreeqcrm"] =               settings.at(47).at(1);
    settings_str["sol_trans_solver"] =            settings.at(48).at(1);
    settings_str["iter_stop_sol"] =               settings.at(49).at(1);
    settings_str["iter_cut_thresh_sol"] =         settings.at(50).at(1);
    settings_str["implicity_sol"] =               settings.at(51).at(1);
    settings_str["mol_diff_rate"] =               settings.at(52).at(1);
    settings_str["network_link_file_path_in"] =   settings.at(53).at(1);
    settings_str["network_junc_file_path_in"] =   settings.at(54).at(1);
    settings_str["mesh_2d_file_path_in"] =        settings.at(55).at(1);
    settings_str["mesh_3d_file_path_in"] =        settings.at(56).at(1);
    settings_str["grid_map_file_path_in"] =       settings.at(57).at(1);
    settings_str["materials_2d_path"] =           settings.at(58).at(1);
    settings_str["materials_3d_path"] =           settings.at(59).at(1);
    settings_str["init_cond_2d_path"] =           settings.at(60).at(1);
    settings_str["init_cond_3d_path"] =           settings.at(61).at(1);
    settings_str["bound_cond_2d_path"] =          settings.at(62).at(1);
    settings_str["bound_cond_3d_path"] =          settings.at(63).at(1);
    settings_str["atmos_forcing_data_path"] =     settings.at(64).at(1);
    settings_str["solute_properties_path"] =      settings.at(65).at(1);
    settings_str["phreeqc_input_path"] =          settings.at(66).at(1);
    settings_str["phreeqc_chem_db_path"] =        settings.at(67).at(1);
    settings_str["phreeqc_spec_map_path"] =       settings.at(68).at(1);
    settings_str["phreeqc_mol_weight_map_path"] = settings.at(69).at(1);
    settings_str["output_csv_path"] =             settings.at(70).at(1);
    settings_str["output_surf_vtk_path"] =        settings.at(71).at(1);
    settings_str["output_subs_vtk_Path"] =        settings.at(72).at(1);

    /*
    std::string runTestsStr = settings.at(1).at(1);
    runTests = (runTestsStr == "1");
    numOfThreads = atoi(settings.at(2).at(1).c_str());
    simTimeMax = atoi(settings.at(3).at(1).c_str());
    timeStep = atoi(settings.at(4).at(1).c_str());
    cellIndex = atoi(settings.at(5).at(1).c_str());
    timePrintThresh = atoi(settings.at(6).at(1).c_str());
    timePrintVtkStart = atoi(settings.at(7).at(1).c_str());
    timePrintVtkEnd = atoi(settings.at(8).at(1).c_str());
    watFlow2dSolver = atoi(settings.at(9).at(1).c_str());
    iterStopWat2d = atoi(settings.at(10).at(1).c_str());
    iterCutThreshWat2d = atof(settings.at(11).at(1).c_str());
    implicityWat2d = atof(settings.at(12).at(1).c_str());
    bisIterCutThreshWat2d = atof(settings.at(13).at(1).c_str());
    bisIterStopWat2d = atoi(settings.at(14).at(1).c_str());
    bisIterCutLeftWat2d = atof(settings.at(15).at(1).c_str());
    bisIterCutRightWat2d = atof(settings.at(16).at(1).c_str());
    watFlowSolver = atoi(settings.at(17).at(1).c_str());
    iterStopWat = atoi(settings.at(18).at(1).c_str());
    iterCutThreshWat = atof(settings.at(19).at(1).c_str());
    implicityWat = atof(settings.at(20).at(1).c_str());
    areaFact = atof(settings.at(21).at(1).c_str());
    petFact = atof(settings.at(22).at(1).c_str());
    heatTransSolver = atoi(settings.at(23).at(1).c_str());
    iterStopHeat = atoi(settings.at(24).at(1).c_str());
    iterCutThreshHeat = atof(settings.at(25).at(1).c_str());
    implicityHeat = atof(settings.at(26).at(1).c_str());
    bisIterCutThreshHeat = atof(settings.at(27).at(1).c_str());
    bisIterStopHeat = atoi(settings.at(28).at(1).c_str());
    bisIterCutLeftHeat = atof(settings.at(29).at(1).c_str());
    bisIterCutRightHeat = atof(settings.at(30).at(1).c_str());
    latHeatFus = atof(settings.at(31).at(1).c_str());
    tempFreez = atof(settings.at(32).at(1).c_str());
    multFreez = atof(settings.at(33).at(1).c_str());
    densAir = atof(settings.at(34).at(1).c_str());
    densIce = atof(settings.at(35).at(1).c_str());
    densWat = atof(settings.at(36).at(1).c_str());
    capAir = atof(settings.at(37).at(1).c_str());
    capIce = atof(settings.at(38).at(1).c_str());
    capWat = atof(settings.at(39).at(1).c_str());
    condAir = atof(settings.at(40).at(1).c_str());
    condIce = atof(settings.at(41).at(1).c_str());
    condWat = atof(settings.at(42).at(1).c_str());
    condMultAir = atof(settings.at(43).at(1).c_str());
    condMultIce = atof(settings.at(44).at(1).c_str());
    condMultWat = atof(settings.at(45).at(1).c_str());
    tempSoilInit = atof(settings.at(46).at(1).c_str());
    std::string usePhreeqcrmStr = settings.at(47).at(1);
    usePhreeqcrm = (usePhreeqcrmStr == "1");
    solTransSolver = atoi(settings.at(48).at(1).c_str());
    iterStopSol = atoi(settings.at(49).at(1).c_str());
    iterCutThreshSol = atof(settings.at(50).at(1).c_str());
    implicitySol = atof(settings.at(51).at(1).c_str());
    molDiffRate = atof(settings.at(52).at(1).c_str());
    mesh2dfilePathIn = settings.at(53).at(1);
    mesh3dfilePathIn = settings.at(54).at(1);
    gridMapfilePathIn = settings.at(55).at(1);
    materials2dPath = settings.at(56).at(1);
    materials3dPath = settings.at(57).at(1);
    initCond2dPath = settings.at(58).at(1);
    initCond3dPath = settings.at(59).at(1);
    boundCond2dPath = settings.at(60).at(1);
    boundCond3dPath = settings.at(61).at(1);
    atmosForcingdataPath = settings.at(62).at(1);
    solutePropertiesPath = settings.at(63).at(1);
    phreeqcInputPath = settings.at(64).at(1);
    phreeqcChemDbPath = settings.at(65).at(1);
    phreeqcSpecMapPath = settings.at(66).at(1);
    phreeqcMolWeightMapPath = settings.at(67).at(1);
    outputCsvPath = settings.at(68).at(1);
    outputSurfVtkPath = settings.at(69).at(1);
    outputSubsVtkPath = settings.at(70).at(1);
    */
}

std::string Settings::get_str(std::string param_name)
{
    if (settings_str.find(param_name) == settings_str.end())
    {
        std::cout << "Parameter " << param_name.c_str() << " not found.\n";
        return "";
    }
    else
    {
        return settings_str[param_name];
    }
}

int Settings::get_int(std::string param_name)
{
    //std::cout << "\n" << settings_str.at(param_name) << "\n";
    if (settings_str.find(param_name) == settings_str.end())
    {
        std::cout << "Parameter \"" << param_name.c_str() << "\" not found.\n";
        return -1;
    }
    else
    {
        return atoi(settings_str[param_name].c_str());
    }
}

double Settings::get_double(std::string param_name)
{
    if (settings_str.find(param_name) == settings_str.end())
    {
        std::cout << "Parameter " << param_name.c_str() << " not found.\n";
        return -1.0;
    }
    else
    {
        return atof(settings_str[param_name].c_str());
    }
}

void Settings::set_value(std::string name, std::string value)
{
    settings_str[name] = value;
}

void Settings::set_value(std::string name, int value)
{
    std::ostringstream sstream;
    sstream << value;
    settings_str[name] = sstream.str();
}

void Settings::set_value(std::string name, double value)
{
    std::ostringstream sstream;
    sstream << value;
    settings_str[name] = sstream.str();
}