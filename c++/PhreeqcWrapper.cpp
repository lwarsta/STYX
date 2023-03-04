#include "PhreeqcWrapper.h"

PhreeqcWrapper::PhreeqcWrapper()
{
    phreeqc_rm = 0;
    num_of_threads = 1;
    num_of_cells_3d = 0;
    num_of_species = 0;
    posWaterMol = 0;
}

PhreeqcWrapper::~PhreeqcWrapper()
{
    delete phreeqc_rm;
}

int PhreeqcWrapper::countWordOccurrences(std::string input, std::string search_word)
{
    // Find a better alternative for checking the number of solutions.
    std::vector<std::string> tokens;
    std::string token;
    for (size_t i = 0; i < input.length(); i++)
    {
        if (token.length() > 0 && (input.at(i) == ';' || input.at(i) == ' ' || input.at(i) == '\n'))
        {
            tokens.push_back(token);
            token = "";
        }
        else if (input.at(i) != ';' && input.at(i) != ' ' && input.at(i) != '\n' && input.at(i) != '\r')
        {
            token += input.at(i);
        }
    }
    int word_count = 0;
    for (size_t i = 0; i < tokens.size(); i++)
    {
        if (search_word == tokens.at(i))
        {
            word_count++;
        }
    }
    return word_count;
}

void PhreeqcWrapper::initialize(
    Settings& settings, 
    std::vector<std::vector<std::string>>& boundCond2d, 
    std::vector<std::vector<std::string>>& boundCond3d, 
    Grid2d& grid2d, 
    Grid3d& grid3d)
{
    try
    {
        // Initialise PHREEQCRM.
        std::cout << "\nInitializing PHREEQC:";
        std::vector<CellGeom3d>* geomCells3d = grid3d.get_geom_cells();
        num_of_cells_3d = (int)geomCells3d->size();
        num_of_threads = 1;
        phreeqc_rm = new PhreeqcRM(num_of_cells_3d, num_of_threads);
        IRM_RESULT status;

        // Create a map to convert species names used in the database used in the simulation.
        std::vector < std::vector<std::string> > specNameMapData = fileIO.load_and_tokenize_file(settings.get_str("phreeqc_spec_map_path"), ',');
        std::map<std::string, std::string> specNameMap;
        for (size_t i = 1; i < specNameMapData.size(); i++) // bypass header
        {
            specNameMap[specNameMapData.at(i).at(0)] = specNameMapData.at(i).at(1);
        }

        // Create a map of species names and molecular weights.
        std::vector < std::vector<std::string> > molWeightMapData = fileIO.load_and_tokenize_file(settings.get_str("phreeqc_mol_weight_map_path"), ',');
        for (size_t i = 1; i < molWeightMapData.size(); i++) // bypass header
        {
            molWeightMap[molWeightMapData.at(i).at(0)] = atof(molWeightMapData.at(i).at(1).c_str());
        }

        // Set PHREEQCRM properties.
        status = phreeqc_rm->SetErrorHandlerMode(1); // OR TRUE?
        status = phreeqc_rm->SetSpeciesSaveOn(true);
        status = phreeqc_rm->SetUnitsSolution(1);
        status = phreeqc_rm->SetUnitsPPassemblage(1);
        status = phreeqc_rm->SetUnitsExchange(1);
        status = phreeqc_rm->SetUnitsSurface(1);
        status = phreeqc_rm->SetUnitsGasPhase(1);
        status = phreeqc_rm->SetUnitsSSassemblage(1);
        status = phreeqc_rm->SetUnitsKinetics(1);
        status = phreeqc_rm->SetSelectedOutputOn(true);

        // Set initial cell volume, porosity, saturation, density, temperature and pressure.
        phreeqcVol.assign(num_of_cells_3d, 0.0);
        phreeqcPor.assign(num_of_cells_3d, 0.0);
        phreeqcSat.assign(num_of_cells_3d, 0.0);
        phreeqcDens.assign(num_of_cells_3d, 0.0);
        phreeqcTemp.assign(num_of_cells_3d, 20.0);
        phreeqcPres.assign(num_of_cells_3d, 2.0);
        std::vector<CellWater3d>* waterCells3d = grid3d.get_water_cells();
        for (int i = 0; i < waterCells3d->size(); i++)
        {
            CellGeom3d* cellGeom = waterCells3d->at(i).getGeom();
            phreeqcVol.at(i) = cellGeom->getVolume() * 1000.0; // m3 -> l TEMPORARILY COMMENTED OUT
            phreeqcPor.at(i) = waterCells3d->at(i).getWatContSat();
            // Saturation varies bewteen 0 ... 1.
            double saturation = 0.0;
            if (waterCells3d->at(i).getWatCont() > 0.0 && waterCells3d->at(i).getWatContSat() > 0.0)
            {
                saturation = waterCells3d->at(i).getWatCont() / waterCells3d->at(i).getWatContSat();
                if (saturation > 1.0)
                {
                    saturation = 1.0;
                }
            }
            phreeqcSat.at(i) = saturation;
        }
        status = phreeqc_rm->SetRepresentativeVolume(phreeqcVol);
        status = phreeqc_rm->SetPorosity(phreeqcPor);
        status = phreeqc_rm->SetSaturation(phreeqcSat);
        //status = phreeqc_rm->SetDensity(cellDens);
        //status = phreeqc_rm->SetTemperature(cellTemp);
        //status = phreeqc_rm->SetPressure(pressure);
        
        // Load chemistry database.
        std::string path_to_file = settings.get_str("phreeqc_chem_db_path");
        status = phreeqc_rm->LoadDatabase(path_to_file.c_str());

        // Define solutions and other simulation properties.
        path_to_file = settings.get_str("phreeqc_input_path");
        std::string chemInit = fileIO.load_file(path_to_file);
        std::cout << "\n" << chemInit.c_str();
        status = phreeqc_rm->RunString(true, true, true, chemInit.c_str());

        // Assess how many solutions the file contains.
        int num_of_solutions = countWordOccurrences(chemInit, "SOLUTION");
        std::cout << "\nNumber of solutions in the chemistry file: " << num_of_solutions;

        // Clear contents of workers and utility - is this needed?
        std::string input = "DELETE; -all";
        status = phreeqc_rm->RunString(true, false, true, input.c_str());

        // Determine component information.
        int ncomps = phreeqc_rm->FindComponents();
        const std::vector<std::string>& components = phreeqc_rm->GetComponents();
        const std::vector < double >& gfw = phreeqc_rm->GetGfw();

        // Run selected output rows separately. The above block must be run before this.
        std::string sel_output = "SELECTED_OUTPUT; -si ";
        const std::vector<std::string>& si = phreeqc_rm->GetSINames();
        for (size_t i = 0; i < si.size(); i++)
        {
            sel_output += si.at(i);
            if (i < si.size() - 1)
            {
                sel_output += " ";
            }
        }
        sel_output += "; END;";
        std::cout << "\nSelected output string: " << sel_output.c_str();
        status = phreeqc_rm->RunString(true, true, true, sel_output.c_str());

        // Determine species information.
        num_of_species = phreeqc_rm->GetSpeciesCount();
        species_names = phreeqc_rm->GetSpeciesNames();
        const std::vector < double >& species_z = phreeqc_rm->GetSpeciesZ();
        const std::vector < double >& species_d = phreeqc_rm->GetSpeciesD25();
        bool species_on = phreeqc_rm->GetSpeciesSaveOn();
        for (size_t i = 0; i < species_names.size(); i++)
        {
            std::cout << "\n" << i << ": " << species_names.at(i);
            if (species_names.at(i) == "H2O")
            {
                posWaterMol = i;
                std::cout << "\nPosition of water molecule in the list: " << posWaterMol;
            }
        }

        // Set array of initial conditions for cells.
        std::vector<int> ic1(num_of_cells_3d * 7, -1);
        for (size_t i = 0; i < num_of_cells_3d; i++)
        {
            // Get boundary condition indices.
            CellGeom3d* geomCell = waterCells3d->at(i).getGeom();
            int boundInd = geomCell->getBoundCondInd();
            // Set boundary condition indices.
            // Solution below was zero!
            ic1.at(0 * num_of_cells_3d + i) = atoi(boundCond3d.at(boundInd + 1).at(2 + 1).c_str()); // Solution
            ic1.at(1 * num_of_cells_3d + i) = atoi(boundCond3d.at(boundInd + 1).at(3 + 1).c_str()); // Equilibrium phases
            ic1.at(2 * num_of_cells_3d + i) = atoi(boundCond3d.at(boundInd + 1).at(4 + 1).c_str()); // Exchange
            ic1.at(3 * num_of_cells_3d + i) = atoi(boundCond3d.at(boundInd + 1).at(5 + 1).c_str()); // Surface - is separate surface needed for each cell?
            ic1.at(4 * num_of_cells_3d + i) = atoi(boundCond3d.at(boundInd + 1).at(6 + 1).c_str()); // Gas phase
            ic1.at(5 * num_of_cells_3d + i) = atoi(boundCond3d.at(boundInd + 1).at(7 + 1).c_str()); // Solid solutions
            ic1.at(6 * num_of_cells_3d + i) = atoi(boundCond3d.at(boundInd + 1).at(8 + 1).c_str()); // Kinetics
        }
        status = phreeqc_rm->InitialPhreeqc2Module(ic1);

        // Get boundary conditions.
        for (int i = 0; i < num_of_solutions; i++)
        {
            std::vector<double> conc;
            std::vector<int> bound_cond_id;
            bound_cond_id.resize(1, i);
            status = phreeqc_rm->InitialPhreeqc2SpeciesConcentrations(conc, bound_cond_id);
            bound_cond_conc.push_back(conc);
        }
        std::cout << "\nBoundary conditions: ";
        for (size_t i = 0; i < bound_cond_conc.size(); i++)
        {
            std::cout << "\n" << i;
            for (size_t j = 0; j < bound_cond_conc.at(i).size(); j++)
            {
                std::cout << "," << bound_cond_conc.at(i).at(j);
            }
        }

        // Initialise cells.
        status = phreeqc_rm->SetTimeStep(0.0);
        status = phreeqc_rm->SetTime(0.0);
        status = phreeqc_rm->RunCells();
        status = phreeqc_rm->GetSpeciesConcentrations(phreeqcConc); // vector is initialised by PHREEQCRM
        std::cout << "\nnum_of_cells: " << num_of_cells_3d;
        std::cout << "\nnum_of_species: " << num_of_species;
        for (size_t i = 0; i < num_of_species; i++)
        {
            size_t index = i * num_of_cells_3d;
            phreeqcConcInit.push_back(phreeqcConc.at(index));
            std::cout << "\n" << i << ": " << phreeqcConc.at(index);
        }
    }
    catch (PhreeqcRMStop)
    {
        std::string e_string = "PHREEQCRM initialization failed with an error in PhreeqcRM.";
        std::cerr << e_string << std::endl;
        //return IRM_FAIL;
    }
    catch (...)
    {
        std::string e_string = "PHREEQCRM initialization failed with an unhandled exception.";
        std::cerr << e_string << std::endl;
        //return IRM_FAIL;
    }
}

double PhreeqcWrapper::get_bound_cond(size_t index)
{ 
    // Unit conversion [mol l-1] * [g mol^-1] = [g l-1] done below. 
    if (bound_cond_conc.size() > 0 &&
        bound_cond_conc.at(bound_cond_conc.size() - 1).size() > index &&
        species_names.size() > index &&
        molWeightMap.find(species_names.at(index)) != molWeightMap.end())
    {
        return bound_cond_conc.at(bound_cond_conc.size() - 1).at(index) * molWeightMap[species_names.at(index)];
    }
    else
    {
        std::cout << "\nRequested index not available in the PHREEQC boundary condition data: " << index;
        return 0.0;
    }
}

int PhreeqcWrapper::run(Settings& settings, Grid2d& grid2d, Grid3d& grid3d)
{
    try
    {
        std::vector<double> concIn(num_of_species * num_of_cells_3d, 0.0); // cellsPerSolute3d
        std::vector<double> concOut(num_of_species * num_of_cells_3d, 0.0); // cellsPerSolute3d
        std::vector<CellWater3d>* waterCells3d = grid3d.get_water_cells();
        std::vector<double> watCont(waterCells3d->size(), 0.0);

        // Extract data from the water cells.
        for (size_t i = 0; i < waterCells3d->size(); i++)
        {
            watCont.at(i) = waterCells3d->at(i).getWatCont();
            double saturation = 0.0;
            if (waterCells3d->at(i).getWatCont() > 0.0 && waterCells3d->at(i).getWatContSat() > 0.0)
            {
                saturation = waterCells3d->at(i).getWatCont() / waterCells3d->at(i).getWatContSat();
                if (saturation > 1.0)
                {
                    saturation = 1.0;
                }
            }

            // Take compression into account.
            if (waterCells3d->at(i).getWatCont() > waterCells3d->at(i).getWatContSat())
            {
                phreeqcPor.at(i) = waterCells3d->at(i).getWatCont();
            }
            else
            {
                phreeqcPor.at(i) = waterCells3d->at(i).getWatContSat();
            }

            phreeqcSat.at(i) = saturation;
        }

        // Extract data from the heat cells.
        std::vector<CellHeat3d>* heatCells3d = grid3d.get_heat_cells();
        for (size_t i = 0; i < heatCells3d->size(); i++)
        {
            phreeqcTemp.at(i) = heatCells3d->at(i).getT();
        }

        // Extract data from the solute cells.
        std::vector<CellSolute3d>* soluteCells3d = grid3d.get_solute_cells();
        for (size_t i = 0; i < num_of_species; i++)
        {
            for (size_t j = 0; j < num_of_cells_3d; j++) //  cellsPerSolute3d
            {
                size_t index = i * num_of_cells_3d + j; // cellsPerSolute3d
                concIn.at(index) = soluteCells3d->at(index).getMass();
                soluteCells3d->at(index).convertMassToConc();
                phreeqcConc.at(index) = 0.001 * soluteCells3d->at(index).getConc() / molWeightMap[species_names.at(i)];
            }
        }

        // Save data to PHREEQCRM.
        IRM_RESULT status;
        status = phreeqc_rm->SetPorosity(phreeqcPor);
        status = phreeqc_rm->SetSaturation(phreeqcSat); // note that this varies between 0.0 ... 1.0
        //status = phreeqc_rm->SetDensity(phreeqcDens);
        status = phreeqc_rm->SetTemperature(phreeqcTemp);
        status = phreeqc_rm->SetPressure(phreeqcPres); // CURRENTLY NOT UPDATED!
        status = phreeqc_rm->SpeciesConcentrations2Module(phreeqcConc);

        // Conduct computations in PHREEQCRM.
        status = phreeqc_rm->SetTimeStep(settings.get_double("time_step"));
        status = phreeqc_rm->SetTime(settings.get_double("sim_time"));
        status = phreeqc_rm->RunCells();

        // Transfer data from PhreeqcRM to the transport model.
        status = phreeqc_rm->GetSpeciesConcentrations(phreeqcConc);
        //status = phreeqc_rm->GetDensity(phreeqcDens);
        const std::vector<double>& phreeqcSolVol = phreeqc_rm->GetSolutionVolume();

        // Save results to the cells.
        for (size_t i = 0; i < num_of_species; i++)
        {
            for (size_t j = 0; j < num_of_cells_3d; j++) // cellsPerSolute3d
            {
                size_t index = i * num_of_cells_3d + j; // cellsPerSolute3d
                CellWater3d* water3d = soluteCells3d->at(index).getWater();
                CellGeom3d* geom3d = soluteCells3d->at(index).getGeom();
                double solVolMult = phreeqcSolVol.at(j) / (1000 * water3d->getWatCont() * geom3d->getVolume()); // mol l-1 / (1000.0 l)
                double concNew = 1000.0 * solVolMult * phreeqcConc.at(index) * molWeightMap[species_names.at(i)]; // mol / l * g / mol = g / l
                soluteCells3d->at(index).setConc(concNew);
                soluteCells3d->at(index).convertConcToMass();
                concOut.at(index) = soluteCells3d->at(index).getMass();
            }
        }
    }
    catch (PhreeqcRMStop)
    {
        std::string e_string = "PHREEQCRM failed with an error in PhreeqcRM.";
        std::cerr << e_string << std::endl;

        // Print cell data.
        for (size_t i = 0; i < phreeqcVol.size(); i++)
        {
            std::cout << phreeqcVol.at(i) << ",";
            std::cout << phreeqcPor.at(i) << ",";
            std::cout << phreeqcSat.at(i) << ",";
            //std::cout << watCont.at(i) << ",";
            std::cout << phreeqcDens.at(i) << ",";
            std::cout << phreeqcTemp.at(i) << ",";
            std::cout << phreeqcPres.at(i);
            std::cout << "\n";
        }

        // Print concentration table.
        std::cout << "\n";
        /*
        for (size_t j = 0; j < cellsPerSolute3d; j++)
        {
            for (size_t i = 0; i < 2 * num_of_species; i++)
            {
                if (i < num_of_species)
                {
                    size_t index = i * cellsPerSolute3d + j;
                    //std::cout << concIn.at(index);
                }
                else if (i >= num_of_species)
                {
                    size_t index = (i - num_of_species) * cellsPerSolute3d + j;
                    //std::cout << concOut.at(index);
                }
                if (i < 2 * num_of_species - 1)
                {
                    std::cout << ",";
                }
            }
            std::cout << "\n";
        }
        */
        return 2; // IRM_FAIL;
    }
    catch (...)
    {
        std::string e_string = "PHREEQCRM failed with an unhandled exception.";
        std::cerr << e_string << std::endl;

        // Print cell data.
        for (size_t i = 0; i < phreeqcVol.size(); i++)
        {
            std::cout << phreeqcVol.at(i) << ",";
            std::cout << phreeqcPor.at(i) << ",";
            std::cout << phreeqcSat.at(i) << ",";
            //std::cout << watCont.at(i) << ",";
            std::cout << phreeqcDens.at(i) << ",";
            std::cout << phreeqcTemp.at(i) << ",";
            std::cout << phreeqcPres.at(i);
            std::cout << "\n";
        }

        // Print concentration table.
        std::cout << "\n";
        /*
        for (size_t j = 0; j < cellsPerSolute3d; j++)
        {
            for (size_t i = 0; i < 2 * num_of_species; i++)
            {
                if (i < num_of_species)
                {
                    size_t index = i * cellsPerSolute3d + j;
                    std::cout << concIn.at(index);
                }
                else if (i >= num_of_species)
                {
                    size_t index = (i - num_of_species) * cellsPerSolute3d + j;
                    std::cout << concOut.at(index);
                }
                if (i < 2 * num_of_species - 1)
                {
                    std::cout << ",";
                }
            }
            std::cout << "\n";
        }
        */
        return 3; // IRM_FAIL;
    }
}