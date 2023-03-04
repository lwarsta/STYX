#ifndef _PHREEQCWRAPPER_H
#define _PHREEQCWRAPPER_H

#include <vector>
#include <string>
#include <map>
#include "FileIO.h"
#include "Settings.h"
#include "PhreeqcRM.h"
#include "Grid2d.h"
#include "Grid3d.h"

class PhreeqcWrapper
{
private:
    int countWordOccurrences(std::string input, std::string search_word);
    FileIO fileIO;
    PhreeqcRM* phreeqc_rm;
    int num_of_threads;
    std::vector<double> phreeqcVol;
    std::vector<double> phreeqcPor;
    std::vector<double> phreeqcSat;
    std::vector<double> phreeqcDens;
    std::vector<double> phreeqcTemp;
    std::vector<double> phreeqcPres;
    std::vector<double> phreeqcConc;
    std::vector<double> phreeqcConcInit;
    int num_of_cells_3d;
    int num_of_species;
    std::vector <std::string> species_names;
    size_t posWaterMol;
    std::map<std::string, double> molWeightMap;
    std::vector <std::vector <double>> bound_cond_conc;
public:
    PhreeqcWrapper();
    ~PhreeqcWrapper();
    void initialize(
        Settings &settings, 
        std::vector<std::vector<std::string>>& boundCond2d, 
        std::vector<std::vector<std::string>> &boundCond3d, 
        Grid2d &grid2d, Grid3d &grid3d);
    int run(Settings& settings, Grid2d& grid2d, Grid3d& grid3d);
    int get_num_of_species() { return num_of_species; }
    std::vector <std::string> get_species_names() { return species_names; }
    double get_bound_cond(size_t index);
    std::vector<double> * get_phreeqc_conc() { return &phreeqcConc; }
    std::map<std::string, double> * get_mol_weight_map() { return &molWeightMap; }
};

#endif