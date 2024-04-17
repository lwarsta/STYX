#ifndef _FRAMEWORK_H
#define _FRAMEWORK_H
#include <sstream>
#include <iomanip>
#include <regex>
#include "Settings.h"
#include "AnalyticalModels.h"
#include "FileIO.h"
#include "AtmosControl.h"
#include "Network.h"
#include "Grid2d.h"
#include "Grid3d.h"
#include "ModelWaterNetDiffBrute.h"
#include "ModelWaterNetExplicit.h"
#include "ModelWater2dDiffBrute.h"
#include "ModelWater3dBrute.h"
#include "ModelWater3dTri.h"
#include "ModelHeat3dBrute.h"
#include "ModelHeat3dTri.h"
#include "ModelSolute3dBrute.h"
#include "ModelSolute3dTri.h"
#include "PhreeqcWrapper.h"

class Framework {
public:
    Framework();
    ~Framework();
    int initialize(std::string pathToSettings);
    int run();
    void saveOutputVTK(std::string path2d, std::string path3d);
    int finalize();
private:
    Settings settings;
    FileIO fileIO;
    AtmosControl atmosControl;
    Network network;
    Grid2d grid2d;
    Grid3d grid3d;
    ModelWaterNetDiffBrute modelWaterNetDiffBrute;
    ModelWaterNetExplicit modelWaterNetExplicit;
    ModelWater2dDiffBrute modelWater2dDiffBrute;
    ModelWater3dBrute modelWater3dBrute;
    ModelWater3dTri modelWater3dTri;
    ModelHeat3dBrute modelHeat3dBrute;
    ModelHeat3dTri modelHeat3dTri;
    ModelSolute3dBrute modelSolute3dBrute;
    ModelSolute3dTri modelSolute3dTri;
    PhreeqcWrapper phreeqcWrapper;
    std::vector < std::vector<std::string>> results;
    std::vector < std::vector<int>> gridMap;
    std::vector <std::string> species_names;
};

#endif
