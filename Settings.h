#ifndef _SETTINGS_H
#define _SETTINGS_H
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <map>

class Settings {
private:
public:
    Settings();
    ~Settings();
	void initialize(std::vector < std::vector<std::string> > settings);
    std::string get_str(std::string param_name);
    int get_int(std::string param_name);
    double get_double(std::string param_name);
    void set_value(std::string name, std::string value);
    void set_value(std::string name, int value);
    void set_value(std::string name, double value);
    // Parameters.
    std::map<std::string, std::string> settings_str;
    /*
    bool runTests;
    int numOfThreads;
    int simTimeMax;
    double timeStep;
    int timePrintThresh;
    int timePrintVtkStart;
    int timePrintVtkEnd;
    int watFlow2dSolver;
    int iterStopWat2d;
    double iterCutThreshWat2d;
    double implicityWat2d;
    double bisIterCutThreshWat2d;
    int bisIterStopWat2d;
    double bisIterCutLeftWat2d;
    double bisIterCutRightWat2d;
    int watFlowSolver;
    int cellIndex;
    int iterStopWat;
    double iterCutThreshWat;
    double implicityWat;
    double areaFact;
    double petFact;
    int heatTransSolver;
    int iterStopHeat;
    double iterCutThreshHeat;
    double implicityHeat;
    double bisIterCutThreshHeat;
    int bisIterStopHeat;
    double bisIterCutLeftHeat;
    double bisIterCutRightHeat;
    double latHeatFus;
    double tempFreez;
    double multFreez;
    double densAir;
    double densIce;
    double densWat;
    double capAir;
    double capIce;
    double capWat;
    double condAir;
    double condIce;
    double condWat;
    double condMultAir;
    double condMultIce;
    double condMultWat;
    double tempSoilInit;
    bool usePhreeqcrm;
    int solTransSolver;
    size_t numOfSolutes;
    int iterStopSol;
    double iterCutThreshSol;
    double implicitySol;
    double molDiffRate ;
    std::string mesh2dfilePathIn ;
    std::string mesh3dfilePathIn;
    std::string gridMapfilePathIn ;
    std::string materials2dPath ;
    std::string materials3dPath ;
    std::string initCond2dPath ;
    std::string initCond3dPath ;
    std::string boundCond2dPath ;
    std::string boundCond3dPath ;
    std::string atmosForcingdataPath ;
    std::string solutePropertiesPath ;
    std::string phreeqcInputPath ;
    std::string phreeqcChemDbPath ;
    std::string phreeqcSpecMapPath ;
    std::string phreeqcMolWeightMapPath;
    std::string outputCsvPath;
    std::string outputSurfVtkPath;
    std::string outputSubsVtkPath;
    */
};

#endif

