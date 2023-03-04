#ifndef _MODELHEAT3DTRI_H
#define _MODELHEAT3DTRI_H
#include <omp.h>
#include "Grid2d.h"
#include "Grid3d.h"

class ModelHeat3dTri
{
public:
    void configure(int iterStopNew, double iterCutThreshNew, 
        double implicityNew, double timeStepNew, double boundTempNew,
        double iterThreshBisNew, int iterStopBisNew, double threshLeftBisNew, 
        double threshRightBisNew, std::vector < std::vector<int>> gridMapNew);
    void run(Grid2d &grid2d, Grid3d &grid);
protected:
    int iterStop;
    double iterCutThresh;
    double implicity;
    double timeStep;
    double boundTemp;
    double iterThreshBis;
    int iterStopBis;
    double threshLeftBis;
    double threshRightBis;
    std::vector < std::vector<int> > gridMap;
    void preprocess(Grid2d &grid2d, Grid3d &grid3d);
    void triDiag(
        std::vector<double>& triA,
        std::vector<double>& triB,
        std::vector<double>& triC,
        std::vector<double>& triD,
        std::vector<double>& triG,
        std::vector<double>& triU,
        std::vector<double>& triV,
        std::vector<double>& triX);
    void iterate(Grid2d &grid2d, Grid3d &grid3d);
    void postprocess(Grid2d &grid2d, Grid3d &grid3d);
    double bisection(CellHeat3d& heat, double left, double right);
};

#endif
