#ifndef _MODELWATER2DDIFFBRUTE_H
#define _MODELWATER2DDIFFBRUTE_H
#include <omp.h>
#include "Grid2d.h"
#include "Grid3d.h"

class ModelWater2dDiffBrute
{
public:
    void configure(int iterStopNew, double iterCutThreshNew, double implicityNew, double timeStepNew, 
                   double iterThreshBisNew, int iterStopBisNew, double threshLeftBisNew, double threshRightBisNew);
    void run(Grid2d &grid2d, Grid3d &grid);
private:
    void preprocess(Grid2d &grid2d, Grid3d &grid3d);
    double calcResidual(CellWater2d& water, double dHw0);
    double bisection(CellWater2d& water, double left, double right);
    void iterate(Grid2d &grid2d, Grid3d &grid3d);
    void postprocess(Grid2d &grid2d, Grid3d &grid3d);
    int iterStop;
    double iterCutThresh;
    double implicity;
    double timeStep;
    double iterThreshBis;
    int iterStopBis;
    double threshLeftBis; // should be zero
    double threshRightBis;
};

#endif
