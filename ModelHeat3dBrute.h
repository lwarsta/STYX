#ifndef _MODELHEAT3DBRUTE_H
#define _MODELHEAT3DBRUTE_H
#include <omp.h>
#include "Grid2d.h"
#include "Grid3d.h"

class ModelHeat3dBrute
{
public:
    void configure(int iterStopNew, double iterCutThreshNew, 
        double implicityNew, double timeStepNew, double boundTempNew, 
        double iterThreshBisNew, int iterStopBisNew, double threshLeftBisNew, 
        double threshRightBisNew);
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
    void preprocess(Grid2d &grid2d, Grid3d &grid3d);
    void iterate(Grid2d &grid2d, Grid3d &grid3d);
    void postprocess(Grid2d &grid2d, Grid3d &grid3d);
    double bisection(CellHeat3d& heat, double left, double right);
};

#endif
