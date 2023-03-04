#ifndef _MODELSOLUTE3DBRUTE_H
#define _MODELSOLUTE3DBRUTE_H
#include <omp.h>
#include "Grid2d.h"
#include "Grid3d.h"

class ModelSolute3dBrute
{
public:
    void configure(int iterStopNew, double iterCutThreshNew, 
        double implicityNew, double timeStepNew, int numOfSolutes);
    void run(Grid2d &grid2d, Grid3d &grid);
protected:
    int iterStop;
    double iterCutThresh;
    double implicity;
    double timeStep;
    size_t numOfSolutes;
    void preprocess(Grid2d &grid2d, Grid3d &grid3d);
    void iterate(Grid2d &grid2d, Grid3d &grid3d);
    void postprocess(Grid2d &grid2d, Grid3d &grid3d);
};

#endif
