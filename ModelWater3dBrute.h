#ifndef _MODELWATER3DBRUTE_H
#define _MODELWATER3DBRUTE_H
#include <omp.h>
#include "Grid2d.h"
#include "Grid3d.h"

class ModelWater3dBrute
{
public:
    void configure(int iterStopNew, double iterCutThreshNew, double implicityNew, double timeStepNew);
    void run(Grid2d &grid2d, Grid3d &grid);
private:
    void preprocess(Grid2d &grid2d, Grid3d &grid3d);
    void iterate(Grid2d &grid2d, Grid3d &grid3d);
    void postprocess(Grid2d &grid2d, Grid3d &grid3d);
    int iterStop;
    double iterCutThresh;
    double implicity;
    double timeStep;
};

#endif
