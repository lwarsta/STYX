#ifndef _MODELWATER3DTRI_H
#define _MODELWATER3DTRI_H
#include <omp.h>
#include "Grid2d.h"
#include "Grid3d.h"

class ModelWater3dTri
{
public:
    void configure(int iterStopNew, double iterCutThreshNew, double implicityNew, double timeStepNew, std::vector < std::vector<int> > gridMapNew);
    void run(Grid2d &grid2d, Grid3d &grid);
private:
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
    int iterStop;
    double iterCutThresh;
    double implicity;
    double timeStep;
    std::vector < std::vector<int> > gridMap;
};

#endif
