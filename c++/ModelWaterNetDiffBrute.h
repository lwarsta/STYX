#ifndef _MODELWATERNETDIFFBRUTE_H
#define _MODELWATERNETDIFFBRUTE_H
#include <omp.h>
#include "Algorithms.h"
#include "Grid2d.h"
#include "Network.h"

class ModelWaterNetDiffBrute
{
public:
    void configure(int iter_stop_new, double iter_cut_thresh_new, double implicity_new, double time_step_new, 
                   double iter_thresh_bis_new, int iter_stop_bis_new, double thresh_left_bis_new, double thresh_right_bis_new);
    void run(Grid2d &grid2d, Network &network);
private:
    void preprocess(Grid2d &grid2d, Network &network);
    double calc_residual(JuncWater& water, double dHw0);
    double bisection(JuncWater& water, double left, double right);
    void iterate(Grid2d &grid2d, Network &network);
    void postprocess(Grid2d &grid2d, Network &network);
    int iter_stop;
    double iter_cut_thresh;
    double implicity;
    double time_step;
    double iter_thresh_bis;
    int iter_stop_bis;
    double thresh_left_bis; // should be zero
    double thresh_right_bis;
};

#endif
