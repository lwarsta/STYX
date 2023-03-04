#ifndef _ANALYTICALMODELS_H
#define _ANALYTICALMODELS_H
#include <iostream>
#include <vector>
#include <cmath>
#include "FileIO.h"
#include "Grid2d.h"
#include "Grid3d.h"
#include "ModelWater3dBrute.h"
#include "ModelWater3dTri.h"
#include "ModelSolute3dBrute.h"
#include "ModelSolute3dTri.h"

class AnalyticalModels {
private:
	double transport1d(double c_in, double v, double d, double kb, double t, 
		double x);
	FileIO fileIO;
	ModelWater3dBrute modelWater3dBrute;
	ModelSolute3dBrute modelSolute3dBrute;
	ModelSolute3dTri modelSolute3dTri;
	void run_1d_flow_comparison();
	double calcWatCont(double wat_cont_sat, double wat_cont_res, 
		double press_head, double press_head_sat, double press_head_res);
	double calcRelCond(double press_head, double press_head_sat, 
		double press_head_res);
	void flow_model(Grid3d& grid3d, double timeStep, double iterCutThresh, 
		int iterStop);
	void run_1d_transport_comparison();
public:
	void run_comparisons();
};

#endif
