#ifndef _CELLHEAT3D_H
#define _CELLHEAT3D_H
#include "CellBase.h"
#include "CellGeom3d.h"
#include "CellWater3d.h"

class CellHeat3d : public CellBase
{
protected:
    CellGeom3d * geom;
    CellWater3d * water;
    std::vector<CellHeat3d*> neighPointers;
	std::vector<double> fluxDiff;
	std::vector < std::vector<double> > fluxConv;
	std::vector<double> interK;
	std::vector<double> gradLimiters;
	double C;
	double K;
	double oldT;
	double T;
	double U;
	double fracSoil;
	double fracWat;
	double fracIce;
	double fracWatTot; // includes both liquid water and ice
	double fracAir;
	double densWat;
	double densIce;
	double densAir;
	double capSoil;
	double capWat;
	double capIce;
	double capAir;
	double condSoil;
	double condWat;
	double condIce;
	double condAir;
	double condMultSoil;
	double condMultWat;
	double condMultIce;
	double condMultAir;
	double latHeatFus;
	double multFreez;
	double tempFreez;
	double gradMax;
	double surfDiffFlux;
public:
    CellHeat3d();
    void assignGeom(CellGeom3d *geomNew){geom = geomNew;}
    CellGeom3d * getGeom(){return geom;}
    void assignWater(CellWater3d *waterNew){water = waterNew;}
    CellWater3d * getWater(){return water;}
    void setNeighPointers(std::vector<CellHeat3d*> neighPointersNew);
    CellHeat3d * getNeigh(int index){return neighPointers.at(index);}
	double calcResidual(double newT);
	double getC() {return C;}
	double getConvFlux(int dir, int inOut) {return fluxConv[dir][inOut];}
	double getDiffFlux(int dir) {return fluxDiff[dir];}
	double getFracAir() {return fracAir;}
	double getFracIce() {return fracIce;}
	double getFracSoil() {return fracSoil;}
	double getFracWater() {return fracWat;}
	double getK() {return K;}
	double getOldT() {return oldT;}
	double getT() {return T;}
	double getU() {return U;}
	void calcC();
	void calcCond();
	void calcConvEnergy(double dt);
	void calcConvFluxes();
	void calcDiffEnergy(double dt);
	void calcDiffFluxes();
	void calcFracWatTot() {fracWatTot = water->getWatCont() + 
		water->getFracIce();}
	void calcInterCond();
	void calcU();
	void calcWatAndIceFracs();
	void setT(double tempNew) {T = tempNew;}
	void swap() { oldT = T; }
	void compGradLimiters();
	void calcSurfDiffFlux();
	double getSurfDiffFlux() {return surfDiffFlux;}
	void setCondAir(double condAirNew) {condAir = condAirNew;}
	void setCondIce(double condIceNew) {condIce = condIceNew;}
	void setCondMultAir(double condMultAirNew) {condMultAir = condMultAirNew;}
	void setCondMultIce(double condMultIceNew) {condMultIce = condMultIceNew;}
	void setCondMultSoil(double condMultSoilNew) {condMultSoil = condMultSoilNew;}
	void setCondMultWat(double condMultWatNew) {condMultWat = condMultWatNew;}
	void setCondSoil(double condSoilNew) {condSoil = condSoilNew;}
	void setCondWat(double condWatNew) {condWat = condWatNew;}
	void setDensAir(double densAirNew) {densAir = densAirNew;}
	void setDensIce(double densIceNew) {densIce = densIceNew;}
	void setDensWat(double densWatNew) {densWat = densWatNew;}
	void setFreezCurvMult(double multFreezNew) {multFreez = multFreezNew;}
	void setHeatCapAir(double capAirNew) {capAir = capAirNew;}
	void setHeatCapIce(double capIceNew) {capIce = capIceNew;}
	void setHeatCapSoil(double capSoilNew) {capSoil = capSoilNew;}
	void setHeatCapWat(double capWatNew) {capWat = capWatNew;}
	void setLatHeatFus(double latHeatFusNew) {latHeatFus = latHeatFusNew;}
	void setTempFreez(double tempFreezNew) {tempFreez = tempFreezNew;}
};

#endif

