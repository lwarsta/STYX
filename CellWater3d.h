#ifndef _CELLWATER3D_H
#define _CELLWATER3D_H
#include "CellBase.h"
#include "CellGeom3d.h"

class CellWater3d : public CellBase
{
protected:
    CellGeom3d * geom;
    std::vector<CellWater3d*> neighPointers;
    std::vector<double> fluxes;
    std::vector<double> condInter;
    std::vector<double> velInter;
    double hydrHead;
    double hydrHeadOld;
    double presHead;
	double watContRes;
	double watContSat;
	double watCont;
	double watContOld;
	double diffWatCap;
	double alpha;
	double m;
	double n;
	double dryWeight;
	double compress;
	double condSat;
	double condUnsat;
	double fluxDrain;
	double drainFactor;
	double dischDrain;
	double fluxInf;
	double infVolume;
	double petCellFact;
	double fluxEvap;
	double evapVolume;
    double fracIce;
public:
    CellWater3d();
    void assignGeom(CellGeom3d *geomNew);
    void setNeighPointers(std::vector<CellWater3d*> neighPointersNew) {
        neighPointers = neighPointersNew;}
    CellGeom3d * getGeom(){return geom;}
    CellWater3d * getNeigh(int index){return neighPointers.at(index);}
    void setAlpha(double alphaNew){alpha = alphaNew;}
    void setN(double nNew);
    void setDryWeight(double dryWeightNew){dryWeight = dryWeightNew;}
    double getDryWeight(){return dryWeight;}
    void setWatContSat(double watContSatNew){watContSat = watContSatNew;}
    double getWatContSat(){return watContSat;}
    void setWatContRes(double watContResNew){watContRes = watContResNew;}
    double getWatContRes() { return watContRes; }
    void setCompress(double compressNew){compress = compressNew;}
    void setCondSat(double condSatNew){condSat = condSatNew;}
    double getCondSat(){return condSat;}
    double getDiffWatCap(){return diffWatCap;}
    double getHydrHead(){return hydrHead;}
    double getHydrHeadOld(){return hydrHeadOld;}
    void setHydrHead(double hydrHeadNew){hydrHead = hydrHeadNew;}
    double getPresHead(){return presHead;}
    double getWatCont(){return watCont;}
    double getFlux(int index){return fluxes.at(index);}
    double getVelInt(int index){return velInter.at(index);}
    void swap(){hydrHeadOld = hydrHead; watContOld = watCont;}
    void calcUnsatCond();
    double getCondUnsat(){return condUnsat;}
    void setCondUnsat(double condUnsatNew) {condUnsat = condUnsatNew;}
    void setDrainFactor(double drainFactorNew){drainFactor = drainFactorNew;}
    void calcPresHead();
    void calcCondInter();
    void calcFluxes();
    void calcFlowVel();
    void calcWatCont();
    void setWatCont(double watContNew) {watCont = watContNew;}
    void calcDiffWatCap();
    double getRelSat();
    double getFluxDrain(){return fluxDrain;}
    void calcFluxDrain();
    void setDrainVol(double dischDrainNew){dischDrain = dischDrainNew;}
    double getDrainVol(){return dischDrain;}
    double getFluxInf(){return fluxInf;}
    void setFluxInf(double fluxInfNew){fluxInf = fluxInfNew;}
    void setInfVol(double infVolumeNew){infVolume = infVolumeNew;}
    double getInfVol(){return infVolume;}
    void setPetFactor(double petCellFactNew){petCellFact = petCellFactNew;}
    void calcFluxEvap(double pet);
    double getFluxEvap(){return fluxEvap;}
    void setEvapVol(double evapVolumeNew){evapVolume = evapVolumeNew;}
    double getEvapVol(){return evapVolume;}
    void setFracIce(double fracIceNew) {fracIce = fracIceNew;}
    double getFracIce() {return fracIce;}
};

#endif
