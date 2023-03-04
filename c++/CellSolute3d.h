#ifndef _CELLSOLUTE3D_H
#define _CELLSOLUTE3D_H
#include "CellBase.h"
#include "CellGeom3d.h"
#include "CellWater3d.h"

class CellSolute3d : public CellBase
{
protected:
    CellGeom3d * geom;
    CellWater3d * water;
    std::vector<CellSolute3d*> neighPointers;
    double conc;
    double concOld;
    double mass;
    double massOld;
    std::vector<double> dispFlux;
    //double dispFlux[6]; // to vector?
    //double disp[3]; // to vector?
    std::vector<double> disp;
    //double dispInter[6]; // to vector?
    std::vector<double> dispInter;
    //double advFlux[6][2]; // to vector?
    std::vector<std::vector<double> > advFlux;
    double alphaCofLong;
    double alphaCofTrans;
    double molDiffRate;
    double decayRate[2]; // to vector?
    int decayTarget[2]; // to vector?
    double adsParam[3]; // to vector?
    double adsMult;
    double ads;
    double decayRateInit; // where is this used?
    double infConc;
    double infMass;
    double drainMass;
    double decayMassOut[2];
    double decayMassIn[2];
public:
    CellSolute3d();
    void assignGeom(CellGeom3d *geomNew){geom = geomNew;}
    CellGeom3d * getGeom(){return geom;}
    void assignWater(CellWater3d *waterNew){water = waterNew;}
    CellWater3d * getWater(){return water;}
    void setNeighPointers(std::vector<CellSolute3d*> neighPointersNew);
    CellSolute3d * getNeigh(int index){return neighPointers.at(index);}
    //void calcAds();
    void calcDisp();
    void calcDispInter();
    double getDisp(int index) {return disp[index];}
    void calcDispFluxes();
    double getDispFlux(int index){return dispFlux[index];}
    void calcAdvFluxes();
    double getAdvFlux(int index0, int index1) {return advFlux[index0][index1];}
    double getDecayRate(int index){if (index >= 0 && index < 2) 
        return decayRate[index]; else return 0.0;}
    double getDecayRateInit(){return decayRateInit;}
    double getFluxDecay(int index){if (index >= 0 && index < 2) 
        return decayRate[index] * concOld * water->getWatCont() * 
        geom->getVolume() * adsMult; else return 0.0;}
    void calcAdsMult();
    double getAdsMult(){return adsMult;}
    double getConcNeigh(int index);
    double getConcNeighOld(int index);
    double getMass() {return mass;}
    void setMass(double massNew) {mass = massNew;}
    double getMassOld() {return massOld;}
    void setConc(double concNew) {conc = concNew;}
    double getConc(){return conc;}
    double getConcOld(){return concOld;}
    void swap(){concOld = conc;}
    void convertMassToConc();
    void convertConcToMass();
    void setDispCofs(double alphaCofLongNew, double alphaCofTransNew) {
        alphaCofLong = alphaCofLongNew; alphaCofTrans = alphaCofTransNew;}
    void setMolDiffRate(double molDiffRateNew){molDiffRate = molDiffRateNew;}
    void setAdsParam(int index, double adsParamNew){
        adsParam[index] = adsParamNew;}
    void setDecayRate(int index, double decayRateNew){ 
        if (index >= 0 && index < 2) decayRate[index] = decayRateNew; }
    void setDecayTarget(int index, int decayTargetNew){ 
        if (index >= 0 && index < 2) decayTarget[index] = decayTargetNew; }
    void setDecayRateInit(double decayRateNew){ 
        decayRateInit = decayRateNew; }
    int getDecayTarget(int index){if (index >= 0 && index < 2)
        return decayTarget[index]; else return -1;}
    void setInfConc(double infConcNew){infConc = infConcNew;}
    double getInfConc(){return infConc;}
    void setInfMass(double infMassNew){infMass = infMassNew;}
    double getInfMass(){return infMass;}
    void setDrainMass(double drainMassNew){drainMass = drainMassNew;}
    double getDrainMass(){return drainMass;}
    void setDecayMassOut(int decayPath, double decayMassOutNew){
        decayMassOut[decayPath] = decayMassOutNew;}
    double getDecayMassOut(int decayPath){return decayMassOut[decayPath];}
    void setDecayMassIn(int decayPath, double decayMassInNew){
        decayMassIn[decayPath] = decayMassInNew;}
    double getDecayMassIn(int decayPath){return decayMassIn[decayPath];}
};

#endif

