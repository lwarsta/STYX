#ifndef _CELLSOLUTE2D_H
#define _CELLSOLUTE2D_H
#include "CellBase.h"
#include "CellGeom2d.h"
#include "CellWater2d.h"

class CellSolute2d : public CellBase
{
protected:
    CellGeom2d * geom;
    CellWater2d * water;
    std::vector<CellSolute2d*> neighPointers;
    double conc;
    double concOld;
    double mass;
    double massOld;
    double decayRate;
    int decayTarget;
    double adsParam[3]; // to vector?
    double depos;
    double deposDryRate;
    double deposWetRate;
public:
    CellSolute2d();
    void assignGeom(CellGeom2d *geomNew){geom = geomNew;}
    CellGeom2d * getGeom(){return geom;}
    void assignWater(CellWater2d *waterNew){water = waterNew;}
    CellWater2d * getWater(){return water;}
    void setNeighPointers(std::vector<CellSolute2d*> neighPointersNew) {
        neighPointers = neighPointersNew;}
    void setConc(double concNew) {conc = concNew;}
    double getConc() {return conc;}
    double getConcOld() {return concOld;}
    void swap() {concOld = conc;}
    void convertMassToConc();
    void convertConcToMass();
    void setMass(double massNew) {mass = massNew;}
    double getMass() {return mass;}
    double getMassOld() {return massOld;}
    void setAdsParam(int param, double adsParamNew) {
        adsParam[param] = adsParamNew;}
    void setDecayRate(double decayRateNew){decayRate = decayRateNew;}
    void setDecayTarget(int decayTargetNew){decayTarget = decayTargetNew;}
    int getDecayTarget(){return decayTarget;}
    double getFluxDecay(){return decayRate * water->getWaterDepth() * 
        geom->getArea();}
    void setFluxDepos(double deposNew){depos = deposNew;}
    double getFluxDepos(){return depos;}
    void setDeposWetRate(double deposWetRateNew){
        deposWetRate = deposWetRateNew;}
    double getDeposWetRate(){return deposWetRate;}
    void setDeposDryRate(double deposDryRateNew) { 
        deposDryRate = deposDryRateNew; }
    double getDeposDryRate() { return deposDryRate; }
};

#endif
