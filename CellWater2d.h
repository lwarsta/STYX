#ifndef _CELLWATER2D_H
#define _CELLWATER2D_H
#include "CellBase.h"
#include "CellGeom2d.h"

struct StormWatItem {
    double watVol;
    double dischTime;
    double distToOutlet;
};

class CellWater2d : public CellBase
{
protected:
    CellGeom2d * geom;
    std::vector<CellWater2d*> neighPointers;
    double waterDepth;
    double waterDepthOld;
    double deprStor;
    double mannN;
    int outletIndex;
    double distToOutlet;
    std::vector<StormWatItem> stormWatItems;
    double avgDisch;
    int sinkIndex;
    double watVolToSinks;
public:
    CellWater2d();
    void assignGeom(CellGeom2d *geomNew){geom = geomNew;}
    CellGeom2d * getGeom(){return geom;}
    CellWater2d * getNeigh(int index) { return neighPointers.at(index); }
    void setNeighPointers(std::vector<CellWater2d*> neighPointersNew) {
        neighPointers = neighPointersNew;}
    void setWaterDepth(double waterDepthNew){waterDepth = waterDepthNew;}
    double getWaterDepth(){return waterDepth;}
    double getWaterDepthOld() { return waterDepthOld; }
    void setMannN(double mannNNew){mannN = mannNNew;}
    double getMannN() { return mannN; }
    void setDeprStor(double deprStorNew){deprStor = deprStorNew;}
    double getDeprStor() { return deprStor; }
    void swap() { waterDepthOld = waterDepth; }
    void setOutletIndex(int outletIndexNew) { outletIndex = outletIndexNew; }
    int getOutletIndex() { return outletIndex; }
    void setDistToOutlet(int distToOutletNew) { 
        distToOutlet = distToOutletNew; }
    double getDistToOutlet() { return distToOutlet; }
    double compFlowToManhole(double timeStep);
    void addWatVolToStorm(double watVolStorm, double distToOutlet, 
        double simTime);
    void removeWatVolFromStorm(double simTime);
    void setAvgDisch(double avgDischNew) { avgDisch = avgDischNew; }
    double getAvgDisch() { return avgDisch; }
    void setSinkIndex(int sinkIndexNew) { sinkIndex = sinkIndexNew; }
    void removeWatBySinks(double timeStep);
    double getWatVolRemovedBySinks() { return watVolToSinks; }
};

#endif
