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
    double up_stor_depth;
    double mannN;
    double deprStor;
    double up_stor_thresh;
    double evap_frac;
    double crop_factor;
    double root_depth;
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
    void set_up_stor_depth(double up_stor_depth_new) { up_stor_depth = up_stor_depth_new; }
    double get_up_stor_depth() { return up_stor_depth; }
    void setMannN(double mannNNew){mannN = mannNNew;}
    double getMannN() { return mannN; }
    void setDeprStor(double deprStorNew){deprStor = deprStorNew;}
    double getDeprStor() { return deprStor; }
    double get_up_stor_thresh() { return up_stor_thresh; }
    void set_up_stor_thresh(double up_stor_thresh_new) { up_stor_thresh = up_stor_thresh_new; }
    double get_evap_frac() { return evap_frac; }
    void set_evap_frac(double evap_frac_new) { evap_frac = evap_frac_new; }
    double get_crop_factor() { return crop_factor; }
    void set_crop_factor(double crop_factor_new) { crop_factor = crop_factor_new; }
    double get_root_depth() { return root_depth; }
    void set_root_depth(double root_depth_new) { root_depth = root_depth_new; }
    void swap() { waterDepthOld = waterDepth; }
    void setOutletIndex(int outletIndexNew) { outletIndex = outletIndexNew; }
    int getOutletIndex() { return outletIndex; }
    void setDistToOutlet(int distToOutletNew) { 
        distToOutlet = distToOutletNew; }
    double getDistToOutlet() { return distToOutlet; }
    void setAvgDisch(double avgDischNew) { avgDisch = avgDischNew; }
    double getAvgDisch() { return avgDisch; }
    void setSinkIndex(int sinkIndexNew) { sinkIndex = sinkIndexNew; }
    void removeWatBySinks(double timeStep);
    double getWatVolRemovedBySinks() { return watVolToSinks; }
};

#endif
