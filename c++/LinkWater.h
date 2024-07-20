#ifndef _LINKWATER_H
#define _LINKWATER_H
#include "Algorithms.h"
#include "LinkGeom.h"
#include "CellBase.h"

class LinkWater : public CellBase
{
protected:
    Algorithms algorithms;
    LinkGeom* geom;
    double mann_n;
    // A common parent could be done for the junction and link classes with these.
    double water_depth;
    double water_depth_old;
public:
    LinkWater();
    void assign_geom(LinkGeom* geom_new) { geom = geom_new; }
    LinkGeom* get_geom() { return geom; }
    void comp_flow_area_and_hydr_rad(double water_depth, double& flow_area, 
                                     double& hydr_rad);
    double calculateWaterDepth(double area, double radius);
    double get_mann_n() { return mann_n; }
    void set_mann_n(double mann_n_new) { mann_n = mann_n_new; }
    // A common parent could be done for the junction and link classes with these.
    void set_water_depth(double water_depth_new) { water_depth = water_depth_new; }
    double get_water_depth() { return water_depth; }
    double get_water_depth_old() { return water_depth_old; }
    void swap() { water_depth_old = water_depth; }
};

#endif