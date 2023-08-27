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
public:
    LinkWater();
    void assign_geom(LinkGeom* geom_new) { geom = geom_new; }
    LinkGeom* get_geom() { return geom; }
    void comp_flow_area_and_hydr_rad(double water_depth, double& flow_area, 
                                     double& hydr_rad);
    double get_mann_n() { return mann_n; }
};

#endif