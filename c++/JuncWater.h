#ifndef _JUNCWATER_H
#define _JUNCWATER_H
#include "LinkJuncBase.h"
#include "JuncGeom.h"

class JuncWater : public LinkJuncBase
{
protected:
    JuncGeom * geom;
    double water_depth;
    double water_depth_old;
public:
    JuncWater();
    void assign_geom(JuncGeom *geom_new){geom = geom_new;}
    JuncGeom * get_geom(){return geom;}
    void set_water_depth(double water_depth_new){water_depth = water_depth_new;}
    double get_water_depth(){return water_depth;}
    double get_water_depth_old() { return water_depth_old; }
    void swap() { water_depth_old = water_depth; }
};

#endif
