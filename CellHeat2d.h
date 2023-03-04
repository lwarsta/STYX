#ifndef _CELLHEAT2D_H
#define _CELLHEAT2D_H
#include "CellBase.h"
#include "CellGeom2d.h"
#include "CellWater2d.h"

class CellHeat2d : public CellBase
{
protected:
    CellGeom2d * geom;
    CellWater2d * water;
    std::vector<CellHeat2d*> neighPointers;
public:
    CellHeat2d();
    void assignGeom(CellGeom2d *geomNew){geom = geomNew;}
    CellGeom2d * getGeom(){return geom;}
    void assignWater(CellWater2d *waterNew){water = waterNew;}
    CellWater2d * getWater(){return water;}
    void setNeighPointers(std::vector<CellHeat2d*> neighPointersNew) {
        neighPointers = neighPointersNew;}
};

#endif
