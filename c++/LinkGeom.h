#ifndef _LINKGEOM_H
#define _LINKGEOM_H
#include "Algorithms.h"
#include "Vertex.h"
#include "LinkJuncGeomBase.h"

class LinkGeom : public LinkJuncGeomBase
{
protected:
	double diameter;
	double area;
	//double length;
	double slope;
public:
	double get_diameter() { return diameter; }
	void set_diameter(double diameter_new) { diameter = diameter_new; }
	double get_slope() { return slope; }
	LinkGeom(int id_new, char geom_type_new, std::vector<int> vert_indices_new,
		     std::vector<Vertex*> vert_pointers_new);
	void comp_geom_properties();
};

#endif
