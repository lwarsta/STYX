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
	double length;
	double length_flat;
	double slope;
public:
	void set_diameter(double diameter_new) { diameter = diameter_new; }
	double get_diameter() { return diameter; }
	double get_area() { return area; }
	double get_length() { return length; }
	double get_length_flat() { return length_flat; }
	double get_slope() { return slope; }
	LinkGeom(int id_new, char geom_type_new, std::vector<int> vert_indices_new,
		     std::vector<Vertex*> vert_pointers_new);
	void comp_geom_properties();
};

#endif
