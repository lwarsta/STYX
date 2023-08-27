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
public:
	double get_diameter() { return diameter; }
	LinkGeom(int id_new, char geom_type_new, std::vector<int> vert_indices_new,
		     std::vector<Vertex*> vert_pointers_new);
};

#endif
