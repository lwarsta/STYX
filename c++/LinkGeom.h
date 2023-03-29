#ifndef _LINKGEOM_H
#define _LINKGEOM_H
#include "Vertex.h"
#include "LinkJuncBase.h"

class LinkGeom : public LinkJuncBase
{
protected:
	double diameter;
	double area;
	// double length; // ??
	double roughness;
	int conn_ind[2]; // -1 = not connected
	int conn_type[2]; // 0 = junction, 1 == cell
	int conn_elev[2];
public:
	LinkGeom(int id_new, double diameter_new, double roughness_new,
		     int conn_ind_00, int conn_ind_01, int conn_type_00,
		     int conn_type_01, double conn_elev_00, double conn_elev_01);
};

#endif
