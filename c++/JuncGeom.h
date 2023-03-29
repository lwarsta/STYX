#ifndef _JUNCGEOM_H
#define _JUNCGEOM_H
#include "Vertex.h"
#include "LinkJuncBase.h"

class JuncGeom : public LinkJuncBase
{
protected:
	Vertex coord;
	double diameter;
	double area;
	bool open;
	double depth;
public:
	JuncGeom(int id_new, double pos_x, double pos_y, double pos_z,
		     double diameter_new, bool open_new, double depth_new);
	void compute_area();
	double get_area(){return area;}
};

#endif
