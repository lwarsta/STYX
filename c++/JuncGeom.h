#ifndef _JUNCGEOM_H
#define _JUNCGEOM_H
#include "Algorithms.h"
#include "Vertex.h"
#include "LinkJuncGeomBase.h"
#include "LinkGeom.h"

class JuncGeom : public LinkJuncGeomBase
{
protected:
	Algorithms algorithms;
	double diameter;
	double area;
	bool open;
	double depth;
	std::vector<int> ids_link_end;
	std::vector<int> ids_link;
	std::vector<LinkGeom*> links;
	std::vector<int> ids_junc;
	std::vector<JuncGeom*> juncs;
public:
	JuncGeom(int idNew, char geomTypeNew, std::vector<int> vertIndicesNew,
		     std::vector<Vertex*> vertPointersNew);
	void save_link(int id_lnk_end, int id_lnk, LinkGeom* link_geom,
		           int id_jnc, JuncGeom* junc_geom);
	std::vector<int> get_ids_lnk_end() { return ids_link_end; }
	std::vector<int> get_ids_lnk() { return ids_link;  }
	std::vector<int> get_ids_junc() { return ids_junc; }
	void compute_area();
	double get_area(){return area;}
	double get_depth() { return depth; }
	double get_diameter() { return diameter; }
};

#endif
