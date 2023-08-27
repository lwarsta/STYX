#include "JuncGeom.h"

// Superclass constructor is also called below.
JuncGeom::JuncGeom(int idNew, char geomTypeNew, std::vector<int> vertIndicesNew,
	std::vector<Vertex*> vertPointersNew) 
	: LinkJuncGeomBase(idNew, geomTypeNew, vertIndicesNew, vertPointersNew) {
	diameter = 0.5;
	area = algorithms.get_pi() * 0.25 * diameter * diameter;
	open = 1;
	
	// Compute junction depth from top and bottom vertices.
	if (vertPointers.size() == 2) {
		depth = vertPointers.at(0)->z - vertPointers.at(1)->z;
	}
	else {
		depth = 0.0;
	}
}

void JuncGeom::save_link(int id_lnk_end, int id_lnk, LinkGeom* link_geom, int id_jnc, JuncGeom* junc_geom)
{
	ids_link_end.push_back(id_lnk_end);
	ids_link.push_back(id_lnk); // is this needed - id can be queried?
	links.push_back(link_geom);
	ids_junc.push_back(id_jnc); // is this needed - id can be queried?
	juncs.push_back(junc_geom);
}

void JuncGeom::compute_area()
{
	
}
