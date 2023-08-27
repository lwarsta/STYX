#include "LinkGeom.h"

LinkGeom::LinkGeom(int id_new, char geom_type_new, std::vector<int> vert_indices_new,
	               std::vector<Vertex*> vert_pointers_new)
	: LinkJuncGeomBase(id_new, geom_type_new, vert_indices_new, vert_pointers_new) {
	// Initialise values - TEMPORARY.
	Algorithms algorithms;
	diameter = 0.1;
	area = algorithms.get_pi() * 0.25 * diameter * diameter;
	Vertex vec_link = algorithms.create_vector(*vertPointers.at(0), *vertPointers.at(1));
	length = algorithms.compute_vector_length(vec_link);
}