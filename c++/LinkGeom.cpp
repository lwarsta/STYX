#include "LinkGeom.h"

LinkGeom::LinkGeom(int id_new, char geom_type_new, std::vector<int> vert_indices_new,
	               std::vector<Vertex*> vert_pointers_new)
	: LinkJuncGeomBase(id_new, geom_type_new, vert_indices_new, vert_pointers_new) {
	diameter = 0.0;
	area = 0.0;
	//length = 0.0;
	slope = 0.0;
}

void LinkGeom::comp_geom_properties() {
	Algorithms algorithms;
	area = algorithms.get_pi() * 0.25 * diameter * diameter;
	//Vertex vec_link = algorithms.create_vector(*vertPointers.at(0), *vertPointers.at(1));
	//length = algorithms.compute_vector_length(vec_link);
	if (vertPointers.size() == 2) {
		Vertex v1 = *vertPointers.at(0);
		Vertex v2 = *vertPointers.at(1);
		double dz = fabs(v1.z - v2.z);
		v1.z = 0.0;
		v2.z = 0.0;
		Vertex vec_link = algorithms.create_vector(v1, v2);
		double length_flat = algorithms.compute_vector_length(vec_link);

		if (length_flat > 0.0)
		{
			slope = dz / length_flat;
		}
	}
}