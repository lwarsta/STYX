#ifndef _ALGORITHMS_H
#define _ALGORITHMS_H
#include <vector>
#include "Vertex.h"

class Algorithms {
private:
	const double pi = 3.14159265359;
public:
	double get_pi() { return pi; }
	Vertex create_vector(Vertex point0, Vertex point1);
	double compute_vector_length(Vertex vector);
	Vertex compute_centre_point(std::vector<Vertex> points);
	double compute_vector_dot_product(Vertex point0, Vertex point1);
	Vertex add_vectors(Vertex point0, Vertex point1);
	Vertex compute_unit_vector(Vertex vector);
	Vertex compute_vector_cross_product(Vertex vector0, Vertex vector1);
};

#endif

