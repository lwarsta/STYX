#include "Algorithms.h"

Vertex Algorithms::create_vector(Vertex point0, Vertex point1)
{
    Vertex vector;
    vector.x = point0.x - point1.x;
    vector.y = point0.y - point1.y;
    vector.z = point0.z - point1.z;

    return vector;
}

double Algorithms::compute_vector_length(Vertex vector)
{
    return sqrt(vector.x * vector.x + vector.y * vector.y + vector.z * vector.z);
}

Vertex Algorithms::compute_centre_point(std::vector<Vertex> points)
{
	Vertex centre_point(0, 0.0, 0.0, 0.0);

	for (size_t i = 0; i < points.size(); i++) {
		centre_point.x += points.at(i).x;
		centre_point.y += points.at(i).y;
		centre_point.z += points.at(i).z;
	}

	if (points.size() > 0) {
		centre_point.x /= points.size();
		centre_point.y /= points.size();
		centre_point.z /= points.size();
	}

	return centre_point;
}

double Algorithms::compute_vector_dot_product(Vertex point0, Vertex point1)
{
	return point0.x * point1.x + point0.y * point1.y + point0.z * point1.z;
}

Vertex Algorithms::add_vectors(Vertex point0, Vertex point1)
{
	point0.x = point0.x + point1.x;
	point0.y = point0.y + point1.y;
	point0.z = point0.z + point1.z;

	return point0;
}

Vertex Algorithms::compute_unit_vector(Vertex vector)
{
	double length = compute_vector_length(vector);

	if (length > 0.0)
	{
		vector.x = vector.x / length;
		vector.y = vector.y / length;
		vector.z = vector.z / length;
		return vector;
	}

	// When a zero length vector is sent in, return the same vector.
	else
	{
		return vector;
	}
}

Vertex Algorithms::compute_vector_cross_product(Vertex vector0,
	Vertex vector1)
{
	Vertex cross_product;
	cross_product.x = vector0.y * vector1.z - vector0.z * vector1.y;
	cross_product.y = vector0.z * vector1.x - vector0.x * vector1.z;
	cross_product.z = vector0.x * vector1.y - vector0.y * vector1.x;

	return cross_product;
}