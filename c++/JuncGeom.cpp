#include "JuncGeom.h"

// Superclass constructor is also called below.
JuncGeom::JuncGeom(int id_new, double pos_x, double pos_y, double pos_z,
	               double diameter_new, bool open_new, double depth_new)
{
	id = id_new;
	coord = Vertex(id_new, pos_x, pos_y, pos_z);
	diameter = diameter_new;
	area = 3.14159265 * 0.25 * diameter * diameter;
	open = open_new;
	depth = depth_new;
}

void JuncGeom::compute_area()
{
	
}
