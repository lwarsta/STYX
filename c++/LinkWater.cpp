#include "LinkWater.h"

LinkWater::LinkWater()
{
    geom = 0;
	mann_n = 0.013; // temporary for concrete pipe 0.013
}

void LinkWater::comp_flow_area_and_hydr_rad(
	double water_depth, 
	double& flow_area,
	double& hydr_rad) {
	if (geom != 0) {
		// Compute trivial cases when water level is below or above the link.
		double radius = 0.5 * geom->get_diameter();

		if (radius <= 0.0 || water_depth <= 0.0) {
			flow_area = 0.0;
			hydr_rad = 0.0;
		}
		else if (water_depth >= 2.0 * radius) {
			flow_area = algorithms.get_pi() * radius * radius;
			hydr_rad = algorithms.get_pi() * 2.0 * radius;
		}
		// Compute flow crossectional area and hydraulic radius
		// when water depth is lower than link diameter but
		// higher than zero.
		else {
			double theta = 2.0 * acos((radius - water_depth) / radius);
			flow_area = 0.5 * radius * radius * (theta - sin(theta));
			double perimeter = radius * theta;
			hydr_rad = flow_area / perimeter;
		}
	}
	else {
		flow_area = 0.0;
		hydr_rad = 0.0;
	}
}