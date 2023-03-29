#include "LinkGeom.h"

LinkGeom::LinkGeom(int id_new, double diameter_new, double roughness_new, 
	               int conn_ind_00, int conn_ind_01, int conn_type_00, 
	               int conn_type_01, double conn_elev_00, double conn_elev_01)
{
	id = id_new;
	diameter = diameter_new;
	area = 3.14159265 * 0.25 * diameter * diameter;
	roughness = roughness_new;
	conn_ind[0] = conn_ind_00;
    conn_ind[1] = conn_ind_01;
	conn_type[0] = conn_type_00;
	conn_type[1] = conn_type_01;
	conn_elev[0] = conn_type_00;
	conn_elev[1] = conn_type_01;
}