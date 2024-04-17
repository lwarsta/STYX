#include "JuncWater.h"

JuncWater::JuncWater() {
    geom = 0;
	water_depth = 0.0;
	water_depth_old = 0.0;
	type = -1;
	outfall_volume = 0.0;
}

void JuncWater::save_link(int id_lnk_end, int id_lnk, LinkWater* link_water, 
	                      int id_jnc, JuncWater* junc_water) {
	ids_link_end.push_back(id_lnk_end);
	ids_link.push_back(id_lnk); // is this needed - id can be queried?
	links.push_back(link_water);
	ids_junc.push_back(id_jnc); // is this needed - id can be queried?
	juncs.push_back(junc_water);
}
