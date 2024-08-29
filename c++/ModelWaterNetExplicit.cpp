#include "ModelWaterNetExplicit.h"

void ModelWaterNetExplicit::configure(int iter_stop_new, double iter_cut_thresh_new, double implicity_new, double time_step_new,
	double iter_thresh_bis_new, int iter_stop_bis_new, double thresh_left_bis_new, double thresh_right_bis_new)
{
	iter_stop = iter_stop_new;
	iter_cut_thresh = iter_cut_thresh_new;
	implicity = implicity_new;
	time_step = time_step_new;
	iter_thresh_bis = iter_thresh_bis_new;
	iter_stop_bis = iter_stop_bis_new;
	thresh_left_bis = thresh_left_bis_new;
	thresh_right_bis = thresh_right_bis_new;
	sub_time_step = 1.0; // Get this from settings.
	num_of_steps = int(time_step / sub_time_step);
}

void ModelWaterNetExplicit::run(Grid2d& grid2d, Network& network)
{
	// Imrpove the sub time step approach taken below.
	for (int i = 0; i < num_of_steps; i++) {
		preprocess(grid2d, network);
		iterate(grid2d, network);
		postprocess(grid2d, network);
	}
}

void ModelWaterNetExplicit::comp_flow_from_cell_to_junc(JuncWater &junc_water, CellWater2d &cell_water)
{
	Algorithms algorithms;
	JuncGeom* junc_geom = junc_water.get_geom();
	CellGeom2d* cell_geom = cell_water.getGeom();
	
	if (junc_geom != 0 && cell_geom != 0) {

		// Decrease cell water depth by depression storage depth.
		double cell_water_depth_tot = cell_water.getWaterDepth();
		double depress_stor = cell_water.getDeprStor();
		double cell_water_depth_act = cell_water_depth_tot - depress_stor;
		double cell_water_depth_dep = 0.0;

		// Compute amount of water stored in the depressions.
		if (cell_water_depth_tot > depress_stor) {
			cell_water_depth_dep = depress_stor;
		}
		else {
			cell_water_depth_dep = cell_water_depth_tot;
		}

		if (cell_water_depth_act < 0.0)
		{
			cell_water_depth_act = 0.0;
		}

		// Compute inflow volume.
		double inflow_vol = 0.0;
		double junc_water_depth = junc_water.get_water_depth();
		double junc_area = junc_geom->get_area();
		double cell_area = cell_geom->getArea();
		double junc_depth = junc_geom->get_depth();
		double cell_mann_n = cell_water.getMannN();

		if (junc_depth > 0.0 && cell_mann_n > 0.0) {
			// Compute flow velocity, flux and volume entering the junction.
			double cell_veloc = pow(cell_water_depth_act, 2.0 / 3.0) / cell_mann_n * sqrt(cell_geom->getAverageSlope());
			//inflow_vol = 2.0 * algorithms.get_pi() * junc_geom->get_diameter() * cell_water_depth_act * cell_veloc * sub_time_step;
			inflow_vol = sqrt(cell_geom->getArea()) * cell_water_depth_act * cell_veloc * sub_time_step;

			// Restrict inflow to existing cell water volume.
			if (inflow_vol > cell_water_depth_act * cell_area) {
				inflow_vol = cell_water_depth_act * cell_area;
			}

			// Restrict inflow to existing space in the junction.
			if (inflow_vol > (junc_depth - junc_water_depth) * junc_area) {
				inflow_vol = (junc_depth - junc_water_depth) * junc_area;
			}
		}
		// When junction depth is zero, set water depth to the same as in the cell.
		else {
			inflow_vol = (cell_water_depth_act - junc_water_depth) * junc_area;
		}

		// Update water depths in the junction and overland cell.
		if (inflow_vol > 0.0 && cell_area > 0.0 && junc_area > 0.0) {
			cell_water_depth_tot -= inflow_vol / cell_area;
			junc_water_depth += inflow_vol / junc_area;
			junc_water.set_water_depth(junc_water_depth);  // Comment out to cut exchange of water between surface
			junc_water.swap();
			cell_water.setWaterDepth(cell_water_depth_tot);  // Comment out to cut exchange of water between surface
			cell_water.swap();
		}
	}
}

void ModelWaterNetExplicit::preprocess(Grid2d& grid2d, Network& network)
{
	std::vector<JuncWater>* water_juncs = network.get_water_juncs();
	std::vector<CellWater2d>* water_cells = grid2d.get_water_cells();
	
	// Transfer water from roof cells to nearby junctions.
	for (size_t i = 0; i < water_cells->size(); i++)
	{
		int outlet_id = water_cells->at(i).getOutletIndex(); // should this information reside in geometry cell?

		if (outlet_id >= 0 && water_juncs->size()) {
			comp_flow_from_cell_to_junc(water_juncs->at(outlet_id), water_cells->at(i));
		}
	}
	
	// Transfer water from overlying cells to junctions.
	for (size_t i = 0; i < water_juncs->size(); i++)
	{
		// Compute water flow from a surface cell into a junction.
		JuncGeom* junc_geom = water_juncs->at(i).get_geom();
		int grid_conn = junc_geom->getGridConnection();

		if (grid_conn >= 0 && grid_conn < water_cells->size()) {
			comp_flow_from_cell_to_junc(water_juncs->at(i), water_cells->at(grid_conn));
		}
	}
}

double ModelWaterNetExplicit::comp_sys_water_volume(Grid2d& grid2d, Network& network)
{
	std::vector<JuncWater>* water_juncs = network.get_water_juncs();
	std::vector<LinkWater>* water_links = network.get_water_links();
	double sys_water_volume = 0.0;

	// Compute water volume in the junctions.
	for (size_t i = 0; i < water_juncs->size(); i++) {
		JuncGeom* geom = water_juncs->at(i).get_geom();
		sys_water_volume += water_juncs->at(i).get_water_depth() * geom->get_area();
	}

	// Compute water volume in the pipes.
	for (size_t i = 0; i < water_links->size(); i++) {
		double water_depth = water_links->at(i).get_water_depth();
		double filled_area = 0.0;
		double hydraulic_rad = 0.0;
		water_links->at(i).comp_flow_area_and_hydr_rad(
			water_depth,
			filled_area,
			hydraulic_rad);
		LinkGeom* geom = water_links->at(i).get_geom();
		sys_water_volume += geom->get_length() * filled_area;
	}

	return sys_water_volume;
}

void ModelWaterNetExplicit::revert_heads(Grid2d& grid2d, Network& network)
{
	std::vector<JuncWater>* water_juncs = network.get_water_juncs();
	std::vector<LinkWater>* water_links = network.get_water_links();

	for (size_t i = 0; i < water_juncs->size(); i++) {
		water_juncs->at(i).set_water_depth(water_juncs->at(i).get_water_depth_old());
	}

	for (size_t i = 0; i < water_links->size(); i++) {
		water_links->at(i).set_water_depth(water_links->at(i).get_water_depth_old());
	}
}

void ModelWaterNetExplicit::iterate(Grid2d& grid2d, Network& network)
{
	std::vector<JuncWater>* water_juncs = network.get_water_juncs();
	std::vector<LinkWater>* water_links = network.get_water_links();
	double time_loc = 0.0;
	double sub_time_step_frac;
	double error_thresh = 0.0001; // Get this constant from settings
	int iterations;
	int iterations_max = 1000; // Get this constant from settings
	double water_depth_min_thresh = -0.01; // Get this constant from settings
	double water_depth_min;
	Algorithms algorithms;
	
	// Compute water flow in the network during a sub time step.
	do { // subtime step fraction loop
		iterations = 0;
		double sys_water_volume_old = comp_sys_water_volume(grid2d, network);
		double volume_error;
		sub_time_step_frac = 1.0;

		if (sub_time_step_frac > sub_time_step) {
			sub_time_step_frac = sub_time_step;
		}

		if (sub_time_step_frac > sub_time_step - time_loc && sub_time_step - time_loc > 0.0) {
			sub_time_step_frac = sub_time_step - time_loc;
		}

		do { // error iteration loop
			// Compute water flow between junctions.
			//#pragma omp parallel for
			for (int i = 0; i < (int)water_juncs->size(); i++)
			{
				// Get current junction properties.
				JuncGeom* geom_junc = water_juncs->at(i).get_geom();

				// COMPUTE HEAD WITHIN THE CLASS/OBJECT.
				std::vector<Vertex*> vrts_junc = geom_junc->getVertPointers();
				double elev_bott_junc = vrts_junc.at(1)->z;
				double depth_water_old_junc = water_juncs->at(i).get_water_depth_old();
				double head_old_junc = elev_bott_junc + depth_water_old_junc;
				double area_junc = geom_junc->get_area();
				std::vector<int> ids_lnk_end = water_juncs->at(i).get_ids_lnk_end(); // this line seems to crash the openmp run
				std::vector<LinkWater*> links_water_neigh = water_juncs->at(i).get_links();
				
				// Get neighbour junctions and links.
				std::vector<JuncWater*> water_juncs_neigh = water_juncs->at(i).get_juncs_neigh();

				for (size_t j = 0; j < water_juncs_neigh.size(); j++)
				{
					// Get neighbour junction properties.
					JuncGeom* geom_junc_neigh = water_juncs_neigh.at(j)->get_geom();
					// COMPUTE HEAD WITHIN THE CLASS/OBJECT.
					std::vector<Vertex*> vrts_junc_neigh = geom_junc_neigh->getVertPointers();
					double elev_bott_junc_neigh = vrts_junc_neigh.at(1)->z;
					double depth_water_old_junc_neigh = water_juncs_neigh.at(j)->get_water_depth_old();
					double head_old_junc_neigh = elev_bott_junc_neigh + depth_water_old_junc_neigh;
					double area_junc_neigh = geom_junc_neigh->get_area();

					// Get link properties.
					LinkGeom* geom_link = links_water_neigh.at(j)->get_geom();
					double diam_link = geom_link->get_diameter();
					double slope_link = geom_link->get_slope();
					double length_flat_link = geom_link->get_length_flat();
					double length_link = geom_link->get_length();
					double area_link = geom_link->get_area();
					// COMPUTE HEAD WITHIN THE CLASS/OBJECT.
					Vertex centre_point_link = geom_link->getCentrePoint();
					double elev_centre_point_link = centre_point_link.z;
					double water_depth_old_link = links_water_neigh.at(j)->get_water_depth_old();
					double head_old_link = elev_centre_point_link + water_depth_old_link;
					double filled_area_old_link = 0.0;
					double hydraulic_rad_old_link = 0.0;
					links_water_neigh.at(j)->comp_flow_area_and_hydr_rad(
						water_depth_old_link,
						filled_area_old_link,
						hydraulic_rad_old_link);
					double filled_volume_old_link = length_link * filled_area_old_link;
					double full_volume_link = length_link * area_link;
					double free_volume_old_link = full_volume_link - filled_volume_old_link;
					int id_lnk_end = ids_lnk_end.at(j);
					int id_lnk_end_neigh = !id_lnk_end;
					std::vector<Vertex*> link_vrts = geom_link->getVertPointers();
					double elev_bott_link = link_vrts.at(id_lnk_end)->z;
					double elev_bott_link_neigh = link_vrts.at(id_lnk_end_neigh)->z;
					double mann_n = links_water_neigh.at(j)->get_mann_n();

					// Check following (> 0.0): length_flat_link, area_link, mann_n ...
					// Check following (>= 0.0): slope_link, hydraulic_rad

					// Unpressurized flow.
					if (free_volume_old_link > 0.0 || head_old_junc < elev_bott_link + diam_link || head_old_junc_neigh < elev_bott_link_neigh + diam_link) {
						// Flow direction from junction to link.
						if (head_old_junc > head_old_link) {
							double depth_water_old_junc_eff = head_old_junc - elev_bott_link;
							if (depth_water_old_junc_eff < 0.0) depth_water_old_junc_eff = 0.0;
							double filled_area = 0.0;
							double hydraulic_rad = 0.0;
							links_water_neigh.at(j)->comp_flow_area_and_hydr_rad(
								depth_water_old_junc_eff,
								filled_area,
								hydraulic_rad);
							double velocity = 1.0 / mann_n * sqrt(slope_link) * pow(hydraulic_rad, 2.0 / 3.0);
							double discharge = velocity * filled_area; // remove ?
							double volume = discharge * sub_time_step_frac;
							double delta_wat_depth = volume / area_junc;
							water_juncs->at(i).set_water_depth(water_juncs->at(i).get_water_depth() - delta_wat_depth);
							double depth_water_link_loc = links_water_neigh.at(j)->get_water_depth();
							double filled_area_link_loc = 0.0;
							double hydraulic_rad_link_loc = 0.0;
							links_water_neigh.at(j)->comp_flow_area_and_hydr_rad(
								depth_water_link_loc,
								filled_area_link_loc,
								hydraulic_rad_link_loc);
							double vol_water_link_loc = filled_area_link_loc * length_link;
							double free_volume_link_loc = full_volume_link - vol_water_link_loc;
							double volume_excess_loc = volume - free_volume_link_loc;

							// Save excess water to the neighbour junction.
							if (volume_excess_loc > 0.0) {
								volume = free_volume_link_loc;
								double delta_wat_depth = volume_excess_loc / geom_junc_neigh->get_area();
								water_juncs_neigh.at(j)->set_water_depth(water_juncs_neigh.at(j)->get_water_depth() + delta_wat_depth);
							}

							// Save remaining water into the link.
							filled_area_link_loc = (vol_water_link_loc + volume) / length_link;
							double depth_link_loc = links_water_neigh.at(j)->calculateWaterDepth(filled_area_link_loc, 0.5 * diam_link);
							links_water_neigh.at(j)->set_water_depth(depth_link_loc);
						}
						// Flow direction from link to junction.
						else {
							double velocity = 1.0 / mann_n * sqrt(slope_link) * pow(hydraulic_rad_old_link, 2.0 / 3.0);
							double discharge = velocity * filled_area_old_link; // remove ?
							double volume = discharge * sub_time_step_frac;
							double delta_wat_depth = volume / area_junc;
							water_juncs->at(i).set_water_depth(water_juncs->at(i).get_water_depth() + delta_wat_depth);
							double depth_water_link_loc = links_water_neigh.at(j)->get_water_depth();
							double filled_area_link_loc = 0.0;
							double hydraulic_rad_link_loc = 0.0;
							links_water_neigh.at(j)->comp_flow_area_and_hydr_rad(
								depth_water_link_loc,
								filled_area_link_loc,
								hydraulic_rad_link_loc);
							double vol_water_link_loc = filled_area_link_loc * length_link - volume;
							filled_area_link_loc = vol_water_link_loc / length_link;
							double depth_link_loc = links_water_neigh.at(j)->calculateWaterDepth(filled_area_link_loc, 0.5 * diam_link);
							links_water_neigh.at(j)->set_water_depth(depth_link_loc);
						}

					}
					// Pressurized flow. Compute flow directly between wells.
					else {
						if (head_old_junc > head_old_junc_neigh && head_old_junc > elev_bott_link + diam_link) {
							double slope = (head_old_junc - head_old_junc_neigh) / length_flat_link;
							double hydraulic_rad = algorithms.get_pi() * diam_link / area_link;
							double velocity = 1.0 / mann_n * sqrt(slope) * pow(hydraulic_rad, 2.0 / 3.0);
							double discharge = velocity * area_link; // remove ?
							double volume = discharge * sub_time_step_frac;
							double delta_wat_depth = volume / area_junc_neigh;
							water_juncs->at(i).set_water_depth(water_juncs->at(i).get_water_depth() - delta_wat_depth);
							//std::cout << "outgoing flow in " << water_juncs->at(i).getId() << "\n";
						}
						else if (head_old_junc <= head_old_junc_neigh && head_old_junc_neigh > elev_bott_link_neigh + diam_link) {
							double slope = (head_old_junc_neigh - head_old_junc) / length_flat_link;
							double hydraulic_rad = algorithms.get_pi() * diam_link / area_link;
							double velocity = 1.0 / mann_n * sqrt(slope) * pow(hydraulic_rad, 2.0 / 3.0);
							double discharge = velocity * area_link; // remove ?
							double volume = discharge * sub_time_step_frac;
							double delta_wat_depth = volume / area_junc;
							water_juncs->at(i).set_water_depth(water_juncs->at(i).get_water_depth() + delta_wat_depth);
							//std::cout << "incoming flow in " << water_juncs->at(i).getId() << "\n";
						}
					}
				}
			}

			// Assess depth and mass balance errors and decrease time step if error(s) exceed threshold value.
			double sys_water_volume = comp_sys_water_volume(grid2d, network);
			volume_error = fabs(sys_water_volume - sys_water_volume_old);
			water_depth_min = 0.0;

			for (size_t i = 0; i < water_juncs->size(); i++) {
				if (water_juncs->at(i).get_water_depth() < water_depth_min) {
					water_depth_min = water_juncs->at(i).get_water_depth();
				}
			}

			for (size_t i = 0; i < water_links->size(); i++) {
				if (water_links->at(i).get_water_depth() < water_depth_min) {
					water_depth_min = water_links->at(i).get_water_depth();
				}
			}

			if (volume_error > error_thresh || water_depth_min < water_depth_min_thresh) {
				sub_time_step_frac *= 0.5;
				revert_heads(grid2d, network);
			}

			iterations++;

		} while ((volume_error > error_thresh || water_depth_min < water_depth_min_thresh) && iterations < iterations_max);

		// Swap new water depths to old in junctions.
		for (size_t i = 0; i < water_juncs->size(); i++) {
			water_juncs->at(i).swap();
		}

		// Swap new water depths to old in links.
		for (size_t i = 0; i < water_links->size(); i++) {
			water_links->at(i).swap();
		}

		time_loc += sub_time_step_frac;

	} while (time_loc < sub_time_step);
}

void ModelWaterNetExplicit::postprocess(Grid2d& grid2d, Network& network)
{

	std::vector<JuncWater>* water_juncs = network.get_water_juncs();
	std::vector<CellWater2d>* water_cells = grid2d.get_water_cells();

	for (size_t i = 0; i < water_juncs->size(); i++)
	{
		// Remove water from outfall junction.
		JuncGeom* junc_geom = water_juncs->at(i).get_geom();
		int junc_type = water_juncs->at(i).get_type();

		if (junc_type == 1 && junc_geom != 0) {
			double junc_water_depth = water_juncs->at(i).get_water_depth();
			water_juncs->at(i).set_water_depth(0.0);
			double outfall_vol = water_juncs->at(i).get_outfall_volume();
			water_juncs->at(i).set_outfall_volume(outfall_vol + junc_water_depth * junc_geom->get_area());
			water_juncs->at(i).swap();
		}
		
		// When water level exceeds junction depth store water on surface cell.
		// Some sort of mass balance problem here. Try to sort this out.
		int grid_conn = junc_geom->getGridConnection();

		if (junc_type == 0 && junc_geom != 0 && grid_conn >= 0 && grid_conn < water_cells->size()) {
			double junc_water_depth = water_juncs->at(i).get_water_depth();
			double junc_area = junc_geom->get_area();
			double junc_depth = junc_geom->get_depth();
			double cell_water_depth = water_cells->at(grid_conn).getWaterDepth();
			CellGeom2d* geom_cell = water_cells->at(grid_conn).getGeom();

			// Update water depths in the junction and overland cell.
			if (geom_cell != 0 && junc_water_depth > junc_depth)
			{
				
				//double cell_area = geom_cell->getArea();
				//double junc_water_volume = (junc_water_depth - junc_depth) * junc_area;
				//cell_water_depth += junc_water_volume / cell_area;
				//junc_water_depth = junc_depth;
				//water_juncs->at(i).set_water_depth(junc_water_depth); // Comment out to cut exchange of water between surface
				//water_juncs->at(i).swap();
				//water_cells->at(grid_conn).setWaterDepth(cell_water_depth); // Comment out to cut exchange of water between surface
				//water_cells->at(grid_conn).swap();
				
				double cell_area = geom_cell->getArea();
				double cell_water_vol = cell_water_depth * cell_area;
				//double excess_water = junc_water_depth * junc_area;
				double excess_water = (junc_water_depth - junc_depth) * junc_area;
				double cell_water_vol_new = cell_water_vol + excess_water;
				//junc_water_depth = 0.0;
				junc_water_depth = junc_depth;
				water_juncs->at(i).set_water_depth(junc_water_depth); // Comment out to cut exchange of water between surface
				water_juncs->at(i).swap();
				water_cells->at(grid_conn).setWaterDepth(cell_water_vol_new / cell_area);
				water_cells->at(grid_conn).swap();
			}
		}
		
	}
}
