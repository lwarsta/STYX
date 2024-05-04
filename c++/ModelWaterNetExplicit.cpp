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
	sub_time_step = 1.0; // Get this to settings
	flow_vel_max = 1.0; // Get this to settings
}

void ModelWaterNetExplicit::run(Grid2d& grid2d, Network& network)
{
	// Imrpove the sub time step approach taken below.
	for (int i = 0; i < (int)time_step; i++) {
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
			inflow_vol = 2.0 * algorithms.get_pi() * junc_geom->get_diameter() * cell_water_depth_act * cell_veloc * sub_time_step;

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

void ModelWaterNetExplicit::iterate(Grid2d& grid2d, Network& network)
{
	std::vector<JuncWater>* water_juncs = network.get_water_juncs();
	double time_loc = 0.0;

	// Compute water volume in the system.
	// Is this the best place for this?
	double sys_water_volume = 0.0;

	for (int i = 0; i < water_juncs->size(); i++) {
		JuncGeom* junc_geom = water_juncs->at(i).get_geom();
		sys_water_volume += water_juncs->at(i).get_water_depth_old() * junc_geom->get_area();
	}

	do {
		std::vector <std::vector<double>> fluxes;
		std::vector <std::vector<double>> velocities;
		std::vector <std::vector<double>> distances;
		fluxes.resize(water_juncs->size());
		velocities.resize(water_juncs->size());
		distances.resize(water_juncs->size());

		//#pragma omp parallel for
		for (int i = 0; i < water_juncs->size(); i++)
		{
			// Calculate water volume in the junction.
			double depth_water = water_juncs->at(i).get_water_depth_old();
			JuncGeom* junc_geom = water_juncs->at(i).get_geom();
			double volume_water = depth_water * junc_geom->get_area();
			double vol_water_change = 0.0;

			// Get neighbours and initialise variables.
			std::vector<JuncWater*> water_juncs_neigh = water_juncs->at(i).get_juncs_neigh();
			std::vector<int> ids_lnk_end = water_juncs->at(i).get_ids_lnk_end();
			std::vector<LinkWater*> links = water_juncs->at(i).get_links();

			// Initialise local inter cell variables.
			fluxes.at(i).assign(water_juncs_neigh.size(), 0.0);
			velocities.at(i).assign(water_juncs_neigh.size(), 0.0);
			distances.at(i).assign(water_juncs_neigh.size(), 0.0);

			// Water depth cannot decrease below zero.
			// This was used to take into account depression storage.
			// Is it still needed?
			if (depth_water < 0.0)
			{
				depth_water = 0.0;
			}

			// Hydraulic head in the junction.
			// Elevation below depicts the bottom of the junction.
			std::vector<Vertex*> junc_vrts = junc_geom->getVertPointers();
			double elevation = junc_vrts.at(1)->z;
			double head = elevation + depth_water;

			for (size_t j = 0; j < water_juncs_neigh.size(); j++)
			{
				// Compute neighbour junction elevation.
				// Elevation below depicts the bottom of the junction.
				JuncWater* junc_water_neigh = water_juncs_neigh.at(j);
				JuncGeom* junc_geom_neigh = junc_water_neigh->get_geom();
				std::vector<Vertex*> junc_vrts_neigh = junc_geom_neigh->getVertPointers();
				double elevation_neigh = junc_vrts_neigh.at(1)->z;

				// Water depth cannot decrease below zero.
				// This was used to take into account depression storage.
				// Is it still needed?
				double depth_water_neigh = junc_water_neigh->get_water_depth_old();

				if (depth_water_neigh < 0.0)
				{
					depth_water_neigh = 0.0;
				}

				double head_neigh = elevation_neigh + depth_water_neigh;

				// Compute hydraulic slope between the junctions.
				Vertex junc_vrt = *junc_vrts.at(1);
				junc_vrt.z = 0.0;
				Vertex junc_vrt_neigh = *junc_vrts_neigh.at(1);
				junc_vrt_neigh.z = 0.0;
				Algorithms algorithms;
				Vertex junc_vec = algorithms.create_vector(junc_vrt, junc_vrt_neigh);
				double distance = algorithms.compute_vector_length(junc_vec);
				double slope = 0.0;
				distances.at(i).at(j) = distance;

				if (distance > 0.0) {
					slope = (head - head_neigh) / distance;
				}

				// Compute flow direction, velocity and flux in a link.
				LinkGeom* link_geom = links.at(j)->get_geom();
				double slope_link = link_geom->get_slope();
				std::vector<Vertex*> link_vrts = link_geom->getVertPointers();
				double filled_area = 0.0;
				double velocity = 0.0;
				double flux = 0.0;
				double mann_n = links.at(j)->get_mann_n();

				if (head > head_neigh) {
					int id_lnk_end = ids_lnk_end.at(j);
					double elevation_link = link_vrts.at(id_lnk_end)->z;

					if (head > elevation_link) {
						double water_depth_lnk = head - elevation_link;

						if (water_depth_lnk > depth_water) {
							water_depth_lnk = depth_water;
						}

						double filled_area = 0.0;
						double hydraulic_rad = 0.0;
						links.at(j)->comp_flow_area_and_hydr_rad(water_depth_lnk, filled_area, hydraulic_rad);
						//double velocity = 1.0 / mann_n * sqrt(slope) * pow(hydraulic_rad, 2.0 / 3.0);
						double velocity = 1.0 / mann_n * sqrt(slope_link) * pow(hydraulic_rad, 2.0 / 3.0);
						
						// Currently flow velocity in a pipe is restricted to a given maximum velocity.
						// This used to ensure stable computation. try to fix this later.
						//if (velocity > flow_vel_max) {
						//	velocity = flow_vel_max;
						//}
						//else if (velocity < -flow_vel_max) {
						//	velocity = -flow_vel_max;
						//}

						velocities.at(i).at(j) = velocity;
						fluxes.at(i).at(j) = filled_area * velocity;
					}
				}
				else if (head <= head_neigh) {
					int id_lnk_end = 0;

					if (ids_lnk_end.at(j) == 0) {
						id_lnk_end = 1;
					}

					double elevation_link = link_vrts.at(id_lnk_end)->z;

					if (head_neigh > elevation_link) {
						double water_depth_lnk = head_neigh - elevation_link;

						if (water_depth_lnk > depth_water_neigh) {
							water_depth_lnk = depth_water_neigh;
						}

						double filled_area = 0.0;
						double hydraulic_rad = 0.0;
						links.at(j)->comp_flow_area_and_hydr_rad(water_depth_lnk, filled_area, hydraulic_rad);
						//double velocity = -1.0 / mann_n * sqrt(-slope) * pow(hydraulic_rad, 2.0 / 3.0);
						double velocity = -1.0 / mann_n * sqrt(slope_link) * pow(hydraulic_rad, 2.0 / 3.0);
						
						// Currently flow velocity in a pipe is restricted to a given maximum velocity.
						// This used to ensure stable computation. try to fix this later.
						//if (velocity > flow_vel_max) {
						//	velocity = flow_vel_max;
						//}
						//else if (velocity < -flow_vel_max) {
						//	velocity = -flow_vel_max;
						//}

						velocities.at(i).at(j) = velocity;
						fluxes.at(i).at(j) = filled_area * velocity;
					}
				}
			}
		}

		// Assess CFL condition and time step length.
		//u* dt / dx <= 1 -> dt = dx / u
		double time_step_sub_new = std::numeric_limits<double>::max();

		for (int i = 0; i < velocities.size(); i++) {
			for (int j = 0; j < velocities.at(i).size(); j++) {
				if (fabs(velocities.at(i).at(j)) > 0.0 &&
					0.075 * distances.at(i).at(j) / fabs(velocities.at(i).at(j)) < time_step_sub_new) {
					time_step_sub_new = 0.075 * distances.at(i).at(j) / fabs(velocities.at(i).at(j));
				}
			}
		}

		if (time_step_sub_new > sub_time_step ) {
			time_step_sub_new = sub_time_step ;
		}

		if (time_step_sub_new > sub_time_step - time_loc && sub_time_step - time_loc > 0.0) {
			time_step_sub_new = sub_time_step - time_loc;
		}
		
		// Compute new water depths in the wells.
		double volume_error = 0.0;
		double error_thresh = 0.000001; // Get this constant from settings
		int iterations_max = 1000; // Get this constant from settings
		int iterations = 0;

		do {
			// Compute change in water volumes.
			double sys_water_volume_new = 0.0;

			for (int i = 0; i < water_juncs->size(); i++) {
				double vol_water_change = 0.0;

				for (int j = 0; j < fluxes.at(i).size(); j++) {
					vol_water_change -= fluxes.at(i).at(j) * time_step_sub_new;
				}

				JuncGeom* junc_geom = water_juncs->at(i).get_geom();
				double depth_water_old = water_juncs->at(i).get_water_depth_old();
				water_juncs->at(i).set_water_depth(depth_water_old + vol_water_change / junc_geom->get_area());
				sys_water_volume_new += water_juncs->at(i).get_water_depth() * junc_geom->get_area();
			}

			// Decrease time step to decrease error.
			volume_error = fabs(sys_water_volume_new - sys_water_volume);

			if (volume_error > error_thresh) {
				time_step_sub_new *= 0.5;
			}

			iterations++;

		} while (volume_error > error_thresh && iterations < iterations_max);

		// Swap new water depths to old water depths.
		for (int i = 0; i < water_juncs->size(); i++) {
			water_juncs->at(i).swap();
		}
				
		time_loc += time_step_sub_new;

	} while (time_loc < time_step);
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
		/*
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
				double cell_area = geom_cell->getArea();
				double junc_water_volume = (junc_water_depth - junc_depth) * junc_area;
				cell_water_depth += junc_water_volume / cell_area;
				junc_water_depth = junc_depth;
				water_juncs->at(i).set_water_depth(junc_water_depth); // Comment out to cut exchange of water between surface
				water_juncs->at(i).swap();
				water_cells->at(grid_conn).setWaterDepth(cell_water_depth); // Comment out to cut exchange of water between surface
				water_cells->at(grid_conn).swap();
			}
		}
		*/
	}
}
