#include "ModelWaterNetDiffBrute.h"

void ModelWaterNetDiffBrute::configure(int iter_stop_new, double iter_cut_thresh_new, double implicity_new, double time_step_new,
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
}

void ModelWaterNetDiffBrute::run(Grid2d &grid2d, Network &network)
{
	preprocess(grid2d, network);
    iterate(grid2d, network);
    postprocess(grid2d, network);
}

void ModelWaterNetDiffBrute::preprocess(Grid2d &grid2d, Network &network)
{
    // Calculate junction properties.
	std::vector<JuncWater>* water_juncs = network.get_water_juncs();
	std::vector<CellWater2d>* water_cells = grid2d.get_water_cells();
	Algorithms algorithms;

	//#pragma omp parallel for
	for (size_t i = 0; i < water_juncs->size(); i++)
	{
		if (water_juncs->at(i).get_water_depth() < 0.0)
		{
			water_juncs->at(i).set_water_depth(0.0);
		}
		/*
		// Compute water volume flowing from a surface cell to a junction.
		JuncGeom* junc_geom = water_juncs->at(i).get_geom();
		int grid_conn = junc_geom->getGridConnection();

		if (junc_geom != 0 && grid_conn >= 0 && grid_conn < water_cells->size()) {
			// Decrease overland water depth by flow threshold depth.
			double cell_water_depth = water_cells->at(grid_conn).getWaterDepth();
			cell_water_depth -= water_cells->at(grid_conn).getDeprStor();

			if (cell_water_depth < 0.0)
			{
				cell_water_depth = 0.0;
			}

			// Compute flow velocity, flux and volume entering the junction.
			double cell_mann_n = water_cells->at(grid_conn).getMannN();
			CellGeom2d* geom_cell = water_cells->at(grid_conn).getGeom();

			if (cell_mann_n > 0.0 && geom_cell != 0) {
				double cell_slope = geom_cell->getAverageSlope();
				double junc_diam = junc_geom->get_diameter();
				double cell_veloc = pow(cell_water_depth, 2.0 / 3.0) / cell_mann_n *
									sqrt(cell_slope);
				double inflow_vol = 2.0 * algorithms.get_pi() * junc_diam *
									cell_water_depth * cell_veloc * time_step;

				// Restrict inflow to existing water volume in the cell.
				if (inflow_vol > cell_water_depth * geom_cell->getArea()) {
					inflow_vol = cell_water_depth * geom_cell->getArea();
				}

				// Restrict inflow to existing space in the junction.
				double junc_water_depth = water_juncs->at(i).get_water_depth();
				double junc_area = junc_geom->get_area();
				double junc_depth = junc_geom->get_depth();
				
				if (inflow_vol > (junc_depth - junc_water_depth) * junc_area) {
					inflow_vol = (junc_depth - junc_water_depth) * junc_area;
				}

				if (inflow_vol < 0.0) {
					inflow_vol = 0.0;
				}
				
				// Update water depths in the junction and overland cell.
				if (geom_cell->getArea() > 0.0 && junc_area > 0.0) {
					cell_water_depth = cell_water_depth - 
										inflow_vol / geom_cell->getArea();
					junc_water_depth = junc_water_depth + inflow_vol / junc_area;
					water_juncs->at(i).set_water_depth(junc_water_depth);
					water_cells->at(grid_conn).setWaterDepth(cell_water_depth);
				}
			}
		}
		*/

		// Compute water flow from a surface cell into a junction.
		JuncGeom* junc_geom = water_juncs->at(i).get_geom();
		int grid_conn = junc_geom->getGridConnection();

		if (junc_geom != 0 && grid_conn >= 0 && grid_conn < water_cells->size()) {
			CellGeom2d* geom_cell = water_cells->at(grid_conn).getGeom();

			if (geom_cell != 0) {
				// Decrease cell water depth by depression storage depth.
				double cell_water_depth = water_cells->at(grid_conn).getWaterDepth();
				cell_water_depth -= water_cells->at(grid_conn).getDeprStor();

				if (cell_water_depth < 0.0)
				{
					cell_water_depth = 0.0;
				}

				// Compute inflow volume.
				double inflow_vol = 0.0;
				double junc_water_depth = water_juncs->at(i).get_water_depth();
				double junc_area = junc_geom->get_area();
				double junc_depth = junc_geom->get_depth();
				double cell_mann_n = water_cells->at(grid_conn).getMannN();
				
				if (junc_depth > 0.0 && cell_mann_n > 0.0) {
					// Compute flow velocity, flux and volume entering the junction.
					double cell_slope = geom_cell->getAverageSlope();
					double junc_diam = junc_geom->get_diameter();
					double cell_veloc = pow(cell_water_depth, 2.0 / 3.0) / cell_mann_n *
						sqrt(cell_slope);
					inflow_vol = 2.0 * algorithms.get_pi() * junc_diam *
									cell_water_depth * cell_veloc * time_step;

					// Restrict inflow to existing cell water volume.
					if (inflow_vol > cell_water_depth * geom_cell->getArea()) {
						inflow_vol = cell_water_depth * geom_cell->getArea();
					}

					// Restrict inflow to existing space in the junction.
					if (inflow_vol > (junc_depth - junc_water_depth) * junc_area) {
						inflow_vol = (junc_depth - junc_water_depth) * junc_area;
					}

					// Do not allow negative inflow.
					if (inflow_vol < 0.0) {
						inflow_vol = 0.0;
					}
				}
				// When junction depth is zero, set water depth to the same as in the cell.
				else {
					inflow_vol = (cell_water_depth - junc_water_depth) * junc_area;
				}

				// Update water depths in the junction and overland cell.
				if (geom_cell->getArea() > 0.0 && junc_area > 0.0) {
					cell_water_depth = cell_water_depth -
						inflow_vol / geom_cell->getArea();
					junc_water_depth = junc_water_depth + inflow_vol / junc_area;
					water_juncs->at(i).set_water_depth(junc_water_depth);
					water_cells->at(grid_conn).setWaterDepth(cell_water_depth);
				}
			}
		}

		// Swap new junction water depth to old water depth variable.
		water_juncs->at(i).swap();
	}
}

double ModelWaterNetDiffBrute::calc_residual(JuncWater& junc_water, double depth_water)
{	
	// Calculate water volume in the junction.
	JuncGeom * junc_geom = junc_water.get_geom();
	double volume_water = (depth_water - junc_water.get_water_depth_old()) *
		                  junc_geom->get_area();
	
	// Get neighbours and initialise variables.
	std::vector<JuncWater*> juncs_water_neigh = junc_water.get_juncs_neigh();
	std::vector<int> ids_lnk_end = junc_water.get_ids_lnk_end();
	std::vector<LinkWater*> links = junc_water.get_links();
	
	// Initialise local inter cell variables.
	std::vector<double> fluxes;
	std::vector<double> velocities;
	fluxes.assign(juncs_water_neigh.size(), 0.0);
	velocities.assign(juncs_water_neigh.size(), 0.0);

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
	
	for (size_t i = 0; i < juncs_water_neigh.size(); i++)
	{
		// Compute neighbour junction elevation.
		// Elevation below depicts the bottom of the junction.
		JuncWater* junc_water_neigh = juncs_water_neigh.at(i);
		JuncGeom * junc_geom_neigh = junc_water_neigh->get_geom();
		std::vector<Vertex*> junc_vrts_neigh = junc_geom_neigh->getVertPointers();
		double elevation_neigh = junc_vrts_neigh.at(1)->z;

		// Water depth cannot decrease below zero.
		// This was used to take into account depression storage.
		// Is it still needed?
		double depth_water_neigh = junc_water_neigh->get_water_depth();
		
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

		if (distance > 0.0) {
			slope = (head - head_neigh) / distance;
		}
		
		// Compute flow direction, velocity and flux in a link.
		LinkGeom* link_geom = links.at(i)->get_geom();
		std::vector<Vertex*> link_vrts = link_geom->getVertPointers();
		double filled_area = 0.0;
		double velocity = 0.0;
		double flux = 0.0;
		double mann_n = links.at(i)->get_mann_n();

		if (head > head_neigh) // the structure below could be improved ...
		{
			// Select link end.
			int id_lnk_end = ids_lnk_end.at(i);
			double elevation_link = link_vrts.at(id_lnk_end)->z;
			
			if (head > elevation_link)
			{
				double water_depth_lnk = head - elevation_link;

				if (water_depth_lnk > depth_water)
				{
					water_depth_lnk = depth_water;
				}

				double filled_area = 0.0;
				double hydraulic_rad = 0.0;
				links.at(i)->comp_flow_area_and_hydr_rad(water_depth_lnk, 
					                                     filled_area,
					                                     hydraulic_rad);
				double velocity = 1.0 / mann_n * sqrt(slope) * 
					              pow(hydraulic_rad, 2.0 / 3.0);
				double flux = filled_area * velocity;
				volume_water += flux * time_step;
				velocities.at(i) = velocity;
				fluxes.at(i) = flux;
			}
		}
		else if (head <= head_neigh) // equality was added?
		{
			// Select link end.
			int id_lnk_end = 0;

			if (ids_lnk_end.at(i) == 0)
			{
				id_lnk_end = 1;
			}
			
			double elevation_link = link_vrts.at(id_lnk_end)->z;
			
			if (head_neigh > elevation_link)
			{
				double water_depth_lnk = head_neigh - elevation_link;

				if (water_depth_lnk > depth_water_neigh)
				{
					water_depth_lnk = depth_water_neigh;
				}

				double filled_area = 0.0;
				double hydraulic_rad = 0.0;
				links.at(i)->comp_flow_area_and_hydr_rad(water_depth_lnk, 
					                                     filled_area,
					                                     hydraulic_rad);
				double velocity = -1.0 / mann_n * sqrt(-slope) * 
					              pow(hydraulic_rad, 2.0 / 3.0);
				double flux = filled_area * velocity;
				volume_water += flux * time_step;
				velocities.at(i) = velocity;
				fluxes.at(i) = flux;
			}
		}
	}

	// Save velocities and fluxes for inspection.
	junc_water.assign_velocities(velocities);
	junc_water.assign_fluxes(fluxes);
	
	return volume_water;
}

double ModelWaterNetDiffBrute::bisection(JuncWater& water, double left, double right)
{
	int iter_count_bis = 0;

	while (fabs(right - left) > 2.0 * iter_thresh_bis && iter_count_bis < iter_stop_bis)
	{
		// Calculate midpoint of domain.
		double midpoint = (right + left) * 0.5;

		// Find f(midpoint).
		if (calc_residual(water, left) * calc_residual(water, midpoint) > 0.0)
		{
			// Throw away left half.
			left = midpoint;
		}
		else
		{
			// Throw away right half.
			right = midpoint;
		}

		iter_count_bis++;
	}

	return (right + left) * 0.5;
}

void ModelWaterNetDiffBrute::iterate(Grid2d &grid2d, Network &network)
{
	// The main iteration loop.
    std::vector<JuncWater>* water_juncs = network.get_water_juncs();
	double iter_dev;
	int iter_count = 0;

	do
    {
		iter_dev = 0.0;

		#pragma omp parallel for
		for (int i = 0; i < water_juncs->size(); i++)
        {
			// Compute new water depth.
            double water_depth_new = bisection(water_juncs->at(i), thresh_left_bis, thresh_right_bis);
            // Save iteration deviation.
			double iterDevNew = fabs(water_juncs->at(i).get_water_depth() - water_depth_new);
			
			// OpenMP double-checked locking.
			if (iterDevNew > iter_dev)
            {
				#pragma omp critical
				if (iterDevNew > iter_dev) iter_dev = iterDevNew;
            }
            
			// Save the new water depth.
			water_juncs->at(i).set_water_depth(water_depth_new);
		}

		iter_count++;

	} while (iter_dev > iter_cut_thresh && iter_count < iter_stop);
}

void ModelWaterNetDiffBrute::postprocess(Grid2d &grid2d, Network &network)
{
	
	std::vector<JuncWater>* water_juncs = network.get_water_juncs();
	std::vector<CellWater2d>* water_cells = grid2d.get_water_cells();

	for (size_t i = 0; i < water_juncs->size(); i++)
	{
		// When water level exceeds junction depth store water on surface cell.
		JuncGeom* junc_geom = water_juncs->at(i).get_geom();
		int grid_conn = junc_geom->getGridConnection();
		
		if (junc_geom != 0 && grid_conn >= 0 && grid_conn < water_cells->size()) {
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
				water_juncs->at(i).set_water_depth(junc_water_depth);
				water_cells->at(grid_conn).setWaterDepth(cell_water_depth);
			}
		}
	}
	/*
	std::cout << "Water volume in the surface cells: " << std::endl;
	double water_vol_cells = 0.0;

	for (size_t i = 0; i < water_cells->size(); i++)
	{
		CellGeom2d* geom_cell = water_cells->at(i).getGeom();
		water_vol_cells += water_cells->at(i).getWaterDepth() * geom_cell->getArea();
	}

	std::cout << water_vol_cells << std::endl;
	
	std::cout << "Printing junction geometry: " << std::endl;

	for (size_t i = 0; i < water_juncs->size(); i++)
	{
		JuncGeom* junc_geom = water_juncs->at(i).get_geom();
		std::vector<Vertex*> junc_vrts = junc_geom->getVertPointers();
		double elevation = junc_vrts.at(1)->z; // bottom
		std::cout << i << "x: " << junc_vrts.at(1)->x << "y: " << junc_vrts.at(1)->y << "z: " << junc_vrts.at(1)->z << std::endl;
	}

	std::cout << "Printing water volumes: " << std::endl;

	for (size_t i = 0; i < water_juncs->size(); i++)
	{
		JuncGeom* junc_geom = water_juncs->at(i).get_geom();
		double junc_area = junc_geom->get_area();
		std::cout << i << ": " << water_juncs->at(i).get_water_depth() * junc_area << std::endl;
	}

	std::cout << "Printing water depths: " << std::endl;

	for (size_t i = 0; i < water_juncs->size(); i++)
	{
		std::cout << i << ": " << water_juncs->at(i).get_water_depth() << std::endl;
	}

	std::cout << "Printing velocities: " << std::endl;

	for (size_t i = 0; i < water_juncs->size(); i++)
	{
		std::cout << i << ": ";
		std::vector<double> velocities = water_juncs->at(i).get_velocities();

		for (size_t j = 0; j < velocities.size(); j++) {
			std::cout << velocities.at(j) << ", ";
		}
		std::cout << std::endl;
	}

	std::cout << "Printing fluxes: " << std::endl;

	for (size_t i = 0; i < water_juncs->size(); i++)
	{
		std::cout << i << ": ";
		std::vector<double> fluxes = water_juncs->at(i).get_fluxes();

		for (size_t j = 0; j < fluxes.size(); j++) {
			std::cout << fluxes.at(j) << ", ";
		}
		std::cout << std::endl;
	}
	*/
}
