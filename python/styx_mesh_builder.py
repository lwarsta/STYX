#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import glob
import configparser
import json
import csv
import math
import rasterio
from functools import reduce
import datetime
import math
import shapely
from shapely.geometry import Point, Polygon, shape, mapping 
import fiona
from rtree import index
import numpy as np
import gmsh
import time
from copy import copy

class Vertex:
    def __init__(self, id, x, y, z):
        """
        """
        self.id = id
        self.x = x
        self.y = y
        self.z = z

class Cell:
    def __init__(self, geom_type, vert_indices, id, grid_connection):
        """
        """
        self.geom_type = geom_type
        self.vert_indices = vert_indices
        self.cp = Vertex(0, 0, 0, 0)
        self.id = id
        self.grid_connection = grid_connection
        self.mask = 0
        self.material = 0
        self.init_cond = id
        self.bound_cond = id
        #self.soiltype = 0        
        self.hydraulic_head = 0.0 # TEMPORARY ?
        self.junction_and_outlet = -1 # TEMPORARY

def avg(lst):
    """
    """
    return reduce(lambda a, b: a + b, lst) / len(lst)

def comp_cell_center_point(cells, vertices):
    """
    """
    for cell in cells:
        vert_x_lst = []
        vert_y_lst = []
        vert_z_lst = []
        for vert_ind in cell.vert_indices:
            vert_x_lst.append( vertices[vert_ind].x )
            vert_y_lst.append( vertices[vert_ind].y )
            vert_z_lst.append( vertices[vert_ind].z )
        cell.cp = Vertex(0, avg(vert_x_lst), avg(vert_y_lst), avg(vert_z_lst))

def load_raster_folder(raster_paths):
    """
    """
    raster_sources = []
    pixel_sizes = []
    for raster_path in raster_paths:
        src = rasterio.open(raster_path, 'r')
        raster_sources.append(src)
        gt = src.transform
        pixel_sizes.append( [abs(gt[0]), abs(gt[4])] )
    return raster_sources, pixel_sizes

def pick_raster_pixel_values(vertices, sources):
    """
    """
    
    # Find the rasters that cover the points.
    src_lst = []
    for vert_ind, vertex in enumerate(vertices):
        src_ind = -1
        for ind, src in enumerate(sources):
            # Last raster is currently selected.
            if vertex.x > src.bounds[0] and \
               vertex.y > src.bounds[1] and \
               vertex.x < src.bounds[2] and \
               vertex.y < src.bounds[3]:
                src_ind = ind
        src_lst.append(src_ind)
    # Get pixel values under the points.
    pixel_values = []
    for vertex, src_ind in zip(vertices, src_lst):
        if src_ind != -1:
            # Is this the best way to do this?
            for val_ind, val in enumerate(sources[src_ind].sample( [(vertex.x, vertex.y)] )):
                if len(val) > 0:
                    pixel_values.append(val[0])
        else:
            # Zero is currently used as a default value.
            pixel_values.append(0)
    return pixel_values

def find_closest_junc(thresh_dist, vert_end, verts, cells):
    """
    """
    dist_sq_min = sys.float_info.max
    cell_id = -1
    for cell in cells:
        ind = cell.vert_indices[0]
        dist_sq = ((verts[ind].x - vert_end.x) * (verts[ind].x - vert_end.x) + 
                  (verts[ind].y - vert_end.y) * (verts[ind].y - vert_end.y))
        if dist_sq < dist_sq_min:
            dist_sq_min = dist_sq
            cell_id = cell.id
    # If the minimum distance to a junction is more than the defined
    # threshold distance, set the cell id to -1 (not found).
    if dist_sq_min > thresh_dist * thresh_dist:
        return -1
    else:
        return cell_id

def main(argv):
    """
    """
    
    # Extract data from the ini file.
    print("-> Extracting parameters from ini file.")
    cfg = configparser.ConfigParser()
    cfg.read(argv[1])
    
    # Generic mesh properties.
    meshing_method = cfg.getint('input', 'meshing_method')
    nz = cfg.getint('input', 'nz')
    lz = cfg.getfloat('input', 'lz')
    z_distr = json.loads(cfg.get('input', 'z_distr'))
    layer_depths_soilt = json.loads(cfg.get('input', 'layer_depths_top_soil'))
    layer_depths_soilb = json.loads(cfg.get('input', 'layer_depths_bot_soil'))
    bedrock_bottom_elev = cfg.getfloat('input', 'bedrock_bottom_elev')
    path_aoi = cfg.get('input', 'path_aoi')
    path_ditches = cfg.get('input', 'path_ditches')
    path_to_dem_folder = cfg.get('input', 'path_to_dem_folder')
    path_to_landuse_folder = cfg.get('input', 'path_to_landuse_folder')
    path_to_soiltype_top_folder = cfg.get('input', 'path_to_soiltype_top_folder')
    path_to_soiltype_bottom_folder = cfg.get('input', 'path_to_soiltype_bottom_folder')
    path_to_mask_folder = cfg.get('input', 'path_to_mask_folder')
    path_to_stormnet_wells = cfg.get('input', 'path_to_stormnet_wells') # CAN BE REMOVED
    path_to_sinks = cfg.get('input', 'path_to_sinks') # CAN BE REMOVED
    path_to_net_junctions = cfg.get('input', 'path_to_net_junctions')
    path_to_net_links = cfg.get('input', 'path_to_net_links')
    # Structured grid.
    nx = cfg.getint('input', 'nx')
    ny = cfg.getint('input', 'ny')
    lx = cfg.getfloat('input', 'lx')
    ly = cfg.getfloat('input', 'ly')
    x_distr = json.loads(cfg.get('input', 'x_distr'))
    y_distr = json.loads(cfg.get('input', 'y_distr'))
    grid_llc = json.loads(cfg.get('input', 'grid_llc'))
    angle_deg = cfg.getfloat('input', 'grid_angle')
    # Unstructured mesh properties.
    lc = cfg.getfloat('input', 'resolution')
    lc_ditch = cfg.getfloat('input', 'resolution_ditch')
    dist_min_ditch = cfg.getfloat('input', 'dist_min_ditch')
    dist_max_ditch = cfg.getfloat('input', 'dist_max_ditch')
    num_pnts_per_curve = cfg.getint('input', 'num_pnts_per_curve')
    # Initial conditions.
    st_init_gw_depth = cfg.getfloat('input', 'st_init_gw_depth')
    st_init_soil_temp = cfg.getfloat('input', 'st_init_soil_temp')
    # Simulation settings file options.
    st_run_tests = cfg.getint('input', 'st_run_tests')
    st_num_of_cpus = cfg.getint('input', 'st_num_of_cpus')
    st_sim_length = cfg.getint('input', 'st_sim_length')
    st_time_step = cfg.getint('input', 'st_time_step')
    st_output_cell_id = cfg.getint('input', 'st_output_cell_id')
    st_result_save_inter = cfg.getint('input', 'st_result_save_inter')
    st_grid_save_start_time = cfg.getint('input', 'st_grid_save_start_time')
    st_grid_save_stop_time = cfg.getint('input', 'st_grid_save_stop_time')
    st_net_flow_solver = cfg.getint('input', 'st_net_flow_solver')
    st_net_flow_num_of_iter = cfg.getint('input', 'st_net_flow_num_of_iter')
    st_net_flow_iter_thresh = cfg.getfloat('input', 'st_net_flow_iter_thresh')
    st_net_flow_iter_implic = cfg.getfloat('input', 'st_net_flow_iter_implic')
    st_net_flow_bis_it_thr = cfg.getfloat('input', 'st_net_flow_bis_it_thr')
    st_net_flow_bis_num_it = cfg.getint('input', 'st_net_flow_bis_num_it')
    st_net_flow_bis_left = cfg.getfloat('input', 'st_net_flow_bis_left')
    st_net_flow_bis_right = cfg.getfloat('input', 'st_net_flow_bis_right')    
    st_sur_flow_solver = cfg.getint('input', 'st_sur_flow_solver')
    st_sur_flow_num_of_iter = cfg.getint('input', 'st_sur_flow_num_of_iter')
    st_sur_flow_iter_thresh = cfg.getfloat('input', 'st_sur_flow_iter_thresh')
    st_sur_flow_iter_implic = cfg.getfloat('input', 'st_sur_flow_iter_implic')
    st_sur_flow_bis_it_thr = cfg.getfloat('input', 'st_sur_flow_bis_it_thr')
    st_sur_flow_bis_num_it = cfg.getint('input', 'st_sur_flow_bis_num_it')
    st_sur_flow_bis_left = cfg.getfloat('input', 'st_sur_flow_bis_left')
    st_sur_flow_bis_right = cfg.getfloat('input', 'st_sur_flow_bis_right')
    st_sub_flow_solver = cfg.getint('input', 'st_sub_flow_solver')
    st_sub_flow_num_of_iter = cfg.getint('input', 'st_sub_flow_num_of_iter')
    st_sub_flow_iter_thresh = cfg.getfloat('input', 'st_sub_flow_iter_thresh')
    st_sub_flow_iter_implic = cfg.getfloat('input', 'st_sub_flow_iter_implic')
    st_precip_corr_factor = cfg.getfloat('input', 'st_precip_corr_factor')
    st_evap_corr_factor = cfg.getfloat('input', 'st_evap_corr_factor')
    st_heat_solver = cfg.getint('input', 'st_heat_solver')
    st_heat_num_of_iter = cfg.getint('input', 'st_heat_num_of_iter')
    st_heat_iter_thresh = cfg.getfloat('input', 'st_heat_iter_thresh')
    st_heat_iter_implic = cfg.getfloat('input', 'st_heat_iter_implic')
    st_heat_bis_iter_thresh = cfg.getfloat('input', 'st_heat_bis_iter_thresh')
    st_heat_bis_num_of_iter = cfg.getint('input', 'st_heat_bis_num_of_iter')
    st_heat_bis_left_limit = cfg.getfloat('input', 'st_heat_bis_left_limit')
    st_heat_bis_right_limit = cfg.getfloat('input', 'st_heat_bis_right_limit')
    st_heat_lat_heat_fus = cfg.getfloat('input', 'st_heat_lat_heat_fus')
    st_heat_freez_point = cfg.getfloat('input', 'st_heat_freez_point')
    st_heat_mult_freez_cur = cfg.getfloat('input', 'st_heat_mult_freez_cur')
    st_heat_dens_air = cfg.getfloat('input', 'st_heat_dens_air')
    st_heat_dens_ice = cfg.getfloat('input', 'st_heat_dens_ice')
    st_heat_dens_water = cfg.getfloat('input', 'st_heat_dens_water')
    st_heat_heat_cap_air = cfg.getfloat('input', 'st_heat_heat_cap_air')
    st_heat_heat_cap_ice = cfg.getfloat('input', 'st_heat_heat_cap_ice')
    st_heat_heat_cap_water = cfg.getfloat('input', 'st_heat_heat_cap_water')
    st_heat_cond_air = cfg.getfloat('input', 'st_heat_cond_air')
    st_heat_cond_ice = cfg.getfloat('input', 'st_heat_cond_ice')
    st_heat_cond_water = cfg.getfloat('input', 'st_heat_cond_water')
    st_heat_cond_mult_air = cfg.getfloat('input', 'st_heat_cond_mult_air')
    st_heat_cond_mult_ice = cfg.getfloat('input', 'st_heat_cond_mult_ice')
    st_heat_cond_mult_water = cfg.getfloat('input', 'st_heat_cond_mult_water')
    st_use_phreeqcrm = cfg.getint('input', 'st_use_phreeqcrm')
    st_sub_trans_solver = cfg.getint('input', 'st_sub_trans_solver')
    st_sub_trans_num_of_iter = cfg.getint('input', 'st_sub_trans_num_of_iter')
    st_sub_trans_iter_thr = cfg.getfloat('input', 'st_sub_trans_iter_thr')
    st_sub_trans_iter_imp = cfg.getfloat('input', 'st_sub_trans_iter_imp')
    st_sub_trans_mol_diff = cfg.getfloat('input', 'st_sub_trans_mol_diff')
    
    st_link_vtk_file = cfg.get('input', 'st_link_vtk_file')
    st_junction_vtk_file = cfg.get('input', 'st_junction_vtk_file')
    
    st_surface_vtk_file = cfg.get('input', 'st_surface_vtk_file')
    st_subsurface_vtk_file = cfg.get('input', 'st_subsurface_vtk_file')
    st_subsurface_grid_map_file = cfg.get('input', 'st_subsurface_grid_map_file')
    st_settings_file = cfg.get('input', 'st_settings_file')
    st_atmos_forc_file = cfg.get('input', 'st_atmos_forc_file')
    st_bound_cond_2d_file = cfg.get('input', 'st_bound_cond_2d_file')
    st_bound_cond_3d_file = cfg.get('input', 'st_bound_cond_3d_file')
    st_init_cond_2d_file = cfg.get('input', 'st_init_cond_2d_file')
    st_init_cond_3d_file = cfg.get('input', 'st_init_cond_3d_file')
    st_materials_2d_file = cfg.get('input', 'st_materials_2d_file')
    st_materials_3d_file = cfg.get('input', 'st_materials_3d_file')
    st_solute_prop_file = cfg.get('input', 'st_solute_prop_file')
    st_phreeqc_sim_descr_file = cfg.get('input', 'st_phreeqc_sim_descr_file')
    st_phreeqc_chem_db_file = cfg.get('input', 'st_phreeqc_chem_db_file')
    st_phreeqc_db_spec_map_file = cfg.get('input', 'st_phreeqc_db_spec_map_file')
    st_phreeqc_mol_weight_map_file = cfg.get('input', 'st_phreeqc_mol_weight_map_file') # too long
    st_csv_output_file = cfg.get('input', 'st_csv_output_file')
    st_link_vtk_output = cfg.get('input', 'st_link_vtk_output')
    st_junction_vtk_output = cfg.get('input', 'st_junction_vtk_output')
    st_surface_vtk_output = cfg.get('input', 'st_surface_vtk_output')
    st_subsurface_vtk_output = cfg.get('input', 'st_subsurface_vtk_output')
    st_path_to_output_folder = cfg.get('output', 'path_to_output_folder')
    st_path_to_project_folder = cfg.get('output', 'path_to_project_folder')
    
    # Set paths to raster files.
    print("-> Finding paths to raster files.")
    paths_to_dem_files = glob.glob(os.path.join(path_to_dem_folder, '*'))
    paths_to_landuse_files = glob.glob(os.path.join(path_to_landuse_folder, '*'))
    paths_to_soiltype_top_files = glob.glob(os.path.join(path_to_soiltype_top_folder, '*'))
    paths_to_soiltype_bottom_files = glob.glob(os.path.join(path_to_soiltype_bottom_folder, '*'))
    paths_to_mask_files = glob.glob(os.path.join(path_to_mask_folder, '*'))
    
    # Load DEM files.
    print("-> Loading digital elevation map raster files.")
    dem_sources, dem_pixel_sizes = load_raster_folder(paths_to_dem_files)
    print ("-> Digital elevation map raster pixel size: {}".format(dem_pixel_sizes))
    
    # Load landuse files.
    print("-> Loading landuse raster files.")
    landuse_sources, landuse_pixel_sizes = load_raster_folder(paths_to_landuse_files)
    print ("-> Land use raster pixel size: {}".format(landuse_pixel_sizes))

    # Load top soil type files.
    print("-> Loading top soil type raster files.")
    soiltype_top_sources, soiltype_top_pixel_sizes = load_raster_folder(paths_to_soiltype_top_files)
    print ("-> Top soil type raster pixel size: {}".format(soiltype_top_pixel_sizes))
    
    # Load bottom soil type files.
    print("-> Loading bottom soil type raster files.")
    soiltype_bottom_sources, soiltype_bottom_pixel_sizes = load_raster_folder(paths_to_soiltype_bottom_files)
    print ("-> Bottom soil type raster pixel size: {}".format(soiltype_bottom_pixel_sizes))

    # Load mask files.
    print("-> Loading mask raster files.")
    mask_sources, mask_pixel_sizes = load_raster_folder(paths_to_mask_files)
    print ("-> Mask raster pixel size: {}".format(mask_pixel_sizes))
    
    # Load stormwater network well geometry.
    print("-> Loading stormwater network well features.")
    stormnet_wells = []
    with fiona.open(path_to_stormnet_wells) as src:
        for feature in src:
            stormnet_wells.append(feature)
            
    # Load sinks.
    print("-> Loading sink features.")
    sinks = []
    with fiona.open(path_to_sinks) as src:
        for feature in src:
            sinks.append(feature)
    
    # Load junctions.
    print("-> Loading junction features.")
    junctions = []
    with fiona.open(path_to_net_junctions) as src:
        for feature in src:
            junctions.append(feature)
            
    # Load links.
    print("-> Loading link features.")
    links = []
    with fiona.open(path_to_net_links) as src:
        for feature in src:
            links.append(feature)
    
    # Create a structured mesh.
    if meshing_method == 1:
        # Compute cell dimensions.
        print("-> Computing cell dimensions.")
        dx_lst = []
        if len(x_distr) == nx and round(sum(x_distr),3) == 1.0:
            for dx_frac in x_distr:
                dx_lst.append(lx * dx_frac)
        else:
            for dx in range(0, nx):
                dx_lst.append(lx / nx)
        dy_lst = []
        if len(y_distr) == ny and round(sum(y_distr),3) == 1.0:
            for dy_frac in y_distr:
                dy_lst.append(ly * dy_frac)
        else:
            for dy in range(0, ny):
                dy_lst.append(ly / ny)
        
        """
        dz_lst = []
        if len(z_distr) == nz and round(sum(z_distr),3) == 1.0:
            for dz_frac in z_distr:
                dz_lst.append(lz * dz_frac)
        else:
            for dz in range(0, nz):
                dz_lst.append(lz / nz)
        """
        
        # Compute grid angle.
        print("-> Computing grid angle.")
        sin_val = math.sin(math.radians(angle_deg))
        cos_val = math.cos(math.radians(angle_deg))
        
        # Create 2d surface vertices.
        print("-> Creating 2d surface grid vertices.")
        vertices_2d = []
        ind = 0
        dy_sum = 0.0
        for j in range(0, ny + 1):
            # Get cell dimension in y direction.
            dy = 0.0
            if j > 0:
                dy = dy_lst[j - 1]
            dx_sum = 0.0
            for i in range(0, nx + 1):
                # Get cell dimension in x direction.
                dx = 0.0
                if i > 0:
                    dx = dx_lst[i - 1]
                # Create vertex and add it to the vertex list.
                x = grid_llc[0] + (dx_sum + dx) * cos_val - (dy_sum + dy) * sin_val
                y = grid_llc[1] + (dx_sum + dx) * sin_val + (dy_sum + dy) * cos_val
                z = 0.0;
                vert = Vertex(ind, x, y, z)
                vertices_2d.append(vert);
                dx_sum += dx
                ind += 1;
            dy_sum += dy
        
        # Extract cell vertex elevation values from digital elevation map rasters.
        # Interesting effect when run after creation of subsurface vertices 
        # - bottom is at the same level everywhere!
        print("-> Extracting vertex elevation values from digital elevation map rasters.")
        pixel_values = pick_raster_pixel_values(vertices_2d, dem_sources)
        for vertex, pixel_value in zip(vertices_2d, pixel_values):
            vertex.z = pixel_value
        
        # Create 2d cells.
        print("-> Creating 2d surface grid cells.")
        cells_2d = [];
        ind = 0;
        for ind_y in range(0, ny):
            for ind_x in range(0, nx):
                vert_indices = []
                vert_indices.append(ind_x + ind_y * (nx + 1) )
                vert_indices.append(ind_x + 1 + ind_y * (nx + 1))
                vert_indices.append(ind_x + 1 + (ind_y + 1) * (nx + 1) )
                vert_indices.append(ind_x + (ind_y + 1) * (nx + 1))
                con_cell_ind = ind_x + ind_y * nx
                cells_2d.append( Cell(9, vert_indices, ind, con_cell_ind) )
                ind += 1;
        
        # Create 3d subsurface vertices.
        layer_depths = []
        layer_depths.extend(layer_depths_soilt)
        layer_depths.extend(layer_depths_soilb)
        layer_depths.append(0.0)
        nz = len(layer_depths)
        dz_cum = 0.0
        vertices_3d = []
        vertices_3d.extend(vertices_2d)
        for layer_depth in layer_depths:
            dz_cum += layer_depth
            for vertex_2d in vertices_2d:
                if layer_depth > 0.0:
                    z = vertex_2d.z - dz_cum
                else:
                    z = bedrock_bottom_elev
                vertex_3d = Vertex(ind, vertex_2d.x, vertex_2d.y, z)
                vertices_3d.append(vertex_3d)
                ind += 1
        
        """
        # Create 3d subsurface vertices.
        print("-> Creating 3d subsurface vertices.")
        vertices_3d = []
        vertices_3d.extend(vertices_2d)
        dz_cum = 0.0
        ind = 0
        for k in range(0, nz):
            dz = dz_lst[k]
            dz_cum += dz
            for vertex_2d in vertices_2d:
                z = vertex_2d.z - dz_cum
                vertex_3d = Vertex(ind, vertex_2d.x, vertex_2d.y, z)
                vertices_3d.append(vertex_3d)
                ind += 1
        """
        # Create 3d cells.
        print("-> Creating 3d subsurface cells.")
        cells_3d = [];
        ind = 0;
        for ind_z in range(0, nz):
            for ind_y in range(0, ny):
                for ind_x in range(0, nx):
                    vert_indices = []
                    vert_indices.append( ind_x + ind_y * (nx + 1) + (ind_z + 1) * (nx + 1) * (ny + 1) )
                    vert_indices.append(ind_x + 1 + ind_y * (nx + 1) + (ind_z + 1) * (nx + 1) * (ny + 1))
                    vert_indices.append(ind_x + 1 + (ind_y + 1) * (nx + 1) + (ind_z + 1) * (nx + 1) * (ny + 1) )
                    vert_indices.append(ind_x + (ind_y + 1) * (nx + 1) + (ind_z + 1) * (nx + 1) * (ny + 1))
                    vert_indices.append( ind_x + ind_y * (nx + 1) + ind_z * (nx + 1) * (ny + 1) )
                    vert_indices.append(ind_x + 1 + ind_y * (nx + 1) + ind_z * (nx + 1) * (ny + 1))
                    vert_indices.append(ind_x + 1 + (ind_y + 1) * (nx + 1) + ind_z * (nx + 1) * (ny + 1) )
                    vert_indices.append(ind_x + (ind_y + 1) * (nx + 1) + ind_z * (nx + 1) * (ny + 1))
                    con_cell_ind = -1
                    if ind_z == 0:
                        con_cell_ind = ind_x + ind_y * nx
                    cells_3d.append( Cell(12, vert_indices, ind, con_cell_ind) )
                    ind += 1;
    # Create an unstructured mesh.
    elif meshing_method == 2:
        # Settings.
        gmsh_model_name = 'generated'

        # Load area of interest geometry.
        print("-> Loading area of interest geometry.")
        polygons_aoi = []
        with fiona.open(path_aoi) as src:
            for feature in src:
                polygons_aoi.append(shape(feature["geometry"]))

        # Load open ditch geometry.
        print("-> Loading open ditch geometry.")
        lines = []
        with fiona.open(path_ditches) as src:
            for idx, feature in enumerate(src):
                geom = shape(feature["geometry"])
                x_lst,y_lst = geom.coords.xy
                first_vert = None
                for x, y in zip(x_lst,y_lst):
                    if first_vert == None:
                        first_vert = [x, y]
                    else:
                        lines.append([first_vert, [x, y]])
                        first_vert = [x, y]

        # Initialize gmsh and the mesh.
        print("-> Initializing gmsh and the mesh.")
        gmsh.initialize()
        gmsh.model.add(gmsh_model_name)

        # Create mesh boundaries of the aoi polygons.
        # Area with holes does currently not work.
        print("-> Creating mesh boundaries of the aoi polygons in gmsh.")
        curve_loop_ids = []
        curve_loop_id = 1
        for polygon in polygons_aoi:
            # Create a list of vertex coordinate pares.
            polygon_x, polygon_y = polygon.exterior.coords.xy
            vertices = []
            for x, y in zip(polygon_x,polygon_y):
                vertices.append([x, y])
            # Remove the last duplicate vertex (same as the first vertex).
            last_point = vertices.pop()
            # Add polygon vertices.
            vtx_ids = []
            for point in vertices:
                vtx_ids.append(gmsh.model.geo.addPoint(point[0], point[1], 0.0, lc))
            # Add polygon lines.
            first_vert = None
            line_ids = []
            for vtx_id in vtx_ids:
                if first_vert == None:
                    first_vert = vtx_id
                else:
                    line_ids.append(gmsh.model.geo.addLine(first_vert, vtx_id))
                    first_vert = vtx_id
            line_ids.append(gmsh.model.geo.addLine(vtx_ids[-1], vtx_ids[0]))
            # Add polygon curve loop (only works with a single polygon).
            gmsh.model.geo.addCurveLoop(line_ids, curve_loop_id)
            curve_loop_ids.append(curve_loop_id)
            curve_loop_id += 1

        # Add aoi polygon surface (does not currently work with holes).
        gmsh.model.geo.addPlaneSurface(curve_loop_ids, 1)

        # Create geometry for the ditches.
        print("-> Creating geometry for the ditches.")
        line_ids_refine = []
        for line in lines:
            # Add line vertices.
            vtx_ids = []
            for point in line:
                vtx_ids.append(gmsh.model.geo.addPoint(point[0], point[1], 0.0, lc))
            # Add lines.
            line_ids_refine.append(gmsh.model.geo.addLine(vtx_ids[0], vtx_ids[1]))
        
        # Synchronize the mesh - unclear why this is done here.
        print("-> Synchronize the mesh in gmsh.")
        gmsh.model.geo.synchronize()

        # Refine the mesh close to the open ditches.
        print("-> Creating rules to refine the mesh close to the open ditches in gmsh.")
        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumbers(1, "CurvesList", line_ids_refine)
        gmsh.model.mesh.field.setNumber(1, "NumPointsPerCurve", num_pnts_per_curve)
        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", lc_ditch)
        gmsh.model.mesh.field.setNumber(2, "SizeMax", lc)
        gmsh.model.mesh.field.setNumber(2, "DistMin", dist_min_ditch)
        gmsh.model.mesh.field.setNumber(2, "DistMax", dist_max_ditch)
        gmsh.model.mesh.field.add("Min", 3)
        gmsh.model.mesh.field.setNumbers(3, "FieldsList", [2])
        gmsh.model.mesh.field.setAsBackgroundMesh(3)
        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

        # Create the mesh.
        print("-> Creating the mesh in gmsh.")
        #gmsh.model.geo.synchronize()
        #gmsh.model.addPhysicalGroup(1, [1, 2, 4], 5)
        #ps = gmsh.model.addPhysicalGroup(2, [1])
        #gmsh.model.setPhysicalName(2, ps, "My surface")
        gmsh.model.mesh.generate(2) # creates a 2d mesh

        # Show the mesh in the UI.
        #if '-nopopup' not in sys.argv:
        #    gmsh.fltk.run()

        # Acquire vertex data from gmsh.
        print("-> Acquiring vertex data from gmsh.")
        node_tags, node_coords, node_params = gmsh.model.mesh.getNodes(2, 1, includeBoundary=True)
        
        # Organize vertex data into list of lists.
        print("-> Organizing vertex data into a list.")
        vert_list = []
        for i in range(len(node_tags)):
            idx = i * 3
            # Do not replace the original node tags because the id data is used later.
            vert_data = [node_tags[i], node_coords[idx], node_coords[idx + 1], node_coords[idx + 2]]
            vert_list.append(vert_data)
        # Sort the list of lists according the vertex index.
        vert_list = sorted(vert_list, key=lambda x: x[0])

        # Acquire element data from gmsh.
        print("-> Acquiring element data from gmsh.")
        elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(2, 1)
        
        # Rename vertex indices in elements so that 
        # the indices refer to positions in the vertex list.
        print("-> Renaming vertex indices in element data.")
        node_id_dic = {}
        for idx, vert in enumerate(vert_list):
            node_id_dic[vert[0]] = idx
        elem_node_tags[0] = [node_id_dic.get(item, item) for item in elem_node_tags[0]]
        
        # Create an element list from gmsh results.
        print("-> Creating a 2d element list.")
        elem_list = []
        for i in range(len(elem_tags[0])):
            idx = i * 3
            elem_data = [elem_tags[0][i], elem_node_tags[0][idx], elem_node_tags[0][idx + 1], elem_node_tags[0][idx + 2]]
            elem_list.append(elem_data)

        # Create 2d surface vertices.
        print("-> Creating 2d surface vertices.")
        vertices_2d = []
        for vertex in vert_list:
            # Add vertex to the vertex list.
            vert = Vertex(vertex[0], vertex[1], vertex[2], vertex[3])
            vertices_2d.append(vert);
        
        # Extract cell vertex elevation values from digital elevation map rasters.
        # Interesting effect when run after creation of subsurface vertices 
        # - bottom is at the same level everywhere!
        print("-> Extracting vertex elevation values from digital elevation map rasters.")
        pixel_values = pick_raster_pixel_values(vertices_2d, dem_sources)
        for vertex, pixel_value in zip(vertices_2d, pixel_values):
            vertex.z = pixel_value
        
        # Create 2d cells.
        print("-> Creating 2d surface grid cells.")
        cells_2d = [];
        ind = 0;
        for element in elem_list:
            # Remove the id number from the list.
            element.pop(0)
            # 5 = VTK triangle
            cells_2d.append( Cell(5, element, ind, ind) )
            ind += 1;
        
        # Create 3d subsurface vertices.
        layer_depths = []
        layer_depths.extend(layer_depths_soilt)
        layer_depths.extend(layer_depths_soilb)
        layer_depths.append(0.0)
        nz = len(layer_depths)
        dz_cum = 0.0
        vertices_3d = []
        vertices_3d.extend(vertices_2d)
        for layer_depth in layer_depths:
            dz_cum += layer_depth
            for vertex_2d in vertices_2d:
                if layer_depth > 0.0:
                    z = vertex_2d.z - dz_cum
                else:
                    z = bedrock_bottom_elev
                vertex_3d = Vertex(ind, vertex_2d.x, vertex_2d.y, z)
                vertices_3d.append(vertex_3d)
                ind += 1
            
        
        """
        # Compute vertical cell dimensions.
        print("-> Computing vertical cell dimensions.")
        dz_lst = []
        if len(z_distr) == nz and round(sum(z_distr),3) == 1.0:
            for dz_frac in z_distr:
                dz_lst.append(lz * dz_frac)
        else:
            for dz in range(0, nz):
                dz_lst.append(lz / nz)
        
        # Create 3d subsurface vertices.
        print("-> Creating 3d subsurface vertices.")
        vertices_3d = []
        vertices_3d.extend(vertices_2d)
        dz_cum = 0.0
        ind = 0
        for k in range(0, nz):
            dz = dz_lst[k]
            dz_cum += dz
            for vertex_2d in vertices_2d:
                z = vertex_2d.z - dz_cum
                vertex_3d = Vertex(ind, vertex_2d.x, vertex_2d.y, z)
                vertices_3d.append(vertex_3d)
                ind += 1
        """
        # Create 3d cells.
        print("-> Creating 3d subsurface cells.")
        cells_3d = [];
        ind = 0;
        for ind_z in range(0, nz):
            for index_2d in range(0, len(cells_2d)):
                # Compute vertex indices.
                vert_indices = []                
                indices_2d = cells_2d[index_2d].vert_indices
                # Bottom indices.
                for vert_ind_2d in indices_2d:
                    vert_indices.append( int(vert_ind_2d + (ind_z + 1) * len(vertices_2d)) ) # for some reason the result is float?
                # Top indices.
                for vert_ind_2d in indices_2d:
                    vert_indices.append( int(vert_ind_2d + ind_z * len(vertices_2d)) ) # for some reason the result is float?
                # Connection index to the surface mesh..
                con_cell_ind = -1
                if ind_z == 0:
                    con_cell_ind = index_2d
                # 13 = VTK wedge
                cells_3d.append( Cell(13, vert_indices, ind, con_cell_ind) )
                ind += 1;
    else:
        vertices_2d = []
        cells_2d = []
        vertices_3d = []
        cells_3d = []
        
    # Compute cell center points.
    print("-> Computing cell center points.")
    comp_cell_center_point(cells_2d, vertices_2d)
    
    # Extract mask values from mask rasters.
    print("-> Extracting mask values from mask rasters.")
    cell_cps = []
    for cell in cells_2d:
        cell_cps.append(cell.cp)
    pixel_values = pick_raster_pixel_values(cell_cps, mask_sources)
    for cell, pixel_value in zip(cells_2d, pixel_values):
        cell.mask = pixel_value

    # Extract cell landuse values from landuse rasters.
    print("-> Extracting cell landuse values from landuse rasters.")
    cell_cps = []
    for cell in cells_2d:
        cell_cps.append(cell.cp)
    pixel_values = pick_raster_pixel_values(cell_cps, landuse_sources)
    for cell, pixel_value in zip(cells_2d, pixel_values):
        if cell.mask == 1:
            cell.material = pixel_value
    
    # Find out how many unique land use classes exist in the surface cells.
    print("-> Finding out how many unique land use classes exist in the surface cells.")
    landuse_classes = {}
    landuse_highest_class = -1
    for cell in cells_2d:
        if cell.material not in landuse_classes:
            landuse_classes[cell.material] = 1
        else:
            landuse_classes[cell.material] += 1
        if cell.material > landuse_highest_class:
            landuse_highest_class = cell.material
    print("-> Number of land use classes: {}".format(landuse_classes))

    # Create shapely polygons from surface grid cells.
    print("-> Creating shapely polygons from surface grid cells.")
    features = []
    for cell in cells_2d:
        point_list = []
        for ind in cell.vert_indices:
            point = Point(vertices_2d[ind].x, vertices_2d[ind].y)
            point_list.append(point)
        feat = Polygon([[p.x, p.y] for p in point_list])
        features.append(feat)
    
    # Save polygons boundaries to rtree structure.
    print("-> Saving polygons boundaries to rtree structure.")
    idx = index.Index()
    for pos, feat in enumerate(features):
        idx.insert(pos, feat.bounds)
    
    # Temporarily here.
    stormnet_data = {}
    cell_outlet_id = {}
    cell_outlet_dist = {}
    # sinks and junctions lists have been loaded
    # Create cells (VTK line) of junctions.
    # Create cells (VTK line) of links
    
    # Material library (links: diameter, roughness, junctions: diameter, roughness)
    # Initial conditions (links: ?, junctions: init_water_elevat)
    # Boundary conditions (links: end_point_type: 0/1, end_point0, end_point1, junctions: ?)
    
    # - Define end point location type (surface cell or junction) of each link and location id (surface cell or well id).
    # - Search locations of junctions in the surface grid (surface cell id).
    # - Define link material (pipe diameter, roughness).
    # - Define junction material (well diameter, open_lid, roughness).
    # - Find junction surface elevation from the digital elevation map.
    
    # Create junction top vertices.
    print("-> Creating junction vertices.")
    verts_junction_top = []
    ind = 0
    for junction in junctions:
        point = shape(junction['geometry'])
        verts_junction_top.append(Vertex(ind, point.x, point.y, 0.0))
        ind += 1
    
    # Extract junction top vertex elevation from the digital elevation map.
    print("-> Extracting vertex elevation values from digital elevation map rasters.")
    pixel_values = pick_raster_pixel_values(verts_junction_top, dem_sources)
    for vertex, pixel_value in zip(verts_junction_top, pixel_values):
        vertex.z = pixel_value
    
    # Add junction bottom vertices.
    print("-> Adding junction bottom vertices.")
    vertices_junction = []
    ind = 0
    for vert_junc, junction in zip(verts_junction_top, junctions):
        vert_junc.id = ind
        vertices_junction.append(vert_junc)
        ind += 1
        vert_junc = copy(vert_junc)
        vert_junc.id = ind
        depth = junction['properties']['depth']
        vert_junc.z = vert_junc.z - depth
        vertices_junction.append(vert_junc)
        ind += 1
    
    # Create vtk junction cells.
    print("-> Creating vtk junction cells.")
    cells_junction = [];
    for ind, junction in enumerate(junctions):
        vert_indices = [2 * ind, 2 * ind + 1]
        # Find the surface cell where the junction is located.
        con_cell_ind = -1
        vert_junc = vertices_junction[ vert_indices[0] ]
        vert_junc_pnt = Point(vert_junc.x, vert_junc.y)
        inters_ids = list(idx.intersection(vert_junc_pnt.bounds))
        for inters_id in inters_ids:
            try:
                if vert_junc_pnt.intersects(features[inters_id]):
                    con_cell_ind = inters_id
            except:
                print("-> Error in intersection computation with outlet "
                      "in system {}".format(network_id))
        # 3 = VTK line
        cells_junction.append( Cell(3, vert_indices, ind, con_cell_ind) )
    
    # Create vtk link vertices.
    print("-> Creating link vertices.")
    vertices_link = []
    ind = 0
    link_to_junc_thresh = 1.0
    link_bound_conds = []
    for link in links:
        line = shape(link['geometry'])
        x_lst,y_lst = line.coords.xy
        # These are currently not used.
        elev0 = link['properties']['elev0']
        elev1 = link['properties']['elev1']
        elevations = [elev0, elev1]
        link_bound_conds_loc = []
        for x, y, elev in zip(x_lst, y_lst, elevations):
            # Initialize connection properties.
            vert_end = Vertex(ind, x, y, elev)
            conn_type = -1
            conn_id = -1
            conn_elev = 0.0
            # Search link location in the junctions.
            if conn_id == -1:
                conn_id = find_closest_junc(link_to_junc_thresh, vert_end, 
                                            vertices_junction, cells_junction)
                if conn_id != -1:
                    conn_type = 0
                    junc_vert_id = cells_junction[conn_id].vert_indices[1]
                    conn_elev = vertices_junction[junc_vert_id].z
            # Search link location in the surface grid.
            if conn_id == -1:
                conn_type = -1
                end_link_geom = Point(x, y)
                inters_ids = list(idx.intersection(end_link_geom.bounds))
                for inters_id in inters_ids:
                    try:
                        if end_link_geom.intersects(features[inters_id]):
                            conn_type = 1
                            conn_id = inters_id
                            conn_elev = cells_2d[inters_id].cp.z
                    except:
                        print("-> Error in intersection computation with outlet "
                              "in system {}".format(network_id))
            # Save link connection properties.
            link_bound_conds_loc.append(conn_type)
            link_bound_conds_loc.append(conn_id)
            vert_end.z = conn_elev
            vertices_link.append(vert_end)
            ind += 1
        link_bound_conds.append(link_bound_conds_loc)
    
    # Create vtk link cells.
    print("-> Creating vtk link cells.")
    cells_link = [];
    for ind, link in enumerate(links):
        vert_indices = [2 * ind, 2 * ind + 1]
        con_cell_ind = -1
        # 3 = VTK line
        cells_link.append( Cell(3, vert_indices, ind, con_cell_ind) )
    
    # Create link mesh output.
    print("-> Creating link mesh output.")
    mesh_link_lst = create_output_data(vertices_link, cells_link, 3)
    
    # Write link mesh to disk.
    print("-> Writing link mesh to disk.")
    path_link_vtk = os.path.join(st_path_to_output_folder, st_link_vtk_file)
    write_output_data_to_disk(path_link_vtk, mesh_link_lst, ' ')
    
    # Create junction mesh output.
    print("-> Creating junction mesh output.")
    mesh_junction_lst = create_output_data(vertices_junction, cells_junction, 3)
    
    # Write junction mesh to disk.
    print("-> Writing junction mesh to disk.")
    path_junction_vtk = os.path.join(st_path_to_output_folder, st_junction_vtk_file)
    write_output_data_to_disk(path_junction_vtk, mesh_junction_lst, ' ')
    
    
    # Link the building roofs to stormwater network junctions.
    # Connect the roof to adjacent street/ground cell if there are no 
    print("-> Link the building roofs to stormwater network junctions.")
    landuse_building = 2
    link_to_junc_thresh = 50.0
    for cell_id, cell in enumerate(cells_2d):
        if cell.material == landuse_building:
            # Find the closest junction in the stormwater drainage network.
            conn_id = find_closest_junc(link_to_junc_thresh, cell.cp, 
                                        vertices_junction, cells_junction)
            #print(conn_id)
            
            """
            distance = sys.float_info.max
            sys_id = -1
            for network_id, network_data in stormnet_data.items():
                # Only search the closed wells to accelerate search.
                for junction in network_data['junct_closed']:
                    geom_junc = shape(junction['geometry'])
                    distance_new = cell_cp.distance(geom_junc)
                    if distance_new < distance:
                        distance = distance_new
                        sys_id = junction['properties']['sys_id']
            # Select first outlet from the outlet list.
            outlet = None
            if sys_id != -1 and len(stormnet_data[sys_id]['outlets']) > 0:
                outlet = stormnet_data[sys_id]['outlets'][0]
            if outlet != None:
                # Process outlet.
                geom_outlet = shape(outlet['geometry'])
                distance = cell_cp.distance(geom_outlet)
                inters_ids = list(idx.intersection(geom_outlet.bounds))
                cell_id_outlet = -1
                for inters_id in inters_ids:
                    try:
                        if geom_outlet.intersects(features[inters_id]):
                            cell_id_outlet = inters_id
                    except:
                        print("-> Error in intersection computation with outlet "
                              "in system {}".format(network_id))
                if cell_id_outlet != -1:
                    cell_outlet_id[cell_id] = cell_id_outlet
                    cell_outlet_dist[cell_id] = distance
                    cells_2d[cell_id].junction_and_outlet = distance
            """
    
    """
    # Sort the stormwater network wells into systems.
    print("-> Sorting the stormwater network wells into systems.")
    stormnet_data = {}
    for feature in stormnet_wells:
        sys_id = feature['properties']['sys_id']
        # Create dictionaries for outlets and junctions.
        if sys_id not in stormnet_data:
            stormnet_data[sys_id] = {}
            stormnet_data[sys_id]['outlets'] = []
            stormnet_data[sys_id]['junct_open'] = []
            stormnet_data[sys_id]['junct_closed'] = []
        # Save outlets.
        if feature['properties']['junc_type'] == 1:
            stormnet_data[sys_id]['outlets'].append(feature)
        # Save junctions with open lid.
        elif feature['properties']['junc_cat'] == 0:
            stormnet_data[sys_id]['junct_open'].append(feature)
        # Save junctions with closed lid.
        elif feature['properties']['junc_cat'] == 1:
            stormnet_data[sys_id]['junct_closed'].append(feature)
    
    # Find locations of junctions and outlets.
    print("-> Finding locations of junctions and outlets.")
    cell_outlet_id = {}
    cell_outlet_dist = {}
    for network_id, network_data in stormnet_data.items():
        # Select first outlet from the outlet list.
        outlet = None
        if len(network_data['outlets']) > 0:
                outlet = network_data['outlets'][0]
        if outlet != None:
            # Process outlet.
            geom_outlet = shape(outlet['geometry'])
            inters_ids = list(idx.intersection(geom_outlet.bounds))
            cell_id_outlet = -1
            for inters_id in inters_ids:
                try:
                    if geom_outlet.intersects(features[inters_id]):
                        cell_id_outlet = inters_id
                except:
                    print("-> Error in intersection computation with outlet "
                          "in system {}".format(network_id))
            # Process junctions.
            for junction in network_data['junct_open']:
                geom_junc = shape(junction['geometry'])
                inters_ids = list(idx.intersection(geom_junc.bounds))
                cell_id_junction = -1
                distance = -1.0
                for inters_id in inters_ids:
                    try:
                        if geom_junc.intersects(features[inters_id]):
                            cell_id_junction = inters_id
                            distance = geom_junc.distance(geom_outlet)
                    except:
                        print("-> Error in intersection computation with junction "
                              "in system {}".format(network_id))
                if cell_id_junction != -1:
                    cell_outlet_id[cell_id_junction] = cell_id_outlet
                    cell_outlet_dist[cell_id_junction] = distance
                    cells_2d[cell_id_junction].junction_and_outlet = distance
        else:
            print("-> System {} does not have an outlet.".format(network_id))
    
    # Find the closest outlet to a cell with builidng landuse type.
    print("-> Finding the closest outlet to a cell with builidng landuse type.")
    landuse_building = 2
    for cell_id, cell in enumerate(cells_2d):
        if cell.material == landuse_building:
            # Find the closest well in the stormwater drainage network.
            cell_cp = Point(cell.cp.x, cell.cp.y)
            distance = sys.float_info.max
            sys_id = -1
            for network_id, network_data in stormnet_data.items():
                # Only search the closed wells to accelerate search.
                for junction in network_data['junct_closed']:
                    geom_junc = shape(junction['geometry'])
                    distance_new = cell_cp.distance(geom_junc)
                    if distance_new < distance:
                        distance = distance_new
                        sys_id = junction['properties']['sys_id']
            # Select first outlet from the outlet list.
            outlet = None
            if sys_id != -1 and len(stormnet_data[sys_id]['outlets']) > 0:
                outlet = stormnet_data[sys_id]['outlets'][0]
            if outlet != None:
                # Process outlet.
                geom_outlet = shape(outlet['geometry'])
                distance = cell_cp.distance(geom_outlet)
                inters_ids = list(idx.intersection(geom_outlet.bounds))
                cell_id_outlet = -1
                for inters_id in inters_ids:
                    try:
                        if geom_outlet.intersects(features[inters_id]):
                            cell_id_outlet = inters_id
                    except:
                        print("-> Error in intersection computation with outlet "
                              "in system {}".format(network_id))
                if cell_id_outlet != -1:
                    cell_outlet_id[cell_id] = cell_id_outlet
                    cell_outlet_dist[cell_id] = distance
                    cells_2d[cell_id].junction_and_outlet = distance
    """
    
    # Find locations of sinks.
    print("-> Finding locations of sinks.")
    sink_id_in_cell = {}
    for sink in sinks:
        sink_id = sink['properties']['id']
        geom_sink = shape(sink['geometry'])
        inters_ids = list(idx.intersection(geom_sink.bounds))
        cell_ids = []
        for inters_id in inters_ids:
            try:
                if geom_sink.intersects(features[inters_id]):
                    cell_ids.append(inters_id)
            except:
                print("-> Error in intersection computation with outlet "
                      "in system {}".format(network_id))
        for cell_id in cell_ids:
            sink_id_in_cell[cell_id] = sink_id
    
    # Create 2d mesh output.
    print("-> Creating 2d mesh output.")
    if meshing_method == 1:
        mesh_2d_lst = create_output_data(vertices_2d, cells_2d, 5)
    elif meshing_method == 2:
        mesh_2d_lst = create_output_data(vertices_2d, cells_2d, 4)
    else:
        mesh_2d_lst = []
   
    # Write 2d mesh to disk.
    print("-> Writing 2d mesh to disk.")
    path_2d_vtk = os.path.join(st_path_to_output_folder, st_surface_vtk_file)
    write_output_data_to_disk(path_2d_vtk, mesh_2d_lst, ' ')
    
    # Compute cell center points.
    print("-> Computing subsurface grid cell center points.")
    comp_cell_center_point(cells_3d, vertices_3d)
    
    # Compute hydraulic head in cells.
    print("-> Computing hydraulic head in subsurface cells.")
    for ind_2d in range(0, len(cells_2d)):
        hydraulic_head = cells_2d[ind_2d].cp.z - st_init_gw_depth
        for ind_z in range(0, nz):
            ind_3d = ind_2d + ind_z * len(cells_2d)
            cells_3d[ind_3d].hydraulic_head = hydraulic_head
    
    # Extract top cell soil type values from top soil type rasters.
    print("-> Extracting top cell soil type values from top soil type rasters.")
    pixel_values_soil_top = pick_raster_pixel_values(cell_cps, soiltype_top_sources)    
    pixel_values_soil_bottom = pick_raster_pixel_values(cell_cps, soiltype_bottom_sources)
    
    
    
    
    # Save soil data into 3d grid.
    # SET MATERIAL 0 as UNDEFINED, 1 AS BEDROCK, THEN CONTINUE NORMALLY WITH INDICES
    #layer_depths.extend(layer_depths_soilt)
    #layer_depths.extend(layer_depths_soilb)
    #layer_depths.append(0.0)
    #nz = len(layer_depths)
    
    print("-> Saving soil data into subsurface grid.")
    for ind_2d in range(0, len(cells_2d)):
        surface_elevation = cells_2d[ind_2d].cp.z
        for ind_z in range(0, nz):
            ind_3d = ind_2d + ind_z * len(cells_2d)
            depth = surface_elevation - cells_3d[ind_3d].cp.z
            if depth <= 1.0:
                cells_3d[ind_3d].material = pixel_values_soil_top[ind_2d]
            else:
                cells_3d[ind_3d].material = pixel_values_soil_bottom[ind_2d]
            # TEMPORARY: SET TOP CELL SOIL TYPE TO IMPERMEABLE FOR CERTAIN LAND USE CLASSES.
            if ind_z == 0 and cells_2d[ind_2d].material in [1, 2, 3, 11]:
                cells_3d[ind_3d].material = 10
            # Set land use of cells outside the mask as undefined.
            if cells_2d[ind_2d].mask == 0:
                cells_3d[ind_3d].material = 0
    
    
    
    
    # Find out how many unique bottom soil classes exist in the subsurface cells.
    print("-> Finding out how many unique bottom soil classes exist in the subsurface cells.")
    soil_classes = {}
    soil_highest_class = -1
    for cell in cells_3d:
        if cell.material not in soil_classes:
            soil_classes[cell.material] = 1
        else:
            soil_classes[cell.material] += 1
        if cell.material > soil_highest_class:
            soil_highest_class = cell.material
    print("-> Number of soil classes: {}".format(soil_classes))

    # Create 3d mesh output.
    print("-> Creating 3d mesh output.")
    if meshing_method == 1:
        mesh_3d_lst = create_output_data(vertices_3d, cells_3d, 9)
    elif meshing_method == 2:
        mesh_3d_lst = create_output_data(vertices_3d, cells_3d, 7)
    else:
        mesh_3d_lst = []
    
    # Write 3d mesh to disk.
    print("-> Writing 3d subsurface mesh to disk.")
    path_3d_vtk = os.path.join(st_path_to_output_folder, st_subsurface_vtk_file)
    write_output_data_to_disk(path_3d_vtk, mesh_3d_lst, ' ')
    
    # Write map of the 3d subsurface grid for tridiagonal matrix algorithm.
    print("-> Writing map of the 3d subsurface grid for tridiagonal matrix algorithm.")
    grid_map_data = [['Cell indices in a column']]    
    for ind_2d in range(0, len(cells_2d)):
        # Bypass inactive areas.
        if cells_2d[ind_2d].material != 0:
            column_indices = []
            for ind_z in range(0, nz):
                ind_3d = ind_2d + ind_z * len(cells_2d)
                column_indices.append(ind_3d)
            grid_map_data.append(column_indices)
    path_grid_map = os.path.join(st_path_to_output_folder, st_subsurface_grid_map_file)
    write_output_data_to_disk(path_grid_map, grid_map_data, ',')

    # Write settings file.
    print("-> Writing settings file.")
    # Parse output paths to project folder.
    path_link_vtk = os.path.join(st_path_to_project_folder, st_link_vtk_file).replace("\\","/")
    path_junction_vtk = os.path.join(st_path_to_project_folder, st_junction_vtk_file).replace("\\","/")
    path_2d_vtk = os.path.join(st_path_to_project_folder, st_surface_vtk_file).replace("\\","/")
    path_3d_vtk = os.path.join(st_path_to_project_folder, st_subsurface_vtk_file).replace("\\","/")
    path_3d_grid_map = os.path.join(st_path_to_project_folder, st_subsurface_grid_map_file).replace("\\","/")
    path_2d_mat_lib = os.path.join(st_path_to_project_folder, st_materials_2d_file).replace("\\","/")
    path_3d_mat_lib = os.path.join(st_path_to_project_folder, st_materials_3d_file).replace("\\","/")
    path_2d_init_cond = os.path.join(st_path_to_project_folder, st_init_cond_2d_file).replace("\\","/")
    path_3d_init_cond = os.path.join(st_path_to_project_folder, st_init_cond_3d_file).replace("\\","/")
    path_2d_bound_cond = os.path.join(st_path_to_project_folder, st_bound_cond_2d_file).replace("\\","/")
    path_3d_bound_cond = os.path.join(st_path_to_project_folder, st_bound_cond_3d_file).replace("\\","/")
    path_atmos_forcing = os.path.join(st_path_to_project_folder, st_atmos_forc_file).replace("\\","/")
    path_solute_prop = os.path.join(st_path_to_project_folder, st_solute_prop_file).replace("\\","/")
    path_phreeqc_sim_descr = os.path.join(st_path_to_project_folder, st_phreeqc_sim_descr_file).replace("\\","/")
    path_phreeqc_chem_db = os.path.join(st_path_to_project_folder, st_phreeqc_chem_db_file).replace("\\","/")
    path_phreeqc_db_spec_map = os.path.join(st_path_to_project_folder, st_phreeqc_db_spec_map_file).replace("\\","/")
    path_phreeqc_mol_weight_map = os.path.join(st_path_to_project_folder, st_phreeqc_mol_weight_map_file).replace("\\","/")
    path_csv_output = os.path.join(st_path_to_project_folder, st_csv_output_file).replace("\\","/")
    path_link_vtk_output = os.path.join(st_path_to_project_folder, st_link_vtk_output).replace("\\","/")
    path_junction_vtk_output = os.path.join(st_path_to_project_folder, st_junction_vtk_output).replace("\\","/")
    path_surface_vtk_output = os.path.join(st_path_to_project_folder, st_surface_vtk_output).replace("\\","/")
    if st_subsurface_vtk_output != '-':
        path_subsurface_vtk_output = os.path.join(st_path_to_project_folder, st_subsurface_vtk_output).replace("\\","/")
    else:
        path_subsurface_vtk_output = '-'
    # Create the settings file in a list format and save it to disk.
    setttings_data = [
        ["Parameter","Value"],
        ["Run tests (0/1)",st_run_tests],
        ["Number of parallel threads (-)",st_num_of_cpus],
        ["Simulation length (s)",st_sim_length],
        ["Time step (s)",st_time_step],
        ["Output cell id (-)",st_output_cell_id],
        ["Results save interval (s)",st_result_save_inter],
        ["Grid save start time (s)",st_grid_save_start_time],
        ["Grid save stop time (s)",st_grid_save_stop_time],
        ["Network water flow solver (0/1)",st_net_flow_solver],
        ["Max. number of iterations for network water flow (-)",st_net_flow_num_of_iter],
        ["Iteration cut threshold for network water flow (m)",st_net_flow_iter_thresh],
        ["Network water flow solution implicity (-)",st_net_flow_iter_implic],
        ["Network water flow bisection iteration threshold (m)",st_net_flow_bis_it_thr],
        ["Network water flow bisection max. number of iterations (-)",st_net_flow_bis_num_it],
        ["Network water flow bisection left depth limit (m)",st_net_flow_bis_left],
        ["Network water flow bisection right depth limit (m)",st_net_flow_bis_right],
        ["Surface water flow solver (0/1)",st_sur_flow_solver],
        ["Max. number of iterations for surface water flow (-)",st_sur_flow_num_of_iter],
        ["Iteration cut threshold for surface water flow (m)",st_sur_flow_iter_thresh],
        ["Surface water flow solution implicity (-)",st_sur_flow_iter_implic],
        ["Surface water flow bisection iteration threshold (m)",st_sur_flow_bis_it_thr],
        ["Surface water flow bisection max. number of iterations (-)",st_sur_flow_bis_num_it],
        ["Surface water flow bisection left depth limit (m)",st_sur_flow_bis_left],
        ["Surface water flow bisection right depth limit (m)",st_sur_flow_bis_right],
        ["Subsurface water flow solver (0/1/2)",st_sub_flow_solver],
        ["Max. number of iterations for subsurface water flow (-)",st_sub_flow_num_of_iter],
        ["Iteration cut threshold for subsurface water flow (m)",st_sub_flow_iter_thresh],
        ["Subsurface Water flow solution implicity (-)",st_sub_flow_iter_implic],
        ["precipitation correction factor (-)",st_precip_corr_factor],
        ["Evapotranspiration correction factor (-)",st_evap_corr_factor],
        ["Heat transport solver (0/1/2)",st_heat_solver],
        ["Max. number of iterations for heat transport (-)",st_heat_num_of_iter],
        ["Iteration cut threshold for heat transport (m) (oC)",st_heat_iter_thresh],
        ["Heat transport solution implicity (-)",st_heat_iter_implic],
        ["Heat transport bisection iteration threshold (oC)",st_heat_bis_iter_thresh],
        ["Heat transport bisection max. number of iterations (-)",st_heat_bis_num_of_iter],
        ["Heat transport bisection left temperature limit (oC)",st_heat_bis_left_limit],
        ["Heat transport bisection right temperature limit (oC)",st_heat_bis_right_limit],
        ["Latent heat of fusion (kJ kg-1)",st_heat_lat_heat_fus],
        ["Freezing temperature (oC)",st_heat_freez_point],
        ["Multiplier of the freezing curve [?]",st_heat_mult_freez_cur],
        ["Density of air in heat transport (kg m-3)",st_heat_dens_air],
        ["Density of ice in heat transport (kg m-3)",st_heat_dens_ice],
        ["Density of water in heat transport (kg m-3)",st_heat_dens_water],
        ["Heat capacity of air (Specific Heat of Moist Air) 1.013 (kJ kg-1 oC-1)",st_heat_heat_cap_air],
        ["Heat capacity of ice (kJ kg-1 oC-1)",st_heat_heat_cap_ice],
        ["Heat capacity of water (kJ kg-1 oC-1)",st_heat_heat_cap_water],
        ["Heat conductivity of Air (kj s-1 m-1 oC-1)",st_heat_cond_air],
        ["Heat conductivity of Ice (kj s-1 m-1 oC-1)",st_heat_cond_ice],
        ["Heat conductivity of Water (kj s-1 m-1 oC-1)",st_heat_cond_water],
        ["Heat conductivity multiplier of air (-)",st_heat_cond_mult_air],
        ["Heat conductivity multiplier of ice (-)",st_heat_cond_mult_ice],
        ["Heat conductivity multiplier of water (-)",st_heat_cond_mult_water],
        ["Soil initial temperature (oC)",st_init_soil_temp],
        ["Use PHREEQCRM",st_use_phreeqcrm],
        ["Solute transport solver (0/1/2)",st_sub_trans_solver],
        ["Max. number of iterations for solute transport (-)",st_sub_trans_num_of_iter],
        ["Iteration cut threshold for solute transport (concentration)",st_sub_trans_iter_thr],
        ["Solute transport solution implicity (-)",st_sub_trans_iter_imp],
        ["Molecular diffusion rate (m2/s)",st_sub_trans_mol_diff],
        ["Network link input path (path to vtk)",path_link_vtk],
        ["Network junction input path (path to vtk)",path_junction_vtk],
        ["Surface grid input path (path to vtk)",path_2d_vtk],
        ["Subsurface grid input path (path to vtk)",path_3d_vtk],
        ["Subsurface grid map input path (path to grid map)",path_3d_grid_map],
        ["Surface surface material library path",path_2d_mat_lib],
        ["Subsurface volumetric material library path",path_3d_mat_lib],
        ["Surface surface initial conditions file path",path_2d_init_cond],
        ["Subsurface volumetric initial conditions file path",path_3d_init_cond],
        ["Surface surface boundary conditions file path",path_2d_bound_cond],
        ["Subsurface volumetric boundary conditions file path",path_3d_bound_cond],
        ["Atmospheric forcing input path",path_atmos_forcing],
        ["Solute properties library path",path_solute_prop],
        ["PHREEQC input data path",path_phreeqc_sim_descr],
        ["PHREEQC chemistry database path",path_phreeqc_chem_db],
        ["PHREEQC database species name map path",path_phreeqc_db_spec_map],
        ["PHREEQC species molecular weight map path",path_phreeqc_mol_weight_map],
        ["Results csv output path",path_csv_output],
        ["Network link VTK output folder path",path_link_vtk_output],
        ["Network junction VTK output folder path",path_junction_vtk_output],
        ["Surface VTK output folder path",path_surface_vtk_output],
        ["Subsurface VTK output folder path",path_subsurface_vtk_output],
    ]
    path_settings = os.path.join(st_path_to_output_folder, st_settings_file)
    write_output_data_to_disk(path_settings, setttings_data, ',')

    # Write surface materials file to disk.
    print("-> Writing surface materials file.")
    materials_2d_data = [
        ['Material','Mannings n [-]','Depression storage [m]'],
    ]
    materials_2d_data.append(['Undefined',0.1,0.001])
    materials_2d_data.append(['Road',0.01,0.001])
    materials_2d_data.append(['Building',0.01,0.001])
    materials_2d_data.append(['OtherImpermeableSurface',0.1,0.001])
    materials_2d_data.append(['Field',0.1,0.002])
    materials_2d_data.append(['LowVegetation',0.1,0.003])
    materials_2d_data.append(['Trees_2-10m',0.1,0.003])
    materials_2d_data.append(['Trees_10-15m',0.1,0.003])
    materials_2d_data.append(['Trees_15-20m',0.1,0.003])
    materials_2d_data.append(['Trees_over_20m',0.1,0.003])
    materials_2d_data.append(['BareSoil',0.1,0.002])
    materials_2d_data.append(['Water',0.01,0.0001])
    path_materials_2d = os.path.join(st_path_to_output_folder, st_materials_2d_file)
    write_output_data_to_disk(path_materials_2d, materials_2d_data, ',')

    # Write subsurface materials file to disk.
    print("-> Writing subsurface materials file.")
    materials_3d_data = [
        ['Material','Saturated hydraulic conductivity [m/s]','Compressibility [-]','Saturated water content [m^3/m^-3]','Residual water content [m^3/m^-3]','Van Genuchten wr parameter alpha (m^-1)','Van Genuchten wr parameter n [-]','Dry weight [kg/m3]','Heat capacity [kJ kg^-1 oC^-1]','Heat conductivity [kj s^-1 m^-1 oC^-1]','Heat conductivity multiplier [-]','Dispersivity trans. [m2]','Dispersivity long. [m2]'],
    ]
    materials_3d_data.append(['Undefined',   0.00000000000,0.0001,0.307,0.0197,13.58,1.8345,1680,0.84,0.0029,1,0.1,0.01])
    materials_3d_data.append(['Bedrock',     0.00000000000,0.0001,0.307,0.0197,13.58,1.8345,1680,0.84,0.0029,1,0.1,0.01])
    materials_3d_data.append(['SandyTill',   1.0000000E-07,0.0001,0.307,0.0197,13.58,1.8345,1680,0.84,0.0029,1,0.1,0.01])
    materials_3d_data.append(['Silt',        1.0000000E-07,0.0001,0.307,0.0197,13.58,1.8345,1680,0.84,0.0029,1,0.1,0.01])
    materials_3d_data.append(['Clay',        1.0000000E-08,0.0001,0.307,0.0197,13.58,1.8345,1680,0.84,0.0029,1,0.1,0.01])
    materials_3d_data.append(['Peat-sara',   1.0000000E-08,0.0001,0.307,0.0197,13.58,1.8345,1680,0.84,0.0029,1,0.1,0.01])
    materials_3d_data.append(['Peat-rahka',  1.0000000E-08,0.0001,0.307,0.0197,13.58,1.8345,1680,0.84,0.0029,1,0.1,0.01])
    materials_3d_data.append(['FillSoil',    1.0000000E-07,0.0001,0.307,0.0197,13.58,1.8345,1680,0.84,0.0029,1,0.1,0.01])
    materials_3d_data.append(['Uncharted',   1.0000000E-07,0.0001,0.307,0.0197,13.58,1.8345,1680,0.84,0.0029,1,0.1,0.01])
    materials_3d_data.append(['Water',       0.00000000000,0.0001,0.307,0.0197,13.58,1.8345,1680,0.84,0.0029,1,0.1,0.01])
    materials_3d_data.append(['impermeable', 0.00000000000,0.0001,0.307,0.0197,13.58,1.8345,1680,0.84,0.0029,1,0.1,0.01])
    path_materials_3d = os.path.join(st_path_to_output_folder, st_materials_3d_file)
    write_output_data_to_disk(path_materials_3d, materials_3d_data, ',')

    # Write surface initial conditions file.
    print("-> Writing surface initial conditions file.")
    init_cond_2d_data = [
        ['Cell id','Water depth (m)'],
    ]
    for cell in cells_2d:
        init_cond_2d_data.append([cell.id,0.0])
    path_init_cond_2d = os.path.join(st_path_to_output_folder, st_init_cond_2d_file)
    write_output_data_to_disk(path_init_cond_2d, init_cond_2d_data, ',')
    
    # Write subsurface initial conditions file.
    print("-> Writing subsurface initial conditions file.")
    init_cond_3d_data = [
        ['Cell id','Hydraulic head (m)'],
    ]
    for cell in cells_3d:
        init_cond_3d_data.append([cell.id,cell.hydraulic_head])
    path_init_cond_3d = os.path.join(st_path_to_output_folder, st_init_cond_3d_file)
    write_output_data_to_disk(path_init_cond_3d, init_cond_3d_data, ',')
    
    # Write surface boundary conditions file.
    print("-> Writing surface boundary conditions file.")
    bound_cond_2d_data = [
        ['Cell id','Outlet id (-)','Distance to outlet (m)','Sink id (-)'],
    ]
    for cell_id, cell in enumerate(cells_2d):
        # Outlet id.
        outlet_id = -1
        if cell_id in cell_outlet_id:
            outlet_id = cell_outlet_id[cell_id]
        # Distance to outlets.
        distance_to_outlet = -1.0
        if cell_id in cell_outlet_dist:
            distance_to_outlet = round(cell_outlet_dist[cell_id],1)
        # Sinks.
        sink_id = -1
        if cell_id in sink_id_in_cell:
            sink_id = sink_id_in_cell[cell_id]
        bound_cond_2d_data.append([cell.id, outlet_id, distance_to_outlet, sink_id])
    path_bound_cond_2d = os.path.join(st_path_to_output_folder, st_bound_cond_2d_file)
    write_output_data_to_disk(path_bound_cond_2d, bound_cond_2d_data, ',')
    
    # Write subsurface boundary conditions file.
    print("-> Writing subsurface boundary conditions file.")
    bound_cond_3d_data = [
        ['Cell id','PET factor (-)','Drain factor (?)','PHREEQC solution','PHREEQC equilibrium phases','PHREEQC exchange','PHREEQC surface','PHREEQC gas phase','PHREEQC solid solutions','PHREEQC kinetics'],
    ]
    for cell in cells_3d:
        bound_cond_3d_data.append([cell.id,0,0,-1,-1,-1,-1,-1,-1,-1])
    path_bound_cond_3d = os.path.join(st_path_to_output_folder, st_bound_cond_3d_file)
    write_output_data_to_disk(path_bound_cond_3d, bound_cond_3d_data, ',')
    
    # Write dummy atmospheric forcing data file for testing purposes.
    print("-> Writing atmospheric forcing data file.")
    atmos_forc_data = [
        ['Datetime','precipitation[mm]','pet [mm]','air temperature [oC]'],
    ]
    datetime_obj = datetime.datetime(2020,1,1,0,0)
    const_air_temp = 20.0
    precip_mult = 10.0 # 1.0
    max_precip_time = 1800
    atmos_forc_time_step = 300
    sim_time = 0
    while sim_time < st_sim_length:
        datetime_str = datetime_obj.strftime('%d.%m.%Y %H:%M')
        precip = 0.0
        if sim_time < 2 * max_precip_time:
            precip = 0.5 * (math.sin(math.pi * (sim_time - 0.5 * max_precip_time) / max_precip_time) + 1)
            precip *= precip_mult
        atmos_forc_data.append([datetime_str,precip,0.0,const_air_temp])
        sim_time += atmos_forc_time_step
        datetime_obj = datetime_obj + datetime.timedelta(minutes=atmos_forc_time_step/60)
    path_atmos_forc = os.path.join(st_path_to_output_folder, st_atmos_forc_file)
    write_output_data_to_disk(path_atmos_forc, atmos_forc_data, ',')
    
    # Write solute properties file.
    print("-> Writing solute properties file.")
    sol_prop_data = [
        ['Solute','Surface dry deposition [g/m2/s]','Surface wet deposition [g/m3]','Adsorption parameter 0','Adsorption parameter 1','Adsorption parameter 2','Decay rate 0','Decay target 0','Decay rate 1','Decay target 1'],
    ]
    path_sol_prop = os.path.join(st_path_to_output_folder, st_solute_prop_file) # this is already done above
    write_output_data_to_disk(path_sol_prop, sol_prop_data, ',')

def create_output_data(vertices, cells, data_mult):
    """
    """
    
    mesh_lst = []
    mesh_lst.append(["#", "vtk", "DataFile", "Version", "2.0"])
    mesh_lst.append(["STYX", "test", "mesh."])
    mesh_lst.append(["ASCII"])
    mesh_lst.append(["DATASET", "UNSTRUCTURED_GRID"])
    mesh_lst.append(["POINTS", str(len(vertices)), "double"])
    for vert in vertices:
        mesh_lst.append([str(vert.x), str(vert.y), str(vert.z)])
    mesh_lst.append(["CELLS", str(len(cells)), str( data_mult * len(cells))])
    for cell in cells:
        row = [str(len(cell.vert_indices))]
        row.extend([str(ind) for ind in cell.vert_indices])
        mesh_lst.append(row)
    mesh_lst.append(["CELL_TYPES", str(len(cells))])
    for cell in cells:
        mesh_lst.append([cell.geom_type])
    mesh_lst.append(["CELL_DATA", str(len(cells))])
    mesh_lst.append(["FIELD", "FieldData", "5"]) # TEMPORARILY CHANGED  FROM 5 TO 6    
    # Save cell id.
    mesh_lst.append(["id", "1", str(len(cells)), "int"])
    for cell in cells:
        mesh_lst.append([str(cell.id)])
    # Save surface-subsurface connection index.
    mesh_lst.append(["grid_connection", "1", str(len(cells)), "int"])
    for cell in cells:
        mesh_lst.append([str(cell.grid_connection)])
    # Save material index.
    mesh_lst.append(["material", "1", str(len(cells)), "int"])
    for cell in cells:
        mesh_lst.append([str(cell.material)])
    # Save initial condition index.
    mesh_lst.append(["init_cond", "1", str(len(cells)), "int"])
    for cell in cells:
        mesh_lst.append([str(cell.init_cond)])
    # Save boundary condition index.
    mesh_lst.append(["bound_cond", "1", str(len(cells)), "int"])
    for cell in cells:
        mesh_lst.append([str(cell.bound_cond)])
    # TEMPORARILY ADDED
    #mesh_lst.append(["junction", "1", str(len(cells)), "float"])
    #for cell in cells:
    #    mesh_lst.append([str(cell.junction_and_outlet)])
    return mesh_lst
        
def write_output_data_to_disk(path, data, delimiter):
    """
    """
    
    with open(path, 'w', newline='', encoding='utf-8') as csv_file:
        writer = csv.writer(csv_file, delimiter=delimiter)
        for row in data:
            writer.writerow(row)
 
def printBanner():
    """
    """
    
    print("----------------------------------------------------------")
    print("  _____ _______   ____   __                      _      _ ")
    print(" /  ___|_   _\\ \\ / /\\ \\ / /                     | |    | |")
    print(" \\ `--.  | |  \\ V /  \\ V /   _ __ ___   ___   __| | ___| |")
    print("  `--. \\ | |   \\ /   /   \\  | '_ ` _ \\ / _ \\ / _` |/ _ \\ |")
    print(" /\\__/ / | |   | |  / /^\\ \\ | | | | | | (_) | (_| |  __/ |")
    print(" \\____/  \\_/   \\_/  \\/   \\/ |_| |_| |_|\\___/ \\__,_|\\___|_|")
    print("----------------------------------------------------------")
    print("                   - STYX mesh builder -                  ")
    print("          2d/3d hydrological modelling framework.         ")
    print("                  Version 0.1 (2021).                     ")
    print("             Developed under the MIT license.             ")
    print("             Administrator: lassi@warsta.net.             ")
    print("----------------------------------------------------------")

# Call the main function and run the program.
if __name__ == "__main__":
    printBanner()
    if len(sys.argv) == 2:
        start_time = time.time()
        main(sys.argv)
        print("Execution time: {} s" .format(time.time() - start_time))
    else:
        print("Error, a path to ini configuration file must be given.")
    print("----------------------------------------------------------")
