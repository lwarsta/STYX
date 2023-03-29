import os
import sys
import time
import vtk
import fiona
import numpy as np
from shapely.geometry import shape, mapping, Polygon
import rasterio
from rasterio.features import rasterize
from typing import Dict, List, Tuple

"""
Command lines:
python "C:/Users/lwlassi/source/repos/lwarsta/STYX/python/vtk_to_vector_and_raster_files.py" "C:/Users/lwlassi/Downloads/styx_cluster/viikinoja_full_coars/surface_grid.vtk"
python "C:/Users/lwlassi/source/repos/lwarsta/STYX/python/vtk_to_vector_and_raster_files.py"  "C:/Users/lwlassi/Downloads/styx_cluster/laajasalo_full_high/results/mesh2d_out_060.vtk"
"""

def main(argv):
    """
    This function takes a path to a VTK file as input, reads the file, extracts the cell data, vertex data,
    and point data from the VTK output object, creates features based on the extracted data, and writes the 
    features to a vector file and data fields to separate raster files.
    
    Args:
    argv (list): Program command line arguments.
    
    Returns:
    dict: A dictionary containing the names of the properties and their data types.
    dict: A dictionary containing the data for each property for each cell.
    
    """
    # Extract the path to the VTK file and the file name.
    path_vtk_file = argv[1]
    path, vtk_file_name = os.path.split(path_vtk_file)
    
    # Print path information for debugging purposes.
    print(path)
    print(vtk_file_name)
    print(path_vtk_file)
    
    # Create a VTK reader object.
    reader = vtk.vtkUnstructuredGridReader()
    
    # Specify the VTK input file name.
    reader.SetFileName(path_vtk_file)
    
    # Read the VTK file.
    reader.Update()
    
    # Get the VTK output object from the reader.
    output = reader.GetOutput() 
    
    # Get cell data from the VTK output object.
    properties, cell_data = extract_cell_data(output)
    
    # Get cell vertex ids from the VTK output object.
    cell_vert_ids = extract_vertex_data(output)
    
    # Get point data from the VTK output object.
    point_data, bounds = extract_point_data(output)
    
    # Create features based on the extracted data.
    # Features with material 0 are currently removed.
    features = create_features(properties, point_data, cell_vert_ids, cell_data)

    # Write the features to a vector file.
    write_vector_file(path, properties, features)
    
    # Write data fields to separate raster files.
    write_raster_files(path, properties, bounds, features)

def extract_cell_data(output):
    """
    Extracts the cell data for each cell in the VTK output.
    
    Args:
    output (vtkOutput): VTK output object.
    
    Returns:
    dict: A dictionary containing the names of the properties and their data types.
    dict: A dictionary containing the data for each property for each cell.
    """
    # Get vtk cell data and number of arrays.
    vtk_cell_data = output.GetCellData()
    vtk_num_arrays = vtk_cell_data.GetNumberOfArrays()

    # Create empty dictionaries to store cell data and properties.
    cell_data = {}
    properties = {}

    # Iterate over each array in vtk cell data.
    for i in range(vtk_num_arrays):
        # Get the i-th array from the cell data object.
        array = vtk_cell_data.GetArray(i)
        # Get the name of the array and data type.
        name = array.GetName()
        dataType = array.GetDataTypeAsString()
        
        # Determine the data type and add to properties dictionary.
        if dataType == 'int':
            properties[name] = dataType
        else: # double, float
            properties[name] = 'float'
        
        # Get the number of tuples (cells) in the array.
        num_tuples = array.GetNumberOfTuples()
        
        # Create a new list to store the values of the current array.
        cell_data[name] = []
        
        # Iterate over each tuple (cell) in the array.
        for j in range(num_tuples):
            # Get the j-th tuple from the array.
            tuple_value = array.GetTuple(j)
            # Convert the tuple to a list and append to the current array list.
            cell_data[name].append(list(tuple_value))

    # Return the properties and cell data as a tuple.
    return properties, cell_data

def extract_vertex_data(output):
    """
    Extracts the vertex data (cell vertex indices) for each cell in the VTK output.

    Args:
    output (vtkOutput): VTK output object.

    Returns:
    dict: A dictionary with cell IDs as keys and a list of vertex indices as values.
    """
    # Get the number of cells in the VTK output
    vtk_num_cells = output.GetNumberOfCells()
    # Create an empty dictionary to store cell vertex indices
    cell_vert_ids = {}
    # Loop through each cell in the VTK output
    for i in range(vtk_num_cells):
        # Get the cell
        cell = output.GetCell(i)
        # Get the number of points in the cell
        num_points = cell.GetNumberOfPoints()
        # Create an empty list to store vertex indices for this cell
        cell_vert_ids[i] = []
        # Loop through each point in the cell and get its index
        for j in range(num_points):
            point_index = cell.GetPointId(j)
            cell_vert_ids[i].append(point_index)
    # Return the dictionary with cell vertex indices
    return cell_vert_ids

def extract_point_data(output):
    """
    Extracts point data and computes geometry bounds from a vtk output.

    Args:
        output (vtkOutput): A vtk output object.

    Returns:
        tuple: A tuple containing:
            - list: A list of point coordinates.
            - list: A list of bounds for the geometry.
    """
    # Create an empty list to hold the point data and initialize the bounds to extreme values.
    point_data = []
    bounds = [sys.float_info.max, sys.float_info.max, sys.float_info.min, sys.float_info.min]

    # Get the points from the vtk output and iterate over them.
    vtk_point_data = output.GetPoints()
    for i in range(vtk_point_data.GetNumberOfPoints()):
        # Get the coordinates of the current point.
        point = vtk_point_data.GetPoint(i)
        
        # Append the point to the list of point data.
        point_data.append(point)
        
        # Compute geometry bounds by comparing current point coordinates to the current bounds.
        if bounds[0] > point[0]:
            bounds[0] = point[0]
        if bounds[2] < point[0]:
            bounds[2] = point[0]
        if bounds[1] > point[1]:
            bounds[1] = point[1]
        if bounds[3] < point[1]:
            bounds[3] = point[1]
    
    # Return the point data and bounds as a tuple.
    return point_data, bounds

def create_features(properties: Dict[str, str], point_data: List[Tuple[float]], cell_vert_ids: Dict[int, List[int]], cell_data: Dict[str, List[Tuple[float]]]) -> List[Dict[str, Dict[str, List[float]]]]:
    """
    Create features from polygon data.

    Args:
        properties (dict): A dictionary containing the names of the properties and their data types.
        point_data (list): A list of tuples containing the (x, y, z) coordinates of each point.
        cell_vert_ids (dict): A dictionary containing the vertex indices for each cell.
        cell_data (dict): A dictionary containing the data for each property for each cell.

    Returns:
        A list of features, where each feature is a dictionary containing the geometry and properties of a polygon.
    """
    
    # Initialize an empty list to store the features.
    features = []

    # Loop through each key and vertice index in the cell_vert_ids dictionary.
    for key, vert_indices in cell_vert_ids.items():
        
        # Set the vertices for the polygon by looping through each vertex index and adding its coordinates to the vertices list.
        vertices = []
        for vert_index in vert_indices:
            vertex = point_data[vert_index]
            vertices.append((vertex[0], vertex[1], vertex[2]))
        
        # Close the polygon by adding the coordinates of the first vertex to the end of the vertices list.
        vertices.append((
            point_data[vert_indices[0]][0], 
            point_data[vert_indices[0]][1], 
            point_data[vert_indices[0]][2]))
        
        # Create a dictionary to store the feature data.
        feature = {}
        
        # Create a Shapely polygon geometry object from the vertices list and add it to the feature dictionary.
        polygon_geom = Polygon(vertices)
        feature['geometry'] = mapping(polygon_geom)
        
        # Loop through each property in the properties dictionary and add it to the feature dictionary.
        feature['properties'] = {}
        for prop_key, prop_type in properties.items():
            data_value = cell_data[prop_key][key][0]
            if prop_type == 'int':
                feature['properties'][prop_key] = int(data_value)
            else:
                feature['properties'][prop_key] = float(data_value)
        
        # Check if the polygon has a 'material' property and if its value is 0.
        # If it does, skip the polygon and continue to the next one.
        if 'material' in feature['properties'] and \
           feature['properties']['material'] == 0:
            continue
        
        # If the polygon does not have a 'material' property or its value is not 0, add the feature to the features list.
        features.append(feature)
        
    # Return the list of features.
    return features

def write_vector_file(path, properties, features):
    """
    Write vector data to a new file with Fiona.

    Args:
        path (str): Path to the directory where the output file will be saved.
        properties (dict): Dictionary of properties for the vector data.
        features (list): List of features to be written to the output file.

    Returns:
        None
    """
    # Define the properties of the vector data
    schema = {'geometry': 'Polygon', 'properties': properties}
    
    # Print properties to console for reference
    print(properties)

    # Set output file parameters
    path_out = os.path.join(path, 'output.gpkg')
    crs_out = 'EPSG:3067'
    driver_out = 'GPKG'
    
    # Use Fiona to open a new output file and write features to it
    with fiona.open(path_out, 'w', crs=crs_out, driver=driver_out, schema=schema) as output:
        output.writerecords(features)

def write_raster_files(path, properties, bounds, features):
    """
    Write raster data to a new file with Rasterio.

    Args:
        path (str): Path to the directory where the output file will be saved.
        properties (dict): Dictionary of properties for the raster data.
        bounds (tuple): Tuple of the form (xmin, ymin, xmax, ymax) representing the extent of the raster data.
        features (list): List of features to be written to the output file.

    Returns:
        None
    """
    # Define some parameters for the output raster files
    res = 4  # Resolution of the output raster (in meters)
    no_data_val = -1  # Value to be used for no data
    compression = 'DEFLATE'  # Compression method for the output raster
    out_crs = 'EPSG:3067'  # Coordinate reference system for the output raster
    out_driver = 'GTiff'  # Driver to be used for the output raster

    # Loop through each property and create a raster file for it
    for prop_key, prop_type in properties.items():
        # Define full path for the raster.
        file_name = 'output_' + prop_key + '.tif'
        path_out = os.path.join(path, file_name)

        # Set data type based on the property type (int or float)
        if prop_type == 'int' or 'int' in prop_type.lower():
            dtype = np.int16
        else:
            dtype = np.float32

        # Create the raster dataset to write the output to
        with rasterio.open(
            path_out,
            'w',
            driver=out_driver,
            height=int((bounds[3] - bounds[1]) / res),
            width=int((bounds[2] - bounds[0]) / res),
            count=1,
            dtype=dtype,
            compress=compression,
            nodata=no_data_val,
            crs=out_crs,
            transform=rasterio.transform.from_bounds(bounds[0], bounds[1], bounds[2], bounds[3], int((bounds[2] - bounds[0]) / res), int((bounds[3] - bounds[1]) / res)),
        ) as dst:
            # Create a numpy array with the same shape and dimensions as the raster dataset
            image = np.full((dst.height, dst.width), no_data_val, dtype=dtype)

            # Define the shapely polygons to be rasterized
            feats_to_rasterio = []
            for feature in features:
                feats_to_rasterio.append((shape(feature['geometry']), feature['properties'][prop_key]))

            # Rasterize the polygon onto the numpy array
            image = rasterize(feats_to_rasterio, out_shape=image.shape, fill=no_data_val, transform=dst.transform)

            # Write the numpy array to the raster dataset
            dst.write(image, 1)

# Define the main entry point of the script
if __name__ == "__main__":
    # Check if the script was called with exactly one argument
    if len(sys.argv) == 2:
        # Record the start time of the script
        start_time = time.time()
        # Call the main function with the script's arguments
        main(sys.argv)
        # Calculate the time it took for the script to execute and print it
        print("Execution time: {} s" .format(time.time() - start_time))
    else:
        # If the script was not called with exactly one argument, print an error message
        print("Error, a path to ini configuration file must be given.")
    # Print a separator to visually distinguish between different script runs
    print("----------------------------------------------------------")
