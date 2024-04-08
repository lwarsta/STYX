import sys
import os
import fiona
from fiona.crs import from_epsg
from shapely.geometry import mapping, Polygon, Point, LineString
import time
import copy

subcatchments_schema = {
    'geometry': 'Polygon',
    'properties': {
        'Name': 'str',
        'RainGage': 'str',
        'Outlet': 'str',
        'Area': 'float',
        '%Imperv': 'float',
        'Width': 'float',
        '%Slope': 'float',
        'CurbLen': 'float',
        'SnowPack': 'str',
        'N-Imperv': 'float',
        'N-Perv': 'float',
        'S-Imperv': 'float',
        'S-Perv': 'float',
        'PctZero': 'float',
        'RouteTo': 'str',
        'PctRouted': 'float',
        'MaxRate': 'float',
        'MinRate': 'float',
        'Decay': 'float',
        'DryTime': 'float',
        'MaxInfil': 'float',
        'Tags': 'str',
    }
}

junctions_schema = {
    'geometry': 'Point',
    'properties': {
        'Name': 'str',
        'Elevation': 'float',
        'MaxDepth': 'float',
        'InitDepth': 'float',
        'SurDepth': 'float',
        'Aponded': 'float',
    }
}

conduits_schema = {
    'geometry': 'LineString',
    'properties': {
        'Name': 'str',
        'FromNode': 'str',
        'ToNode': 'str',
        'Length': 'float',
        'Roughness': 'float',
        'InOffset': 'float',
        'OutOffset': 'float',
        'InitFlow': 'float',
        'MaxFlow': 'float',
        'Shape': 'str',
        'Geom1': 'float',
        'Geom2': 'float',
        'Geom3': 'float',
        'Geom4': 'float',
        'Barrels': 'int',
        'Culvert': 'str'
    }
}

outfalls_schema = {
    'geometry': 'Point',
    'properties': {
        'Name': 'str',
        'Elevation': 'float',
        'Type': 'str',
        'StageData': 'str',
        'Gated': 'str',
        'RouteTo': 'str'
    }
}

raingages_schema = {
    'geometry': 'Point',
    'properties': {
        'Name': 'str',
        'Format': 'str',
        'Interval': 'str',
        'SCF': 'str',
        'Source': 'str'
    }
}

links_schema_styx = {
    'geometry': 'LineString',
    'properties': {
        'id': 'int',
        'elev0': 'float',
        'elev1': 'float',
        'diameter': 'float',
        'roughness': 'float',
    }
}

junctions_schema_styx = {
    'geometry': 'Point',
    'properties': {
        'id': 'int',
        'diameter': 'float',
        'depth': 'float',
        'lid_open': 'bool',
        'type': 'int',
    }
}

def extract_subcatchments(file_path):
    """
    Extracts subcatchment data from an SWMM .inp file.

    Parameters:
    file_path (str): Path to the SWMM .inp file.

    Returns:
    list: A list of dictionaries containing subcatchment data and geometry.
    """
    subcatchments = {}
    sections = ["[SUBCATCHMENTS]", "[SUBAREAS]", "[INFILTRATION]", "[TAGS]", "[Polygons]"]
    current_section = None
    vertices = {}
    
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(';') or not line:
                continue
            
            if line in sections:
                current_section = line
                continue
            
            if line.startswith('[') and line not in sections:
                current_section = ""
            
            if current_section == "[SUBCATCHMENTS]":
                parts = line.split()
                subcatchments[parts[0]] = {
                    'RainGage': parts[1],
                    'Outlet': parts[2],
                    'Area': float(parts[3]),
                    '%Imperv': float(parts[4]),
                    'Width': float(parts[5]),
                    '%Slope': float(parts[6]),
                    'CurbLen': float(parts[7]),
                    'SnowPack': parts[8] if len(parts) > 8 else None,  # Adjust index as necessary
                    'Tags': '',
                    'N-Imperv': 0,
                    'N-Perv': 0,
                    'S-Imperv': 0,
                    'S-Perv': 0,
                    'PctZero': 0,
                    'RouteTo': '',
                    'PctRouted': 0,
                    'MaxRate': 0,
                    'MinRate': 0,
                    'Decay': 0,
                    'DryTime': 0,
                    'MaxInfil': 0,
                }
            
            elif current_section == "[SUBAREAS]":
                parts = line.split()
                sub = subcatchments.get(parts[0])
                if sub:
                    sub['N-Imperv'] = float(parts[1])
                    sub['N-Perv'] = float(parts[2])
                    sub['S-Imperv'] = float(parts[3])
                    sub['S-Perv'] = float(parts[4])
                    sub['PctZero'] = float(parts[5])
                    sub['RouteTo'] = parts[6]
                    sub['PctRouted'] = float(parts[7]) if len(parts) > 7 else None
            
            elif current_section == "[INFILTRATION]":
                parts = line.split()
                sub = subcatchments.get(parts[0])
                if sub:
                    sub['MaxRate'] = float(parts[1])
                    sub['MinRate'] = float(parts[2])
                    sub['Decay'] = float(parts[3])
                    sub['DryTime'] = float(parts[4])
                    sub['MaxInfil'] = float(parts[5])
            
            elif current_section == "[TAGS]":
                parts = line.split()
                sub = subcatchments.get(parts[1])  # Assuming second part is the subcatchment name
                if sub:
                    sub['Tags'] = parts[2]  # Assuming third part is the tag
                    
            elif current_section == "[Polygons]":
                parts = line.split()
                subcatchment_name, x_coord, y_coord = parts[0], float(parts[1]), float(parts[2])
                if subcatchment_name in vertices:
                    vertices[subcatchment_name].append((x_coord, y_coord))
                else:
                    vertices[subcatchment_name] = [(x_coord, y_coord)]
                    
    # Construct polygons from vertices and add to subcatchments
    for subcatchment_name, points in vertices.items():
        if subcatchment_name in subcatchments:
            polygon = Polygon(points)
            subcatchments[subcatchment_name]['geometry'] = polygon

    # Convert the subcatchments dictionary to a list of dictionaries for uniformity with other extraction functions
    subcatchments_list = []
    for key, value in subcatchments.items():
        value['Name'] = key
        subcatchments_list.append(value)

    return subcatchments_list

def extract_junctions(file_path):
    """
    Extracts junction data from an SWMM .inp file.

    Parameters:
    file_path (str): Path to the SWMM .inp file.

    Returns:
    list: A list of dictionaries containing junction data and geometry.
    """
    junctions = {}
    coordinates = {}
    sections = ["[JUNCTIONS]", "[COORDINATES]"]
    current_section = None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(';') or not line:
                continue

            if line in sections:
                current_section = line
                continue

            if line.startswith('[') and line not in sections:
                current_section = ""

            if current_section == "[JUNCTIONS]":
                parts = line.split()
                junctions[parts[0]] = {
                    'Elevation': float(parts[1]),
                    'MaxDepth': float(parts[2]),
                    'InitDepth': float(parts[3]),
                    'SurDepth': float(parts[4]),
                    'Aponded': float(parts[5]),
                }

            elif current_section == "[COORDINATES]":
                parts = line.split()
                node_name, x_coord, y_coord = parts[0], float(parts[1]), float(parts[2])
                coordinates[node_name] = Point(x_coord, y_coord)

    # Add geometry to junctions
    for name, junction in junctions.items():
        if name in coordinates:
            junction['geometry'] = coordinates[name]
        
    # Convert the junctions dictionary to a list of dictionaries for uniformity with other extraction functions
    junctions_list = [dict(Name=name, **data) for name, data in junctions.items()]
    
    return junctions_list

def extract_conduits(file_path):
    """
    Extracts conduit data from an SWMM .inp file and creates LineString geometries.

    Parameters:
    file_path (str): Path to the SWMM .inp file.

    Returns:
    list: A list of dictionaries containing conduit data and LineString geometries.
    """
    conduits = {}
    nodes = {}
    sections = ["[CONDUITS]", "[XSECTIONS]", "[COORDINATES]"]
    current_section = None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(';') or not line:
                continue
            
            if line in sections:
                current_section = line
                continue
            
            if line.startswith('[') and line not in sections:
                current_section = ""
                continue
            
            if current_section == "[CONDUITS]":
                parts = line.split()
                conduits[parts[0]] = {
                    'FromNode': parts[1],
                    'ToNode': parts[2],
                    'Length': float(parts[3]),
                    'Roughness': float(parts[4]),
                    'InOffset': float(parts[5]),
                    'OutOffset': float(parts[6]),
                    'InitFlow': float(parts[7]) if len(parts) > 7 else 0,
                    'MaxFlow': float(parts[8]) if len(parts) > 8 else 0,
                }

            elif current_section == "[XSECTIONS]":
                parts = line.split()
                if parts[0] in conduits:
                    conduits[parts[0]].update({
                        'Shape': parts[1],
                        'Geom1': float(parts[2]),
                        'Geom2': float(parts[3]) if len(parts) > 3 else 0,
                        'Geom3': float(parts[4]) if len(parts) > 4 else 0,
                        'Geom4': float(parts[5]) if len(parts) > 5 else 0,
                        'Barrels': int(parts[6]) if len(parts) > 6 else 1,
                        'Culvert': parts[7] if len(parts) > 7 else '',
                    })

            elif current_section == "[COORDINATES]":
                parts = line.split()
                nodes[parts[0]] = {
                    'X-Coord': float(parts[1]),
                    'Y-Coord': float(parts[2])
                }

    # Create LineString geometries for conduits
    for key, value in conduits.items():
        from_node = value['FromNode']
        to_node = value['ToNode']
        if from_node in nodes and to_node in nodes:
            from_coord = (nodes[from_node]['X-Coord'], nodes[from_node]['Y-Coord'])
            to_coord = (nodes[to_node]['X-Coord'], nodes[to_node]['Y-Coord'])
            value['geometry'] = LineString([from_coord, to_coord])
    
    # Convert the conduits dictionary to a list of dictionaries for uniformity with other extraction functions
    conduits_list = []
    for key, value in conduits.items():
        value['Name'] = key
        conduits_list.append(value)

    return conduits_list

def extract_outfalls(file_path):
    """
    Extracts outfall data from an SWMM .inp file, considering fixed-width fields.

    Parameters:
    file_path (str): Path to the SWMM .inp file.

    Returns:
    list: A list of dictionaries containing outfall data and geometry.
    """
    outfalls = {}
    coordinates = {}
    sections = ["[OUTFALLS]", "[COORDINATES]"]
    current_section = None
    field_widths = [16, 10, 10, 16, 8, 16]  # Fixed widths for each field

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(';') or not line:
                continue

            if line in sections:
                current_section = line
                continue

            if line.startswith('[') and line not in sections:
                current_section = None
                continue

            if current_section == "[OUTFALLS]":
                # Extract fields considering fixed widths
                parts = [line[sum(field_widths[:i]):sum(field_widths[:i+1])].strip() for i in range(len(field_widths))]
                name, elevation, type_, stage_data, gated, route_to = parts

                outfalls[name] = {
                    'Elevation': float(elevation) if elevation else None,
                    'Type': type_,
                    'StageData': stage_data if stage_data else None,
                    'Gated': gated,
                    'RouteTo': route_to if route_to else None,
                }

            elif current_section == "[COORDINATES]":
                parts = line.split()
                node_name, x_coord, y_coord = parts[0], float(parts[1]), float(parts[2])
                coordinates[node_name] = Point(x_coord, y_coord)

    # Add geometry to outfalls
    for name, outfall in outfalls.items():
        if name in coordinates:
            outfall['geometry'] = coordinates[name]

    # Convert the outfalls dictionary to a list of dictionaries for uniformity with other extraction functions
    outfalls_list = [dict(Name=name, **data) for name, data in outfalls.items()]

    return outfalls_list

def extract_raingages(file_path):
    """
    Extracts raingage data from an SWMM .inp file, including their coordinates and attributes.

    Parameters:
    file_path (str): Path to the SWMM .inp file.

    Returns:
    list: A list of dictionaries containing raingage data and point geometries.
    """
    raingages = {}
    symbols = {}
    sections = ["[RAINGAGES]", "[SYMBOLS]"]
    current_section = None

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith(';') or not line:
                continue

            if line in sections:
                current_section = line
                continue

            if line.startswith('[') and line not in sections:
                current_section = None

            if current_section == "[RAINGAGES]":
                parts = line.split()
                raingages[parts[0]] = {
                    'Name': parts[0],
                    'Format': parts[1],
                    'Interval': parts[2],
                    'SCF': parts[3],
                    'Source': parts[4],
                }

            elif current_section == "[SYMBOLS]":
                parts = line.split()
                gage_name, x_coord, y_coord = parts[0], float(parts[1]), float(parts[2])
                symbols[gage_name] = Point(x_coord, y_coord)

    # Merge coordinates from SYMBOLS into RAINGAGES data
    for gage_name, details in raingages.items():
        if gage_name in symbols:
            details['geometry'] = symbols[gage_name]

    # Convert the raingages dictionary to a list for consistency
    raingages_list = [dict(**data) for data in raingages.values()]
    
    return raingages_list

def construct_styx_links(conduits_data):
    """
    Constructs STYX model links from conduit data extracted from an SWMM .inp file for saving as GPKG, including geometry.

    Parameters:
    conduits_data (list): A list of dictionaries containing conduit data.

    Returns:
    list: A list of GeoJSON-like mappings formatted for saving as GPKG.
    """
    links_styx = []
    for conduit in conduits_data:
        link = {
            'id': 0, #conduit['Name'],
            'elev0': -1,  # Placeholder, update if elevation data becomes available
            'elev1': -1,  # Placeholder, update if elevation data becomes available
            'diameter': conduit['Geom1'],
            'roughness': conduit['Roughness'],
            'geometry': conduit['geometry'],
        }
        links_styx.append(link)

    return links_styx

def construct_styx_junctions(junctions_data, outfalls_data):
    """
    Constructs STYX model junctions from junction and outfall data extracted from an SWMM .inp file for saving as GPKG, including geometry.

    Parameters:
    junctions_data (list): A list of dictionaries containing junction data.
    outfalls_data (list): A list of dictionaries containing outfall data.

    Returns:
    list: A list of GeoJSON-like mappings formatted for saving as GPKG.
    """

    junctions_styx = []
    for junction in junctions_data:
        junction_styx = {
            'id': 0, #junction['Name'],
            'diameter': 0.5,  # Default value, adjust as needed
            'depth': junction['MaxDepth'],
            'lid_open': 1,  # Assuming open lid for all junctions, adjust as needed
            'type': 0,  # 0 for regular junction
            'geometry': junction['geometry'],
        }
        junctions_styx.append(junction_styx)

    for outfall in outfalls_data:
        outfall_styx = {
            'id': 0, #outfall['Name'],
            'diameter': 0.5,  # Default value, adjust as needed
            'depth': 0.0,  # outfall['Elevation']
            'lid_open': True,  # Assuming open lid for all outfalls, adjust as needed
            'type': 1,  # 1 for outfall
            'geometry': outfall['geometry'],
        }
        junctions_styx.append(outfall_styx)
    
    return junctions_styx


def write_to_gpkg(output_file, layer_name, schema, data, epsg_code, overwrite=False):
    """
    Writes given data to a specified layer in a GPKG file using Fiona and Shapely, with options to overwrite and set CRS.

    Parameters:
    output_file (str): Path to the GPKG output file.
    layer_name (str): Name of the layer to be written.
    schema (dict): Fiona schema dictionary defining the structure of the layer.
    data (list): List of dictionaries containing the data to be written. Each dictionary must include a 'geometry' key with a Shapely geometry object.
    epsg_code (int): EPSG code for the coordinate reference system.
    overwrite (bool, optional): If True, the GPKG file will be overwritten. Defaults to False.
    """
    # If overwrite is True and file exists, remove the file
    #if overwrite and os.path.exists(output_file):
    #    os.remove(output_file)

    # Setting the mode to 'w' if overwrite is True or file does not exist, otherwise 'a'
    mode = 'w' if overwrite or not os.path.exists(output_file) else 'a'
    
    with fiona.open(output_file, mode=mode, driver='GPKG', layer=layer_name, schema=schema, crs=from_epsg(epsg_code)) as layer:
        for feature in data:
            geom = feature.pop('geometry')  # Remove the geometry from the feature dictionary
            layer.write({
                'geometry': mapping(geom),  # Convert Shapely geometry to GeoJSON format
                'properties': feature  # Use the remaining dictionary as the feature properties
            })

def main(inp_file, crs, output_file):
    """
    Main function to extract SWMM .inp file data and write it to GPKG file.

    Parameters:
    inp_file (str): Path to the SWMM .inp file.
    crs (str): Coordinate reference system for the output file.
    output_file (str): Path to the output GPKG file.
    """
    start_time = time.time()

    subcatchments_data = extract_subcatchments(inp_file)
    junctions_data = extract_junctions(inp_file)
    conduits_data = extract_conduits(inp_file)
    outfalls_data = extract_outfalls(inp_file)
    raingages_data = extract_raingages(inp_file)
    styx_links_data = construct_styx_links(conduits_data)
    styx_junctions_data = construct_styx_junctions(junctions_data, outfalls_data)
    
    # Write data to GPKG using fiona
    write_to_gpkg(output_file, 'subcatchments', subcatchments_schema, subcatchments_data, crs, overwrite=True)
    write_to_gpkg(output_file, 'junctions', junctions_schema, junctions_data, crs, overwrite=True)
    write_to_gpkg(output_file, 'conduits', conduits_schema, conduits_data, crs, overwrite=True)
    write_to_gpkg(output_file, 'outfalls', outfalls_schema, outfalls_data, crs, overwrite=True)
    write_to_gpkg(output_file, 'raingages', raingages_schema, raingages_data, crs, overwrite=True)
    write_to_gpkg(output_file, 'junctions_styx', junctions_schema_styx, styx_junctions_data, crs, overwrite=True)
    write_to_gpkg(output_file, 'links_styx', links_schema_styx, styx_links_data, crs, overwrite=True)
    
    end_time = time.time()
    print(f"Data extraction and GPKG file creation completed in {end_time - start_time} seconds.")

if __name__ == '__main__':
    # Command line arguments
    if len(sys.argv) < 4:
        print("Usage: python script.py <inp_file> <crs> <output_file>")
        sys.exit(1)

    inp_file_path = sys.argv[1]
    crs = sys.argv[2]
    output_file_path = sys.argv[3]
    
    main(inp_file_path, crs, output_file_path)
