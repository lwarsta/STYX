import sys
from typing import List, Dict, Any, Optional
import fiona
from shapely.geometry import shape, Point, mapping
import argparse
import json

def process_feature(feature: Dict[str, Any], material_value: int, water_depth_threshold: float) -> Optional[Dict[str, Any]]:
    """
    Process a single feature and create a point if it meets the criteria.

    Args:
        feature (Dict[str, Any]): The input feature from the GIS file.
        material_value (int): The material value to match.
        water_depth_threshold (float): The minimum water depth threshold.

    Returns:
        Optional[Dict[str, Any]]: A new feature with a point geometry if criteria are met, None otherwise.
    """
    properties = feature['properties']
    geometry = shape(feature['geometry'])
    
    if (properties['material'] == material_value and 
        properties['water_depth'] > water_depth_threshold):
        
        center_point = geometry.centroid
        
        new_feature = {
            'geometry': mapping(center_point),
            'properties': properties
        }
        return new_feature
    return None

def main(input_file: str, output_file: str, parameter_sets: List[Dict[str, Any]]) -> None:
    """
    Main function to process the input file and create an output file with points.

    Args:
        input_file (str): Path to the input GPKG file.
        output_file (str): Path to the output GPKG file.
        parameter_sets (List[Dict[str, Any]]): List of parameter sets to check against.

    Returns:
        None
    """
    with fiona.open(input_file, 'r') as source:
        output_schema = {
            'geometry': 'Point',
            'properties': source.schema['properties']
        }
        
        with fiona.open(output_file, 'w', driver='GPKG', crs=source.crs, schema=output_schema) as sink:
            for feature in source:
                for params in parameter_sets:
                    new_feature = process_feature(feature, params['material'], params['water_depth'])
                    if new_feature:
                        sink.write(new_feature)
                        break  # Stop checking other parameter sets if a point is created

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process GIS vector file and create points based on multiple criteria.')
    parser.add_argument('input_file', help='Path to input GPKG file')
    parser.add_argument('output_file', help='Path to output GPKG file')
    parser.add_argument('parameter_sets', help='JSON string of parameter sets')
    
    args = parser.parse_args()
    
    try:
        parameter_sets = json.loads(args.parameter_sets)
    except json.JSONDecodeError:
        print("Error: Invalid JSON format for parameter sets")
        sys.exit(1)
    
    main(args.input_file, args.output_file, parameter_sets)