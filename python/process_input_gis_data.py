import os
import sys
import shutil
from osgeo import gdal, ogr, osr
import fiona
from shapely.geometry import shape, Polygon, MultiPolygon
from fiona.transform import transform_geom
import math

# Define mapping of land use codes
mapping_land_use = {
    0: 0,    # Empty
    111: 1,  # Paved road
    112: 10, # Non-paved road
    120: 2,  # Building
    121: 2,  # Building
    130: 1,  # Other impermeable surface
    211: 4,  # Field
    212: 5,  # Other low vegetation
    221: 6,  # Trees 2-10 m
    222: 7,  # Trees 10-15 m
    223: 8,  # Trees 15-20 m
    224: 9,  # Trees >20 m
    310: 3,  # Rock outcrop
    410: 10, # Bare soil
    510: 11, # Water
    520: 11, # Sea
}

# Define mapping of soil type codes to local codes used by the model.
mapping_soils = {
    195602: 1,	  # Kartoittamaton (0)
    195512: 2,	  # Saraturve (Ct) RT
    195411: 3,	  # hieno Hieta (HHt) RT
    195315: 4,	  # karkea Hieta (KHt) RT
    195314: 5,	  # Hiekka (Hk) RT
    195313: 6,	  # Sora (Sr) RT
    1953142: 7,	  # hieno Hiekka (HHk) RT
    195215: 8, 	  # Hienoainesmoreeni (HMr) RT
    195412: 9, 	  # Hiesu (Hs) RT
    195111: 10,	  # Kalliomaa (Ka) RT
    195214: 11,	  # Hiekkamoreeni (Mr) RT
    195312: 12,	  # Kiviä (Ki) RT
    195112: 13,	  # Rakka (RaKa) RT
    195511: 14,	  # Lieju (Lj) RT
    19541321: 15, # Liejusavi (LjSa) RT
    19541221: 16, # Liejuhiesu (LjHs) RT
    19531521: 17, # liejuinen Hieta (karkea) (LjHt) RT
    19541121: 18, # liejuinen hieno Hieta (LjHHt) RT
    19531421: 19, # liejuinen Hiekka (LjHk) RT
    195311: 20,	  # Lohkareita (Lo) RT
    195213: 21,	  # Soramoreeni (SrMr) RT
    195413: 22,	  # Savi (Sa) RT
    195113: 23,	  # Rapakallio (RpKa) RT
    195513: 24,	  # Rahkaturve (St) RT
    195601: 25,	  # Täytemaa (Ta)
    195603: 26,	  # Vesi (Ve)
    195514: 27,	  # Turvetuotantoalue (Tu) RT
}

# Define mapping of local top soil type codes according to land use.
remapping_top_soils = {
    1: 29,  # Paved road
    10: 5,  # Non-paved road
    2: 29,  # Building
    2: 29,  # Building
    3: 29,  # Other impermeable surface
    4: 28,  # Field
    5: 28,  # Other low vegetation
    6: 28,  # Trees 2-10 m
    7: 28,  # Trees 10-15 m
    8: 28,  # Trees 15-20 m
    9: 28,  # Trees >20 m
    3: 10,  # Rock outcrop
    10: 5,  # Bare soil
}

def get_adjusted_extents(aoi_layer):
    """
    Extracts and adjusts the extents of the AOI layer to be divisible by 32.

    Parameters:
    aoi_layer (str): Path to the AOI vector file that defines the area of interest.

    Returns:
    tuple: A tuple containing the adjusted extents (xmin, xmax, ymin, ymax).
    
    Raises:
    Exception: Prints an error message if there are issues during the reading or processing of the AOI layer.
    """
    try:
        aoi_ds = ogr.Open(aoi_layer)
        aoi_layer = aoi_ds.GetLayer()
        aoi_feature = aoi_layer.GetNextFeature()
        aoi_geom = aoi_feature.GetGeometryRef()
        xmin, xmax, ymin, ymax = aoi_geom.GetEnvelope()

        # Adjust extents to be divisible by 32
        xmin = math.floor(xmin / 32) * 32 - 16
        xmax = math.ceil(xmax / 32) * 32 + 16
        ymin = math.floor(ymin / 32) * 32 - 16
        ymax = math.ceil(ymax / 32) * 32 + 16

        return xmin, xmax, ymin, ymax
    except Exception as e:
        print(f"Error processing AOI layer: {e}")
        return None

def cut_to_aoi(input_file, extents, output_file): # aoi_layer
    """
    Clips the input file to the adjusted extents of the Area of Interest (AOI) layer and saves it to an output file. 
    The extents of the AOI are adjusted to be divisible by 32 in both the x and y directions.

    Parameters:
    input_file (str): Path to the input vector file that will be clipped. This should be a file path to a vector 
                      data format supported by OGR, like a GeoPackage (.gpkg), Shapefile (.shp), etc.
    aoi_layer (str): Path to the AOI vector file that defines the area to which the input file will be clipped. 
                     This should be a vector file containing a single polygon geometry that outlines the AOI.
    output_file (str): Path where the clipped vector data will be saved. If the file already exists, it will be overwritten. 
                       The format of the output file is GeoPackage (.gpkg).

    Outputs:
    None: The function does not return any values but writes the clipped vector data to the specified output file.

    Raises:
    Exception: Prints an error message if there are issues during the processing, such as reading input files, 
               creating output files, or clipping operations.
    """
    try:
        print(f"Cutting {input_file} to AOI and saving to {output_file}")
        input_ds = ogr.Open(input_file)
        input_layer = input_ds.GetLayer()

        # Get the extents of the area.
        xmin, xmax, ymin, ymax = extents

        # Create an envelope geometry for the adjusted extents
        ring = ogr.Geometry(ogr.wkbLinearRing)
        ring.AddPoint(xmin, ymin)
        ring.AddPoint(xmin, ymax)
        ring.AddPoint(xmax, ymax)
        ring.AddPoint(xmax, ymin)
        ring.AddPoint(xmin, ymin)
        envelope = ogr.Geometry(ogr.wkbPolygon)
        envelope.AddGeometry(ring)

        driver = ogr.GetDriverByName('GPKG')
        if os.path.exists(output_file):
            driver.DeleteDataSource(output_file)
        output_ds = driver.CreateDataSource(output_file)
        output_layer = output_ds.CreateLayer('clipped_layer', geom_type=ogr.wkbMultiPolygon, srs=input_layer.GetSpatialRef())

        input_layer_def = input_layer.GetLayerDefn()
        for i in range(input_layer_def.GetFieldCount()):
            field_def = input_layer_def.GetFieldDefn(i)
            output_layer.CreateField(field_def)

        input_layer.SetSpatialFilter(envelope)
        for feature in input_layer:
            if feature.GetGeometryRef().Intersects(envelope):
                clipped_geom = feature.GetGeometryRef().Intersection(envelope)
                feature.SetGeometry(clipped_geom)
                output_layer.CreateFeature(feature)

        output_ds = None
    except Exception as e:
        print(f"Error processing {input_file}: {e}")

def merge_layers(output_folder, output_gpkg, combined_layer_name):
    """
    Merge all layers from GeoPackage files found in the specified output folder into a single layer in a new GeoPackage file.
    All polygons are upgraded to multipolygons.

    Args:
    output_folder (str): Folder where the GeoPackage files are stored.
    output_gpkg (str): Path to the output GeoPackage file where the combined layer will be stored.
    combined_layer_name (str): Name of the combined layer in the output GeoPackage.

    Returns:
    None
    """
    gpkg_files = [os.path.join(output_folder, f) for f in os.listdir(output_folder) if f.endswith('.gpkg')]
    features = []
    crs = None
    schema = None

    for gpkg_file in gpkg_files:
        with fiona.open(gpkg_file) as src:
            if src:
                if not crs:
                    crs = src.crs  # Capture the CRS from the first file
                if not schema:
                    # Use the schema of the first non-empty layer
                    schema = src.schema
                    schema['geometry'] = 'MultiPolygon'  # Ensure the geometry type is MultiPolygon

                for feature in src:
                    geom = shape(feature['geometry'])
                    # Ensure all geometries are converted to MultiPolygon
                    if isinstance(geom, Polygon):
                        geom = MultiPolygon([geom])
                    elif isinstance(geom, MultiPolygon):
                        geom = MultiPolygon([shape(poly) for poly in geom.geoms])
                    else:
                        print(f"Unsupported geometry type: {geom.type} in feature {feature['id']}")

                    # Update the feature's geometry
                    feature['geometry'] = transform_geom(src.crs, src.crs, geom.__geo_interface__)
                    features.append(feature)

    # Write all features into the combined layer of a new GeoPackage
    with fiona.open(output_gpkg, 'w', driver='GPKG', schema=schema, layer=combined_layer_name, crs=crs) as dst:
        for feature in features:
            try:
                dst.write(feature)
            except Exception as e:
                print(f"Error writing feature {feature['id']}: {e}")

def rasterize_vector(input_gpkg, output_tif, extents, attribute_name, pixel_size=2):
    """
    Converts vector data from a GeoPackage into a raster TIFF format based on a specified attribute.
    If the attribute is not present in the vector data, a default value of 1 is burned into the raster.

    Parameters:
    input_gpkg (str): Path to the input GeoPackage containing the vector layer to be rasterized.
    output_tif (str): Path where the rasterized output will be saved as a TIFF file.
    extents (tuple): A tuple containing the extents (xmin, xmax, ymin, ymax) to use for rasterization.
    attribute_name (str): The attribute field in the vector data that will be used to assign raster cell values.
                          If the attribute does not exist, a default value of 1 is used for all raster cells covered by polygons.
    pixel_size (int): The size of each pixel in the output raster, representing both the width and height of each pixel in the same units as the input data.

    Outputs:
    None: The function writes the rasterized data to the specified TIFF file and does not return a value.

    Raises:
    Exception: If there are any errors during the file handling or rasterization process.
    """
    ds = ogr.Open(input_gpkg)
    layer = ds.GetLayer()
    x_min, x_max, y_min, y_max = extents
    x_res = int((x_max - x_min) / pixel_size)
    y_res = int((y_max - y_min) / pixel_size)

    # Get spatial reference from vector
    spatial_ref = layer.GetSpatialRef()
    wkt = spatial_ref.ExportToWkt()

    target_ds = gdal.GetDriverByName('GTiff').Create(output_tif, x_res, y_res, 1, gdal.GDT_UInt32) #  gdal.GDT_Int16
    target_ds.SetGeoTransform((x_min, pixel_size, 0, y_max, 0, -pixel_size))
    target_ds.SetProjection(wkt)  # Set the projection

    band = target_ds.GetRasterBand(1)
    band.Fill(0)  # initialize raster with zeros

    # Check if the attribute exists in the layer
    layer_def = layer.GetLayerDefn()
    if layer_def.GetFieldIndex(attribute_name) != -1:
        rasterize_options = [f"ATTRIBUTE={attribute_name}"]
    else:
        rasterize_options = ["ALL_TOUCHED=TRUE", "BURN_VALUE=1"]
    
    gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[1], options=rasterize_options)
    target_ds = None

def remap_values(input_tif, output_tif, mapping):
    """
    Remaps the pixel values of a raster file according to a specified mapping dictionary. This function
    reads an input TIFF file, applies the remapping to transform pixel values, and writes the result to an
    output TIFF file.

    Parameters:
    input_tif (str): Path to the input raster TIFF file whose pixel values are to be remapped.
    output_tif (str): Path for saving the output raster TIFF file with the remapped pixel values.
    mapping (dict): A dictionary where each key-value pair represents a mapping from the source value (key)
                    to the target value (value) that should replace it in the output raster.

    Outputs:
    None: The function does not return any values but writes the remapped raster to the specified output file.

    Raises:
    Exception: Errors during the processing, such as issues opening the input raster, applying the remapping,
               or writing to the output raster, will raise an exception with details of the failure.

    This function ensures that the spatial reference and geotransform properties of the output raster are identical 
    to those of the input raster, maintaining the geographic context and alignment of the data.
    """
    ds = gdal.Open(input_tif)
    band = ds.GetRasterBand(1)
    arr = band.ReadAsArray()

    # Create an empty array for output with the same shape as input
    out_arr = arr.copy()

    # Apply mapping
    for source, target in mapping.items():
        out_arr[arr == source] = target

    # Create output raster
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(output_tif, ds.RasterXSize, ds.RasterYSize, 1, gdal.GDT_Int16)
    out_ds.SetGeoTransform(ds.GetGeoTransform())
    
    # Set the projection from the input
    out_ds.SetProjection(ds.GetProjection())

    out_band = out_ds.GetRasterBand(1)
    out_band.WriteArray(out_arr)
    out_band.FlushCache()

    # Clean up
    out_band = None
    out_ds = None
    ds = None

def crop_raster(input_raster, output_raster, extents):
    """
    Crops the input raster based on the specified extents and saves the cropped raster to a new file.

    Parameters:
    input_raster (str): Path to the input raster file.
    output_raster (str): Path where the cropped raster will be saved.
    extents (list): A list of four floats that define the extents of the crop area in the format [xmin, xmax, ymin, ymax].

    Outputs:
    None: The function does not return any values but writes the cropped raster to the specified output file.

    Raises:
    Exception: Raises an exception if there are issues opening the input raster, creating the output raster, or during the cropping operation.
    """
    # Open input raster
    ds = gdal.Open(input_raster)
    if ds is None:
        raise Exception(f"Could not open input raster file: {input_raster}")

    # Get GeoTransform and calculate new geotransform for the output
    gt = ds.GetGeoTransform()
    x_min, x_max, y_min, y_max = extents
    x_off = int((x_min - gt[0]) / gt[1])  # cols to skip
    y_off = int((y_max - gt[3]) / gt[5])  # rows to skip (note gt[5] is negative)
    x_size = int((x_max - x_min) / gt[1])  # cols to include
    y_size = int((y_min - y_max) / gt[5])  # rows to include (note gt[5] is negative)

    # Create output raster
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(output_raster, x_size, y_size, ds.RasterCount, ds.GetRasterBand(1).DataType)
    if out_ds is None:
        raise Exception(f"Could not create output raster file: {output_raster}")

    # Set new geotransform and projection
    new_gt = [
        gt[0] + x_off * gt[1],
        gt[1],
        0.0,
        gt[3] + y_off * gt[5],
        0.0,
        gt[5]
    ]
    out_ds.SetGeoTransform(new_gt)
    out_ds.SetProjection(ds.GetProjection())

    # Rasterize the output raster
    for i in range(ds.RasterCount):
        in_band = ds.GetRasterBand(i + 1)
        out_band = out_ds.GetRasterBand(i + 1)
        data = in_band.ReadAsArray(x_off, y_off, x_size, y_size)
        out_band.WriteArray(data)
        out_band.FlushCache()

    # Cleanup
    out_ds = None
    ds = None
    print("Raster cropping completed successfully.")

def adjust_raster_values(input_path_raster, input_path_raster_adjustment, output_path_raster, value_map):
    """
    Adjust raster values based on a mapping from a land use raster.

    Args:
        input_path_raster (str): Path to the input raster file representing soil types.
        input_path_raster_adjustment (str): Path to the land use raster used for adjustments.
        output_path_raster (str): Path where the adjusted raster will be saved.
        value_map (dict): A dictionary mapping land use values (keys) to soil type values (values).

    Returns:
        None: Outputs a new raster with adjusted values.
    """

    # Open the input soil type raster
    src = gdal.Open(input_path_raster)
    soil_band = src.GetRasterBand(1)
    soil_data = soil_band.ReadAsArray()

    # Open the adjustment raster
    src_adjust = gdal.Open(input_path_raster_adjustment)
    adjust_band = src_adjust.GetRasterBand(1)
    adjustment_data = adjust_band.ReadAsArray()

    # Ensure both rasters have the same shape
    if adjustment_data.shape != soil_data.shape:
        print("Error: Rasters do not match in size.")
        return

    # Modify the soil raster based on the adjustment raster and the value map
    for land_use_value, soil_type_value in value_map.items():
        mask = (adjustment_data == land_use_value)
        soil_data[mask] = soil_type_value

    # Create the output raster
    driver = gdal.GetDriverByName('GTiff')
    out_ds = driver.Create(output_path_raster, src.RasterXSize, src.RasterYSize, 1, gdal.GDT_Int32)
    out_band = out_ds.GetRasterBand(1)

    # Set GeoTransform and Projection
    out_ds.SetGeoTransform(src.GetGeoTransform())
    out_ds.SetProjection(src.GetProjection())

    # Write the modified data and close the datasets
    out_band.WriteArray(soil_data)
    out_band.FlushCache()
    out_band = None
    out_ds = None
    src = None
    src_adjust = None

def copy_file(source_path, destination_path):
    """
    Copies a file from a source path to a destination path and creates the destination folder if it does not exist.

    Args:
        source_path (str): Path to the source file to be copied.
        destination_path (str): Path where the copied file should be placed.

    Returns:
        None: The function does not return a value but performs the file copy operation and ensures the destination directory exists.

    Raises:
        OSError: An error occurred creating the directory or copying the file.
    """
    # Check if the destination directory exists, if not create it
    if not os.path.exists(os.path.dirname(destination_path)):
        os.makedirs(os.path.dirname(destination_path))
    
    # Copy the file
    shutil.copy(source_path, destination_path)

def main(
        path_input_file_aoi, 
        path_input_folder_land_use, 
        path_input_file_soil_types, 
        path_input_file_dem, 
        path_input_file_bedrock, 
        path_output_folder_tmp,
        path_output_final_mask,
        path_output_final_dem,
        path_output_final_bedrock_elev,
        path_output_final_land_use,
        path_output_final_soil_type_top,
        path_output_final_soil_type_bottom):
    """
    """
    
    # Create the temporary work directory if it does not exist.
    if not os.path.exists(os.path.dirname(path_output_folder_tmp)):
        os.makedirs(os.path.dirname(path_output_folder_tmp))
    
    # Create output paths.
    path_output_file_dem = os.path.join(path_output_folder_tmp, 'dem.tif')
    path_output_file_bedrock_elev = os.path.join(path_output_folder_tmp, 'bedrock_elevation.tif')
    path_output_file_land_use_combined = os.path.join(path_output_folder_tmp, 'land_use_combined.gpkg')
    path_output_file_land_use = os.path.join(path_output_folder_tmp, 'land_use.tif')
    path_output_file_soil_types_top = os.path.join(path_output_folder_tmp, 'soil_types_top.tif')
    path_output_file_soil_types_bottom = os.path.join(path_output_folder_tmp, 'soil_types_bottom.tif')
    path_output_file_mask = os.path.join(path_output_folder_tmp, 'mask.tif')
    path_output_file_land_use_remapped = os.path.join(path_output_folder_tmp, 'land_use_remapped.tif')
    path_output_file_soil_types_top_remapped = os.path.join(path_output_folder_tmp, 'soil_types_top_remapped.tif')
    path_output_file_soil_types_bottom_remapped = os.path.join(path_output_folder_tmp, 'soil_types_bottom_remapped.tif')
    path_output_file_soil_types_top_remapped_land_use = os.path.join(path_output_folder_tmp, 'soil_types_top_remapped_with_landuse.tif')
        
    # Get adjusted extents.
    extents = get_adjusted_extents(path_input_file_aoi)
    
    # Rasterize a mask layer.
    rasterize_vector(path_input_file_aoi, path_output_file_mask, extents, 'xxx')
    
    # Crop the elevation model.
    crop_raster(path_input_file_dem, path_output_file_dem, extents)

    # Crop the bedrock elevation model.
    crop_raster(path_input_file_bedrock, path_output_file_bedrock_elev, extents)
    
    # Process each land use file
    for file_name in os.listdir(path_input_folder_land_use):
        print(file_name)
        if file_name.endswith('.gpkg'):
            input_file = os.path.join(path_input_folder_land_use, file_name)
            output_file = os.path.join(path_output_folder_tmp, 'clipped_' + file_name)
            cut_to_aoi(input_file, extents, output_file)
    
    # Combine clipped land use files
    merge_layers(path_output_folder_tmp, path_output_file_land_use_combined, 'combined_land_use')
    
    # Rasterize the combined land use file
    rasterize_vector(path_output_file_land_use_combined, path_output_file_land_use, extents, 'koodi')

    # Remap raster land use values
    remap_values(path_output_file_land_use, path_output_file_land_use_remapped, mapping_land_use)
    
    # Rasterize soil types.
    rasterize_vector(path_input_file_soil_types, path_output_file_soil_types_top, extents, 'PINTAMAALAJI_KOODI')
    rasterize_vector(path_input_file_soil_types, path_output_file_soil_types_bottom, extents, 'POHJAMAALAJI_KOODI')
    
    # Remap raster soil type values
    remap_values(path_output_file_soil_types_top, path_output_file_soil_types_top_remapped, mapping_soils)
    remap_values(path_output_file_soil_types_bottom, path_output_file_soil_types_bottom_remapped, mapping_soils)
    
    # Remap top soil type values to take into account impermeable areas.
    adjust_raster_values(path_output_file_soil_types_top_remapped, 
                         path_output_file_land_use_remapped, 
                         path_output_file_soil_types_top_remapped_land_use, 
                         remapping_top_soils)
    
    # Copy the files into final directories.
    copy_file(path_output_file_mask, path_output_final_mask)
    copy_file(path_output_file_dem, path_output_final_dem)
    copy_file(path_output_file_bedrock_elev, path_output_final_bedrock_elev)
    copy_file(path_output_file_land_use_remapped, path_output_final_land_use)
    copy_file(path_output_file_soil_types_top_remapped_land_use, path_output_final_soil_type_top)
    copy_file(path_output_file_soil_types_bottom_remapped, path_output_final_soil_type_bottom)

if __name__ == '__main__':
    # Command line arguments
    if len(sys.argv) != 13:
        print("Usage: python script.py "
              "<path_input_file_aoi> "
              "<path_input_folder_land_use> "
              "<path_input_file_soil_types> "
              "<path_input_file_dem> "
              "<path_input_file_bedrock> "
              "<path_output_folder_tmp> "
              "<path_output_file_mask> "
              "<path_output_file_dem> "
              "<path_output_file_bedrock_elev> "
              "<path_output_file_land_use> "
              "<path_output_file_soil_type_top> "
              "<path_output_file_soil_type_bottom>")
        sys.exit(1)
    
    # Extract input and output paths.
    path_input_file_aoi = sys.argv[1]
    path_input_folder_land_use = sys.argv[2]
    path_input_file_soil_types = sys.argv[3]
    path_input_file_dem = sys.argv[4]
    path_input_file_bedrock = sys.argv[5]
    path_output_folder_tmp = sys.argv[6]
    path_output_final_mask = sys.argv[7]
    path_output_final_dem = sys.argv[8]
    path_output_final_bedrock_elev = sys.argv[9]
    path_output_final_land_use = sys.argv[10]
    path_output_final_soil_type_top = sys.argv[11]
    path_output_final_soil_type_bottom = sys.argv[12]
    
    main(
        path_input_file_aoi, 
        path_input_folder_land_use, 
        path_input_file_soil_types, 
        path_input_file_dem, 
        path_input_file_bedrock, 
        path_output_folder_tmp,
        path_output_final_mask,
        path_output_final_dem,
        path_output_final_bedrock_elev,
        path_output_final_land_use,
        path_output_final_soil_type_top,
        path_output_final_soil_type_bottom
    )
