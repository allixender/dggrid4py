from pathlib import Path
import copy
import math

import numpy as np
import pandas as pd

import geopandas as gpd
from shapely.geometry import Point, Polygon
from shapely.ops import transform

from dggrid4py import igeo7


def dggrid_get_res(dggrid_instance, dggrid_dggs="ISEA7H", max_res=16):

    isea7h_res = dggrid_instance.grid_stats_table(dggrid_dggs, max_res)
    isea7h_res = isea7h_res.rename(
        columns={
            "Resolution": f"{dggrid_dggs}_resolution",
            "Area (km^2)": "average_hexagon_area_km2",
            "CLS (km)": "cls_km",
        }
    )
    isea7h_res = isea7h_res.rename(
        columns={col: col.lower().replace(" ", "_") for col in isea7h_res.columns}
    ).set_index(f"{dggrid_dggs.lower()}_resolution")
    isea7h_res["average_hexagon_area_m2"] = np.float32(
        isea7h_res["average_hexagon_area_km2"] * 1000000
    )
    isea7h_res["cls_m"] = np.float32(isea7h_res["cls_km"] * 1000)
    return isea7h_res


def dggrid_igeo7_grid_cell_centroids_from_cellids(series, dggrid_instance, address_type='Z7_STRING'):
    resolution = igeo7.get_z7string_resolution(series[0])
    
    gdf = dggrid_instance.grid_cell_centroids_from_cellids(series,
                                                dggs_type='IGEO7',
                                                resolution=resolution,
                                                input_address_type=address_type,
                                                output_address_type=address_type)
    gdf.crs = 4326
    return gdf


def dggrid_igeo7_grid_cell_polygons_from_cellids(series, dggrid_instance, address_type='Z7_STRING'):
    resolution = igeo7.get_z7string_resolution(series[0])
    
    gdf = dggrid_instance.grid_cell_polygons_from_cellids(series,
                                                dggs_type='IGEO7',
                                                resolution=resolution,
                                                input_address_type=address_type,
                                                output_address_type=address_type)
    gdf.crs = 4326
    return gdf



def dggrid_igeo7_q2di_from_cellids(series, dggrid_instance, address_type='Z7_STRING'):
    resolution = igeo7.get_z7string_resolution(series[0])
    
    q2di_df = dggrid_instance.address_transform(series,
                                                'IGEO7',
                                                resolution,
                                                input_address_type=address_type,
                                                output_address_type='Q2DI')
    
    q2di_df[["Q", "I", "J"]] = q2di_df['Q2DI'].apply(lambda s: pd.Series( s.split(" ") ))
    for c in ["Q", "I", "J"]:
        q2di_df[c] = q2di_df[c].astype(np.int64)
        
    return q2di_df


def to_parent_series(series):
    parents = series.apply(lambda z: igeo7.get_z7string_local_pos(z)[0])
    return parents.values


def z7_base_pentagons():
    base_pentagons = [str(b).zfill(2) for b in range(0, 12)]
    return base_pentagons


def z7_get_base_pentagon(z7_str):
    return z7_str[:2]


def z7_is_pentagon(z7_str):
    if z7_str in z7_base_pentagons():
        return True
        
    resolution = igeo7.get_z7string_resolution(z7_str)
    base_pent = z7_get_base_pentagon(z7_str)
    
    if z7_str == base_pent + str(0).zfill(resolution):
        return True
    return False



def z7_k1_ring_neighbours(z7_str, dggrid_instance, cls_m, stricter_clip=True):

    import pyproj

    resolution = igeo7.get_z7string_resolution(z7_str)
    parent, local_pos, is_center = igeo7.get_z7string_local_pos(z7_str)

    if is_center:
        # print("should be straightforward return")
        if z7_is_pentagon(z7_str):
            return np.array([parent + str(n) for n in [1, 3,4,5,6]])
        else:
            return np.array([parent + str(n) for n in [1, 2, 3,4,5,6]])
        
    the_one = dggrid_instance.grid_cell_centroids_from_cellids([z7_str],
                                                dggs_type='IGEO7',
                                                resolution=resolution,
                                                input_address_type='Z7_STRING',
                                                output_address_type='Z7_STRING').iloc[0]['geometry']

    if not (-180 <= the_one.x <= 180) and not (-90 <= the_one.y <= 90):
        raise ValueError(f"Not a valid WGS84 geom: {str(the_one.wkt)}")
        
    local_proj_str = f"+proj=laea +lat_0={the_one.y} +lon_0={the_one.x}"
    local_projection = pyproj.Transformer.from_crs('EPSG:4236', local_proj_str, always_xy=True).transform
    lamb_geom = transform(local_projection, the_one)

    neighbour_field_local = lamb_geom.buffer(cls_m)
    neighbour_field_local_clip = None
    if stricter_clip:
        neighbour_field_local_clip = neighbour_field_local.buffer(cls_m / 6)

    inverse = pyproj.Transformer.from_crs(local_proj_str, 'EPSG:4236', always_xy=True).transform
    neighbour_field = transform(inverse, neighbour_field_local)
    if stricter_clip:
        neighbour_field_clip = transform(inverse, neighbour_field_local_clip)
    
    k_ring_group = dggrid_instance.grid_cell_centroids_for_extent('IGEO7',
                                                                  resolution,
                                                                  clip_geom=neighbour_field,
                                                                  output_address_type='Z7_STRING')
    k_ring_group.crs = 4326
    # we couldtechnically already drop the geometry? But maybe need for the orientations
    if stricter_clip:
        k_ring_group = k_ring_group[k_ring_group.within(neighbour_field_clip)]

    # we should drop the centre, right?!
    has_in = k_ring_group.loc[k_ring_group['name'] == z7_str]
    if len(has_in.index) > 0:
        k_ring_group = k_ring_group.drop(has_in.index[0])
    # ideally we order/index them according Z7 local_pos?
    
    # k_ring_group['local_pos'] = k_ring_group['name'].apply(lambda cellid: cellid[-1:])
    return k_ring_group['name'].values


def suggest_window_blocks_per_chunk(rs_src, mem_use_mb):
    print("########### suggest window blocks for mem use ##########")
    shapes = []
    for bs in rs_src.block_shapes:
        print("block shapes:" + str(bs))
        shapes.append(bs)
    block_shape = shapes[0]

    mem_per_block = block_shape[0] * block_shape[1] * 64
    mem_use_byte = mem_use_mb * 1024 * 1024
    blocks_for_mem_byte = mem_use_byte / mem_per_block
    base_num_blocks = math.sqrt(blocks_for_mem_byte)
    proposed_squared_window_blocks = int(math.floor(base_num_blocks))
    baselength = base_num_blocks * block_shape[0]
    
    print(
        "suggested window_blocks per chunk -> {} (max number pixels of chuck side {} for {} mb in-mem use)".format(
            proposed_squared_window_blocks, baselength, mem_use_mb
        )
    )
    return proposed_squared_window_blocks


def extract_windows_with_bounds(raster_path, window_blocks_per_chunk=None, mem_use_mb=500):

    from pyproj import CRS

    import rasterio
    from rasterio.windows import Window
    from affine import Affine

    """
    Extracts windows from a raster with their corresponding geographic bounds.
    
    Args:
        raster_path: Path to the GeoTIFF file
        window_blocks_per_chunk: number of blocks per dimension in each window (optional)
        mem_use_mb: Memory usage constraint in MB (used if chunk_size not provided)
    
    Yields:
        tuple: (window, bounds, data_array) where:
               - window is the rasterio Window object
               - bounds is (minx, miny, maxx, maxy) in geographic coordinates
               - data_array is the array data for that window
    """
    with rasterio.open(raster_path) as src:
        # Determine window_blocks_per_chunk if not provided
        block_height, block_width = src.block_shapes[0]
        
        if window_blocks_per_chunk is None:
            window_blocks_per_chunk = suggest_window_blocks_per_chunk(src, mem_use_mb)

        # Calculate total number of blocks in each dimension
        num_block_rows = math.ceil(src.height / block_height)
        num_block_cols = math.ceil(src.width / block_width)
        
        # Process in chunks (multiple blocks)
        for row_chunk in range(0, num_block_rows, window_blocks_per_chunk):
            for col_chunk in range(0, num_block_cols, window_blocks_per_chunk):
                # Calculate window boundaries in pixels
                start_row = row_chunk * block_height
                start_col = col_chunk * block_width
                
                # Make sure we don't go past the raster boundaries
                end_row = min(start_row + (window_blocks_per_chunk * block_height), src.height)
                end_col = min(start_col + (window_blocks_per_chunk * block_width), src.width)
                
                # Create the window
                window = Window(col_off=start_col, row_off=start_row, 
                               width=end_col - start_col, height=end_row - start_row)
                
                # Get the geographic bounds of this window
                window_transform = rasterio.windows.transform(window, src.transform)
                minx, miny = rasterio.transform.xy(window_transform, end_row - start_row, 0)
                maxx, maxy = rasterio.transform.xy(window_transform, 0, end_col - start_col)
                bounds = (minx, miny, maxx, maxy)
                
                # Read the data for this window
                data = src.read(1, window=window, masked=True)
                transform = copy.deepcopy(src.transform)
                # Yield the window, bounds, and data
                yield window, bounds, data, transform


def __haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = (
        math.sin(dlat / 2) ** 2
        + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    )
    c = 2 * math.asin(math.sqrt(a))
    r = 6371  # Radius of earth in kilometers. Use 3956 for miles
    return c * r * 1000


def get_crs_info(crs_wkt):

    from pyproj import CRS

    crs = CRS.from_wkt(crs_wkt)
    is_projected = crs.is_projected
    is_geographic = crs.is_geographic

    crs_info = {
        "type": None,
        "unit_name": crs.axis_info[0].unit_name,
        "unit_conversion_factor": crs.axis_info[0].unit_conversion_factor
    }
    
    if is_projected:
        crs_info["type"] = "Projected"
    
    elif is_geographic:
        crs_info["type"] = "Geographic"

    else:
        raise ValueError("Neither projected nor geographic CRS !?")

    return crs_info


def projected_distance(east1, north1, east2, north2):
    deast = east2 - east1
    dnorth = north2 - north1
    return abs(dnorth)


def get_raster_pixel_edge_len(rs_src, adjust_latitudes, pix_size_factor):
    xtransform = rs_src.transform
    xheight = rs_src.height
    crs_wkt = rs_src.meta['crs'].wkt

    crs_info = get_crs_info(crs_wkt)

    a_hor_step = xtransform[0]
    a_vert_step = xtransform[4]

    ax1 = xtransform[2]  # origin
    ay1 = xtransform[5]  # origin
    ax2 = xtransform[2] + a_hor_step  # upper left pixel x2
    ay2 = xtransform[5] + a_vert_step  # upper left pixel y2

    pixel_edge_len = None
    
    if "metre" in crs_info["unit_name"].lower():
        pixel_edge_len = projected_distance(ax1, ay1, ax1, ay2)
        
    elif "degree" in crs_info["unit_name"].lower():
        pixel_edge_len = __haversine(ax1, ay1, ax1, ay2)

        if adjust_latitudes:
            widths = []
            lats = []
            ty = ay2
            for i in range(xheight - 1):
                hor_cell_side = __haversine(ax1, ty, ax2, ty)
                widths.append(hor_cell_side)
                lats.append(ty)
                ty = ty + a_vert_step
    
            dfb = np.array(widths)
            pixel_edge_len = dfb.std() + dfb.min()
            
    else:
        raise ValueError(f"Unclear how to calculate pixel side length with this crs: {crs_info} ")

    return pixel_edge_len


def propose_dggs_level_for_pixel_length(dggrid_instance, pixel_edge_len, pix_size_factor, dggrid_dggs="ISEA7H", max_res=16):
    cls_m = 0
    average_hexagon_area_m2 = 0
    dggrid_res = dggrid_get_res(dggrid_instance, dggrid_dggs, max_res)
    resolution = 1
        
    for idx, row in dggrid_res.iterrows():
        if row["cls_m"] < pixel_edge_len / pix_size_factor:
            cls_m = row["cls_m"]
            average_hexagon_area_m2 = row["average_hexagon_area_m2"]
            resolution = idx
            break

    print(
        f"pixel_edge_len: {pixel_edge_len} - resolution: {resolution} - {dggrid_dggs} cls_m: {cls_m} m - avg cell_size: {average_hexagon_area_m2} m2"
    )

    return {'resolution': resolution, 'cls_m': cls_m, 'average_hexagon_area_m2': average_hexagon_area_m2}


def create_geopoints_for_window(full_transform, window, data_array, crs_ref="EPSG:4326"):

    import rasterio
    """
    Extract array values within a polygon from a window array
    
    Args:
        array: The window data array
        full_transform: Affine transform of the full raster
        window: The Window object representing this chunk
        polygon: A shapely Polygon or GeoDataFrame
    
    Returns:
        GeoDataFrame with points and values that fall within the polygon
    """
    
    # Get window transform
    window_transform = rasterio.windows.transform(window, full_transform)
    
    # Create indices for the window
    rows, cols = np.indices((window.height, window.width))
    
    # Create points for all pixel centers in the window
    points = []
    row_indices = []
    col_indices = []
    data_vals = []
    
    for row in range(window.height):
        for col in range(window.width):
            # Get geographic coordinates (using center of pixel)
            x, y = window_transform * (col + 0.5, row + 0.5)
            points.append(Point(x, y))
            row_indices.append(row)
            col_indices.append(col)
            pixel_value = None
            if not data_array.mask[row, col]:
                # Access the valid pixel value
                pixel_value = data_array[row, col]
            data_vals.append(pixel_value)
            
            
    
    # Create GeoDataFrame of all points
    pixel_points = gpd.GeoDataFrame({
        'row': row_indices,
        'col': col_indices,
        'data': data_vals,
        'geometry': points
    }, crs=crs_ref) 

    return pixel_points
