"""
IGEO7 convenience wrappers around :class:`dggrid4py.DGGRIDv8`.

All wrappers apply the current IGEO7 best practice:

- DGGRIDv8 only, with the Z7 hierarchical index (``HIERNDX`` / ``Z7``)
- ``dggs_vert0_lon = 11.20`` (DGGRID's default is 11.25)
- WGS84 geometries are converted to the authalic sphere before they are passed to DGGRID,
  and DGGRID output geometries are converted back to WGS84 (geodetic latitude)
"""
from pathlib import Path
import copy
import math
import warnings

import numpy as np
import pandas as pd

import geopandas as gpd
from shapely.geometry import Point, Polygon
from shapely.ops import transform

from dggrid4py import igeo7
from dggrid4py.auxlat import geoseries_to_authalic, geoseries_to_geodetic
from dggrid4py.dggrid_runner import DGGRIDv8

IGEO7_VERT0_LON = 11.20

_legacy_address_types = {'Z7_STRING': 'DIGIT_STRING', 'Z7': 'INT64'}


def dggrid_get_res(dggrid_instance, dggrid_dggs="ISEA7H", max_res=16):

    # IGEO7 has the same cells as ISEA7H, and DGGRID 8.42 fails on OUTPUT_STATS with the IGEO7 preset
    stats_dggs = "ISEA7H" if dggrid_dggs == "IGEO7" else dggrid_dggs
    isea7h_res = dggrid_instance.grid_stats_table(stats_dggs, max_res)
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


def igeo7_meta_config(hier_ndx_form='DIGIT_STRING', **overrides):
    """
    DGGRIDv8 meta configuration for IGEO7 with the Z7 index.

    Args:
        hier_ndx_form (str): ``DIGIT_STRING`` (e.g. ``'003456231'``) or ``INT64``
            (DGGRID writes this as the 16 character Z7 hex string, e.g. ``'0042097fffffffff'``)
        **overrides: any other DGGRID parameters, e.g. ``dggs_vert0_lon`` to change the default 11.20

    Returns:
        dict: keyword arguments for the DGGRIDv8 functions
    """
    hier_ndx_form = _normalise_hier_ndx_form(hier_ndx_form)
    meta = {
        "input_address_type": "HIERNDX",
        "input_hier_ndx_system": "Z7",
        "input_hier_ndx_form": hier_ndx_form,
        "output_address_type": "HIERNDX",
        "output_cell_label_type": "OUTPUT_ADDRESS_TYPE",
        "output_hier_ndx_system": "Z7",
        "output_hier_ndx_form": hier_ndx_form,
        "dggs_vert0_lon": IGEO7_VERT0_LON,
    }
    meta.update(overrides)
    return meta


def _normalise_hier_ndx_form(hier_ndx_form):
    if hier_ndx_form in _legacy_address_types:
        new_form = _legacy_address_types[hier_ndx_form]
        warnings.warn(
            f"address_type '{hier_ndx_form}' is the DGGRIDv7 form, use hier_ndx_form='{new_form}' with DGGRIDv8",
            DeprecationWarning,
            stacklevel=3,
        )
        return new_form
    if hier_ndx_form not in ('DIGIT_STRING', 'INT64'):
        raise ValueError(f"hier_ndx_form must be 'DIGIT_STRING' or 'INT64', got {hier_ndx_form!r}")
    return hier_ndx_form


def _require_dggrid_v8(dggrid_instance):
    if not isinstance(dggrid_instance, DGGRIDv8):
        raise TypeError(
            "IGEO7 wrappers require a DGGRIDv8 instance (DGGRID >= 8.41 with the Z7 hierarchical index), "
            f"got {type(dggrid_instance).__name__}"
        )


def _cell_id_list(cell_ids):
    # positional access that works for lists, numpy arrays and pandas Series with any index
    cell_id_list = [str(c) for c in np.asarray(cell_ids).tolist()]
    if len(cell_id_list) == 0:
        raise ValueError("no cell ids given")
    return cell_id_list


def z7_resolution(cell_id, hier_ndx_form='DIGIT_STRING'):
    """
    Resolution of a Z7 cell id in ``DIGIT_STRING`` or ``INT64`` (hex string) form.
    """
    if hier_ndx_form == 'INT64':
        return igeo7.get_z7hex_resolution(str(cell_id))
    return igeo7.get_z7string_resolution(str(cell_id))


def _from_cellids(func_name, cell_ids, dggrid_instance, hier_ndx_form, to_geodetic, meta_overrides):
    _require_dggrid_v8(dggrid_instance)
    hier_ndx_form = _normalise_hier_ndx_form(hier_ndx_form)
    cell_id_list = _cell_id_list(cell_ids)
    resolution = z7_resolution(cell_id_list[0], hier_ndx_form)
    gdf = getattr(dggrid_instance, func_name)(
        cell_id_list,
        dggs_type='IGEO7',
        resolution=resolution,
        **igeo7_meta_config(hier_ndx_form, **meta_overrides),
    )
    if to_geodetic:
        gdf['geometry'] = geoseries_to_geodetic(gdf.geometry)
    gdf = gdf.set_crs(4326, allow_override=True)
    return gdf


def dggrid_igeo7_grid_cell_centroids_from_cellids(series, dggrid_instance, hier_ndx_form='DIGIT_STRING', to_geodetic=True, **meta_overrides):
    """
    Cell centroids for a list of IGEO7 Z7 cell ids (all of the same resolution).

    Returns a GeoDataFrame in WGS84 (geodetic latitude), unless ``to_geodetic=False``
    (then the coordinates stay on DGGRID's authalic sphere).
    """
    return _from_cellids('grid_cell_centroids_from_cellids', series, dggrid_instance, hier_ndx_form, to_geodetic, meta_overrides)


def dggrid_igeo7_grid_cell_polygons_from_cellids(series, dggrid_instance, hier_ndx_form='DIGIT_STRING', to_geodetic=True, **meta_overrides):
    """
    Cell polygons for a list of IGEO7 Z7 cell ids (all of the same resolution).

    Returns a GeoDataFrame in WGS84 (geodetic latitude), unless ``to_geodetic=False``
    (then the coordinates stay on DGGRID's authalic sphere).
    """
    return _from_cellids('grid_cell_polygons_from_cellids', series, dggrid_instance, hier_ndx_form, to_geodetic, meta_overrides)


def dggrid_igeo7_grid_cell_polygons_for_extent(clip_geom, resolution, dggrid_instance, hier_ndx_form='DIGIT_STRING', to_geodetic=True, **meta_overrides):
    """
    IGEO7 cell polygons at ``resolution`` covering ``clip_geom`` (a shapely geometry in WGS84).

    The clip geometry is converted to the authalic sphere before it is passed to DGGRID, and the
    resulting cell polygons are converted back to WGS84 (unless ``to_geodetic=False``).
    """
    _require_dggrid_v8(dggrid_instance)
    hier_ndx_form = _normalise_hier_ndx_form(hier_ndx_form)
    clip_authalic = geoseries_to_authalic(gpd.GeoSeries([clip_geom])).iloc[0]
    gdf = dggrid_instance.grid_cell_polygons_for_extent(
        'IGEO7',
        resolution,
        clip_geom=clip_authalic,
        **igeo7_meta_config(hier_ndx_form, **meta_overrides),
    )
    if to_geodetic:
        gdf['geometry'] = geoseries_to_geodetic(gdf.geometry)
    gdf = gdf.set_crs(4326, allow_override=True)
    return gdf


def dggrid_igeo7_cells_for_geo_points(geodf_points_wgs84, resolution, dggrid_instance, hier_ndx_form='DIGIT_STRING', column='name', **meta_overrides):
    """
    Assign IGEO7 Z7 cell ids at ``resolution`` to WGS84 points.

    The point latitudes are converted to the authalic sphere before they are passed to DGGRID.
    Returns a copy of ``geodf_points_wgs84`` (geometry unchanged, in WGS84) with the cell ids in ``column``.
    Use :func:`dggrid_igeo7_grid_cell_polygons_from_cellids` to get the cell polygons for these ids.
    """
    _require_dggrid_v8(dggrid_instance)
    hier_ndx_form = _normalise_hier_ndx_form(hier_ndx_form)
    meta = igeo7_meta_config(hier_ndx_form, **meta_overrides)
    points_authalic = gpd.GeoDataFrame(geometry=geoseries_to_authalic(geodf_points_wgs84.geometry.reset_index(drop=True)), crs=4326)
    cells = dggrid_instance.cells_for_geo_points(
        points_authalic,
        True,
        'IGEO7',
        resolution,
        **meta,
    )
    result = geodf_points_wgs84.copy()
    result[column] = cells['name'].astype(str).values
    return result


def dggrid_igeo7_q2di_from_cellids(series, dggrid_instance, hier_ndx_form='DIGIT_STRING', **meta_overrides):
    """
    Q2DI (quad, i, j) addresses for a list of IGEO7 Z7 cell ids (all of the same resolution).
    """
    _require_dggrid_v8(dggrid_instance)
    hier_ndx_form = _normalise_hier_ndx_form(hier_ndx_form)
    cell_id_list = _cell_id_list(series)
    resolution = z7_resolution(cell_id_list[0], hier_ndx_form)
    meta = {k: v for k, v in igeo7_meta_config(hier_ndx_form, **meta_overrides).items() if not k.startswith('output_')}

    q2di_df = dggrid_instance.address_transform(cell_id_list,
                                                'IGEO7',
                                                resolution,
                                                output_address_type='Q2DI',
                                                **meta)

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



def z7_k1_ring_neighbours(z7_str, dggrid_instance, cls_m, stricter_clip=True, **meta_overrides):
    """
    Z7 ids (DIGIT_STRING form) of the direct neighbours of ``z7_str``.

    ``cls_m`` is the characteristic length scale of the resolution in metres (see :func:`dggrid_get_res`).
    All geometry work happens on DGGRID's authalic sphere, so no WGS84 conversion is needed here.
    """
    import pyproj

    _require_dggrid_v8(dggrid_instance)
    z7_str = str(z7_str)
    resolution = igeo7.get_z7string_resolution(z7_str)
    parent, local_pos, is_center = igeo7.get_z7string_local_pos(z7_str)

    if is_center:
        # centre child: the neighbours are its siblings
        if z7_is_pentagon(z7_str):
            return np.array([parent + str(n) for n in [1, 3,4,5,6]])
        else:
            return np.array([parent + str(n) for n in [1, 2, 3,4,5,6]])

    the_one = dggrid_igeo7_grid_cell_centroids_from_cellids(
        [z7_str], dggrid_instance, to_geodetic=False, **meta_overrides
    ).iloc[0]['geometry']

    if not (-180 <= the_one.x <= 180) or not (-90 <= the_one.y <= 90):
        raise ValueError(f"Not a valid lon/lat geom: {str(the_one.wkt)}")

    local_proj_str = f"+proj=laea +lat_0={the_one.y} +lon_0={the_one.x}"
    local_projection = pyproj.Transformer.from_crs('EPSG:4326', local_proj_str, always_xy=True).transform
    lamb_geom = transform(local_projection, the_one)

    neighbour_field_local = lamb_geom.buffer(cls_m)
    neighbour_field_local_clip = None
    if stricter_clip:
        neighbour_field_local_clip = neighbour_field_local.buffer(cls_m / 6)

    inverse = pyproj.Transformer.from_crs(local_proj_str, 'EPSG:4326', always_xy=True).transform
    neighbour_field = transform(inverse, neighbour_field_local)
    if stricter_clip:
        neighbour_field_clip = transform(inverse, neighbour_field_local_clip)

    # neighbour_field is already on the authalic sphere, pass it to DGGRID as is
    k_ring_group = dggrid_instance.grid_cell_centroids_for_extent('IGEO7',
                                                                  resolution,
                                                                  clip_geom=neighbour_field,
                                                                  **igeo7_meta_config('DIGIT_STRING', **meta_overrides))
    k_ring_group = k_ring_group.set_crs(4326, allow_override=True)
    k_ring_group['name'] = k_ring_group['name'].astype(str)
    if stricter_clip:
        k_ring_group = k_ring_group[k_ring_group.within(neighbour_field_clip)]

    # drop the centre cell itself
    k_ring_group = k_ring_group[k_ring_group['name'] != z7_str]
    return np.array(k_ring_group['name'].tolist())


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
