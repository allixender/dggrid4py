# -*- coding: utf-8 -*-

from pathlib import Path
import uuid
import shutil
import os
import sys
import subprocess
import json
import traceback
import chardet
import numpy as np
import pandas as pd

import fiona
from fiona.crs import from_epsg

import gdal
import geopandas as gpd
import rasterio

import shapely
from shapely.geometry import Polygon, box, shape

from dggrid_runner import DGGRIDv7, Dggs, dgselect


def create_grid_cells(dggs_type, resolution, mixed_aperture_level=None, clip_geom=None, capture_logs=False, silent=False):

    tmp_id = uuid.uuid4()
    tmp_dir = '/tmp/grids' 
    dggrid_instance = DGGRIDv7(executable='../src/apps/dggrid/dggrid', working_dir=tmp_dir, capture_logs=capture_logs, silent=silent )
    dggs = dgselect(dggs_type = dggs_type, res= resolution, mixed_aperture_level=mixed_aperture_level)
    
    subset_conf = { 'update_frequency': 100000, 'clip_subset_type': 'WHOLE_EARTH' }

    if not clip_geom is None and clip_geom.area > 0:

        clip_gdf = gpd.GeoDataFrame(pd.DataFrame({'id' : [1], 'geometry': [clip_geom]}), geometry='geometry', crs=from_epsg(4326))
        clip_gdf.to_file(Path(tmp_dir) / f"temp_clip_{tmp_id}.geojson", driver='GeoJSON' )

        subset_conf.update({
            'clip_subset_type': 'GDAL',
            'clip_region_files': str( (Path(tmp_dir) / f"temp_clip_{tmp_id}.geojson").resolve()),
            })

    output_conf = {
        'cell_output_type': 'GDAL',
        'cell_output_gdal_format' : 'GeoJSON',
        'cell_output_file_name': f"temp_{dggs_type}_{resolution}_out_{tmp_id}.geojson",
        }
    
    dggs_ops = dggrid_instance.grid_gen(dggs, subset_conf, output_conf )

    gdf = gpd.read_file( Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.geojson", driver='GeoJSON' )

    try:
        os.remove( str( Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.geojson") )
        os.remove( str( Path(tmp_dir) / f"temp_clip_{tmp_id}.geojson") )
    except Exception:
        pass

    return gdf


def grid_stats_table(dggs_type, resolution, mixed_aperture_level=None, clip_geom=None, capture_logs=False, silent=False):
    tmp_id = uuid.uuid4()
    tmp_dir = '/tmp/grids' 
    dggrid_instance = DGGRIDv7(executable='../src/apps/dggrid/dggrid', working_dir=tmp_dir, capture_logs=capture_logs, silent=silent )
    dggs = dgselect(dggs_type = dggs_type, res= resolution, mixed_aperture_level=mixed_aperture_level)

    dggs_ops = dggrid_instance.grid_stats(dggs)
    df = pd.DataFrame(dggs_ops['output_conf']['stats_output'])

    df.rename(columns={0: 'Resolution', 1: "Cells", 2:"Area (km^2)", 3: "CLS (km)"}, inplace=True)
    df['Resolution'] = df['Resolution'].astype(int)
    df['Cells'] = df['Cells'].astype(np.int64)
    return df



def get_cell_ids():
    pass


def proto_helpers(dggrid_instance):

    """
    - get_all_indexes/cell_ids for dggs at resolution / get_eaggr_indexes_at_level
    - geometry_from_cellid / poly_outline or point/centre
    - fill extent/subset with cells at resolution
    - get parent_for_cell_id at coarser resolution
    - get children_for_cell_id at finer resolution
    
    - sample raster values into s2 dggs cells
    - sample vector values into s2 dggs cells
    """


if __name__ == '__main__':

    est_bound = shapely.geometry.box(20.37,57.52, 28.2,60.0 )
    
    # gdf1 = create_grid_cells('ISEA4T', 3, clip_geom=None, silent=False)
    # print(gdf1.head())

    # gdf2 = create_grid_cells('ISEA7H', 5, clip_geom=est_bound, silent=True)
    # print(gdf2.head())
    # gdf2.to_file('/tmp/grids/est_shape_isea7h_5.shp')

    gdf3 = create_grid_cells('ISEA7H', 8, clip_geom=est_bound, silent=True)
    print(gdf3.head())
    # gdf3.to_file('/tmp/grids/est_shape_isea7h_8.shp')

    df1 = grid_stats_table('ISEA7H', 8, clip_geom=est_bound, silent=False)
    print(df1.head(8))



    