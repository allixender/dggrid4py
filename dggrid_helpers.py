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


def grid_cell_polygons_for_extent(dggs_type, resolution, mixed_aperture_level=None, clip_geom=None, capture_logs=False, silent=True):
    """
    if extent aka clip_geom is empty/None then the whole_earth is default area
    """
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


def grid_stats_table(dggs_type, resolution, mixed_aperture_level=None, capture_logs=False, silent=True):
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


def grid_cellids_for_extent(dggs_type, resolution, mixed_aperture_level=None, clip_geom=None, capture_logs=False, silent=True):
    """
    if extent aka clip_geom is empty/None then the whole_earth is default area
    """
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
        'point_output_type': 'TEXT',
        'point_output_file_name': f"temp_{dggs_type}_{resolution}_out_{tmp_id}",
        }
    
    dggs_ops = dggrid_instance.grid_gen(dggs, subset_conf, output_conf )

    df = pd.read_csv( Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.txt" , header=None)
    df = df.dropna()

    try:
        os.remove( str( Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.txt") )
        os.remove( str( Path(tmp_dir) / f"temp_clip_{tmp_id}.geojson") )
    except Exception:
        pass

    return df


def grid_cell_polygons_from_cellids(cell_id_list, dggs_type, resolution, mixed_aperture_level=None, capture_logs=True, silent=False):
    """
    if cell_id_list is empty/None then the whole_earth is default area
    """

    tmp_id = uuid.uuid4()
    tmp_dir = '/tmp/grids' 
    dggrid_instance = DGGRIDv7(executable='../src/apps/dggrid/dggrid', working_dir=tmp_dir, capture_logs=capture_logs, silent=silent )
    dggs = dgselect(dggs_type = dggs_type, res= resolution, mixed_aperture_level=mixed_aperture_level)
    
    subset_conf = { 'update_frequency': 100000, 'clip_subset_type': 'WHOLE_EARTH' }
    seq_df = None

    if not cell_id_list is None and len(cell_id_list) > 0:

        seq_df = pd.DataFrame({ 'seqs': cell_id_list})
        seq_df.to_csv( str( (Path(tmp_dir) / f"temp_clip_{tmp_id}.txt").resolve()) , header=False, index=False, columns=['seqs'], sep=' ')

        subset_conf.update({
            'clip_subset_type': 'SEQNUMS',
            'clip_region_files': str( (Path(tmp_dir) / f"temp_clip_{tmp_id}.txt").resolve()),
            })

    output_conf = {
        'cell_output_type': 'GDAL',
        'cell_output_gdal_format' : 'GeoJSON',
        'cell_output_file_name': f"temp_{dggs_type}_{resolution}_out_{tmp_id}.geojson",
        }
    
    dggs_ops = dggrid_instance.grid_gen(dggs, subset_conf, output_conf )

    gdf = gpd.read_file( Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.geojson", driver='GeoJSON' )

    if not cell_id_list is None and len(cell_id_list) > 0 and not seq_df is None:
        seq_df['cell_exists'] = True
        seq_df.set_index('seqs', inplace=True)
        gdf = gdf.join( seq_df, how='inner', on='Name')
        gdf = gdf.loc[gdf['cell_exists']].drop(columns=['cell_exists'])

    try:
        os.remove( str( Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.geojson") )
        os.remove( str( Path(tmp_dir) / f"temp_clip_{tmp_id}.txt") )
    except Exception:
        pass

    return gdf


def cellids_for_geo_points():
    pass


def proto_helpers(dggrid_instance):

    """
    - grid_cell_polygons_for_extent(): fill extent/subset with cells at resolution (clip or world)
    - grid_cell_polygons_from_cellids(): geometry_from_cellid for dggs at resolution (from id list)
    - grid_cellids_for_extent(): get_all_indexes/cell_ids for dggs at resolution (clip or world)
    - poly_outline for point/centre at resolution?
    - get parent_for_cell_id at coarser resolution
    - get children_for_cell_id at finer resolution
    
    - sample raster values into s2 dggs cells
    - sample vector values into s2 dggs cells
    """


if __name__ == '__main__':

    est_bound = shapely.geometry.box(20.37,57.52, 28.2,60.0 )
    
    # gdf1 = grid_cell_polygons_for_extent('ISEA4T', 3, clip_geom=None, silent=False)
    # print(gdf1.head())

    # gdf2 = grid_cell_polygons_for_extent('ISEA7H', 5, clip_geom=est_bound, silent=True)
    # print(gdf2.head())
    # gdf2.to_file('/tmp/grids/est_shape_isea7h_5.shp')

    # gdf3 = grid_cell_polygons_for_extent('ISEA7H', 8, clip_geom=est_bound)
    # print(gdf3.head())
    # gdf3.to_file('/tmp/grids/est_shape_isea7h_8.shp')

    # df1 = grid_stats_table('ISEA7H', 8, silent=False)
    # print(df1.head(8))
    # df1.to_csv('/tmp/grids/eisea7h_8_stats.csv', index=False)

    df2 = grid_cellids_for_extent('ISEA7H', 5, clip_geom=est_bound, capture_logs=False, silent=True)
    print(df2)
    df2.to_csv('/tmp/grids/est_isea7h_5_seqnums.csv', index=False, header=None)

    cell_list_est = df2[0].values
    print(cell_list_est)
    print(cell_list_est.shape)

    gdf4 = grid_cell_polygons_from_cellids(cell_list_est, 'ISEA7H', 5, capture_logs=False, silent=True)
    print(gdf4.head())
    gdf4.to_file('/tmp/grids/seqnums_isea7h_3.shp')


