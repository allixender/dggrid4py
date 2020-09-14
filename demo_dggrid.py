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

from dggrid_runner import DGGRIDv7, Dggs, dgselect

    
def example_generate(dggrid_instance):

    try:
        dgselect( dggs_type='CUSTOM',
                projection='ISEA', 
                aperture=3,
                topology='HEXAGON',
                res=2)
    except ValueError as ex:
        print(ex)
        pass
    
    dggs = dgselect(dggs_type = 'ISEA3H',
                res= 2)
    
    subset_conf = {
        'clip_subset_type' : 'WHOLE_EARTH',
        # 'update_frequency' : 100000,
        # 'geodetic_densify' : 0.0
        }

    output_conf = {
        'cell_output_type' : 'SHAPEFILE',
        'cell_output_file_name' : 'dggrid_isea3h_grid',
        # 'densification' : 5
        }

    dggs_ops = dggrid_instance.grid_gen(dggs, subset_conf, output_conf )

    return dggs_ops
    

def example_aigenerate(dggrid_instance, example_src):
    """
    # specify the operation
    dggrid_operation GENERATE_GRID

    # specify a ISEA3H; override the default resolution
    dggs_type ISEA3H
    dggs_res_spec 6

    # control grid generation
    clip_subset_type SHAPEFILE
    clip_region_files ./inputfiles/orbuff.shp
    update_frequency 10000000

    # specify the output
    cell_output_type AIGEN
    cell_output_file_name ./outputfiles/orCells
    densification 5

    point_output_type AIGEN
    point_output_file_name ./outputfiles/orPts
    """

    dggs = dgselect(dggs_type = 'ISEA3H', res= 6)
    
    subset_conf = {
        'clip_subset_type': 'SHAPEFILE',
        'clip_region_files': str(Path(example_src / "aigenerate/inputfiles/orbuff.shp").resolve()),
        'update_frequency': 10000000
        }

    output_conf = {
        'cell_output_type': 'AIGEN',
        'cell_output_file_name': "orCells",
        'densification': 5,

        'point_output_type': 'AIGEN',
        'point_output_file_name': "orPts"
        }

    dggs_ops = dggrid_instance.grid_gen(dggs, subset_conf, output_conf )

    return dggs_ops


def example_gdal(dggrid_instance, example_src):
    """
    # specify the operation
    dggrid_operation GENERATE_GRID

    # specify the DGG
    dggs_type ISEA7H
    dggs_res_spec 9

    # control grid generation
    clip_subset_type GDAL
    # gdal uses the file extension to determine the file format
    clip_region_files inputfiles/corvallis.shp

    # specify the output
    cell_output_type GDAL
    # for output the format needs to be explicit
    cell_output_gdal_format KML
    cell_output_file_name ./outputfiles/corvallisCells.kml

    point_output_type GDAL
    point_output_gdal_format GeoJSON
    point_output_file_name ./outputfiles/corvallisPts.geojson
    """

    dggs = dgselect(dggs_type = 'ISEA7H', res= 9)
    
    subset_conf = {
        'clip_subset_type': 'GDAL',
        'clip_region_files': str(Path(example_src / "gdalExample/inputfiles/corvallis.shp").resolve()),
        'update_frequency': 10000000
        }

    output_conf = {
        'cell_output_type': 'GDAL',
        'cell_output_gdal_format' : 'KML',
        'cell_output_file_name': "corvallisCells.kml",

        'point_output_type': 'GDAL',
        'point_output_gdal_format': 'GeoJSON',
        'point_output_file_name': "corvallisPts.geojson"
        }

    dggs_ops = dggrid_instance.grid_gen(dggs, subset_conf, output_conf )

    return dggs_ops


if __name__ == '__main__':

    dggrid = DGGRIDv7(executable='../src/apps/dggrid/dggrid', working_dir='/tmp/grids')

    example_src = Path('../examples')

    print( dggrid.is_runnable() == True )

    result_info = example_generate(dggrid)
    print(result_info)

    result_info = example_aigenerate(dggrid, example_src)
    print(result_info)

    result_info = example_gdal(dggrid, example_src)
    print(result_info)
    






