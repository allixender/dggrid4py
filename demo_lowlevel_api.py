# -*- coding: utf-8 -*-

from pathlib import Path

import pandas as pd

import geopandas as gpd
from shapely.geometry import Polygon, box, shape

from dggrid4py import DGGRIDv7, Dggs, dgselect, dggs_types


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
    dggs_ops = dggrid_instance.dgapi_grid_gen(dggs, subset_conf, output_conf )
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

    dggs_ops = dggrid_instance.dgapi_grid_gen(dggs, subset_conf, output_conf )

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

    dggs_ops = dggrid_instance.dgapi_grid_gen(dggs, subset_conf, output_conf )

    return dggs_ops


def generate_defs_1_5(dggrid_instance):

    for gridtype in filter(lambda x: x.startswith('CUSTOM') == False, dggs_types):
        for res in range(1,6):
            print(f"dgselect(dggs_type = '{gridtype}', res= {res})")
            fname = f"DGGRID_{gridtype}_res_{res}"

            if '43' in gridtype and res >= 4:
                mixed_aperture_level = 4

            dggs = dgselect(dggs_type = gridtype, res= res)

            subset_conf = {
                'clip_subset_type': 'WHOLE_EARTH',
                'update_frequency': 100000,
                'geodetic_densify': 0.0
                }

            output_conf = {
                'cell_output_type': 'GEOJSON',
                'cell_output_file_name': fname,
                # 'densification': 6
                }

            dggs_ops = dggrid_instance.dgapi_grid_gen(dggs, subset_conf, output_conf )

            print(dggs_ops)



def example_table_stats(dggrid_instance):
    """
    # specify the operation
    dggrid_operation OUTPUT_STATS

    # specify the DGG
    dggs_type ISEA43H
    dggs_num_aperture_4_res 5
    dggs_res_spec 15
    """

    dggs = dgselect(dggs_type = 'ISEA43H', res= 10, mixed_aperture_level= 5)
    dggs_ops = dggrid_instance.dgapi_grid_stats(dggs)

    return dggs_ops


def example_transform_geo_to_seqnum(dggrid_instance, example_src):
    """
    # specify the operation
    dggrid_operation TRANSFORM_POINTS

    # specify the DGG
    dggs_type ISEA3H
    dggs_res_spec 9

    # specify bin controls
    input_file_name inputfiles/20k.txt
    input_address_type GEO
    input_delimiter " "

    # specify the output
    output_file_name outputfiles/cities3h9.txt
    output_address_type SEQNUM
    output_delimiter ","
    """

    dggs = dgselect(dggs_type = 'ISEA3H', res= 9)

    subset_conf = {
        'input_file_name':  str(Path(example_src / "transform/inputfiles/20k.txt").resolve()),
        'input_address_type': 'GEO',
        'input_delimiter': "\" \""
        }

    output_conf = {
        'output_file_name': 'cities3h9.txt',
        'output_address_type': 'SEQNUM',
        'output_delimiter': "\",\""
        }

    dggs_ops = dggrid_instance.dgapi_grid_transform(dggs, subset_conf, output_conf)

    return dggs_ops


def example_transform_reverse_seqnum_to_any(dggrid_instance, out_type):

    dggs = dgselect(dggs_type = 'ISEA3H', res= 9)

    subset_conf = {
        'input_file_name':  "/tmp/grids/cities3h9.txt",
        'input_address_type': 'SEQNUM',
        'input_delimiter': "\",\""
        }

    output_conf = {
        'output_file_name': f"cities3h9.out.{out_type}",
        'output_address_type': out_type,
        'output_delimiter': "\" \""
        }

    dggs_ops = dggrid_instance.dgapi_grid_transform(dggs, subset_conf, output_conf)

    return dggs_ops


if __name__ == '__main__':

    dggrid = DGGRIDv7(executable='../src/apps/dggrid/dggrid', working_dir='/tmp/grids', capture_logs=True, silent=False)

    example_src = Path('../examples')

    print( dggrid.is_runnable() == True )

    # ------------- GENERATE_GRID
    result_info = example_generate(dggrid)
    print(result_info)

    result_info = example_aigenerate(dggrid, example_src)
    print(result_info)

    result_info = example_gdal(dggrid, example_src)
    print(result_info)

    # generate_defs_1_5(dggrid)

    # ----------- OUTPUT_STATS ----
    result_info = example_table_stats(dggrid)
    print(result_info)

    import pandas as pd
    print(pd.DataFrame(result_info['output_conf']['stats_output']))

    # ---------- TRANSFORM ---------------
    result_info = example_transform_geo_to_seqnum(dggrid, example_src)
    print(result_info)

    for out_type in [
            'GEO', # geodetic coordinates -123.36 43.22 20300 Roseburg
            'Q2DI', # quad number and (i, j) coordinates on that quad
            'SEQNUM', # DGGS index - linear address (1 to size-of-DGG), not supported for parameter input_address_type if dggs_aperture_type is SEQUENCE
            'INTERLEAVE', # digit-interleaved form of Q2DI, only supported for parameter output_address_type; only available for hexagonal aperture 3 and 4 grids
            'PLANE', # (x, y) coordinates on unfolded ISEA plane,  only supported for parameter output_address_type;
            'Q2DD', # quad number and (x, y) coordinates on that quad
            'PROJTRI', # PROJTRI - triangle number and (x, y) coordinates within that triangle on the ISEA plane
            'VERTEX2DD', # vertex number, triangle number, and (x, y) coordinates on ISEA plane
            'AIGEN'
        ]:
        result_info = example_transform_reverse_seqnum_to_any(dggrid, out_type)
        print(result_info)
