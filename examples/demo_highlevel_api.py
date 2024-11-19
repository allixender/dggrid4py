# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 - Alexander Kmoch
# Licenced under GNU AFFERO GENERAL PUBLIC LICENSE. Please consult the LICENCE 
# file for details.
#
# Author: Alexander Kmoch (alexander.kmoch@ut.ee)
# Date: 07-12-2022 
#

import geopandas as gpd
import shapely

import dggrid4py
from dggrid4py import DGGRIDv7, dggs_types


def highlevel_grid_gen_and_transform(dggrid_instance):

    """
    - grid_cell_polygons_for_extent(): fill extent/subset with cells at resolution (clip or world)
    - grid_cell_polygons_from_cellids(): geometry_from_cellid for dggs at resolution (from id list)
    - grid_cellids_for_extent(): get_all_indexes/cell_ids for dggs at resolution (clip or world)
    - cells_for_geo_points(): poly_outline for point/centre at resolution
    """

    est_bound = shapely.geometry.box(20.2,57.00, 28.4,60.0 )

    gdf1 = dggrid_instance.grid_cell_polygons_for_extent('ISEA4T', 5, clip_geom=est_bound)
    print(gdf1.head(3))
    # gdf1.to_file('/tmp/est_shape_isea4t_5.shp')

    gdf2 = dggrid_instance.grid_cell_polygons_for_extent('ISEA7H', 5, clip_geom=est_bound)
    print(gdf2.head(3))
    # gdf2.to_file('/tmp/est_shape_isea7h_5.shp')

    gdf2_a = dggrid_instance.grid_cell_polygons_for_extent('ISEA7H', 6, clip_geom=est_bound)
    print(gdf2_a.head(3))
    # gdf2_a.to_file('/tmp/est_shape_isea7h_6.shp')

    gdf3 = dggrid_instance.grid_cell_polygons_for_extent('ISEA7H', 8, clip_geom=est_bound, output_address_type='Z7')
    # print(gdf3.head())
    # gdf3.to_file('/tmp/est_shape_isea7h_8.shp')

    gdf3_a = dggrid_instance.grid_cell_polygons_for_extent('ISEA3H', 9, clip_geom=est_bound, output_address_type='Z3_STRING')
    print(gdf3.head(3))
    # gdf3_a.to_file('/tmp/est_shape_isea7h_9.shp')

    gdf_centroids = dggrid_instance.grid_cell_centroids_for_extent(dggs_type='ISEA7H', resolution=4, mixed_aperture_level=None, clip_geom=None)
    
    df1 = dggrid_instance.grid_stats_table('ISEA7H', 20)
    print(df1.head(8))
    # df1.to_csv('/tmp/eisea7h_8_stats.csv', index=False)

    df2 = dggrid_instance.grid_cellids_for_extent('ISEA7H', 5, clip_geom=est_bound, output_address_type='SEQNUM')
    print(df2.head(3))
    # df2.to_csv('/tmp/est_isea7h_5_gridgen_from_seqnums.csv', index=False, header=None)

    cell_list_est = df2[0].values
    # # print(cell_list_est)
    # # print(cell_list_est.shape)

    gdf4 = dggrid_instance.grid_cell_polygons_from_cellids(cell_list_est, 'ISEA7H', 5)
    print(gdf4.head(3))
    # gdf4.to_file('/tmp/from_seqnums_isea7h_5.shp')

    gdf4 = dggrid_instance.grid_cell_polygons_from_cellids(cell_list_est, 'ISEA7H', 5, clip_subset_type='SEQNUMS', input_address_type='SEQNUM')
    print(gdf4.head(3))
    # gdf4.to_file('/tmp/from_seqnums_isea7h_5.shp')

    gdf4['centroid_geo'] = gdf4['geometry'].centroid

    tgeo = gdf4.copy()
    tgeo = tgeo.rename(columns={'geometry': 'old_geo', 'centroid_geo': 'geometry' }).drop(columns=['old_geo'])
    geodf_points_wgs84 = gpd.GeoDataFrame(tgeo, geometry='geometry', crs=4326)

    gdf5 = dggrid_instance.cells_for_geo_points(geodf_points_wgs84, False, 'ISEA7H', 5)
    print(gdf5.head(3))
    # gdf5.to_file('/tmp/polycells_from_points_isea7h_5.shp')

    gdf6 = dggrid_instance.cells_for_geo_points(geodf_points_wgs84, True, 'ISEA7H', 5)
    print(gdf6.head(3))
    # gdf6.to_file('/tmp/geopoint_cellids_from_points_isea7h_5.shp')

    # v0.2.6 API update split at dateline for cartesian GIS tools
    gdf7 = dggrid_instance.grid_cell_polygons_for_extent('ISEA7H', 3, split_dateline=True)
    print(gdf7.head(3))
    # gdf7.to_file('/tmp/global_isea7h_3_interrupted.shp')

    gdf_z1 = dggrid_instance.grid_cell_polygons_for_extent('IGEO7', 5, clip_geom=est_bound, output_address_type='Z7_STRING')
    print(gdf_z1.head(3))

    df_z1 = dggrid_instance.guess_zstr_resolution(gdf_z1['name'].values, 'IGEO7', input_address_type='Z7_STRING')
    print(df_z1.head(3))

    df_q2di = dggrid_instance.address_transform(gdf_z1['name'].values, 'IGEO7', 5, input_address_type='Z7_STRING', output_address_type='Q2DI')
    print(df_q2di.head(3))

    df_tri = dggrid_instance.address_transform(gdf_z1['name'].values, 'IGEO7', 5, input_address_type='Z7_STRING', output_address_type='PROJTRI')
    print(df_tri.head(3))

    children = dggrid_instance.grid_cell_polygons_from_cellids(
        cell_id_list=['00012502340'],    # the input/parent cell id
        dggs_type='IGEO7',               # dggs type
        resolution=11,                   # target resolution of children
        clip_subset_type='COARSE_CELLS', # new parameter
        clip_cell_res=9,                 # resolution of parent cell
        input_address_type='Z7_STRING',  # address_type
        output_address_type='Z7_STRING'  # address_type
    )
    print(children.head(3))



def highlevel_grid_stats(dggrid_instance):

    # generate grid_stats for all predefined DGGS types (except CUSTOM, PLANETGRID will not complete either at higher resolutions)
    for gridtype in filter(lambda x: x.startswith('CUSTOM') == False, dggs_types):

        mixed_aperture_level=None
        # for mixed aperture eg ISEA34H define the level of switch, see dggrid manual
        if '43' in gridtype:
            mixed_aperture_level = 7

        try:
            print(f"{gridtype} - {15}")
            df = dggrid_instance.grid_stats_table(dggs_type=gridtype, resolution=15, mixed_aperture_level=mixed_aperture_level, )
            print(df.head(10))
            # df.to_csv(f"/tmp/{gridtype}_{15}_stats.csv", index=False)
        except ValueError as ex:
            print(ex)
            pass


if __name__ == '__main__':

    print(dggrid4py.__version__)

    # DGGRID from https://github.com/sahrk/DGGRID
    # with a /tmp dir, e.g. on Linux/Mac
    dggrid_path = '/usr/local/bin/dggrid'

    import os

    if not os.environ['DGGRID_PATH'] is None:
        dggrid_path = os.environ['DGGRID_PATH']
    
    debug_mode = False
    if not os.environ['DGGRID_DEBUG'] is None:
        debug_mode = True if str(os.environ['DGGRID_DEBUG']).lower() == 'true' else False

    dggrid = DGGRIDv7(executable=dggrid_path, working_dir='/tmp', capture_logs=True, silent=False)

    highlevel_grid_gen_and_transform(dggrid)

    # highlevel_grid_stats(dggrid)

