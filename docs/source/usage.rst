Usage
=====

.. _installation:

Installation
------------

To use dggrid4py, first install it using pip:

.. code-block:: console

   (.venv) $ pip install dggrid4py


You need the ddgrid tool compiled available on the system.

Besides some lowlevel access influence the dggrid operations’ metafile
creation, a few highlevel functions are integrated to work with the more
comfortable geopython libraries, like shapely and geopandas

-  grid_cell_polygons_for_extent(): fill extent/subset with cells at
   resolution (clip or world)
-  grid_cell_polygons_from_cellids(): geometry_from_cellid for dggs at
   resolution (from id list)
-  grid_cellids_for_extent(): get_all_indexes/cell_ids for dggs at
   resolution (clip or world)
-  cells_for_geo_points(): poly_outline for point/centre at resolution
-  address_transform():  conversion betwenn cell_id address types, like SEQNUM, Z7, or Q2DI

.. code:: python

   import geopandas
   import shapely

   from dggrid4py import DGGRIDv7

   # create an inital instance that knows where the dggrid tool lives, configure temp workspace and log/stdout output
   dggrid_instance = DGGRIDv7(executable='<path_to>/dggrid', working_dir='.', capture_logs=False, silent=False, tmp_geo_out_legacy=False, debug=False)

   # global ISEA4T grid at resolution 5 into GeoDataFrame to Shapefile
   gdf1 = dggrid_instance.grid_cell_polygons_for_extent('ISEA4T', 5)
   print(gdf1.head())
   gdf1.to_file('isea4t_5.shp')

   gdf_centroids = dggrid_instance.grid_cell_centroids_for_extent(dggs_type='ISEA7H', resolution=4, mixed_aperture_level=None, clip_geom=None)
   print(gdf_centroids.head())

   # clip extent
   clip_bound = shapely.geometry.box(20.2,57.00, 28.4,60.0 )

   # ISEA7H grid at resolution 9, for extent of provided WGS84 rectangle into GeoDataFrame to Shapefile
   gdf3 = dggrid_instance.grid_cell_polygons_for_extent('ISEA7H', 9, clip_geom=est_bound)
   print(gdf3.head())
   gdf3.to_file('grids/est_shape_isea7h_9.shp')

   # generate cell and areal statistics for a ISEA7H grids from resolution 0 to 8 (return a pandas DataFrame)
   df1 = dggrid_instance.grid_stats_table('ISEA7H', 8)
   print(df1.head(8))
   df1.to_csv('isea7h_8_stats.csv', index=False)

   # generate the DGGS grid cells that would cover a GeoDataFrame of points, return Polygons with cell IDs as GeoDataFrame
   gdf4 = dggrid_instance.cells_for_geo_points(geodf_points_wgs84, False, 'ISEA7H', 5)
   print(gdf4.head())
   gdf4.to_file('polycells_from_points_isea7h_5.shp')

   # generate the DGGS grid cells that would cover a GeoDataFrame of points, return cell IDs added as column to the points GDF
   gdf5 = dggrid_instance.cells_for_geo_points(geodf_points_wgs84=geodf_points_wgs84, cell_ids_only=True, dggs_type='ISEA4H', resolution=8)
   print(gdf5.head())
   gdf5.to_file('geopoint_cellids_from_points_isea4h_8.shp')

   # generate DGGS grid cell polygons based on 'cell_id_list' (a list or np.array of provided cell_ids)
   gdf6 = dggrid_instance.grid_cell_polygons_from_cellids(cell_id_list=[1, 4, 8], 'ISEA7H', 5)
   print(gdf6.head())
   gdf6.to_file('from_seqnums_isea7h_5.shp')

   # v0.2.6 API update split at dateline for cartesian GIS tools
   gdf7 = dggrid_instance.grid_cell_polygons_for_extent('ISEA7H', 3, split_dateline=True)
   gdf7.to_file('global_isea7h_3_interrupted.shp')

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


Portable DGGRID binary
----------------------

If you don't have a special local distribution of the dggrid-tool or if you didn't install with conda-forge, you can use a provided portable:

.. code:: python

   from dggrid4py import tool

   dggrid_exec = tool.get_portable_executable(".")
   dggrid_instance_portable = DGGRIDv7(executable=dggrid_exec, working_dir='.', capture_logs=False, silent=True, has_gdal=False, tmp_geo_out_legacy=True, debug=False)



TODO
----

Contributions are welcome.

-  get parent_for_cell_id at coarser resolution

-  get children_for_cell_id at finer resolution

With the IGEO7/Z7 index system parent/child relations are available.
