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

.. code:: python

   import geopandas
   import shapely

   from dggrid4py import DGGRIDv7

   # create an inital instance that knows where the dggrid tool lives, configure temp workspace and log/stdout output
   dggrid_instance = DGGRIDv7(executable='<path_to>/dggrid', working_dir='.', capture_logs=False, silent=False)

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

TODO
----

Contributions are welcome:

-  sample vector values into dggs cells (aka binning)

-  sample raster values into dggs cells

-  get parent_for_cell_id at coarser resolution

-  get children_for_cell_id at finer resolution

