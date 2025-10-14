dggrid4py.DGGRIDv8
==================

The DGGRIDv8 class has been introduced to provide custom parameters to the DGGRID version 8.x series.

For "default" or "classical" behavior of dggrid4py parameter handling, the DGGRIDv7 class should be used.

How to pass additional parameters to DGGRIDv8:

.. image:: config_v8_extra.png

.. currentmodule:: dggrid4py

.. autoclass:: DGGRIDv8

   
   .. automethod:: __init__

   
   .. rubric:: Methods

   .. autosummary::

      ~DGGRIDv8.__init__
      ~DGGRIDv8.is_runnable
      ~DGGRIDv8.check_gdal_support
      ~DGGRIDv8.post_process_split_dateline
      ~DGGRIDv8.run

      ~DGGRIDv8.grid_cell_polygons_for_extent
      ~DGGRIDv8.grid_cell_centroids_for_extent
      ~DGGRIDv8.grid_cell_polygons_from_cellids
      ~DGGRIDv8.grid_cell_centroids_from_cellids
      ~DGGRIDv8.grid_cellids_for_extent
      ~DGGRIDv8.cells_for_geo_points
      ~DGGRIDv8.grid_stats_table

      ~DGGRIDv8.dgapi_grid_gen
      ~DGGRIDv8.dgapi_grid_stats
      ~DGGRIDv8.dgapi_grid_transform
      ~DGGRIDv8.dgapi_point_value_binning
      ~DGGRIDv8.dgapi_pres_binning

      ~DGGRIDv8.cells_for_geo_points
      ~DGGRIDv8.address_transform
   
   
   
