dggrid4py.DGGRIDv7
==================

.. deprecated:: 0.6.0
   The DGGRIDv7 class is deprecated and will be removed in a future version. Use :doc:`dggrid4py.DGGRIDv8`
   with DGGRID 8.42 or newer; it accepts the former address type names (e.g. ``Z7_STRING``) with a
   ``DeprecationWarning``, see :doc:`../IGEO7`.

.. currentmodule:: dggrid4py

.. autoclass:: DGGRIDv7

   
   .. automethod:: __init__

   
   .. rubric:: Methods

   .. autosummary::
   
      ~DGGRIDv7.__init__
      ~DGGRIDv7.is_runnable
      ~DGGRIDv7.check_gdal_support
      ~DGGRIDv7.post_process_split_dateline
      ~DGGRIDv7.run

      ~DGGRIDv7.grid_cell_polygons_for_extent
      ~DGGRIDv7.grid_cell_centroids_for_extent
      ~DGGRIDv7.grid_cell_polygons_from_cellids
      ~DGGRIDv7.grid_cell_centroids_from_cellids
      ~DGGRIDv7.grid_cellids_for_extent
      ~DGGRIDv7.cells_for_geo_points
      ~DGGRIDv7.grid_stats_table

      ~DGGRIDv7.dgapi_grid_gen
      ~DGGRIDv7.dgapi_grid_stats
      ~DGGRIDv7.dgapi_grid_transform
      ~DGGRIDv7.dgapi_point_value_binning
      ~DGGRIDv7.dgapi_pres_binning

      ~DGGRIDv7.cells_for_geo_points
      ~DGGRIDv7.address_transform