dggrid4py.igeo7_ext
===================

This extension module provides the IGEO7 wrappers around :class:`dggrid4py.DGGRIDv8`. The wrappers apply the IGEO7 definition
(``dggs_vert0_lon = 11.20``, Z7 hierarchical index, conversion between WGS84 and authalic latitudes), see :doc:`../IGEO7` for
the background and examples. The module also contains experimental helper functions for raster data.

The :mod:`dggrid4py.igeo7` module provides the functions on the Z7 index itself.

.. automodule:: dggrid4py.igeo7_ext
   :members:

   .. rubric:: Functions

   .. autosummary::
   
      igeo7_meta_config
      dggrid_igeo7_grid_cell_polygons_for_extent
      dggrid_igeo7_grid_cell_polygons_from_cellids
      dggrid_igeo7_grid_cell_centroids_from_cellids
      dggrid_igeo7_cells_for_geo_points
      dggrid_igeo7_q2di_from_cellids
      dggrid_get_res
      z7_resolution
      z7_k1_ring_neighbours
      to_parent_series
      z7_base_pentagons
      z7_get_base_pentagon
      z7_is_pentagon
      suggest_window_blocks_per_chunk
      extract_windows_with_bounds
      get_crs_info
      projected_distance
      get_raster_pixel_edge_len
      propose_dggs_level_for_pixel_length
      create_geopoints_for_window
