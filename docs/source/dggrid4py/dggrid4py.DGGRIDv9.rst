dggrid4py.DGGRIDv9
==================

The DGGRIDv9 class is the class for DGGRID 9.1 or newer, and it raises a ``ValueError`` with an older DGGRID. In contrast to the :doc:`dggrid4py.DGGRIDv8` class, it uses the DGGS types of DGGRID 9.1 as they are defined there, which differs in the following points:

- only the parameters that are given and that the operation needs are written to the metafile, while everything else comes from the preset of DGGRID (``DGGRIDv8`` always writes the projection, aperture, topology and orientation);
- ``IGEO7`` runs as ``IGEO7v2``, which applies ``dggs_vert0_lon = 11.20`` and takes WGS84 coordinates in input and output, so no conversion of the latitudes is needed (see :ref:`igeo7_dggrid91`);
- the cell identifiers have the form of the preset by default, which is ``HIERNDX`` with the ``INT64`` form for a preset with a hierarchical index and ``SEQNUM`` otherwise;
- a parameter that DGGRID 9.1 does not know, or that dggrid4py sets by itself (e.g. the file names), raises a ``ValueError``, while ``DGGRIDv8`` ignores it without an error;
- the coordinate reference system of the returned GeoDataFrame is the one that DGGRID declares, which is ``EPSG:4326`` for an ellipsoidal preset and DGGRID's authalic sphere for a spherical one (``DGGRIDv8`` sets ``EPSG:4326`` for all).

The names of the DGGS types follow DGGRID: ``S`` stands for the spherical form on the authalic sphere (e.g. ``ISEA7HS``, the former name ``ISEA7H`` is an alias), ``L`` for the ellipsoidal form with WGS84 coordinates in input and output (``ISEA7HL``), and a final ``Z`` for a hierarchical index (``ISEA7HLZ``, which is the same as ``IGEO7v2``). ``dggrid4py.dggs_types_v9`` lists all of them. ``H3`` is not included, because the available DGGRID builds do not support it. The full list of parameters is given in the `DGGRID manual <https://dggrid.readthedocs.io/v9.1/dggrid_manual.html>`__.

.. code:: python

   import shutil
   import tempfile

   import shapely

   from dggrid4py import DGGRIDv9

   dggrid_instance = DGGRIDv9(executable=shutil.which("dggrid"), working_dir=tempfile.mkdtemp(), capture_logs=False, silent=True)

   # ISEA7H on the WGS84 ellipsoid, EPSG:4326 in and out
   extent = shapely.box(26.664593, 58.348705, 26.785607, 58.422495)
   cells = dggrid_instance.grid_cell_polygons_for_extent("ISEA7HL", 9, clip_geom=extent)

   # new parameters of DGGRID 9.1 are passed as keyword arguments
   cells = dggrid_instance.grid_cell_polygons_for_extent("ISEA7HL", 9, clip_geom=extent, polyclip_algo="CLIPPER1")

The IVEA types use the Icosahedral Vertex-oriented great-circle Equal Area (IVEA) projection instead of ISEA, in the same spherical and ellipsoidal forms (e.g. ``IVEA7HS``, ``IVEA7HL`` and ``IVEA7HLZ``).

The Fuller types need attention when moving from ``DGGRIDv8``, because DGGRID 9.1 changes their default orientation to DYMAXION. All cell identifiers therefore differ from the ones of ``DGGRIDv8``, which keeps the former orientation, while ``dggs_orient_preset="ISEA"`` returns the former identifiers. The low-level functions (``dgapi_*``) take a ``Dggs`` object, which is created with ``dgselect(dggs_type, minimal=True, res=...)`` for the same behaviour as above.

.. currentmodule:: dggrid4py

.. autoclass:: DGGRIDv9
   :members:
   :inherited-members:

   
   .. automethod:: __init__

   
   .. rubric:: Methods

   .. autosummary::

      ~DGGRIDv9.__init__
      ~DGGRIDv9.is_runnable
      ~DGGRIDv9.check_gdal_support
      ~DGGRIDv9.post_process_split_dateline
      ~DGGRIDv9.run

      ~DGGRIDv9.grid_cell_polygons_for_extent
      ~DGGRIDv9.grid_cell_centroids_for_extent
      ~DGGRIDv9.grid_cell_polygons_from_cellids
      ~DGGRIDv9.grid_cell_centroids_from_cellids
      ~DGGRIDv9.grid_cellids_for_extent
      ~DGGRIDv9.cells_for_geo_points
      ~DGGRIDv9.grid_stats_table

      ~DGGRIDv9.dgapi_grid_gen
      ~DGGRIDv9.dgapi_grid_stats
      ~DGGRIDv9.dgapi_grid_transform
      ~DGGRIDv9.dgapi_point_value_binning
      ~DGGRIDv9.dgapi_pres_binning

      ~DGGRIDv9.cells_for_geo_points
      ~DGGRIDv9.address_transform
      ~DGGRIDv9.cell_parents
      ~DGGRIDv9.cell_children
      ~DGGRIDv9.cell_neighbours
