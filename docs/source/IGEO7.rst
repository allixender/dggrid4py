IGEO7
=====

Overview: A Hierarchically Indexed Hexagonal Equal-Area Discrete Global Grid System
-----------------------------------------------------------------------------------

Hexagonal Discrete Global Grid Systems (DGGS) offer significant advantages for spatial analysis due to their uniform cell shapes and efficient indexing. Among the three central place apertures (3, 4, and 7), aperture 7 subdivisions exhibit very desirable properties, including the preservation of hexagonal symmetry and the formation of unambiguous indexing hierarchies.

Interest in hierarchically indexed aperture 7 hexagonal DGGS has recently increased due to the popularity of the H3 DGGS. But there are currently no open-source equal-area aperture 7 hexagonal DGGS available that provide indexing capabilities similar to H3.

We present **IGEO7**, a novel pure aperture 7 hexagonal DGGS, and **Z7**, its associated hierarchical integer indexing system. In contrast to H3, where cell sizes vary by up to ±50% across the globe, IGEO7 uses cells of equal area, making it a true equal-area DGGS.

IGEO7 and Z7 are implemented in the open-source software DGGRID. We also present a use case for on-demand suitability modeling to demonstrate a practical application of this new DGGS.


.. _igeo7_definition:

Orientation and ellipsoid: which IGEO7 you get
----------------------------------------------

IGEO7 is available in DGGRID since version 8.41 as the ``IGEO7`` DGGS type, which is an ISEA7H grid with the Z7 index. In the DGGRID 8 series this type inherits the orientation of ISEA7H (``dggs_vert0_lon = 11.25``), and DGGRID takes all coordinates as longitude and latitude on the authalic sphere.

Recent developments in the OGC API DGGS standard draft make reference to refinement ratio 7 DGGS with the ISEA projection. In order to align IGEO7 with this work, the IGEO7 used in dggrid4py, on `igeo7.org <https://igeo7.org>`_ and in the related data stores applies two adjustments to the DGGRID preset:

- the base icosahedron is rotated by 0.05 degrees (``dggs_vert0_lon = 11.20`` instead of ``11.25``), which places its vertices better over water bodies;
- WGS84 coordinates are converted to authalic latitudes before they are passed to DGGRID, and the DGGRID output is converted back, so that the grid refers to the WGS84 ellipsoid and not to a spherical approximation (see :mod:`dggrid4py.auxlat`, which is based on ``pygeodesy``).

Both adjustments are part of the grid definition and not an option, because each of them changes the cell identifiers. The following table shows the Z7 cell at resolution 9 for two locations (longitude, latitude in WGS84) and the four possible combinations, calculated with DGGRID 8.44:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Configuration
     - Lisbon (-9.1393, 38.7223)
     - Tartu (26.7220, 58.3776)
   * - 11.20 with authalic conversion (IGEO7)
     - ``00641565463``
     - ``00010224545``
   * - 11.20 without conversion
     - ``00641565515``
     - ``00010202462``
   * - 11.25 with authalic conversion
     - ``00641542231``
     - ``00010261335``
   * - 11.25 without conversion (DGGRID 8 preset)
     - ``00641543636``
     - ``00010261152``

In practice this means that IGEO7 as described in the original publication (the last row) is not the same grid as the ellipsoid-adjusted IGEO7 (the first row), although both can be generated with dggrid4py. We hope that the required changes become available in DGGRID itself, so that no further confusion can arise. Until then, implementers should state explicitly which orientation and which latitude conversion they use, and the cells in the first row can serve as a check.


Recommended: the wrappers in ``dggrid4py.igeo7_ext``
----------------------------------------------------

The wrappers in :mod:`dggrid4py.igeo7_ext` apply the IGEO7 definition from above, so that WGS84 data goes in and WGS84 data comes out:

- they only accept a ``DGGRIDv8`` instance (a ``DGGRIDv7`` instance raises ``TypeError``);
- they use the Z7 hierarchical index and ``dggs_vert0_lon = 11.20``, see ``igeo7_meta_config()``;
- WGS84 inputs (clip geometries, points) are converted to the authalic sphere before they go to DGGRID, and output geometries are converted back to WGS84.

The example below needs DGGRID 8.42 or newer, e.g. version 8.44 from conda-forge, which is installed with the pixi quickstart in :ref:`installation`. It was tested with DGGRID 8.42, 8.44 and the 9.0b pre-release.

.. code:: python

    import shutil
    import tempfile

    import geopandas as gpd
    import shapely

    from dggrid4py import DGGRIDv8, igeo7_ext

    # inside a pixi/conda environment, dggrid is on the PATH; otherwise pass the path to the DGGRID executable
    dggrid_instance = DGGRIDv8(shutil.which("dggrid"), working_dir=tempfile.mkdtemp(), capture_logs=False, silent=True)

    # Tartu bbox in wgs84
    extent = shapely.box(26.664593, 58.348705, 26.785607, 58.422495)

    # cell polygons for an extent, in and out in wgs84
    cells = igeo7_ext.dggrid_igeo7_grid_cell_polygons_for_extent(extent, 9, dggrid_instance)

    # centroids and polygons for a list of Z7 cell ids (all of the same resolution)
    centroids = igeo7_ext.dggrid_igeo7_grid_cell_centroids_from_cellids(cells["name"], dggrid_instance)
    polygons = igeo7_ext.dggrid_igeo7_grid_cell_polygons_from_cellids(cells["name"], dggrid_instance)

    # Z7 cell ids for wgs84 points, returns a copy of the points with a 'name' column
    points = gpd.GeoDataFrame(
        {"site": ["Lisbon", "Tartu"]},
        geometry=gpd.points_from_xy([-9.1393, 26.7220], [38.7223, 58.3776]),
        crs=4326,
    )
    points = igeo7_ext.dggrid_igeo7_cells_for_geo_points(points, 9, dggrid_instance)
    print(points["name"].tolist())   # ['00641565463', '00010224545']

    # Q2DI addresses and direct neighbours
    q2di = igeo7_ext.dggrid_igeo7_q2di_from_cellids(cells["name"], dggrid_instance)
    cls_m = igeo7_ext.dggrid_get_res(dggrid_instance, "IGEO7", 9).loc[9, "cls_m"]
    neighbours = igeo7_ext.z7_k1_ring_neighbours(cells["name"].iloc[0], dggrid_instance, cls_m)

``hier_ndx_form='INT64'`` switches all wrappers to the 64 bit Z7 form. DGGRID writes this form as a 16 character hex string (e.g. ``'0042529bffffffff'``), which :mod:`dggrid4py.igeo7` can decode (``igeo7.z7hex_to_z7string('0042529bffffffff')`` returns ``'00010224515'``). Any other DGGRID parameter can be passed as keyword argument and overrides the defaults, e.g. ``dggs_vert0_lon=11.25`` for the orientation of the DGGRID preset.


Explicit configuration with ``DGGRIDv8``
----------------------------------------

The functions of the ``DGGRIDv8`` class can also be used directly, for example if a function has no wrapper yet. In this case the same two adjustments have to be applied by hand. ``igeo7_ext.igeo7_meta_config()`` returns the DGGRID parameters for the Z7 index and the orientation as a dictionary:

.. code:: python

    from dggrid4py.auxlat import geoseries_to_authalic, geoseries_to_geodetic

    meta_config = igeo7_ext.igeo7_meta_config()
    # {
    #     "input_address_type": "HIERNDX",           # cell ids are passed in as a hierarchical index,
    #     "input_hier_ndx_system": "Z7",             # of the Z7 system,
    #     "input_hier_ndx_form": "DIGIT_STRING",     # in the textual form (e.g. '00010224545')
    #     "output_address_type": "HIERNDX",          # and the same for the cell ids that come back
    #     "output_cell_label_type": "OUTPUT_ADDRESS_TYPE",
    #     "output_hier_ndx_system": "Z7",
    #     "output_hier_ndx_form": "DIGIT_STRING",
    #     "dggs_vert0_lon": 11.2,                    # the DGGRID default is 11.25
    # }

``input_hier_ndx_form`` and ``output_hier_ndx_form`` accept either ``DIGIT_STRING`` or ``INT64`` (``igeo7_meta_config('INT64')``). The dictionary is passed to the dggrid4py functions as keyword arguments. The WGS84 input is converted with ``geoseries_to_authalic`` before the call, and the output geometries are converted back with ``geoseries_to_geodetic``:

.. code:: python

    # Tartu bbox in wgs84, around 50 km^2
    extent = gpd.GeoSeries([shapely.box(26.664593, 58.348705, 26.785607, 58.422495)], crs=4326)

    # convert the extent from wgs84 to the authalic sphere
    extent_authalic = geoseries_to_authalic(extent)

    # generate IGEO7 cells with Z7 cell ids for the extent
    igeo7_cells = dggrid_instance.grid_cell_polygons_for_extent("IGEO7", 12, clip_geom=extent_authalic.iloc[0], **meta_config)

    # convert the cell geometries back to wgs84
    igeo7_cells["geometry"] = geoseries_to_geodetic(igeo7_cells.geometry)

For points, only the cell identifiers come back, so there is nothing to convert after the call. ``cells_for_geo_points`` returns a copy of the points it was given, which are the authalic ones here; the cell identifiers are in the same order as the original points.

.. code:: python

    points = gpd.GeoDataFrame(
        {"site": ["Lisbon", "Tartu"]},
        geometry=gpd.points_from_xy([-9.1393, 26.7220], [38.7223, 58.3776]),
        crs=4326,
    )
    points_authalic = points.set_geometry(geoseries_to_authalic(points.geometry))
    cell_ids = dggrid_instance.cells_for_geo_points(points_authalic, True, "IGEO7", 9, **meta_config)
    points["name"] = cell_ids["name"].values   # ['00641565463', '00010224545']

Skipping the authalic conversion matters. In our tests, the WGS84 centroids of resolution 9 cells in Tartu (around 58.4°N) fell into a different cell for every single point when they were passed to DGGRID without conversion.


Parents, children and neighbours
--------------------------------

The Z7 index encodes the hierarchy in the identifier itself. In the ``DIGIT_STRING`` form, the first two digits are the base cell and every further digit adds one resolution, so that the resolution and the parent of a cell can be derived without DGGRID:

.. code:: python

    from dggrid4py import igeo7

    cell_id = "00010224545"                        # Tartu at resolution 9
    igeo7.get_z7string_resolution(cell_id)         # 9
    parent, local_pos, is_center = igeo7.get_z7string_local_pos(cell_id)
    # ('0001022454', '5', False)

For the opposite direction, DGGRID offers ``clip_subset_type='COARSE_CELLS'``, which generates the cells of a finer resolution for a list of coarser cells. This is a spatial clip and not a lookup of index children. It returns all finer cells that intersect the coarse cell, which includes cells of the neighbouring parents. In addition, the index descendants of an aperture 7 cell do not nest exactly into the hexagon of their ancestor, so that a clip over several resolutions does not contain all of them (for the cell below, 331 of the 343 descendants at resolution 8). For one resolution step, the index children can be selected by their prefix:

.. code:: python

    coarse = "0001022"                             # the cell around Tartu at resolution 5
    clipped = dggrid_instance.grid_cell_polygons_from_cellids(
        [coarse], "IGEO7", 6, clip_subset_type="COARSE_CELLS", clip_cell_res=5, **meta_config
    )
    children = clipped[clipped["name"].str.startswith(coarse)]
    len(clipped), len(children)                    # (13, 7)

The direct neighbours of a cell are available through ``igeo7_ext.z7_k1_ring_neighbours`` (see the wrapper example above). A function that enumerates the index children over several resolutions is not part of dggrid4py yet.


Migrating from ``Z7_STRING`` and ``DGGRIDv7``
---------------------------------------------

Earlier examples for IGEO7 used the ``DGGRIDv7`` class with ``output_address_type='Z7_STRING'``. These address type names are deprecated in DGGRID 8.44 and removed in DGGRID 9, where the hierarchical indexes are selected with ``HIERNDX`` and two additional parameters:

.. list-table::
   :header-rows: 1
   :widths: 25 25 25 25

   * - ``DGGRIDv7`` address type
     - ``DGGRIDv8`` address type
     - ``hier_ndx_system``
     - ``hier_ndx_form``
   * - ``Z7_STRING``
     - ``HIERNDX``
     - ``Z7``
     - ``DIGIT_STRING``
   * - ``Z7``
     - ``HIERNDX``
     - ``Z7``
     - ``INT64``
   * - ``Z3_STRING``, ``Z3``
     - ``HIERNDX``
     - ``Z3``
     - ``DIGIT_STRING``, ``INT64``
   * - ``ZORDER_STRING``, ``ZORDER``
     - ``HIERNDX``
     - ``ZORDER``
     - ``DIGIT_STRING``, ``INT64``

Since version 0.6.0, the ``DGGRIDv8`` class accepts the old names, maps them according to this table and issues a ``DeprecationWarning``. An address type that is not known at all raises a ``ValueError``. Up to version 0.5.3, ``DGGRIDv8`` ignored an address type that it did not know, without an error, and returned sequence numbers or the cell of another identifier. Results that were produced with ``Z7_STRING`` and ``DGGRIDv8`` in these versions should therefore be checked.

The ``DGGRIDv7`` class keeps the old names and should only be used with a DGGRID 7 executable. Code that moves to ``DGGRIDv8`` continues to work, but it uses the DGGRID preset (11.25, no authalic conversion) as before, as long as the orientation and the conversion are not added as described above.


API Reference
-------------

The IGEO7 functionality is spread over three modules of the dggrid4py Python package.

.. seealso::

   - :mod:`dggrid4py.igeo7_ext` - IGEO7 wrappers around ``DGGRIDv8``, ``igeo7_meta_config()`` and neighbour functions
   - :mod:`dggrid4py.igeo7` - functions on the Z7 index itself (resolution, parent, conversion between the forms)
   - :mod:`dggrid4py.auxlat` - conversion between geodetic (WGS84) and authalic latitudes

Publication
-----------

For more details, see the published paper:

.. seealso::

   Kmoch, A., Sahr, K., Chan, W. T., and Uuemaa, E. (2025). IGEO7: A new hierarchically indexed hexagonal equal-area discrete global grid system. *AGILE GIScience Series*, 6, 32. https://doi.org/10.5194/agile-giss-6-32-2025

   Full text available at: https://agile-giss.copernicus.org/articles/6/32/2025/


If you use IGEO7 in your research, please cite:

.. code-block:: bibtex

   @article{kmoch2025igeo7,
     author = {Kmoch, A. and Sahr, K. and Chan, W. T. and Uuemaa, E.},
     title = {IGEO7: A new hierarchically indexed hexagonal equal-area discrete global grid system},
     journal = {AGILE GIScience Series},
     volume = {6},
     pages = {32},
     year = {2025},
     doi = {10.5194/agile-giss-6-32-2025},
     url = {https://doi.org/10.5194/agile-giss-6-32-2025}
   }
