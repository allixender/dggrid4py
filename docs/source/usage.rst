Usage
=====

.. _installation:

Installation
------------

Quickstart with pixi (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The easiest way to get started is a `pixi <https://pixi.sh>`_ environment. It installs
DGGRID (with GDAL support) and the geospatial Python stack from conda-forge, and
dggrid4py from PyPI, so you don't need to compile anything.

.. code-block:: console

   $ pixi init my-dggs-project -c conda-forge
   $ cd my-dggs-project
   $ pixi add python=3.12 dggrid geopandas pyogrio
   $ pixi add --pypi dggrid4py

``dggrid4py`` pulls in ``pygeodesy`` from PyPI (it is not available on conda-forge).
Check that the ``dggrid`` executable is available in the environment:

.. code-block:: console

   $ pixi run which dggrid
   $ pixi run dggrid        # prints "usage: dggrid metaFileName"

Inside the pixi environment, DGGRID lives on the ``PATH``, so you can point dggrid4py at it with ``shutil.which``:

.. code:: python

   import shutil
   import tempfile
   from dggrid4py import DGGRIDv8

   dggrid_instance = DGGRIDv8(executable=shutil.which("dggrid"), working_dir=tempfile.mkdtemp(), capture_logs=False, silent=True)

Run your scripts with ``pixi run python my_script.py`` or open a shell in the environment with ``pixi shell``.

For IGEO7 with the Z7 index, continue with :ref:`igeo7_usage`.

Developing dggrid4py with pixi
""""""""""""""""""""""""""""""

To run the dggrid4py test suite against the conda-forge DGGRID, create a pixi environment
outside the repository and point ``DGGRID_PATH`` at its ``dggrid``:

.. code-block:: console

   $ pixi init dggrid4py-dev -c conda-forge && cd dggrid4py-dev
   $ pixi add python=3.12 dggrid pytest geopandas shapely pandas numpy pyogrio pip
   $ pixi run python -m pip install pygeodesy
   $ cd <path_to>/dggrid4py
   $ DGGRID_PATH=<path_to>/dggrid4py-dev/.pixi/envs/default/bin/dggrid \
       PYTHONPATH=. <path_to>/dggrid4py-dev/.pixi/envs/default/bin/python -m pytest tests

Installation with pip
^^^^^^^^^^^^^^^^^^^^^

To use dggrid4py, first install it using pip:

.. code-block:: console

   (.venv) $ pip install dggrid4py


dggrid4py runs the ``dggrid`` command line tool, which has to be available on the system.

You can install DGGRID from conda-forge:

.. code-block:: console

   (.venv) $ conda install -c conda-forge dggrid

Or compile from source: https://github.com/sahrk/DGGRID


Portable DGGRID binary
----------------------

.. warning::

   The precompiled portable DGGRID binaries are still **experimental**. They are built
   without GDAL support, and the current download is a pre-release of DGGRID 9. For a
   reliable setup, use the pixi or conda-forge installation described above.

If DGGRID can neither be installed from conda-forge nor compiled, dggrid4py can download a portable binary
from the `DGGRID_portables <https://github.com/allixender/DGGRID_portables>`_ releases. These binaries are
available for Linux, macOS and Windows (x86_64 and arm64 each) and have no further dependencies.

.. code:: python

   import tempfile

   from dggrid4py import DGGRIDv8, tool

   # downloads and unpacks the binary for the current platform into the given folder (only once)
   dggrid_exec = tool.get_portable_executable("dggrid_portable")
   dggrid_instance_portable = DGGRIDv8(executable=dggrid_exec, working_dir=tempfile.mkdtemp(), capture_logs=False, silent=True, has_gdal=False)

``get_portable_executable`` verifies the downloaded archive against the SHA256 checksums of the release, unpacks it
and returns the path of the executable. A binary that is already in the folder is used again, and it is also
returned if there is no network connection. The following points should be considered:

- the binaries are built without GDAL, thus the instance has to be created with ``has_gdal=False``; dggrid4py
  then exchanges Shapefiles with DGGRID, and the returned GeoDataFrames have the same columns and cell identifiers
  as with a GDAL build;
- the default release ``edge`` is a rolling pre-release that follows the DGGRID development version (currently
  DGGRID 9.0b), and it is downloaded again after it was rebuilt; a tagged release can be selected with
  ``release="<tag>"`` once one is published;
- DGGRID 9 removes the address types ``Z7``, ``Z7_STRING`` and their Z3 and ZORDER counterparts, as well as the
  ``DGGRIDv7`` way of clipping by sequence numbers, so the portable binary has to be used with the ``DGGRIDv8`` class;
- on Windows arm64, coordinates can differ by up to about 1e-5 degrees from the other platforms (see the release notes).


.. _basic_usage:

Basic Usage
-----------

Besides some low-level access to influence the metafile creation of the DGGRID operations, a few high-level
functions are integrated to work with the more comfortable geopython libraries, like shapely and geopandas:

-  grid_cell_polygons_for_extent(): fill extent/subset with cells at
   resolution (clip or world);
-  grid_cell_centroids_for_extent(): the same for the cell centroids;
-  grid_cell_polygons_from_cellids(): geometry_from_cellid for dggs at
   resolution (from id list);
-  grid_cellids_for_extent(): get_all_indexes/cell_ids for dggs at
   resolution (clip or world);
-  cells_for_geo_points(): poly_outline for point/centre at resolution;
-  address_transform(): conversion between cell_id address types, like SEQNUM, Q2DI or the hierarchical indexes;
-  grid_stats_table(): number of cells, cell area and spacing per resolution.

The functions are methods of a DGGRID instance, which knows where the ``dggrid`` executable lives. Use the
``DGGRIDv8`` class with DGGRID 8.42 or newer (it also runs the 9.0b pre-release). The following example uses
the classical DGGS types of DGGRID with their default orientation.

.. deprecated:: 0.6.0
   The ``DGGRIDv7`` class for DGGRID 7 is deprecated and issues a ``DeprecationWarning`` when an instance is
   created. It will be removed in a future version. Existing code can replace ``DGGRIDv7`` with ``DGGRIDv8``;
   the former hierarchical address types (e.g. ``Z7_STRING``) are still accepted there, see :doc:`IGEO7`.

.. code:: python

   import shutil
   import tempfile

   import geopandas
   import shapely

   from dggrid4py import DGGRIDv8

   # create an initial instance that knows where the dggrid tool lives, configure temp workspace and log/stdout output
   dggrid_instance = DGGRIDv8(executable=shutil.which("dggrid"), working_dir=tempfile.mkdtemp(), capture_logs=False, silent=True)

   # global ISEA4T grid at resolution 5 into GeoDataFrame to Shapefile
   gdf1 = dggrid_instance.grid_cell_polygons_for_extent('ISEA4T', 5)
   print(gdf1.head())
   gdf1.to_file('isea4t_5.shp')

   # cell centroids of the global ISEA7H grid at resolution 4
   gdf_centroids = dggrid_instance.grid_cell_centroids_for_extent(dggs_type='ISEA7H', resolution=4, mixed_aperture_level=None, clip_geom=None)
   print(gdf_centroids.head())

   # clip extent
   clip_bound = shapely.geometry.box(20.2, 57.00, 28.4, 60.0)

   # ISEA7H grid at resolution 9, for extent of provided WGS84 rectangle into GeoDataFrame to Shapefile
   gdf3 = dggrid_instance.grid_cell_polygons_for_extent('ISEA7H', 9, clip_geom=clip_bound)
   print(gdf3.head())
   gdf3.to_file('est_shape_isea7h_9.shp')

   # generate cell and areal statistics for a ISEA7H grids from resolution 0 to 8 (return a pandas DataFrame)
   df1 = dggrid_instance.grid_stats_table('ISEA7H', 8)
   print(df1.head(8))
   df1.to_csv('isea7h_8_stats.csv', index=False)

   # generate the DGGS grid cells that would cover a GeoDataFrame of points, return Polygons with cell IDs as GeoDataFrame
   points = [shapely.Point(20.5, 57.5), shapely.Point(23.5, 58.5)]
   geodf_points_wgs84 = geopandas.GeoDataFrame({'site': ['A', 'B']}, geometry=points, crs='EPSG:4326')
   gdf4 = dggrid_instance.cells_for_geo_points(geodf_points_wgs84, False, 'ISEA7H', 5)
   print(gdf4.head())
   gdf4.to_file('polycells_from_points_isea7h_5.shp')

   # generate the DGGS grid cells that would cover a GeoDataFrame of points, return a copy of the points with the cell IDs in the column 'name'
   gdf5 = dggrid_instance.cells_for_geo_points(geodf_points_wgs84=geodf_points_wgs84, cell_ids_only=True, dggs_type='ISEA4H', resolution=8)
   print(gdf5.head())
   gdf5.to_file('geopoint_cellids_from_points_isea4h_8.shp')

   # generate DGGS grid cell polygons based on 'cell_id_list' (a list or np.array of provided cell_ids)
   gdf6 = dggrid_instance.grid_cell_polygons_from_cellids(cell_id_list=[1, 4, 8], dggs_type='ISEA7H', resolution=5)
   print(gdf6.head())
   gdf6.to_file('from_seqnums_isea7h_5.shp')

   # split cells at the dateline for cartesian GIS tools
   gdf7 = dggrid_instance.grid_cell_polygons_for_extent('ISEA7H', 3, split_dateline=True)
   gdf7.to_file('global_isea7h_3_interrupted.shp')

   # convert cell IDs between address types, here from sequence numbers to Q2DI (quad number and i, j coordinates)
   df_q2di = dggrid_instance.address_transform([1, 4, 8], 'ISEA7H', 5, input_address_type='SEQNUM', output_address_type='Q2DI')
   print(df_q2di.head(3))

   # cells at resolution 6 that intersect a coarser cell at resolution 5 (a spatial clip, it includes overlapping neighbours)
   gdf8 = dggrid_instance.grid_cell_polygons_from_cellids([100], 'ISEA7H', 6, clip_subset_type='COARSE_CELLS', clip_cell_res=5)
   print(len(gdf8))

All functions that return geometries return a GeoDataFrame with the cell identifier in the column ``name`` and
with ``EPSG:4326`` as coordinate reference system. For the classical DGGS types, DGGRID takes and returns longitude
and latitude in degrees on its sphere, and dggrid4py passes the coordinates through without a datum conversion.
``cells_for_geo_points`` does not modify the GeoDataFrame it is given. It returns a copy with the columns ``lon``,
``lat`` and ``name`` (the cell identifier), or with ``cell_ids_only=False`` one cell polygon per point, with the
cell identifier in the column ``zone`` and the columns of the points attached.

Further DGGRID parameters can be passed to all functions as keyword arguments, e.g. ``dggs_vert0_lon``,
``dggs_vert0_lat`` and ``dggs_vert0_azimuth`` for the orientation of the icosahedron. Only the parameters that
are given are changed, the others keep the DGGRID default. An input or output address type that the DGGRID
class does not know raises a ``ValueError``.


.. _igeo7_usage:

IGEO7 Usage
-----------

IGEO7 is the ISEA7H grid with the Z7 hierarchical index. In contrast to the classical DGGS types above, the
IGEO7 used in dggrid4py is defined with a different orientation (``dggs_vert0_lon = 11.20``) and on the WGS84
ellipsoid, which requires a conversion to authalic latitudes around each DGGRID call. The wrappers in
:mod:`dggrid4py.igeo7_ext` apply both for you:

.. code:: python

   from dggrid4py import igeo7_ext

   # Tartu bbox in wgs84, IGEO7 cell polygons at resolution 9 with Z7 cell ids in the column 'name'
   extent = shapely.box(26.664593, 58.348705, 26.785607, 58.422495)
   cells = igeo7_ext.dggrid_igeo7_grid_cell_polygons_for_extent(extent, 9, dggrid_instance)

The background on orientation and ellipsoid, the full examples, the explicit configuration with ``DGGRIDv8``
and the migration from the former ``Z7_STRING`` address type are described in :doc:`IGEO7`.

TODO
----

Contributions are welcome.

-  get parent_for_cell_id at coarser resolution

-  get children_for_cell_id at finer resolution

With the IGEO7/Z7 index system, the parent of a cell can be derived from the cell identifier. A function that
enumerates the index children is still missing, see :doc:`IGEO7`.
