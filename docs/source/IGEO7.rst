IGEO7
=====

Overview: A Hierarchically Indexed Hexagonal Equal-Area Discrete Global Grid System
-----------------------------------------------------------------------------------

Hexagonal Discrete Global Grid Systems (DGGS) offer significant advantages for spatial analysis due to their uniform cell shapes and efficient indexing. Among the three central place apertures (3, 4, and 7), aperture 7 subdivisions exhibit very desirable properties, including the preservation of hexagonal symmetry and the formation of unambiguous indexing hierarchies.

Interest in hierarchically indexed aperture 7 hexagonal DGGS has recently increased due to the popularity of the H3 DGGS. But there are currently no open-source equal-area aperture 7 hexagonal DGGS available, that provide similar indexing capabilities like H3.

We present **IGEO7**, a novel pure aperture 7 hexagonal DGGS, and **Z7**, its associated hierarchical integer indexing system. In contrast to H3, where cell sizes vary by up to ±50% across the globe, IGEO7 uses cells of equal area, making it a true equal-area DGGS.

IGEO7 and Z7 are implemented in the open-source software DGGRID. We also present a use case for on-demand suitability modeling to demonstrate a practical application of this new DGGS.


Practical information and sphere vs ellipsoid
---------------------------------------------

The original IGEO7 implementation is available in the DGGRID software since version 8.41 and is a ISEA7H type DGGS with the new Z7 indexing system.

Recent developments in the OGC API DGGS standard draft make reference to refinement ratio 7 DGGRS with ISEA projection.
In order to have DGGH and ZIRS compliant with the OGC DGGS standard,
we aim to enable a few minor adjustments to IGEO7 through the use of dggrid4py:

- to apply authalic conversion to the WGS84 ellipsoid instead of the spherical approximation, using `pygeodesy` (based off geographiclib), see :mod:`dggrid4py.auxlat`
- a rotation of the base icosahedron of 0.05 degrees to align the vertices better with water bodies through a specific parameter in DGGRID: ``dggs_vert0_lon = 11.20`` instead of the DGGRID default of ``11.25``, see the ``meta_config`` in the example below

In practice this means that IGEO7 as described in the original publication is not the same as the ellipsoid-adjusted IGEO7 version, though both can be generated through dggrid4py.

We hope to make the required changes available in DGGRID in the future, so that now further confusions can arise. Until then, there might 
be two slightly different IGEO7 implementations in use and implementers shall be explicit. 

Example of how to generate cells in IGEO7 DGGRS with the Z7 indexing system using dggrid4py
-------------------------------------------------------------------------------------------
In this example, we demonstrate how to use ``dggrid4py`` to generate cells in **IGEO7** DGGRS with the **Z7** indexing system for an input extent in WGS84. The ``DGGRID`` version we use in this example is ``8.43``; any newer version works too, e.g. ``8.44`` from conda-forge, which you get with the pixi quickstart in :ref:`installation`.

First, we instantiate a DGGRIDv8 object from dggrid4py.

.. code:: python

    from dggrid4py import DGGRIDv8
    from dggrid4py.auxlat import geoseries_to_authalic, geoseries_to_geodetic
    import shutil
    import tempfile
    import shapely
    from geopandas import GeoSeries

    # inside a pixi/conda environment, dggrid is on the PATH; otherwise pass the path to the DGGRID executable
    dggrid_instance = DGGRIDv8(shutil.which("dggrid"), working_dir=tempfile.mkdtemp())

Then we create a ``meta_config`` dictionary for use by DGGRIDv8's functions. This ``meta_config`` specifies parameters used by ``DGGRID``, such as which indexing system to use and the position of the initial vertex, etc. Users can override the dggrid4py default settings or introduce additional parameters to DGGRID using this dictionary.

.. code:: python

    meta_config = {
        #input cell ids representation
        "input_address_type": 'HIERNDX', # hierarchical index
        "input_hier_ndx_system": 'Z7', # Z7 hierarchy
        "input_hier_ndx_form": 'DIGIT_STRING', # Z7 textual representation (e.g. '003456231')
        # output cell ids representation
        "output_address_type": 'HIERNDX',
        "output_cell_label_type": 'OUTPUT_ADDRESS_TYPE',
        "output_hier_ndx_system": 'Z7',
        "output_hier_ndx_form": 'DIGIT_STRING',
        # initial vertex longitude, 0.05 degrees off the DGGRID default of 11.25
        "dggs_vert0_lon": 11.20
    }

``input_hier_ndx_form`` and ``output_hier_ndx_form`` accept either ``DIGIT_STRING`` (Z7 textual form, e.g. ``'003456231'``) or ``INT64`` (Z7 integer form).

Then we can use the ``grid_cell_polygons_for_extent`` function to generate **IGEO7** cells using the **Z7** indexing system for an extent. However, as mentioned above, ``DGGRID`` (as of current version 8.43) uses an authalic sphere as the Earth's reference model, so passing geopoints or extents in WGS84 using an ellipsoid as the reference model causes discrepancies. Therefore, we need to convert the input coordinates from WGS84 to authalic for input to the functions, and converting the output back from authalic to WGS84. Users can perform the conversion using  ``geoseries_to_authalic`` and ``geoseries_to_geodetic`` from ``dggrid4py.auxlat``.

.. code:: python

    # Tartu bbox in wgs84, around 50 km^2
    extent = GeoSeries([shapely.box(26.664593, 58.348705, 26.785607, 58.422495)], crs="wgs84")

    # convert the extent in wgs84 to authalic sphere
    extent = geoseries_to_authalic(extent)

    # generate IGEO7 cells with Z7 zone ID for the extent, by passing the meta_config dictionary to specify 
    # the parameters
    igeo7_cells_df = dggrid_instance.grid_cell_polygons_for_extent(dggs_type="IGEO7", resolution=12, clip_geom=extent.geometry[0], **meta_config)

    # convert the cell geometries back to wgs84
    igeo7_cells_df['geometry'] = geoseries_to_geodetic(igeo7_cells_df['geometry'])


Convenience wrappers in ``dggrid4py.igeo7_ext``
-----------------------------------------------

The wrappers in :mod:`dggrid4py.igeo7_ext` apply all of the above for you:

- they only accept a ``DGGRIDv8`` instance (a ``DGGRIDv7`` instance raises ``TypeError``)
- they use the Z7 hierarchical index and ``dggs_vert0_lon = 11.20``, see ``igeo7_meta_config()``
- WGS84 inputs (clip geometries, points) are converted to the authalic sphere before they go to DGGRID,
  and output geometries are converted back to WGS84

.. code:: python

    import shutil
    import tempfile
    import shapely
    import geopandas as gpd
    from dggrid4py import DGGRIDv8, igeo7_ext

    dggrid_instance = DGGRIDv8(shutil.which("dggrid"), working_dir=tempfile.mkdtemp())

    # Tartu bbox in wgs84
    extent = shapely.box(26.664593, 58.348705, 26.785607, 58.422495)

    # cell polygons for an extent, in and out in wgs84
    cells = igeo7_ext.dggrid_igeo7_grid_cell_polygons_for_extent(extent, 9, dggrid_instance)

    # centroids and polygons for a list of Z7 cell ids (all of the same resolution)
    centroids = igeo7_ext.dggrid_igeo7_grid_cell_centroids_from_cellids(cells["name"], dggrid_instance)
    polygons = igeo7_ext.dggrid_igeo7_grid_cell_polygons_from_cellids(cells["name"], dggrid_instance)

    # Z7 cell ids for wgs84 points, returns a copy of the points with a 'name' column
    points = gpd.GeoDataFrame(geometry=[shapely.Point(26.72, 58.38)], crs=4326)
    points = igeo7_ext.dggrid_igeo7_cells_for_geo_points(points, 9, dggrid_instance)

    # Q2DI addresses and direct neighbours
    q2di = igeo7_ext.dggrid_igeo7_q2di_from_cellids(cells["name"], dggrid_instance)
    cls_m = igeo7_ext.dggrid_get_res(dggrid_instance, "IGEO7", 9).loc[9, "cls_m"]
    neighbours = igeo7_ext.z7_k1_ring_neighbours(cells["name"].iloc[0], dggrid_instance, cls_m)

``hier_ndx_form='INT64'`` switches all wrappers to the 64 bit Z7 form (DGGRID writes it as a 16 character hex string,
e.g. ``'0042529bffffffff'``, which :mod:`dggrid4py.igeo7` can decode). The former ``address_type='Z7_STRING'``
argument still works, with a ``DeprecationWarning``, and maps to ``DIGIT_STRING``.

Skipping the authalic conversion matters: in our tests, the WGS84 centroids of resolution 9 cells in Tartu
(around 58.4°N) fell into a different cell for every single point when passed to DGGRID without conversion.


API Reference
-------------

The IGEO7 implementation is available through the dggrid4py Python package.

.. seealso::

   :mod:`dggrid4py.igeo7` - IGEO7 module API documentation

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
