dggrid4py.auxlat
================

Authalic to WGS84 latitude conversions
--------------------------------------

This is experimental and based on `pygeodesy <https://pygeodesy.readthedocs.io/en/latest/>`_.

For example, use like below to convert a clipping polygon from WGS84 to authalic latitudes for DGGRID input.

.. code:: python

   import shapely

   from dggrid4py.auxlat import geodetic_to_authalic, authalic_to_geodetic

   Polygon([ (lon, geodetic_to_authalic(lat)) for (lon, lat) in clip_poly.exterior.coords])


And then use  to convert DGGRID output geometries (centroids or grid_Cell_polygons) from authalic to WGS84 latitudes.

.. code:: python

   import geopandas as gpd
   from dggrid4py.auxlat import geoseries_to_authalic, geoseries_to_geodetic
   
   geoseries_to_geodetic(gdf.geometry)


.. automodule:: dggrid4py.auxlat

   .. rubric:: Functions

   .. autosummary::
   
      authalic_to_geodetic
      geodetic_to_authalic
      geoseries_to_authalic
      geoseries_to_geodetic