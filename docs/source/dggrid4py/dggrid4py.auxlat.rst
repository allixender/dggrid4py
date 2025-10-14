dggrid4py.auxlat
================

Authalic to WGS84 latitude conversions
--------------------------------------

This is experimental and based on `pygeodesy <https://pygeodesy.readthedocs.io/en/latest/>`_.

E.g. use `Polygon([ (lon, func(lat)) for (lon, lat) in polygon.exterior.coords])`to convert a clipping polygon from WGS84 to authalic latitudes for DGGRID input.


E.g. use `geoseries_to_geodetic(gdf.geometry)` to convert DGGRID output geometries (centroids or grid_Cell_polygons) from authalic to WGS84 latitudes.


.. automodule:: dggrid4py.auxlat

   .. rubric:: Functions

   .. autosummary::
   
      authalic_to_geodetic
      geodetic_to_authalic
      geoseries_to_authalic
      geoseries_to_geodetic