from pygeodesy.ellipsoids import Ellipsoids
from shapely.geometry import Point, Polygon

wgs84 = Ellipsoids.WGS84

def authalic_to_geodetic(lat_authalic):
    return wgs84.auxAuthalic(lat_authalic, inverse=True)

def geodetic_to_authalic(lat_geodetic):
    return wgs84.auxAuthalic(lat_geodetic, inverse=False)

def _apply_to_simple_shapely_polygon(polygon, func):
    # cannot be a complex polygon, e.g. with holes, just a simple one like the hexagons or bbox
    # multigeometries should be treated separately
    # it is assumed that coordinates are in lat/lon (either on authalic sphere or ellipsoid/wgs84)
    return Polygon([ (lon, func(lat)) for (lon, lat) in polygon.exterior.coords])

def _apply_to_shapely_point(point, func):
    # points is a list of shapely Point objects
    return Point(point.x, func(point.y))

def _apply_to_shapely_points(points, func):
    # points is a list of shapely Point objects
    return [ Point(point.x, func(point.y)) for point in points]


def geoseries_to_authalic(geoseries):
    # geoseries is a geopandas GeoSeries of shapely geometries
    return geoseries.apply(lambda geom: geom if geom.is_empty else
                          (_apply_to_shapely_point(geom, geodetic_to_authalic) if isinstance(geom, Point) else
                           _apply_to_simple_shapely_polygon(geom, geodetic_to_authalic)))

def geoseries_to_geodetic(geoseries):
    # geoseries is a geopandas GeoSeries of shapely geometries
    return geoseries.apply(lambda geom: geom if geom.is_empty else
                          (_apply_to_shapely_point(geom, authalic_to_geodetic) if isinstance(geom, Point) else
                           _apply_to_simple_shapely_polygon(geom, authalic_to_geodetic)))
