import numpy as np
import shapely
import shapely.ops
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


def transform_xy(func, geom):
    """
    Apply ``func(x, y)``, which takes and returns arrays, to the coordinates of a shapely geometry. A z coordinate stays.
    """
    if not hasattr(shapely, 'transform'):
        # shapely < 2
        def _with_z(x, y, z=None):
            new_x, new_y = func(x, y)
            return (new_x, new_y) if z is None else (new_x, new_y, z)
        return shapely.ops.transform(_with_z, geom)

    def _coords(coords):
        new_x, new_y = func(coords[:, 0], coords[:, 1])
        return np.column_stack([new_x, new_y, *coords.T[2:]])
    return shapely.transform(geom, _coords, include_z=geom.has_z)


def _apply_to_geometry(geom, func):
    # any shapely geometry (points, lines, polygons with holes, multi-geometries, collections),
    # coordinates in lon/lat, only the latitude is converted
    if geom is None or geom.is_empty:
        return geom
    vfunc = np.vectorize(lambda lat: float(func(lat)), otypes=[float])

    def _lat(x, y):
        return x, vfunc(y) if np.ndim(y) else float(func(y))

    return transform_xy(_lat, geom)


def geoseries_to_authalic(geoseries):
    """
    Convert the latitudes of a GeoSeries in WGS84 (geodetic) to the authalic sphere used by DGGRID.
    """
    return geoseries.apply(lambda geom: _apply_to_geometry(geom, geodetic_to_authalic))


def geoseries_to_geodetic(geoseries):
    """
    Convert the latitudes of a GeoSeries on the authalic sphere (DGGRID output) back to WGS84 (geodetic).
    """
    return geoseries.apply(lambda geom: _apply_to_geometry(geom, authalic_to_geodetic))
