#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil
import tempfile

import pytest
import shapely
import geopandas as gpd

from dggrid4py import DGGRIDv7, DGGRIDv8, igeo7_ext
from dggrid4py.auxlat import geoseries_to_authalic, geoseries_to_geodetic

# Tartu bbox in WGS84
TARTU = shapely.box(26.664593, 58.348705, 26.785607, 58.422495)


def _dggrid():
    path = os.getenv("DGGRID_PATH") or shutil.which("dggrid")
    if not path or not os.path.isfile(path):
        pytest.skip("DGGRID executable not available")
    return DGGRIDv8(executable=path, working_dir=tempfile.mkdtemp(), capture_logs=True, silent=True)


def test_igeo7_meta_config_and_v8_only():
    meta = igeo7_ext.igeo7_meta_config()
    assert meta["dggs_vert0_lon"] == 11.20
    assert meta["input_address_type"] == meta["output_address_type"] == "HIERNDX"
    assert meta["input_hier_ndx_system"] == meta["output_hier_ndx_system"] == "Z7"
    assert meta["output_hier_ndx_form"] == "DIGIT_STRING"
    assert igeo7_ext.igeo7_meta_config("INT64", dggs_vert0_lon=11.25)["dggs_vert0_lon"] == 11.25

    # DGGRIDv7 Z7_STRING form is mapped with a deprecation warning, unknown forms are rejected
    with pytest.warns(DeprecationWarning):
        assert igeo7_ext.igeo7_meta_config("Z7_STRING")["input_hier_ndx_form"] == "DIGIT_STRING"
    with pytest.raises(ValueError):
        igeo7_ext.igeo7_meta_config("Z3")

    # the wrappers refuse DGGRIDv7 before DGGRID runs
    with pytest.raises(TypeError):
        igeo7_ext.dggrid_igeo7_grid_cell_polygons_from_cellids(["000102022"], DGGRIDv7(executable="dggrid"))


def test_authalic_roundtrip_multigeometries():
    geoms = gpd.GeoSeries([
        TARTU,
        shapely.Point(26.7, 58.4),
        shapely.MultiPolygon([shapely.box(0, 0, 1, 1), shapely.box(2, 45, 3, 46)]),
        shapely.box(0, 0, 10, 10).difference(shapely.box(2, 4, 3, 5)),  # with hole
        shapely.Polygon(),
    ])
    authalic = geoseries_to_authalic(geoms)
    # authalic latitude is smaller than geodetic in the northern hemisphere, longitude unchanged
    assert authalic.iloc[1].x == 26.7
    assert authalic.iloc[1].y < 58.4
    assert authalic.iloc[2].geom_type == "MultiPolygon"
    assert len(authalic.iloc[3].interiors) == 1
    assert authalic.iloc[4].is_empty
    back = geoseries_to_geodetic(authalic)
    for orig, rt in zip(geoms.iloc[:4], back.iloc[:4]):
        assert orig.hausdorff_distance(rt) < 1e-9


def test_igeo7_extent_points_roundtrip():
    dggrid = _dggrid()
    cells = igeo7_ext.dggrid_igeo7_grid_cell_polygons_for_extent(TARTU, 9, dggrid)
    assert len(cells) > 0
    assert cells.crs.to_epsg() == 4326
    ids = cells["name"].tolist()
    assert all(igeo7_ext.z7_resolution(i) == 9 for i in ids)

    # WGS84 centroids fall inside their WGS84 polygons
    centroids = igeo7_ext.dggrid_igeo7_grid_cell_centroids_from_cellids(ids, dggrid).set_index("name")
    polygons = igeo7_ext.dggrid_igeo7_grid_cell_polygons_from_cellids(ids, dggrid).set_index("name")
    assert all(polygons.loc[i].geometry.contains(centroids.loc[i].geometry) for i in ids)

    # WGS84 centroids map back to the same cell ids
    points = gpd.GeoDataFrame({"src": centroids.index}, geometry=centroids.geometry.values, crs=4326)
    result = igeo7_ext.dggrid_igeo7_cells_for_geo_points(points, 9, dggrid)
    assert (result["src"] == result["name"]).all()
    assert result.geometry.equals(points.geometry)

    # without the authalic conversion the WGS84 points land in other cells
    naive = dggrid.cells_for_geo_points(points[["geometry"]].copy(), True, "IGEO7", 9, **igeo7_ext.igeo7_meta_config())
    assert (naive["name"].astype(str).values != result["src"].values).any()


def test_igeo7_k1_ring_neighbours():
    dggrid = _dggrid()
    cells = igeo7_ext.dggrid_igeo7_grid_cell_polygons_for_extent(TARTU, 9, dggrid)
    cell = next(c for c in cells["name"] if not c.endswith("0"))
    cls_m = igeo7_ext.dggrid_get_res(dggrid, "IGEO7", 9).loc[9, "cls_m"]

    neighbours = igeo7_ext.z7_k1_ring_neighbours(cell, dggrid, cls_m)
    assert len(neighbours) == 6
    assert cell not in neighbours
    assert all(igeo7_ext.z7_resolution(n) == 9 for n in neighbours)

    # the neighbour polygons touch the cell polygon
    polys = igeo7_ext.dggrid_igeo7_grid_cell_polygons_from_cellids([cell, *neighbours], dggrid).set_index("name")
    centre = polys.loc[cell].geometry
    assert all(centre.buffer(1e-7).intersects(polys.loc[n].geometry) for n in neighbours)
