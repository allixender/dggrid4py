#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The wrappers have to return the same IGEO7 for DGGRIDv8 and DGGRIDv9 on DGGRID 9.1.
"""
import os
import shutil
import tempfile

import geopandas as gpd
import numpy as np
import pytest
import shapely

from dggrid4py import DGGRIDv8, DGGRIDv9, igeo7_ext
from dggrid4py.dggrid_runner import _read_dggrid_version

TARTU_BOX = shapely.box(26.6, 58.3, 26.85, 58.45)
SITES = {"Lisbon": (-9.1393, 38.7223), "Tartu": (26.7220, 58.3776)}
REFERENCE = {5: ["0064156", "0001022"], 6: ["00641565", "00010224"], 9: ["00641565463", "00010224545"]}


def _instances():
    path = os.getenv("DGGRID_PATH") or shutil.which("dggrid")
    if not path or not os.path.isfile(path):
        pytest.skip("DGGRID executable not available")
    version = _read_dggrid_version(path)
    if version is None or version < (9, 1):
        pytest.skip(f"needs DGGRID 9.1 or newer, local DGGRID {version}")
    kwargs = dict(executable=path, capture_logs=True, silent=True)
    return DGGRIDv8(working_dir=tempfile.mkdtemp(), **kwargs), DGGRIDv9(working_dir=tempfile.mkdtemp(), **kwargs)


def _sites():
    return gpd.GeoDataFrame(
        {"site": list(SITES)},
        geometry=gpd.points_from_xy(*zip(*SITES.values())),
        crs=4326,
    )


@pytest.mark.parametrize("resolution", [5, 6, 9])
def test_reference_cells_with_both_classes(resolution):
    for dggrid in _instances():
        cells = igeo7_ext.dggrid_igeo7_cells_for_geo_points(_sites(), resolution, dggrid)
        assert cells["name"].tolist() == REFERENCE[resolution]


@pytest.mark.parametrize("resolution", [9, 12])
def test_sample_of_points_same_ids(resolution):
    v8, v9 = _instances()
    rng = np.random.default_rng(42)
    lon = rng.uniform(-180, 180, 200)
    lat = np.degrees(np.arcsin(rng.uniform(-1, 1, 200)))
    points = gpd.GeoDataFrame(geometry=gpd.points_from_xy(lon, lat), crs=4326)
    a = igeo7_ext.dggrid_igeo7_cells_for_geo_points(points, resolution, v8)["name"].tolist()
    b = igeo7_ext.dggrid_igeo7_cells_for_geo_points(points, resolution, v9)["name"].tolist()
    assert a == b


@pytest.mark.parametrize("hier_ndx_form", ["DIGIT_STRING", "INT64"])
@pytest.mark.parametrize("to_geodetic", [True, False])
def test_extent_and_from_cellids_same_cells_and_geometry(hier_ndx_form, to_geodetic):
    v8, v9 = _instances()
    a = igeo7_ext.dggrid_igeo7_grid_cell_polygons_for_extent(TARTU_BOX, 8, v8, hier_ndx_form, to_geodetic)
    b = igeo7_ext.dggrid_igeo7_grid_cell_polygons_for_extent(TARTU_BOX, 8, v9, hier_ndx_form, to_geodetic)
    assert len(a) > 0
    a, b = a.sort_values("name").reset_index(drop=True), b.sort_values("name").reset_index(drop=True)
    assert a["name"].tolist() == b["name"].tolist()
    # even resolution, there is no difference in the cell boundaries between the two paths
    assert a.geometry.hausdorff_distance(b.geometry).max() < 1e-5
    assert a.crs.to_epsg() == b.crs.to_epsg() == 4326

    ids = a["name"].tolist()[:5]
    c = igeo7_ext.dggrid_igeo7_grid_cell_centroids_from_cellids(ids, v8, hier_ndx_form, to_geodetic).sort_values("name").reset_index(drop=True)
    d = igeo7_ext.dggrid_igeo7_grid_cell_centroids_from_cellids(ids, v9, hier_ndx_form, to_geodetic).sort_values("name").reset_index(drop=True)
    assert c["name"].tolist() == d["name"].tolist()
    assert c.geometry.distance(d.geometry).max() < 1e-5


def test_to_geodetic_false_is_on_the_sphere():
    _, v9 = _instances()
    wgs84 = igeo7_ext.dggrid_igeo7_grid_cell_centroids_from_cellids(["00010224545"], v9).geometry.iloc[0]
    sphere = igeo7_ext.dggrid_igeo7_grid_cell_centroids_from_cellids(["00010224545"], v9, to_geodetic=False).geometry.iloc[0]
    # around 0.1 degrees between geodetic and authalic latitude at 58 degrees north
    assert 0.05 < wgs84.y - sphere.y < 0.2
    assert wgs84.x == pytest.approx(sphere.x, abs=1e-6)


def test_q2di_same():
    v8, v9 = _instances()
    ids = ["00641565463", "00010224545"]
    a = igeo7_ext.dggrid_igeo7_q2di_from_cellids(ids, v8)
    b = igeo7_ext.dggrid_igeo7_q2di_from_cellids(ids, v9)
    assert a[["Q", "I", "J"]].values.tolist() == b[["Q", "I", "J"]].values.tolist()


@pytest.mark.parametrize("cell_id", ["00010224545", "0001022454"])
def test_neighbours_same(cell_id):
    v8, v9 = _instances()
    resolution = len(cell_id) - 2
    cls_m = igeo7_ext.dggrid_get_res(v8, "IGEO7", resolution).loc[resolution, "cls_m"]
    a = igeo7_ext.z7_k1_ring_neighbours(cell_id, v8, cls_m)
    b = igeo7_ext.z7_k1_ring_neighbours(cell_id, v9, cls_m)
    assert len(a) == 6
    assert sorted(a) == sorted(b)


def test_meta_config_per_class():
    v8 = DGGRIDv8(executable="dggrid", working_dir=tempfile.mkdtemp())
    v9 = DGGRIDv9(executable="dggrid", working_dir=tempfile.mkdtemp())
    assert igeo7_ext.igeo7_meta_config(dggrid_instance=v8) == igeo7_ext.igeo7_meta_config()
    assert igeo7_ext.igeo7_meta_config(dggrid_instance=v8)["dggs_vert0_lon"] == 11.20
    assert igeo7_ext.igeo7_meta_config("INT64", v9) == {"input_hier_ndx_form": "INT64", "output_hier_ndx_form": "INT64"}
    assert igeo7_ext.igeo7_meta_config("INT64", v9, polyclip_algo="CLIPPER1")["polyclip_algo"] == "CLIPPER1"


def test_neighbours_v9_without_cls_m_v8_needs_it():
    v8, v9 = _instances()
    assert len(igeo7_ext.z7_k1_ring_neighbours("00010224545", v9)) == 6
    with pytest.raises(ValueError, match="needs cls_m"):
        igeo7_ext.z7_k1_ring_neighbours("00010224545", v8)
