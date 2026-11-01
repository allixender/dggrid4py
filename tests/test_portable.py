#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
The portable DGGRID binary (no GDAL, Shapefile in and out) has to return the same cells
as the local DGGRID build (usually with GDAL).
"""
import os
import shutil
import tempfile

import geopandas as gpd
import pytest
import shapely

from dggrid4py import DGGRIDv8, DGGRIDv9, igeo7_ext, tool

clip_bound = shapely.geometry.box(27.2, 57.5, 29.3, 59.2)

IGEO7 = {
    "input_address_type": "HIERNDX",
    "input_hier_ndx_system": "Z7",
    "input_hier_ndx_form": "DIGIT_STRING",
    "output_address_type": "HIERNDX",
    "output_cell_label_type": "OUTPUT_ADDRESS_TYPE",
    "output_hier_ndx_system": "Z7",
    "output_hier_ndx_form": "DIGIT_STRING",
    "dggs_vert0_lon": 11.20,
}


# the lines of portable binaries that have to work with this dggrid4py version
@pytest.fixture(scope="module", params=["stable", "edge"])
def portable_dggrid(request, tmp_path_factory):
    try:
        executable = tool.get_portable_executable(tmp_path_factory.mktemp("portable"), line=request.param)
    except (OSError, ValueError) as e:
        pytest.skip(f"portable DGGRID '{request.param}' not available: {e}")
    return DGGRIDv8(executable=executable, working_dir=tempfile.mkdtemp(), capture_logs=True, silent=True, has_gdal=False)


@pytest.fixture(scope="module")
def local_dggrid():
    path = os.getenv("DGGRID_PATH") or shutil.which("dggrid")
    if not path or not os.path.isfile(path):
        pytest.skip("DGGRID executable not available")
    return DGGRIDv8(executable=path, working_dir=tempfile.mkdtemp(), capture_logs=True, silent=True)


@pytest.fixture(scope="module")
def cellids100(local_dggrid):
    cells = local_dggrid.grid_cell_polygons_for_extent("IGEO7", 3, clip_geom=clip_bound, **IGEO7)
    return cells["name"][:100].tolist()


def _names(gdf):
    return sorted(gdf["name"])


def test_grid_cell_polygons_for_extent(portable_dggrid, local_dggrid):
    portable = portable_dggrid.grid_cell_polygons_for_extent("IGEO7", 3, clip_geom=clip_bound, **IGEO7)
    local = local_dggrid.grid_cell_polygons_for_extent("IGEO7", 3, clip_geom=clip_bound, **IGEO7)
    assert len(local) > 0
    assert _names(portable) == _names(local)
    assert portable.crs == local.crs == "EPSG:4326"


def test_grid_cell_centroids_for_extent(portable_dggrid, local_dggrid):
    portable = portable_dggrid.grid_cell_centroids_for_extent("IGEO7", 3, clip_geom=clip_bound, **IGEO7)
    local = local_dggrid.grid_cell_centroids_for_extent("IGEO7", 3, clip_geom=clip_bound, **IGEO7)
    assert _names(portable) == _names(local)


def test_grid_cell_polygons_from_cellids(portable_dggrid, local_dggrid, cellids100):
    portable = portable_dggrid.grid_cell_polygons_from_cellids(cellids100, "IGEO7", 3, **IGEO7)
    local = local_dggrid.grid_cell_polygons_from_cellids(cellids100, "IGEO7", 3, **IGEO7)
    assert _names(portable) == _names(local) == sorted(cellids100)


def test_grid_cell_centroids_from_cellids(portable_dggrid, local_dggrid, cellids100):
    portable = portable_dggrid.grid_cell_centroids_from_cellids(cellids100, "IGEO7", 3, **IGEO7)
    local = local_dggrid.grid_cell_centroids_from_cellids(cellids100, "IGEO7", 3, **IGEO7)
    assert _names(portable) == _names(local) == sorted(cellids100)


def _same_dggrid_version(portable_dggrid, local_dggrid):
    # COARSE_CELLS returns slightly different neighbouring cells in DGGRID 8.44 and 9.1, with GDAL and without
    if portable_dggrid.dggrid_version != local_dggrid.dggrid_version:
        pytest.skip(f"portable DGGRID {portable_dggrid.dggrid_version}, local DGGRID {local_dggrid.dggrid_version}")


def test_igeo7_reference_cells(portable_dggrid):
    points = gpd.GeoDataFrame(geometry=gpd.points_from_xy([-9.1393, 26.7220], [38.7223, 58.3776]), crs=4326)
    cells = igeo7_ext.dggrid_igeo7_cells_for_geo_points(points, 9, portable_dggrid)
    assert cells["name"].tolist() == ["00641565463", "00010224545"]


def test_grid_cell_polygons_from_cellids_coarse_cells(portable_dggrid, local_dggrid, cellids100):
    _same_dggrid_version(portable_dggrid, local_dggrid)
    portable = portable_dggrid.grid_cell_polygons_from_cellids(cellids100, "IGEO7", 5, clip_subset_type="COARSE_CELLS", clip_cell_res=3, **IGEO7)
    local = local_dggrid.grid_cell_polygons_from_cellids(cellids100, "IGEO7", 5, clip_subset_type="COARSE_CELLS", clip_cell_res=3, **IGEO7)
    assert _names(portable) == _names(local)


def test_grid_cell_centroids_from_cellids_coarse_cells(portable_dggrid, local_dggrid, cellids100):
    _same_dggrid_version(portable_dggrid, local_dggrid)
    portable = portable_dggrid.grid_cell_centroids_from_cellids(cellids100, "IGEO7", 5, clip_subset_type="COARSE_CELLS", clip_cell_res=3, **IGEO7)
    local = local_dggrid.grid_cell_centroids_from_cellids(cellids100, "IGEO7", 5, clip_subset_type="COARSE_CELLS", clip_cell_res=3, **IGEO7)
    assert _names(portable) == _names(local)


def test_grid_cellids_for_extent(portable_dggrid, local_dggrid):
    portable = portable_dggrid.grid_cellids_for_extent("IGEO7", 3, clip_geom=clip_bound, **IGEO7)
    local = local_dggrid.grid_cellids_for_extent("IGEO7", 3, clip_geom=clip_bound, **IGEO7)
    assert sorted(portable[0]) == sorted(local[0])


def test_seqnum_cells(portable_dggrid, local_dggrid):
    portable = portable_dggrid.grid_cell_polygons_from_cellids([1, 4, 8], "ISEA7H", 5)
    local = local_dggrid.grid_cell_polygons_from_cellids([1, 4, 8], "ISEA7H", 5)
    assert _names(portable) == _names(local) == [1, 4, 8]


def test_dggrid_v9_on_portable_reference_cells(portable_dggrid):
    if (portable_dggrid.dggrid_version or (0, 0)) < (9, 1):
        pytest.skip(f"needs DGGRID 9.1 or newer, portable DGGRID {portable_dggrid.dggrid_version}")
    dggrid = DGGRIDv9(executable=portable_dggrid.executable, working_dir=tempfile.mkdtemp(), capture_logs=True, silent=True, has_gdal=False)
    sites = gpd.GeoDataFrame(geometry=gpd.points_from_xy([-9.1393, 26.7220], [38.7223, 58.3776]), crs=4326)
    cells = dggrid.cells_for_geo_points(sites, True, "IGEO7", 9, output_hier_ndx_form="DIGIT_STRING")
    assert cells["name"].tolist() == ["00641565463", "00010224545"]
    box = shapely.box(26.6, 58.3, 26.85, 58.45)
    polygons = dggrid.grid_cell_polygons_for_extent("IGEO7", 6, clip_geom=box, output_hier_ndx_form="DIGIT_STRING")
    assert polygons.crs.to_epsg() == 4326
    assert "00010224" in set(polygons["name"])
