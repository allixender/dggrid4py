#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pytest

import os
import shutil
import tempfile

import geopandas as gpd
import numpy as np
import shapely

from dggrid4py import DGGRIDv8
from dggrid4py.dggrid_runner import DGGRIDv9, Dggs, dg_grid_meta, dgselect, dggs_types_v9, _read_dggrid_version

LISBON = (-9.1393, 38.7223)
TARTU = (26.7220, 58.3776)
# IGEO7 reference cells, Z7 digit strings at resolution 5, 6 and 9
REFERENCE = {
    LISBON: {5: "0064156", 6: "00641565", 9: "00641565463"},
    TARTU: {5: "0001022", 6: "00010224", 9: "00010224545"},
}
TARTU_BOX = shapely.box(26.6, 58.3, 26.85, 58.45)
DIGITS = {
    "output_hier_ndx_form": "DIGIT_STRING",
    "input_hier_ndx_form": "DIGIT_STRING",
}


def _dggrid_path():
    path = os.getenv("DGGRID_PATH") or shutil.which("dggrid")
    if not path or not os.path.isfile(path):
        pytest.skip("DGGRID executable not available")
    return path


def _dggrid91():
    path = _dggrid_path()
    version = _read_dggrid_version(path)
    if version is None or version < (9, 1):
        pytest.skip(f"needs DGGRID 9.1 or newer, local DGGRID {version}")
    return DGGRIDv9(executable=path, working_dir=tempfile.mkdtemp(), capture_logs=True, silent=True)


def _points():
    return gpd.GeoDataFrame(
        {"city": ["Lisbon", "Tartu"]},
        geometry=gpd.points_from_xy([LISBON[0], TARTU[0]], [LISBON[1], TARTU[1]]),
        crs=4326,
    )


def test_dgselect_minimal_only_given_parameters():
    assert dgselect("IGEO7v2", minimal=True, res=9).metafile == ["dggs_type IGEO7v2", "dggs_res_spec 9"]


def test_dgselect_minimal_mixed_aperture_level_and_precision():
    meta = dgselect("ISEA43HS", minimal=True, res=5, mixed_aperture_level=2, precision=9).metafile
    assert sorted(meta) == ["dggs_num_aperture_4_res 2", "dggs_res_spec 5", "dggs_type ISEA43HS", "precision 9"]


def test_dgselect_minimal_unknown_type():
    with pytest.raises(ValueError, match="unknown dggs_type: 'ISEA5H'"):
        dgselect("ISEA5H", minimal=True, res=5)


def test_dgselect_default_keeps_explicit_parameters():
    meta = dgselect("ISEA7H", res=5).metafile
    assert "dggs_vert0_lon 11.25" in meta
    assert "dggs_proj ISEA" in meta


def test_extra_params_go_to_the_metafile():
    dggs = dgselect("ISEA7HL", minimal=True, res=5)
    dggs.extra_params["polyclip_algo"] = "CLIPPER1"
    assert "polyclip_algo CLIPPER1" in dg_grid_meta(dggs)
    assert "extra_params" not in dggs.to_dict()


def test_extra_params_orientation_is_checked():
    dggs = Dggs(dggs_type="ISEA7HL", projection=None, aperture=None, topology=None, precision=None,
                azimuth_deg=None, pole_lat_deg=None, pole_lon_deg=None)
    dggs.extra_params["dggs_orient_specify_type"] = "NOT_A_TYPE"
    with pytest.raises(ValueError, match="dggs_orient_specify_type"):
        dg_grid_meta(dggs)


def test_all_v9_types_in_list_are_unique():
    assert len(dggs_types_v9) == len(set(dggs_types_v9))


def test_v9_on_older_dggrid_clear_error(monkeypatch):
    dggrid = DGGRIDv9(executable="dggrid", working_dir=tempfile.mkdtemp(), silent=True)
    monkeypatch.setattr(DGGRIDv9, "dggrid_version", (8, 44))
    with pytest.raises(ValueError, match="needs DGGRID 9.1 or newer.*version 8.44"):
        dggrid.run(["dggrid_operation OUTPUT_STATS"])


def test_v9_unknown_version_tries_anyway(monkeypatch):
    dggrid = DGGRIDv9(executable="/does/not/exist/dggrid", working_dir=tempfile.mkdtemp(), silent=True)
    monkeypatch.setattr(DGGRIDv9, "dggrid_version", None)
    assert dggrid.run(["dggrid_operation OUTPUT_STATS"]) == -1


def test_v9_metafile_has_no_implicit_parameters():
    dggrid = _dggrid91()
    dggrid.grid_stats_table("ISEA7HL", 3)
    assert sorted(dggrid.last_ops_meta) == ["dggrid_operation OUTPUT_STATS", "dggs_res_spec 3", "dggs_type ISEA7HL"]


def test_v9_new_parameters_go_to_the_metafile():
    dggrid = _dggrid91()
    dggrid.grid_cellids_for_extent("ISEA7HL", 6, clip_geom=TARTU_BOX, polyclip_algo="CLIPPER1")
    assert "polyclip_algo CLIPPER1" in dggrid.last_ops_meta


@pytest.mark.parametrize("param", ["not_a_parameter", "output_file_name", "clip_region_files", "dggrid_operation"])
def test_v9_unknown_or_managed_parameter_raises(param):
    dggrid = DGGRIDv9(executable="dggrid", working_dir=tempfile.mkdtemp(), silent=True)
    with pytest.raises(ValueError, match=f"'{param}' can't be used with DGGRIDv9"):
        dggrid._select_dggs("ISEA7HL", 5, None, {param: "x"})


@pytest.mark.parametrize("dggs_type", sorted(dggs_types_v9))
def test_v9_every_preset_generates_cells_and_stats(dggs_type):
    dggrid = _dggrid91()
    resolution = 3
    stats = dggrid.grid_stats_table(dggs_type, resolution)
    # the resolutions of SUPERFUND are its own, DGGRID lists more of them in the statistics
    if dggs_type != "SUPERFUND":
        assert len(stats) == resolution + 1
    cells = dggrid.grid_cell_polygons_for_extent(dggs_type, resolution)
    assert len(cells) == stats["Cells"].iloc[-1]
    assert cells["name"].is_unique


@pytest.mark.parametrize("resolution", [5, 6, 9])
def test_v9_igeo7_reference_cells_in_digit_string_form(resolution):
    dggrid = _dggrid91()
    cells = dggrid.cells_for_geo_points(_points(), True, "IGEO7", resolution, output_hier_ndx_form="DIGIT_STRING")
    assert cells["name"].tolist() == [REFERENCE[LISBON][resolution], REFERENCE[TARTU][resolution]]


def test_v9_igeo7_default_form_is_int64_hex_and_round_trips():
    dggrid = _dggrid91()
    ids = dggrid.cells_for_geo_points(_points(), True, "IGEO7", 9)["name"].tolist()
    assert all(len(cell_id) == 16 for cell_id in ids)
    # a Z7 id can't be taken for a sequence number
    polygons = dggrid.grid_cell_polygons_from_cellids(ids, "IGEO7", 9)
    assert sorted(polygons["name"]) == sorted(ids)
    digits = dggrid.address_transform(ids, "IGEO7", 9, output_hier_ndx_form="DIGIT_STRING", output_cell_label_type="OUTPUT_ADDRESS_TYPE")
    assert digits["HIERNDX"].tolist() == [REFERENCE[LISBON][9], REFERENCE[TARTU][9]]


def test_v9_digit_string_input():
    dggrid = _dggrid91()
    digits = [REFERENCE[LISBON][9], REFERENCE[TARTU][9]]
    polygons = dggrid.grid_cell_polygons_from_cellids(digits, "IGEO7", 9, **DIGITS)
    assert sorted(polygons["name"]) == sorted(digits)
    for x, y in (LISBON, TARTU):
        assert polygons.contains(shapely.Point(x, y)).sum() == 1
    centroids = dggrid.grid_cell_centroids_from_cellids(digits, "IGEO7", 9, **DIGITS)
    assert sorted(centroids["name"]) == sorted(digits)


def test_v9_cells_for_geo_points_polygons():
    dggrid = _dggrid91()
    cells = dggrid.cells_for_geo_points(_points(), False, "IGEO7", 9, output_hier_ndx_form="DIGIT_STRING")
    assert cells["zone"].tolist() == [REFERENCE[LISBON][9], REFERENCE[TARTU][9]]
    assert cells["city"].tolist() == ["Lisbon", "Tartu"]


def test_v9_non_hierarchical_preset_stays_seqnum():
    dggrid = _dggrid91()
    polygons = dggrid.grid_cell_polygons_from_cellids([1, 4, 8], "ISEA7HL", 3)
    assert sorted(polygons["name"]) == [1, 4, 8]
    assert dggrid.address_transform([1, 4], "ISEA7HL", 3, output_address_type="Q2DI").columns.tolist() == ["SEQNUM", "Q2DI"]


def test_v9_coarse_cells():
    dggrid = _dggrid91()
    cells = dggrid.grid_cell_polygons_from_cellids(
        ["000102245"], "IGEO7", 9, clip_subset_type="COARSE_CELLS", clip_cell_res=7, **DIGITS,
    )
    assert sum(str(name).startswith("000102245") for name in cells["name"]) == 49


def test_v9_extent_same_cells_as_v8_recipe():
    path = _dggrid91().executable
    v8 = DGGRIDv8(executable=path, working_dir=tempfile.mkdtemp(), capture_logs=True, silent=True)
    from dggrid4py import igeo7_ext
    wrapper = igeo7_ext.dggrid_igeo7_grid_cell_polygons_for_extent(TARTU_BOX, 8, v8)
    cells = _dggrid91().grid_cell_polygons_for_extent("IGEO7", 8, clip_geom=TARTU_BOX, output_hier_ndx_form="DIGIT_STRING")
    assert set(cells["name"]) == set(wrapper["name"])


@pytest.mark.parametrize("func_name", ["grid_cell_polygons_for_extent", "grid_cell_centroids_for_extent"])
def test_v9_extent_hierarchical_form_is_not_dropped(func_name):
    dggrid = _dggrid91()
    cells = getattr(dggrid, func_name)("IGEO7", 6, clip_geom=TARTU_BOX, output_hier_ndx_form="DIGIT_STRING")
    assert all(name.startswith("000") and len(name) == 8 for name in cells["name"])
    assert "output_hier_ndx_form DIGIT_STRING" in dggrid.last_ops_meta


def test_v9_cellids_for_extent_hierarchical_form():
    dggrid = _dggrid91()
    cells = dggrid.grid_cellids_for_extent("IGEO7", 6, clip_geom=TARTU_BOX, output_hier_ndx_form="DIGIT_STRING")
    assert all(len(name) == 8 for name in cells[0])


def test_v9_igeo7_runs_as_igeo7v2_without_log_warning():
    dggrid = _dggrid91()
    result = dggrid.grid_cellids_for_extent("IGEO7", 6, clip_geom=TARTU_BOX)
    assert "dggs_type IGEO7v2" in dggrid.last_ops_meta
    assert "alias" not in dggrid.last_run_logs
    assert len(result) > 0


@pytest.mark.parametrize("dggs_type, is_epsg4326", [("IGEO7", True), ("ISEA7HL", True), ("ISEA7H", False), ("IGEO7v1", False)])
def test_v9_keeps_crs_declared_by_dggrid(dggs_type, is_epsg4326):
    dggrid = _dggrid91()
    cells = dggrid.grid_cell_polygons_for_extent(dggs_type, 2)
    centroids = dggrid.grid_cell_centroids_for_extent(dggs_type, 2)
    assert cells.crs is not None and centroids.crs is not None
    assert (cells.crs.to_epsg() == 4326) == is_epsg4326
    assert (centroids.crs.to_epsg() == 4326) == is_epsg4326


def test_v8_still_forces_epsg4326():
    path = _dggrid91().executable
    v8 = DGGRIDv8(executable=path, working_dir=tempfile.mkdtemp(), capture_logs=True, silent=True)
    assert v8.grid_cell_polygons_for_extent("ISEA7H", 2).crs.to_epsg() == 4326


def test_v9_fuller_is_dymaxion_and_isea_orientation_gives_v8_ids():
    dggrid = _dggrid91()
    v8 = DGGRIDv8(executable=dggrid.executable, working_dir=tempfile.mkdtemp(), capture_logs=True, silent=True)
    v8_ids = v8.cells_for_geo_points(_points(), True, "FULLER7H", 7, output_address_type="SEQNUM")["name"].tolist()
    v9_ids = dggrid.cells_for_geo_points(_points(), True, "FULLER7H", 7)["name"].tolist()
    v9_isea = dggrid.cells_for_geo_points(_points(), True, "FULLER7H", 7, dggs_orient_preset="ISEA")["name"].tolist()
    assert v9_ids != v8_ids
    assert v9_isea == v8_ids


def test_v9_vert0_from_caller_still_used():
    dggrid = _dggrid91()
    dggrid.grid_cellids_for_extent("ISEA7HS", 2, dggs_vert0_lon=11.20)
    assert "dggs_vert0_lon 11.2" in dggrid.last_ops_meta
    assert "dggs_orient_specify_type SPECIFIED" in dggrid.last_ops_meta


def test_v9_bool_parameter_upper_case():
    dggrid = _dggrid91()
    dggrid.grid_cellids_for_extent("ISEA7HL", 6, clip_geom=TARTU_BOX, clip_using_holes=True)
    assert "clip_using_holes TRUE" in dggrid.last_ops_meta


@pytest.mark.parametrize("datum, expected", [
    ("WGS84_AUTHALIC_SPHERE", {"sphere_radius_type": "AUTHALIC_SPHERE"}),
    ("CUSTOM_SPHERE", {"sphere_radius_type": "CUSTOM_SPHERE"}),
    ("WGS84_MEAN_SPHERE", {"sphere_radius_type": "CUSTOM_SPHERE", "custom_sphere_radius": 6371.00877141505983}),
])
def test_v9_proj_datum_deprecated(datum, expected):
    dggrid = DGGRIDv9(executable="dggrid", working_dir=tempfile.mkdtemp(), silent=True)
    with pytest.warns(DeprecationWarning, match="use sphere_radius_type"):
        dggs = dggrid._select_dggs("ISEA7HS", 5, None, {"proj_datum": datum})
    assert dggs.extra_params == expected


def test_v9_proj_datum_invalid():
    dggrid = DGGRIDv9(executable="dggrid", working_dir=tempfile.mkdtemp(), silent=True)
    with pytest.warns(DeprecationWarning), pytest.raises(ValueError, match="proj_datum"):
        dggrid._select_dggs("ISEA7HS", 5, None, {"proj_datum": "MARS"})


def test_v9_proj_datum_runs():
    dggrid = _dggrid91()
    with pytest.warns(DeprecationWarning):
        dggrid.grid_cellids_for_extent("ISEA7HS", 2, proj_datum="WGS84_AUTHALIC_SPHERE")
    assert "sphere_radius_type AUTHALIC_SPHERE" in dggrid.last_ops_meta
