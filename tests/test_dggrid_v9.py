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
