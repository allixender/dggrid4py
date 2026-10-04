#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil
import tempfile
import warnings

import pytest
import shapely
import geopandas as gpd

from dggrid4py import DGGRIDv7, DGGRIDv8
from dggrid4py.auxlat import geoseries_to_authalic
from dggrid4py.dggrid_runner import AnyDGGRID, dg_grid_meta, dgselect, get_geo_out, specify_clip_settings

# IGEO7 reference cells: dggs_vert0_lon 11.20, authalic latitudes, Z7 digit strings at resolution 5, 6 and 9
LISBON = (-9.1393, 38.7223)
TARTU = (26.7220, 58.3776)
REFERENCE = {
    LISBON: {5: "0064156", 6: "00641565", 9: "00641565463"},
    TARTU: {5: "0001022", 6: "00010224", 9: "00010224545"},
}
TARTU_BOX = shapely.box(26.6, 58.3, 26.85, 58.45)

Z7 = {
    "input_address_type": "HIERNDX",
    "input_hier_ndx_system": "Z7",
    "input_hier_ndx_form": "DIGIT_STRING",
    "output_address_type": "HIERNDX",
    "output_cell_label_type": "OUTPUT_ADDRESS_TYPE",
    "output_hier_ndx_system": "Z7",
    "output_hier_ndx_form": "DIGIT_STRING",
}
VERT0_FULL = {"dggs_vert0_lon": 11.20, "dggs_vert0_lat": 58.28252559, "dggs_vert0_azimuth": 0.0}


def _dggrid_path():
    path = os.getenv("DGGRID_PATH") or shutil.which("dggrid")
    if not path or not os.path.isfile(path):
        pytest.skip("DGGRID executable not available")
    return path


def _dggrid(cls=DGGRIDv8):
    return cls(executable=_dggrid_path(), working_dir=tempfile.mkdtemp(), capture_logs=True, silent=True)


def _points(authalic=True):
    points = gpd.GeoDataFrame(
        {"city": ["Lisbon", "Tartu"]},
        geometry=gpd.points_from_xy([LISBON[0], TARTU[0]], [LISBON[1], TARTU[1]]),
        crs=4326,
    )
    if authalic:
        points["geometry"] = geoseries_to_authalic(points.geometry)
    return points


# address types: unknown ones raise, the DGGRIDv7 hierarchical index names are mapped on DGGRIDv8

def test_resolve_address_type_known_and_none():
    dggrid = DGGRIDv8(executable="dggrid")
    conf_extra = {}
    assert dggrid.resolve_address_type("output", None, conf_extra) is None
    assert dggrid.resolve_address_type("output", "SEQNUM", conf_extra) == "SEQNUM"
    assert dggrid.resolve_address_type("input", "HIERNDX", conf_extra) == "HIERNDX"
    assert conf_extra == {}


@pytest.mark.parametrize(
    "legacy, system, form",
    [
        ("Z7_STRING", "Z7", "DIGIT_STRING"),
        ("Z7", "Z7", "INT64"),
        ("Z3_STRING", "Z3", "DIGIT_STRING"),
        ("ZORDER", "ZORDER", "INT64"),
    ],
)
def test_resolve_address_type_maps_legacy_names_on_v8(legacy, system, form):
    dggrid = DGGRIDv8(executable="dggrid")

    conf_extra = {}
    with pytest.warns(DeprecationWarning, match=legacy):
        assert dggrid.resolve_address_type("output", legacy, conf_extra) == "HIERNDX"
    assert conf_extra == {
        "output_hier_ndx_system": system,
        "output_hier_ndx_form": form,
        "output_cell_label_type": "OUTPUT_ADDRESS_TYPE",
    }

    conf_extra = {}
    with pytest.warns(DeprecationWarning, match=legacy):
        assert dggrid.resolve_address_type("input", legacy, conf_extra) == "HIERNDX"
    assert conf_extra == {"input_hier_ndx_system": system, "input_hier_ndx_form": form}


def test_resolve_address_type_from_conf_extra_keeps_explicit_fields():
    dggrid = DGGRIDv8(executable="dggrid")
    conf_extra = {"input_address_type": "Z7_STRING", "input_hier_ndx_form": "INT64"}
    with pytest.warns(DeprecationWarning):
        assert dggrid.resolve_address_type("input", None, conf_extra) == "HIERNDX"
    assert conf_extra == {"input_address_type": "HIERNDX", "input_hier_ndx_system": "Z7", "input_hier_ndx_form": "INT64"}


def test_resolve_address_type_v7_keeps_legacy_names():
    dggrid = DGGRIDv7(executable="dggrid")
    conf_extra = {}
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert dggrid.resolve_address_type("output", "Z7_STRING", conf_extra) == "Z7_STRING"
    assert conf_extra == {}
    with pytest.raises(ValueError, match="use DGGRIDv8"):
        dggrid.resolve_address_type("output", "HIERNDX", conf_extra)


@pytest.mark.parametrize("cls", [DGGRIDv7, DGGRIDv8])
@pytest.mark.parametrize(
    "call",
    [
        lambda d: d.grid_cell_polygons_for_extent("IGEO7", 5, clip_geom=TARTU_BOX, output_address_type="NONSENSE"),
        lambda d: d.grid_cell_centroids_for_extent("IGEO7", 5, clip_geom=TARTU_BOX, output_address_type="NONSENSE"),
        lambda d: d.grid_cellids_for_extent("IGEO7", 5, clip_geom=TARTU_BOX, output_address_type="NONSENSE"),
        lambda d: d.grid_cell_polygons_from_cellids(["0001022"], "IGEO7", 5, input_address_type="NONSENSE"),
        lambda d: d.grid_cell_centroids_from_cellids(["0001022"], "IGEO7", 5, output_address_type="NONSENSE"),
        lambda d: d.cells_for_geo_points(_points(), True, "IGEO7", 5, output_address_type="NONSENSE"),
        lambda d: d.address_transform(["0001022"], "IGEO7", 5, input_address_type="NONSENSE"),
        lambda d: d.grid_cell_polygons_for_extent("IGEO7", 5, clip_geom=TARTU_BOX, input_address_type="NONSENSE"),
    ],
)
def test_unknown_address_type_raises_before_dggrid_runs(cls, call, monkeypatch):
    dggrid = cls(executable="dggrid", working_dir=tempfile.mkdtemp())

    def no_run(__metafile):
        raise AssertionError("DGGRID must not run with an unknown address type")

    monkeypatch.setattr(dggrid, "run", no_run)
    with pytest.raises(ValueError, match="unknown (in|out)put_address_type"):
        call(dggrid)


# metafile parameters

def test_clip_settings_seqnum_input():
    tmp_dir = tempfile.mkdtemp()
    # clip_subset_type SEQNUMS is removed in DGGRID 9, it stays for the DGGRIDv7 class only
    settings, _ = specify_clip_settings("WHOLE_EARTH", tmp_dir, "a", input_address_type="SEQNUM", cell_id_list=[1, 4, 8])
    assert settings["clip_subset_type"] == "INPUT_ADDRESS_TYPE"
    assert settings["input_address_type"] == "SEQNUM"

    settings, _ = specify_clip_settings("WHOLE_EARTH", tmp_dir, "b", input_address_type="SEQNUM", cell_id_list=[1, 4, 8], seqnums_clip=True)
    assert settings["clip_subset_type"] == "SEQNUMS"
    assert "input_address_type" not in settings

    # COARSE_CELLS is no longer dropped for SEQNUM input
    settings, _ = specify_clip_settings("COARSE_CELLS", tmp_dir, "c", input_address_type="SEQNUM", cell_id_list=[100], clip_cell_res=5)
    assert settings["clip_subset_type"] == "COARSE_CELLS"
    assert settings["clip_cell_addresses"] == "100"
    assert settings["clip_cell_res"] == 5


def test_centroids_from_cellids_pass_resolution_to_clip_settings(monkeypatch):
    dggrid = DGGRIDv8(executable="dggrid", working_dir=tempfile.mkdtemp())
    metafile = []

    def mock_run(__metafile):
        metafile[:] = __metafile
        return -1

    monkeypatch.setattr(dggrid, "run", mock_run)
    with pytest.raises(ValueError):
        dggrid.grid_cell_centroids_from_cellids(
            ["023255620345"], "IGEO7", 18, clip_subset_type="COARSE_CELLS", clip_cell_res=10, **Z7,
        )
    assert "clipper_scale_factor 100000000" in metafile


@pytest.mark.parametrize("dggs_type", ["ISEA43H", "FULLER43H"])
def test_mixed_aperture_metafile(dggs_type):
    metafile = dg_grid_meta(dgselect(dggs_type, res=3, mixed_aperture_level=2))
    assert "dggs_aperture_type MIXED43" in metafile
    assert "dggs_num_aperture_4_res 2" in metafile
    assert not any(line.startswith("dggs_aperture ") for line in metafile)


def test_geo_out_without_gdal():
    # has_gdal=False alone is enough, FlatGeobuf needs a DGGRID with GDAL
    assert get_geo_out(legacy=False, has_gdal=False)["ext"] == "shp"
    assert DGGRIDv8(executable="dggrid", has_gdal=False).tmp_geo_out["ext"] == "shp"
    assert DGGRIDv8(executable="dggrid", has_gdal=False, tmp_geo_out_legacy=True).tmp_geo_out["ext"] == "shp"


@pytest.mark.parametrize("column", ["name", "Name", "global_id"])
def test_read_geo_out_cell_id_column_and_crs(column):
    path = os.path.join(tempfile.mkdtemp(), "out.shp")
    points = gpd.GeoDataFrame({column: ["1", "2"]}, geometry=[shapely.Point(1, 2), shapely.Point(3, 4)])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # written without a CRS on purpose, as DGGRID with GDAL does
        points.to_file(path, engine="pyogrio")
    gdf = DGGRIDv8(executable="dggrid").read_geo_out(path)
    assert list(gdf.columns) == ["name", "geometry"]
    assert gdf.crs == "EPSG:4326"


def test_any_dggrid_extra_fields():
    dggrid = AnyDGGRID(executable="dggrid")
    assert dggrid.check_output_extra_fields({"output_hier_ndx_form": "DIGIT_STRING"}) == {"output_hier_ndx_form": "DIGIT_STRING"}
    assert dggrid.check_input_extra_fields({"input_hier_ndx_system": "Z7"}) == {"input_hier_ndx_system": "Z7"}


# with DGGRID

@pytest.mark.parametrize("resolution", [5, 6, 9])
def test_igeo7_reference_cells(resolution):
    dggrid = _dggrid()
    expected = [REFERENCE[LISBON][resolution], REFERENCE[TARTU][resolution]]

    lon_only = dggrid.cells_for_geo_points(_points(), True, "IGEO7", resolution, **Z7, dggs_vert0_lon=11.20)
    assert lon_only["name"].tolist() == expected

    all_three = dggrid.cells_for_geo_points(_points(), True, "IGEO7", resolution, **Z7, **VERT0_FULL)
    assert all_three["name"].tolist() == expected


def test_cells_for_geo_points_does_not_modify_input():
    dggrid = _dggrid()
    points = _points()
    expected = [REFERENCE[LISBON][9], REFERENCE[TARTU][9]]

    first = dggrid.cells_for_geo_points(points, True, "IGEO7", 9, **Z7, dggs_vert0_lon=11.20)
    assert list(points.columns) == ["city", "geometry"]
    assert first is not points
    assert list(first.columns) == ["city", "geometry", "lon", "lat", "name"]

    # a second call with the same GeoDataFrame used to send 'lon lon lat lat' to DGGRID
    second = dggrid.cells_for_geo_points(points, True, "IGEO7", 9, **Z7, dggs_vert0_lon=11.20)
    assert first["name"].tolist() == second["name"].tolist() == expected

    # also when the result (with its 'lon', 'lat' and 'name' columns) is passed in again
    third = dggrid.cells_for_geo_points(first, True, "IGEO7", 9, **Z7, dggs_vert0_lon=11.20)
    assert third["name"].tolist() == expected


def test_cells_for_geo_points_polygons_keep_hierarchical_index():
    dggrid = _dggrid()
    points = _points()
    # a third point in the same cell as Tartu: still one row per point
    points.loc[2] = ["Tartu 2", shapely.Point(points.geometry.iloc[1].x + 0.0001, points.geometry.iloc[1].y + 0.0001)]

    cells = dggrid.cells_for_geo_points(points, False, "IGEO7", 5, **Z7, dggs_vert0_lon=11.20)
    assert cells["zone"].tolist() == [REFERENCE[LISBON][5], REFERENCE[TARTU][5], REFERENCE[TARTU][5]]
    assert cells["city"].tolist() == ["Lisbon", "Tartu", "Tartu 2"]
    assert all(cell.contains(point) for cell, point in zip(cells.geometry, points.geometry))


def test_legacy_z7_string_on_v8_equals_hierndx():
    dggrid = _dggrid()

    expected = sorted(dggrid.grid_cell_polygons_for_extent("IGEO7", 6, clip_geom=TARTU_BOX, **Z7)["name"])
    assert len(expected) > 0 and all(len(cell_id) == 8 and cell_id.isdigit() for cell_id in expected)

    with pytest.warns(DeprecationWarning, match="Z7_STRING"):
        polygons = dggrid.grid_cell_polygons_for_extent("IGEO7", 6, clip_geom=TARTU_BOX, output_address_type="Z7_STRING")
    assert sorted(polygons["name"]) == expected

    with pytest.warns(DeprecationWarning, match="Z7_STRING"):
        centroids = dggrid.grid_cell_centroids_for_extent("IGEO7", 6, clip_geom=TARTU_BOX, output_address_type="Z7_STRING")
    assert sorted(centroids["name"]) == expected

    with pytest.warns(DeprecationWarning, match="Z7_STRING"):
        cell_ids = dggrid.grid_cellids_for_extent("IGEO7", 6, clip_geom=TARTU_BOX, output_address_type="Z7_STRING")
    assert sorted(cell_ids[0]) == expected

    with pytest.warns(DeprecationWarning, match="Z7_STRING"):
        points = dggrid.cells_for_geo_points(_points(), True, "IGEO7", 9, output_address_type="Z7_STRING", dggs_vert0_lon=11.20)
    assert points["name"].tolist() == [REFERENCE[LISBON][9], REFERENCE[TARTU][9]]

    # the input side used to return another cell (the Z7 string was read as a sequence number)
    tartu = REFERENCE[TARTU][9]
    with pytest.warns(DeprecationWarning, match="Z7_STRING"):
        from_ids = dggrid.grid_cell_polygons_from_cellids(
            [tartu], "IGEO7", 9, input_address_type="Z7_STRING", output_address_type="Z7_STRING", dggs_vert0_lon=11.20,
        )
    assert from_ids["name"].tolist() == [tartu]

    with pytest.warns(DeprecationWarning, match="Z7_STRING"):
        q2di = dggrid.address_transform([tartu], "IGEO7", 9, input_address_type="Z7_STRING", output_address_type="Q2DI", dggs_vert0_lon=11.20)
    assert q2di["Z7_STRING"].tolist() == [tartu]
    assert len(q2di["Q2DI"].iloc[0].split(" ")) == 3


@pytest.mark.parametrize("func_name", ["grid_cell_polygons_from_cellids", "grid_cell_centroids_from_cellids"])
def test_from_cellids_seqnum(func_name):
    dggrid = _dggrid()
    gdf = getattr(dggrid, func_name)([1, 4, 8], "ISEA7H", 5)
    assert sorted(gdf["name"].tolist()) == [1, 4, 8]
    assert gdf.crs == "EPSG:4326"
    assert "clip_subset_type INPUT_ADDRESS_TYPE" in dggrid.last_ops_meta
    assert "input_address_type SEQNUM" in dggrid.last_ops_meta


def test_coarse_cells_with_seqnum_input():
    dggrid = _dggrid()
    parent = dggrid.grid_cell_polygons_from_cellids([100], "ISEA7H", 5).geometry.iloc[0]
    cells = dggrid.grid_cell_polygons_from_cellids([100], "ISEA7H", 6, clip_subset_type="COARSE_CELLS", clip_cell_res=5)
    # a spatial clip: the 7 cells centred in the parent, plus the neighbours that overlap it
    assert sum(shapely.centroid(cell).within(parent) for cell in cells.geometry) == 7
    assert len(cells) > 7
    assert all(cells.geometry.intersects(parent))


def test_grid_stats_table_igeo7():
    dggrid = _dggrid()
    igeo7 = dggrid.grid_stats_table("IGEO7", 5)
    assert igeo7["Cells"].tolist() == [12, 72, 492, 3432, 24012, 168072]
    assert igeo7.equals(dggrid.grid_stats_table("ISEA7H", 5))


def test_mixed_aperture_grid():
    dggrid = _dggrid()
    assert dggrid.grid_stats_table("ISEA43H", 3, mixed_aperture_level=2)["Cells"].tolist() == [12, 42, 162, 482]
    assert len(dggrid.grid_cell_polygons_for_extent("ISEA43H", 3, mixed_aperture_level=2)) == 482


def test_extent_outputs_have_name_and_crs():
    dggrid = _dggrid()
    polygons = dggrid.grid_cell_polygons_for_extent("ISEA7H", 5, clip_geom=TARTU_BOX)
    centroids = dggrid.grid_cell_centroids_for_extent("ISEA7H", 5, clip_geom=TARTU_BOX)
    for gdf in (polygons, centroids):
        assert list(gdf.columns) == ["name", "geometry"]
        assert gdf.crs == "EPSG:4326"
    assert sorted(polygons["name"]) == sorted(centroids["name"])
