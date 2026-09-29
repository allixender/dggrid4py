#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil
import tempfile

import pytest
import shapely
import geopandas as gpd

from dggrid4py import DGGRIDv8
from dggrid4py.dggrid_runner import specify_orient_type_args


def test_specify_orient_type_args_default_empty():
    assert specify_orient_type_args() == {}


def test_specify_orient_type_args_omits_unset_vert0():
    # issue #46: unset values must not end up as "dggs_vert0_lat None" in the metafile
    assert specify_orient_type_args(dggs_vert0_lon=11.20) == {
        "dggs_orient_specify_type": "SPECIFIED",
        "dggs_vert0_lon": 11.20,
    }


def test_specify_orient_type_args_full_specified_keeps_values():
    conf = specify_orient_type_args(
        dggs_vert0_lon="11.20",
        dggs_vert0_lat="58.282525588538994675786",
        dggs_vert0_azimuth=0.0,
    )
    assert conf == {
        "dggs_orient_specify_type": "SPECIFIED",
        "dggs_vert0_lon": "11.20",
        "dggs_vert0_lat": "58.282525588538994675786",
        "dggs_vert0_azimuth": 0.0,
    }


def test_specify_orient_type_args_specified_without_values():
    assert specify_orient_type_args("SPECIFIED") == {"dggs_orient_specify_type": "SPECIFIED"}


def test_specify_orient_type_args_random():
    assert specify_orient_type_args(dggs_orient_specify_type="RANDOM") == {
        "dggs_orient_specify_type": "RANDOM",
        "dggs_orient_rand_seed": 42,
    }
    assert specify_orient_type_args(orient_type="RANDOM", dggs_orient_rand_seed=7) == {
        "dggs_orient_specify_type": "RANDOM",
        "dggs_orient_rand_seed": 7,
    }


def test_specify_orient_type_args_region_center():
    assert specify_orient_type_args(
        dggs_orient_specify_type="REGION_CENTER", region_center_lon=26.7, region_center_lat=58.4
    ) == {
        "dggs_orient_specify_type": "REGION_CENTER",
        "region_center_lon": 26.7,
        "region_center_lat": 58.4,
    }


@pytest.mark.parametrize(
    "kwargs",
    [
        {"dggs_vert0_lon": 181},
        {"dggs_vert0_lat": "-90.5"},
        {"dggs_vert0_azimuth": 361.0},
        {"dggs_vert0_lon": "abc"},
        {"dggs_vert0_lon": float("nan")},
        {"dggs_orient_specify_type": "UNKNOWN"},
        {"dggs_orient_specify_type": "RANDOM", "dggs_vert0_lon": 11.20},
        {"dggs_orient_specify_type": "REGION_CENTER", "region_center_lat": 91},
    ],
)
def test_specify_orient_type_args_invalid(kwargs):
    with pytest.raises(ValueError):
        specify_orient_type_args(**kwargs)


def _dggrid_executable():
    path = os.getenv("DGGRID_PATH") or shutil.which("dggrid")
    if not path or not os.path.isfile(path):
        pytest.skip("DGGRID executable not available")
    return path


@pytest.mark.parametrize(
    "orient",
    [
        {},
        {"dggs_vert0_lon": 11.20},
        {"dggs_vert0_lon": 11.20, "dggs_vert0_lat": 58.28252559, "dggs_vert0_azimuth": 0.0},
    ],
)
def test_cells_for_geo_points_roundtrip_vert0(orient):
    # issue #46: centroids of generated cells must map back to the same cell ids
    dggrid = DGGRIDv8(executable=_dggrid_executable(), working_dir=tempfile.mkdtemp(), capture_logs=True, silent=True)
    z7 = {
        "output_address_type": "HIERNDX",
        "output_cell_label_type": "OUTPUT_ADDRESS_TYPE",
        "output_hier_ndx_system": "Z7",
        "output_hier_ndx_form": "DIGIT_STRING",
    }
    clip = shapely.box(26.66, 58.34, 26.79, 58.43)
    centroids = dggrid.grid_cell_centroids_for_extent("IGEO7", 9, clip_geom=clip, **z7, **orient)
    assert len(centroids) > 0

    points = gpd.GeoDataFrame({"src": centroids["name"].astype(str).values}, geometry=centroids.geometry.values, crs=4326)
    result = dggrid.cells_for_geo_points(
        points, True, "IGEO7", 9,
        output_address_type="HIERNDX", output_hier_ndx_system="Z7", output_hier_ndx_form="DIGIT_STRING",
        **orient,
    )
    assert (result["src"] == result["name"].astype(str)).all()
