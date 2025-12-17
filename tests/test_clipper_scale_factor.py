#!/usr/bin/env python
# -*- coding: utf-8 -*-
import decimal
import inspect
import os
import tempfile

import pytest
import shapely
import geopandas as gpd
from geopandas.testing import assert_geodataframe_equal

from dggrid4py import DGGRIDv8, Dggs

dggrid_path = os.getenv("DGGRID_PATH")
dggrid = DGGRIDv8(executable=dggrid_path, debug=True)


def mock_dggrid_run(__metafile):
    return 0


def test_grid_cell_polygons_for_extent(monkeypatch):
    metafile = []

    def mock_dggrid_grid_gen_run(__metafile):
        metafile[:] = __metafile
        return -1  # cause grid_gen to early-exit in error

    monkeypatch.setattr(dggrid, "run", mock_dggrid_grid_gen_run)

    clip_bound = shapely.geometry.box(27.2, 57.5, 29.3, 59.2)
    with pytest.raises(ValueError):  # catch and ignore (early-abort "run error")
        dggrid.grid_cell_polygons_for_extent(
            dggs_type="IGEO7",
            resolution=17,
            densification=5,
            geodetic_densify=0.01,
            clip_geom=clip_bound,
            # use string to preserve precision and training zeros explicitly
            dggs_vert0_azimuth=0.0,
            dggs_vert0_lat="58.282525588538994675786",  # default: 58.28252559
            dggs_vert0_lon="11.20",   # default: 11.25
        )

    # pre-check temp file paths to ignore in check of specific values
    meta_args = dict([line.split(" ") for line in metafile])
    assert meta_args["clip_region_files"].startswith("/tmp/dggrid")
    assert meta_args["cell_output_file_name"].startswith("/tmp/dggrid")
    meta_args.pop("clip_region_files")
    meta_args.pop("cell_output_file_name")
    metafile_patched = [f"{key} {val}" for key, val in meta_args.items()]
    assert set(metafile_patched) == {
        "dggrid_operation GENERATE_GRID",
        "dggs_type IGEO7",
        "dggs_proj ISEA",
        "dggs_aperture 7",
        "dggs_topology HEXAGON",
        # NOTE: 'dggs_res_spec' to be inferred from Dggs() attribute, not passing it explicitly to 'specify_resolution'
        "dggs_res_spec 17",
        "precision 7",
        "densification 5",
        "geodetic_densify 0.01",
        "clip_subset_type GDAL",
        "clipper_scale_factor 10000000",
        # "clip_region_files /tmp/dggrid/...",
        # "cell_output_file_name /tmp/dggrid/...",
        "cell_output_type GDAL",
        "cell_output_gdal_format FlatGeobuf",
        # following set explicitly by input parameters
        "dggs_orient_specify_type SPECIFIED",
        "dggs_vert0_azimuth 0.0",
        "dggs_vert0_lat 58.282525588538994675786",
        "dggs_vert0_lon 11.20",
        # WARNING: following technically not set by Dggs(), though it probably should for 'IGEO7' ?
        # "output_cell_label_type OUTPUT_ADDRESS_TYPE",
        # "output_address_type HIERNDX",
        # "output_hier_ndx_system Z7",
        # "output_hier_ndx_form DIGIT_STRING",
        "point_output_type NONE"
    }

def test_grid_cell_polygons_from_cellids(monkeypatch):
    metafile = []

    def mock_dggrid_grid_gen_run(__metafile):
        metafile[:] = __metafile
        return -1  # cause grid_gen to early-exit in error

    monkeypatch.setattr(dggrid, "run", mock_dggrid_grid_gen_run)

    with pytest.raises(ValueError):  # catch and ignore (early-abort "run error")
        dggrid.grid_cell_polygons_from_cellids(
            dggs_type="IGEO7",
            resolution=18,
            cell_id_list=["023255620345"],
            # clip_cell_densification=5,
            clip_subset_type="COARSE_CELLS",  # required for densification to take effect
            clip_cell_res=10,
            input_address_type="HIERNDX",
            input_hier_ndx_forms="DIGIT_STRING",
            input_hier_ndx_systems="Z7",
            output_cell_label_type="OUTPUT_ADDRESS_TYPE",
            output_address_type="HIERNDX",
            output_hier_ndx_forms="DIGIT_STRING",
            output_hier_ndx_systems="Z7",
            # use string to preserve precision and training zeros explicitly
            dggs_vert0_azimuth=0.0,
            dggs_vert0_lat="58.282525588538994675786",  # default: 58.28252559
            dggs_vert0_lon="11.20",   # default: 11.25
        )

    # pre-check temp file paths to ignore in check of specific values
    meta_args = dict([line.split(" ") for line in metafile])
    assert meta_args["clip_region_files"].startswith("/tmp/dggrid")
    assert meta_args["cell_output_file_name"].startswith("/tmp/dggrid")
    meta_args.pop("clip_region_files")
    meta_args.pop("cell_output_file_name")
    metafile_patched = [f"{key} {val}" for key, val in meta_args.items()]

    assert set(metafile_patched) == {
        "dggrid_operation GENERATE_GRID",
        "dggs_type IGEO7",
        "dggs_proj ISEA",
        "dggs_aperture 7",
        "dggs_topology HEXAGON",
        # NOTE: 'dggs_res_spec' to be inferred from Dggs() attribute, not passing it explicitly to 'specify_resolution'
        "dggs_res_spec 18",
        "precision 7",
        # "clip_region_files /tmp/dggrid/...",
        # "cell_output_file_name /tmp/dggrid/...",
        "cell_output_type GDAL",
        "cell_output_gdal_format FlatGeobuf",
        "clip_cell_addresses 023255620345",
        "clip_subset_type COARSE_CELLS",
        "clip_cell_res 10",
        "clipper_scale_factor 100000000",
        # following set explicitly by input parameters
        "dggs_orient_specify_type SPECIFIED",
        "dggs_vert0_azimuth 0.0",
        "dggs_vert0_lat 58.282525588538994675786",
        "dggs_vert0_lon 11.20",
        "input_address_type HIERNDX",
        "output_cell_label_type OUTPUT_ADDRESS_TYPE",
        "output_address_type HIERNDX",
        # WARNING: following technically not set by Dggs(), though it probably should for 'IGEO7' ?
        # "output_hier_ndx_system Z7",
        # "output_hier_ndx_form DIGIT_STRING",
        "point_output_type NONE"
    }

