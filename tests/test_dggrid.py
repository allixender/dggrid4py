#!/usr/bin/env python
# -*- coding: utf-8 -*-
import decimal
import os

import pytest
import shapely
import geopandas as gpd
from geopandas.testing import assert_geodataframe_equal

from dggrid4py import DGGRIDv8, Dggs

dggrid_path = os.getenv("DGGRID_PATH")
dggrid = DGGRIDv8(executable=dggrid_path)


def mock_dggrid_run(__metafile):
    return 0


def test_dggs_param_handling():
    dggs = Dggs(
        dggs_type="ISEA3H",
        resolution=4,
        precision=12,
        pole_lon_deg=11.20,
        pole_lat_deg=58.282525588538994675786,
        azimuth_deg=0.0,
    )
    assert dggs.dggs_type == "ISEA3H"
    assert dggs.resolution == 4
    assert dggs.precision == 12
    assert dggs.pole_lon_deg == 11.20
    assert dggs.pole_lat_deg == 58.282525588538994675786
    assert dggs.azimuth_deg == 0.0

    assert dggs.get_par("resolution") == 4
    assert dggs.get_par("precision") == 12
    assert dggs.get_par("pole_lon_deg") == 11.20
    assert dggs.get_par("pole_lat_deg") == 58.282525588538994675786
    assert dggs.get_par("azimuth_deg") == 0.0

    assert dggs.get_par("dggs_res_spec") == 4
    assert dggs.get_par("dggs_vert0_lon") == 11.20
    assert dggs.get_par("dggs_vert0_lat") == 58.282525588538994675786
    assert dggs.get_par("dggs_vert0_azimuth") == 0.0


def test_dgapi_grid_gen_params(monkeypatch):
    monkeypatch.setattr(dggrid, "run", mock_dggrid_run)

    dggs = Dggs(
        dggs_type="IGEO7",
        aperture=7,
        resolution=4,
        precision=12,
        densification=5,
        pole_lon_deg="11.20",
        pole_lat_deg=decimal.Decimal("58.282525588538994675786"),
        azimuth_deg=0.0,
        geodetic_densify=0.0,
    )
    result = dggrid.dgapi_grid_gen(
        dggs,
        {
            "clip_subset_type": "GDAL",
            "clip_region_files": "tests/data/clip_shapefile.shp",
            "dggs_orient_specify_type": "SPECIFIED",
        },
        {
            "output_cell_label_type": "OUTPUT_ADDRESS_TYPE",
            "output_address_type": "HIERNDX",
            "output_hier_ndx_system": "Z7",
            "output_hier_ndx_form": "DIGIT_STRING",
            "collection_output_gdal_format": "GeoJSON",
            "collection_output_file_name": "test-collection.geojson",
            "cell_output_type": "GDAL_COLLECTION",
            "cell_output_gdal_format": "GeoJSON",
            "cell_output_file_name": "test-cells.geojson",
            "point_output_type": "GDAL_COLLECTION",
            "point_output_gdal_format": "GeoJSON",
            "point_output_file_name": "test-points.geojson",
            "neighbor_output_type": "GDAL_COLLECTION",
            "neighbor_output_file_name": "test-neighbors.geojson",
            "children_output_type": "GDAL_COLLECTION",
            "children_output_file_name": "test-neighbors.geojson",
        }
    )

    assert set(result["metafile"]) == {
        "dggrid_operation GENERATE_GRID",
        "dggs_type IGEO7",
        "dggs_proj ISEA",
        "dggs_aperture 7",
        "dggs_topology HEXAGON",
        "dggs_res_spec 4",
        "precision 12",
        "densification 5",
        "clip_subset_type GDAL",
        "clip_region_files tests/data/clip_shapefile.shp",
        "geodetic_densify 0.0",
        "dggs_orient_specify_type SPECIFIED",
        "dggs_vert0_lon 11.20",
        "dggs_vert0_lat 58.282525588538994675786",
        "dggs_vert0_azimuth 0.0",
        "output_cell_label_type OUTPUT_ADDRESS_TYPE",
        "output_address_type HIERNDX",
        "output_hier_ndx_system Z7",
        "output_hier_ndx_form DIGIT_STRING",
        "collection_output_gdal_format GeoJSON",
        "collection_output_file_name test-collection.geojson",
        "cell_output_type GDAL_COLLECTION",
        "cell_output_gdal_format GeoJSON",
        "cell_output_file_name test-cells.geojson",
        "point_output_type GDAL_COLLECTION",
        "point_output_gdal_format GeoJSON",
        "point_output_file_name test-points.geojson",
        "neighbor_output_type GDAL_COLLECTION",
        "neighbor_output_file_name test-neighbors.geojson",
        "children_output_type GDAL_COLLECTION",
        "children_output_file_name test-neighbors.geojson",
    }

    assert result["output_conf"] == {
        "output_cell_label_type": "OUTPUT_ADDRESS_TYPE",
        "output_address_type": "HIERNDX",
        "output_hier_ndx_system": "Z7",
        "output_hier_ndx_form": "DIGIT_STRING",
        "collection_output_gdal_format": "GeoJSON",
        "collection_output_file_name": "test-collection.geojson",
        "cell_output_type": "GDAL_COLLECTION",
        "cell_output_gdal_format": "GeoJSON",
        "cell_output_file_name": "test-cells.geojson",
        "point_output_type": "GDAL_COLLECTION",
        "point_output_gdal_format": "GeoJSON",
        "point_output_file_name": "test-points.geojson",
        "neighbor_output_type": "GDAL_COLLECTION",
        "neighbor_output_file_name": "test-neighbors.geojson",
        "children_output_type": "GDAL_COLLECTION",
        "children_output_file_name": "test-neighbors.geojson",
    }


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
            resolution=3,
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
        "dggs_res_spec 3",
        "precision 7",
        "densification 5",
        "geodetic_densify 0.01",
        "clip_subset_type GDAL",
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
            resolution=3,
            cell_id_list=["023"],
            clip_cell_densification=5,
            clip_subset_type="COARSE_CELLS",  # required for densification to take effect
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
        "dggs_res_spec 3",
        "precision 7",
        # "clip_region_files /tmp/dggrid/...",
        # "cell_output_file_name /tmp/dggrid/...",
        "cell_output_type GDAL",
        "cell_output_gdal_format FlatGeobuf",
        "clip_cell_addresses 023",
        "clip_cell_densification 5",
        "clip_subset_type COARSE_CELLS",
        "clip_cell_res 1",
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


def test_cells_for_geo_points(monkeypatch):
    dgapi_grid_transform = dggrid.dgapi_grid_transform
    dgapi_grid_gen = dggrid.dgapi_grid_gen

    dgapi_grid_transform_metafile = []
    dgapi_grid_transform_out_conf = {}

    dgapi_grid_gen_metafile = []
    dgapi_grid_gen_out_conf = {}

    def mock_dgapi_grid_transform(*args, **kwargs):
        _res = dgapi_grid_transform(*args, **kwargs)
        dgapi_grid_transform_metafile[:] = _res["metafile"]
        dgapi_grid_transform_out_conf.update(_res["output_conf"])
        return _res

    def mock_dgapi_grid_gen(*args, **kwargs):
        _res = dgapi_grid_gen(*args, **kwargs)
        dgapi_grid_gen_metafile[:] = _res["metafile"]
        dgapi_grid_gen_out_conf.update(_res["output_conf"])
        return _res

    monkeypatch.setattr(dggrid, "dgapi_grid_transform", mock_dgapi_grid_transform)
    monkeypatch.setattr(dggrid, "dgapi_grid_gen", mock_dgapi_grid_gen)

    points = [shapely.Point(20.5, 57.5), shapely.Point(21.0, 58.0)]
    geodf_points_wgs84 = gpd.GeoDataFrame({'name': ['A', 'B']}, geometry=points, crs='EPSG:4326')

    result = dggrid.cells_for_geo_points(
        geodf_points_wgs84,
        cell_ids_only=False,
        dggs_type="ISEA7H",
        resolution=5,
        # use string to preserve precision and training zeros explicitly
        dggs_vert0_azimuth=0.0,
        dggs_vert0_lat="58.282525588538994675786",  # default: 58.28252559
        dggs_vert0_lon="11.20",   # default: 11.25
    )

    # pre-check temp file paths to ignore in check of specific values
    meta_args = dict([line.split(" ", 1) for line in dgapi_grid_transform_metafile])
    dgapi_grid_transform_input_file_name = meta_args.pop("input_file_name")
    dgapi_grid_transform_output_file_name = meta_args.pop("output_file_name")
    assert dgapi_grid_transform_input_file_name.startswith("/tmp/dggrid")
    assert dgapi_grid_transform_output_file_name.startswith("/tmp/dggrid")
    dgapi_grid_transform_metafile_patched = [f"{key} {val}" for key, val in meta_args.items()]

    meta_args = dict([line.split(" ", 1) for line in dgapi_grid_gen_metafile])
    dgapi_grid_gen_cell_output_file_name = meta_args.pop("cell_output_file_name")
    dgapi_grid_gen_clip_region_files = meta_args.pop("clip_region_files")
    assert dgapi_grid_gen_cell_output_file_name.startswith("/tmp/dggrid")
    assert dgapi_grid_gen_clip_region_files.startswith("/tmp/dggrid")
    dgapi_grid_gen_metafile_patched = [f"{key} {val}" for key, val in meta_args.items()]

    assert dgapi_grid_transform_output_file_name != dgapi_grid_gen_cell_output_file_name

    assert set(dgapi_grid_transform_metafile_patched) == {
        "dggrid_operation TRANSFORM_POINTS",
        "dggs_type ISEA7H",
        "dggs_proj ISEA",
        "dggs_aperture 7",
        "dggs_topology HEXAGON",
        "dggs_res_spec 5",
        "precision 7",
        # following set explicitly by input parameters
        "dggs_orient_specify_type SPECIFIED",
        "dggs_vert0_azimuth 0.0",
        "dggs_vert0_lat 58.282525588538994675786",
        "dggs_vert0_lon 11.20",
        # following enforced by function to align with input GeoDataFrame
        "input_address_type GEO",
        "input_delimiter \" \"",
        "output_address_type SEQNUM",
        "output_delimiter \",\"",
        # pre-checked temp file locations
        # "input_file_name": "/tmp/dggrid/...",
        # "output_file_name": "/tmp/dggrid/...",
    }

    assert dgapi_grid_transform_out_conf == {
        "output_address_type": "SEQNUM",
        "output_delimiter": "\",\"",
        "output_file_name": dgapi_grid_transform_output_file_name,
    }

    assert set(dgapi_grid_gen_metafile_patched) == {
        "dggrid_operation GENERATE_GRID",
        "dggs_type ISEA7H",
        "dggs_proj ISEA",
        "dggs_aperture 7",
        "dggs_topology HEXAGON",
        "dggs_res_spec 5",
        "precision 7",
        # following set explicitly by input parameters
        "dggs_orient_specify_type SPECIFIED",
        "dggs_vert0_azimuth 0.0",
        "dggs_vert0_lat 58.282525588538994675786",
        "dggs_vert0_lon 11.20",
        # following enforced by function to align with input of previous transform step
        "clip_subset_type SEQNUMS",
        "cell_output_type GDAL",
        "cell_output_gdal_format FlatGeobuf",
        "point_output_type NONE",
        # pre-checked temp file locations
        # "cell_output_file_name": "/tmp/dggrid/...",
    }

    assert dgapi_grid_gen_out_conf == {
        "cell_output_type": "GDAL",
        "cell_output_gdal_format": "FlatGeobuf",
        "cell_output_file_name": dgapi_grid_gen_cell_output_file_name,
        "point_output_type": "NONE",
    }

    expect = gpd.GeoDataFrame.from_features(
        [
            {
                "type": "Feature",
                "properties": {"zone": "51695", "name": "A", "lon": 20.5, "lat": 57.5},
                "geometry": {"type": "Polygon", "coordinates": [[
                    [21.189511635794823, 58.2893639588515],
                    [20.950622414770574, 58.00135094074029],
                    [21.232988751914608, 57.69437147048094],
                    [21.74777706813212, 57.674593541688374],
                    [21.990765226755794, 57.96161815328658],
                    [21.71493417563983, 58.2694113005297],
                    [21.189511635794823, 58.2893639588515],
                ]]}
            },
            {
                "type": "Feature",
                "properties": {"zone": "51548", "name": "B", "lon": 21.0, "lat": 58.0},
                "geometry": {"type": "Polygon", "coordinates": [[
                    [20.430061092474467, 58.01819524792648],
                    [20.20251122281021, 57.727966938664856],
                    [20.491335795187357, 57.42178746002092],
                    [21.001317147943144, 57.405038598901655],
                    [21.232988751914608, 57.69437147048094],
                    [20.950622414770574, 58.00135094074029],
                    [20.430061092474467, 58.01819524792648],
                ]]}
            },
        ],
        columns=result.columns,  # ensure ordering matches to allow compare
    )
    assert_geodataframe_equal(result, expect)
