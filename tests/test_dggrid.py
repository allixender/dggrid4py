#!/usr/bin/env python
# -*- coding: utf-8 -*-
import decimal

from dggrid4py import DGGRIDv8, Dggs

dggrid = DGGRIDv8()

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
