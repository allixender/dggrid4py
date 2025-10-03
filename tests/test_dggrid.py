#!/usr/bin/env python
# -*- coding: utf-8 -*-

from dggrid4py import DGGRIDv8, Dggs

dggrid = DGGRIDv8()

def mock_dggrid_run(__metafile):
    return 0


def test_dgapi_grid_gen_params(monkeypatch):
    monkeypatch.setattr(dggrid, "run", mock_dggrid_run)

    dggs = Dggs(
        dggs_type="IGEO7",
        resolution=4,
        precision=12,
        pole_lon_deg=11.20,
        pole_lat_deg=58.282525588538994675786,
        azimuth_deg=0.0,
    )
    result = dggrid.dgapi_grid_gen(
        dggs,
        {
            "clip_subset_type": "GDAL",
            "clip_region_files": "tests/data/clip_shapefile.shp",
            "geodetic_densify": 0.0,
            "dggs_orient_specify_type": "SPECIFIED",
            "densification": 5,
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
        "dggrid_operation": "GENERATE_GRID",
        "dggs_type": "IGEO7",
        "clip_subset_type": "GDAL",
        "clip_region_files": "tests/data/clip_shapefile.shp",
        "geodetic_densify": 0.0,
        "dggs_orient_specify_type": "SPECIFIED",
        "dggs_vert0_lon": 11.20,
        "dggs_vert0_lat": 58.282525588538994675786,
        "dggs_vert0_azimuth": 0.0,
        "densification": 5,
        "precision": 12,
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
