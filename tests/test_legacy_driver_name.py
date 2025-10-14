from dggrid4py import DGGRIDv8
from dggrid4py import tool
import os
from pathlib import Path
import shapely
import pytest
import shutil
import tempfile

clip_bound = shapely.geometry.box(27.2, 57.5, 29.3, 59.2)

# TODO: make this a python tempfile dir
working_dir_legacy = tempfile.TemporaryDirectory()
working_dir = tempfile.TemporaryDirectory()

dggrid_executable_legacy = tool.get_portable_executable('.')
dggrid_executable=os.environ['DGGRID_PATH'] if 'DGGRID_PATH' in os.environ else shutil.which('dggrid83')

portable_dggrid = DGGRIDv8(executable=dggrid_executable_legacy, working_dir=working_dir_legacy.name, capture_logs=False, silent=True, has_gdal=False,
                           tmp_geo_out_legacy=True, debug=False)
gdal_dggrid = DGGRIDv8(executable=dggrid_executable, working_dir=working_dir.name, capture_logs=False, silent=True, has_gdal=True,
                       tmp_geo_out_legacy=False, debug=False)

cellids100 = gdal_dggrid.grid_cell_polygons_for_extent("IGEO7", 3, clip_geom=clip_bound, output_address_type='Z7_STRING')

cellids100 = cellids100['name'][:100].tolist()

# test grid_cell_polygons_for_extent


def test_grid_cell_polygons_for_extent():
    legacy_test = portable_dggrid.grid_cell_polygons_for_extent("IGEO7", 3, clip_geom=clip_bound, output_address_type='Z7_STRING')
    test = gdal_dggrid.grid_cell_polygons_for_extent("IGEO7", 3, clip_geom=clip_bound, output_address_type='Z7_STRING')
    test = set(test.sort_values('name')['name'])
    legacy_test = set(legacy_test.sort_values('global_id')['global_id'])
    assert len(test - legacy_test) == 0


def test_grid_cell_centroid_for_extent():
    legacy_test = portable_dggrid.grid_cell_centroids_for_extent("IGEO7", 3, clip_geom=clip_bound, output_address_type='Z7_STRING')
    test = gdal_dggrid.grid_cell_centroids_for_extent("IGEO7", 3, clip_geom=clip_bound, output_address_type='Z7_STRING')
    test = set(test.sort_values('name')['name'])
    legacy_test = set(legacy_test.sort_values('global_id')['global_id'])
    assert len(test - legacy_test) == 0


def test_grid_cell_polygons_from_cellids():
    legacy_test = portable_dggrid.grid_cell_polygons_from_cellids(cellids100, "IGEO7", 3, input_address_type='Z7_STRING', output_address_type='Z7_STRING')
    test = gdal_dggrid.grid_cell_polygons_from_cellids(cellids100, "IGEO7", 3, input_address_type='Z7_STRING', output_address_type='Z7_STRING')
    test = set(test.sort_values('name')['name'])
    legacy_test = set(legacy_test.sort_values('global_id')['global_id'])
    assert len(test - legacy_test) == 0


def test_grid_cell_centroids_from_cellids():
    legacy_test = portable_dggrid.grid_cell_centroids_from_cellids(cellids100, "IGEO7", 3, input_address_type='Z7_STRING', output_address_type='Z7_STRING')
    test = gdal_dggrid.grid_cell_centroids_from_cellids(cellids100, "IGEO7", 3, input_address_type='Z7_STRING', output_address_type='Z7_STRING')
    test = set(test.sort_values('name')['name'])
    legacy_test = set(legacy_test.sort_values('global_id')['global_id'])
    assert len(test - legacy_test) == 0


def test_grid_cell_polygons_from_cellids_coarse_cells():
    legacy_test = portable_dggrid.grid_cell_polygons_from_cellids(cellids100, "IGEO7", 8, clip_subset_type='COARSE_CELLS',
                                                                  clip_cell_res=3, input_address_type='Z7_STRING', output_address_type='Z7_STRING')
    test = gdal_dggrid.grid_cell_polygons_from_cellids(cellids100, "IGEO7", 8, clip_subset_type='COARSE_CELLS', clip_cell_res=3,
                                                       input_address_type='Z7_STRING', output_address_type='Z7_STRING')

    test = set(test.sort_values('name')['name'])
    legacy_test = set(legacy_test.sort_values('global_id')['global_id'])
    assert len(test - legacy_test) == 0

def test_grid_cell_centroid_from_cellids_coarse_cells():
    legacy_test = portable_dggrid.grid_cell_centroids_from_cellids(cellids100, "IGEO7", 8, clip_subset_type='COARSE_CELLS',
                                                                   clip_cell_res=3, input_address_type='Z7_STRING', output_address_type='Z7_STRING')
    test = gdal_dggrid.grid_cell_centroids_from_cellids(cellids100, "IGEO7", 8, clip_subset_type='COARSE_CELLS', clip_cell_res=3,
                                                        input_address_type='Z7_STRING', output_address_type='Z7_STRING')
    test = set(test.sort_values('name')['name'])
    legacy_test = set(legacy_test.sort_values('global_id')['global_id'])
    assert len(test - legacy_test) == 0


def test_grid_cellids_for_extent():
    legacy_test = portable_dggrid.grid_cellids_for_extent("IGEO7", 9, clip_geom=clip_bound, output_address_type='Z7_STRING')
    test = gdal_dggrid.grid_cellids_for_extent("IGEO7", 9, clip_geom=clip_bound, output_address_type='Z7_STRING')
    assert all(test[0] == legacy_test[0])
