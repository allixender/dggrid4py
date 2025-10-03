# dggrid4py - a Python library to run highlevel functions of DGGRID

[![PyPI version](https://badge.fury.io/py/dggrid4py.svg)](https://badge.fury.io/py/dggrid4py) [![DOI](https://zenodo.org/badge/295495597.svg)](https://zenodo.org/badge/latestdoi/295495597) [![Documentation Status](https://readthedocs.org/projects/dggrid4py/badge/?version=latest)](https://dggrid4py.readthedocs.io/en/latest/?badge=latest) [GitHub](https://github.com/allixender/dggrid4py/)

[![Population Gridded](day-04-hexa.png)](https://twitter.com/allixender/status/1324055326111485959)

GNU AFFERO GENERAL PUBLIC LICENSE

[DGGRID](https://www.discreteglobalgrids.org/software/) is a free software program for creating and manipulating Discrete Global Grids created and maintained by Kevin Sahr. DGGRID version 8.34 was released 13. November 2024

- [DGGRID Version 8.3 on GitHub](https://github.com/sahrk/DGGRID)
- [DGGRID User Manual](https://github.com/sahrk/DGGRID/blob/d08e10d761f7bedd72a253ab1057458f339de51e/dggridManualV81b.pdf)

You need the `dggrid` tool compiled available on the system.

Besides some low-level access influence the dggrid operations' metafile creation, a few high-level functions are integrated to work with the more comfortable geopython libraries, like shapely and geopandas

- grid_cell_polygons_for_extent(): fill extent/subset with cells at resolution (clip or world)
- grid_cell_polygons_from_cellids(): geometry_from_cellid for dggs at resolution (from id list)
- grid_cellids_for_extent(): get_all_indexes/cell_ids for dggs at resolution (clip or world)
- cells_for_geo_points(): poly_outline for point/centre at resolution
- address_transform():  conversion between cell_id address types, like SEQNUM, Z7, or Q2DI


```python
import geopandas
import shapely

from dggrid4py import DGGRIDv7


# create an initial instance that knows where the dggrid tool lives, configure temp workspace and log/stdout output
# if you have 
dggrid_instance = DGGRIDv7(executable='<path_to>/dggrid', working_dir='.', capture_logs=False, silent=False, tmp_geo_out_legacy=False, debug=False)


# global ISEA4T grid at resolution 5 into GeoDataFrame to Shapefile
gdf1 = dggrid_instance.grid_cell_polygons_for_extent('ISEA4T', 5)
print(gdf1.head())
gdf1.to_file('isea4t_5.shp')

gdf_centroids = dggrid_instance.grid_cell_centroids_for_extent(dggs_type='ISEA7H', resolution=4, mixed_aperture_level=None, clip_geom=None)

# clip extent
clip_bound = shapely.geometry.box(20.2,57.00, 28.4,60.0)

# ISEA7H grid at resolution 9, for extent of provided WGS84 rectangle into GeoDataFrame to Shapefile
gdf3 = dggrid_instance.grid_cell_polygons_for_extent('ISEA7H', 9, clip_geom=clip_bound)
print(gdf3.head())
gdf3.to_file('grids/est_shape_isea7h_9.shp')

# generate cell and areal statistics for a ISEA7H grids from resolution 0 to 8 (return a pandas DataFrame)
df1 = dggrid_instance.grid_stats_table('ISEA7H', 8)
print(df1.head(8))
df1.to_csv('isea7h_8_stats.csv', index=False)

# generate the DGGS grid cells that would cover a GeoDataFrame of points, return Polygons with cell IDs as GeoDataFrame
points = [shapely.Point(20.5, 57.5), shapely.Point(21.0, 58.0)]
geodf_points_wgs84 = geopandas.GeoDataFrame({'name': ['A', 'B']}, geometry=points, crs='EPSG:4326')
gdf4 = dggrid_instance.cells_for_geo_points(geodf_points_wgs84, False, 'ISEA7H', 5)
print(gdf4.head())
gdf4.to_file('polycells_from_points_isea7h_5.shp')

# generate the DGGS grid cells that would cover a GeoDataFrame of points, return cell IDs added as column to the points GDF
gdf5 = dggrid_instance.cells_for_geo_points(geodf_points_wgs84=geodf_points_wgs84, cell_ids_only=True, dggs_type='ISEA4H', resolution=8)
print(gdf5.head())
gdf5.to_file('geopoint_cellids_from_points_isea4h_8.shp')

# generate DGGS grid cell polygons based on 'cell_id_list' (a list or np.array of provided cell_ids)
gdf6 = dggrid_instance.grid_cell_polygons_from_cellids(cell_id_list=[1, 4, 8], dggs_type='ISEA7H', resolution=5)
print(gdf6.head())
gdf6.to_file('from_seqnums_isea7h_5.shp')

# v0.2.6 API update split at dateline for cartesian GIS tools
gdf7 = dggrid_instance.grid_cell_polygons_for_extent('ISEA7H', 3, split_dateline=True)
gdf7.to_file('global_isea7h_3_interrupted.shp')

gdf_z1 = dggrid_instance.grid_cell_polygons_for_extent('IGEO7', 5, clip_geom=clip_bound, output_address_type='Z7_STRING')
print(gdf_z1.head(3))

df_z1 = dggrid_instance.guess_zstr_resolution(gdf_z1['name'].values, 'IGEO7', input_address_type='Z7_STRING')
print(df_z1.head(3))

df_q2di = dggrid_instance.address_transform(gdf_z1['name'].values, 'IGEO7', 5, input_address_type='Z7_STRING', output_address_type='Q2DI')
print(df_q2di.head(3))

df_tri = dggrid_instance.address_transform(gdf_z1['name'].values, 'IGEO7', 5, input_address_type='Z7_STRING', output_address_type='PROJTRI')
print(df_tri.head(3))

```


### Portable DGGRID binary

if you don't have a special local distribution of the dggrid-tool or if you didn't install with conda-forge, you can use a provided portable:

```python

from dggrid4py import DGGRIDv7, tool

dggrid_exec = tool.get_portable_executable(".")
dggrid_instance_portable = DGGRIDv7(executable=dggrid_exec, working_dir='.', capture_logs=False, silent=True, has_gdal=False, tmp_geo_out_legacy=True, debug=False)

```

## TODO:

- get parent_for_cell_id at coarser resolution
- get children_for_cell_id at finer resolution

Remark: This is now possible with the IGEO7/Z7 system.

## Related work:

Originally insprired by [dggridR](https://github.com/r-barnes/dggridR), Richard Barnes’ R interface to DGGRID. However, dggridR is directly linked via Rcpp to DGGRID and calls native C/C++ functions.

After some unsuccessful trials with ctypes, cython, CFFI, pybind11 or cppyy (rather due to lack of experience) I found [am2222/pydggrid](https://github.com/am2222/pydggrid) ([on PyPI](https://pypi.org/project/pydggrid/)) which made apparently some initial scaffolding for the transform operation with [pybind11](https://pybind11.readthedocs.io/en/master/) including some sophisticated conda packaging for Windows. This might be worth following up. Interestingly, its todos include "Adding GDAL export Geometry Support" and "Support GridGeneration using DGGRID" which this dggrid4py module supports with integration of GeoPandas.


## Bundling for different operating systems

Having to compile DGGRID for Windows can be a bit challenging. We are
working on an updated conda package. Currently DGGRID v8.3 is available on conda-forge:

[![Latest version on conda-forge](https://anaconda.org/conda-forge/dggrid/badges/version.svg)](https://anaconda.org/conda-forge/dggrid)

## greater context DGGS in Earth Sciences and GIS

Some reading to be excited about: [discourse.pangeo.io](https://discourse.pangeo.io/t/discrete-global-grid-systems-dggs-use-with-pangeo/2274)
