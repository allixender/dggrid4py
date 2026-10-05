# dggrid4py - a Python library to run highlevel functions of DGGRID

[![PyPI version](https://badge.fury.io/py/dggrid4py.svg)](https://badge.fury.io/py/dggrid4py) [![DOI](https://zenodo.org/badge/295495597.svg)](https://zenodo.org/badge/latestdoi/295495597) [![Documentation Status](https://readthedocs.org/projects/dggrid4py/badge/?version=latest)](https://dggrid4py.readthedocs.io/en/latest/?badge=latest) [GitHub](https://github.com/allixender/dggrid4py/)

[![Population Gridded](day-04-hexa.png)](https://twitter.com/allixender/status/1324055326111485959)

GNU AFFERO GENERAL PUBLIC LICENSE

[DGGRID](https://www.discreteglobalgrids.org/software/) is a free software program for creating and manipulating Discrete Global Grids created and maintained by Kevin Sahr. DGGRID version 8.44 was released 1. December 2025.

- [DGGRID on GitHub](https://github.com/sahrk/DGGRID)
- [DGGRID User Manual](https://github.com/sahrk/DGGRID/blob/d08e10d761f7bedd72a253ab1057458f339de51e/dggridManualV81b.pdf)

dggrid4py runs the `dggrid` command line tool, which has to be available on the system.

### Quickstart with pixi (recommended)

The easiest way to get started is a [pixi](https://pixi.sh) environment. DGGRID (with GDAL) and the geospatial stack come from conda-forge, dggrid4py (and `pygeodesy`) from PyPI:

```bash
pixi init my-dggs-project -c conda-forge
cd my-dggs-project
pixi add python=3.12 dggrid geopandas pyogrio
pixi add --pypi dggrid4py
pixi run dggrid   # prints "usage: dggrid metaFileName"
```

Inside the environment, `dggrid` is on the `PATH`:

```python
import shutil, tempfile
from dggrid4py import DGGRIDv8

dggrid_instance = DGGRIDv8(executable=shutil.which("dggrid"), working_dir=tempfile.mkdtemp(), capture_logs=False, silent=True)
```

### Basic usage

Besides some low-level access to influence the metafile creation of the DGGRID operations, a few high-level functions are integrated to work with the more comfortable geopython libraries, like shapely and geopandas:

- grid_cell_polygons_for_extent(): fill extent/subset with cells at resolution (clip or world);
- grid_cell_centroids_for_extent(): the same for the cell centroids;
- grid_cell_polygons_from_cellids(): geometry_from_cellid for dggs at resolution (from id list);
- grid_cellids_for_extent(): get_all_indexes/cell_ids for dggs at resolution (clip or world);
- cells_for_geo_points(): poly_outline for point/centre at resolution;
- address_transform(): conversion between cell_id address types, like SEQNUM, Q2DI or the hierarchical indexes;
- grid_stats_table(): number of cells, cell area and spacing per resolution.

Use the `DGGRIDv8` class with DGGRID 8.42 or newer (it also runs the 9.0b pre-release). The `DGGRIDv7` class for DGGRID 7 is deprecated since version 0.6.0 and will be removed in a future version.

```python
import shutil
import tempfile

import geopandas
import shapely

from dggrid4py import DGGRIDv8

# create an initial instance that knows where the dggrid tool lives, configure temp workspace and log/stdout output
dggrid_instance = DGGRIDv8(executable=shutil.which("dggrid"), working_dir=tempfile.mkdtemp(), capture_logs=False, silent=True)

# global ISEA4T grid at resolution 5 into GeoDataFrame to Shapefile
gdf1 = dggrid_instance.grid_cell_polygons_for_extent('ISEA4T', 5)
print(gdf1.head())
gdf1.to_file('isea4t_5.shp')

# cell centroids of the global ISEA7H grid at resolution 4
gdf_centroids = dggrid_instance.grid_cell_centroids_for_extent(dggs_type='ISEA7H', resolution=4, mixed_aperture_level=None, clip_geom=None)
print(gdf_centroids.head())

# clip extent
clip_bound = shapely.geometry.box(20.2, 57.00, 28.4, 60.0)

# ISEA7H grid at resolution 9, for extent of provided WGS84 rectangle into GeoDataFrame to Shapefile
gdf3 = dggrid_instance.grid_cell_polygons_for_extent('ISEA7H', 9, clip_geom=clip_bound)
print(gdf3.head())
gdf3.to_file('est_shape_isea7h_9.shp')

# generate cell and areal statistics for a ISEA7H grids from resolution 0 to 8 (return a pandas DataFrame)
df1 = dggrid_instance.grid_stats_table('ISEA7H', 8)
print(df1.head(8))
df1.to_csv('isea7h_8_stats.csv', index=False)

# generate the DGGS grid cells that would cover a GeoDataFrame of points, return Polygons with cell IDs as GeoDataFrame
points = [shapely.Point(20.5, 57.5), shapely.Point(23.5, 58.5)]
geodf_points_wgs84 = geopandas.GeoDataFrame({'site': ['A', 'B']}, geometry=points, crs='EPSG:4326')
gdf4 = dggrid_instance.cells_for_geo_points(geodf_points_wgs84, False, 'ISEA7H', 5)
print(gdf4.head())
gdf4.to_file('polycells_from_points_isea7h_5.shp')

# generate the DGGS grid cells that would cover a GeoDataFrame of points, return a copy of the points with the cell IDs in the column 'name'
gdf5 = dggrid_instance.cells_for_geo_points(geodf_points_wgs84=geodf_points_wgs84, cell_ids_only=True, dggs_type='ISEA4H', resolution=8)
print(gdf5.head())
gdf5.to_file('geopoint_cellids_from_points_isea4h_8.shp')

# generate DGGS grid cell polygons based on 'cell_id_list' (a list or np.array of provided cell_ids)
gdf6 = dggrid_instance.grid_cell_polygons_from_cellids(cell_id_list=[1, 4, 8], dggs_type='ISEA7H', resolution=5)
print(gdf6.head())
gdf6.to_file('from_seqnums_isea7h_5.shp')

# split cells at the dateline for cartesian GIS tools
gdf7 = dggrid_instance.grid_cell_polygons_for_extent('ISEA7H', 3, split_dateline=True)
gdf7.to_file('global_isea7h_3_interrupted.shp')

# convert cell IDs between address types, here from sequence numbers to Q2DI (quad number and i, j coordinates)
df_q2di = dggrid_instance.address_transform([1, 4, 8], 'ISEA7H', 5, input_address_type='SEQNUM', output_address_type='Q2DI')
print(df_q2di.head(3))

# cells at resolution 6 that intersect a coarser cell at resolution 5 (a spatial clip, it includes overlapping neighbours)
gdf8 = dggrid_instance.grid_cell_polygons_from_cellids([100], 'ISEA7H', 6, clip_subset_type='COARSE_CELLS', clip_cell_res=5)
print(len(gdf8))
```

All functions that return geometries return a GeoDataFrame with the cell identifier in the column `name` and with `EPSG:4326` as coordinate reference system. An input or output address type that the DGGRID class does not know raises a `ValueError`.

### IGEO7 and the Z7 index

IGEO7 is the ISEA7H grid with the Z7 hierarchical index. The IGEO7 used in dggrid4py is defined with `dggs_vert0_lon = 11.20` (the DGGRID default is 11.25) and on the WGS84 ellipsoid, which requires a conversion to authalic latitudes around each DGGRID call. The wrappers in `dggrid4py.igeo7_ext` apply both, so that WGS84 data goes in and WGS84 data comes out:

```python
from dggrid4py import igeo7_ext

# Tartu bbox in wgs84, IGEO7 cell polygons at resolution 9 with Z7 cell ids in the column 'name'
extent = shapely.box(26.664593, 58.348705, 26.785607, 58.422495)
cells = igeo7_ext.dggrid_igeo7_grid_cell_polygons_for_extent(extent, 9, dggrid_instance)

# Z7 cell ids for wgs84 points
sites = geopandas.GeoDataFrame({'site': ['Lisbon', 'Tartu']}, geometry=geopandas.points_from_xy([-9.1393, 26.7220], [38.7223, 58.3776]), crs=4326)
sites = igeo7_ext.dggrid_igeo7_cells_for_geo_points(sites, 9, dggrid_instance)
print(sites['name'].tolist())   # ['00641565463', '00010224545']
```

The former address type `Z7_STRING` of the `DGGRIDv7` class is deprecated in DGGRID 8.44 and removed in DGGRID 9. The `DGGRIDv8` class still accepts it with a `DeprecationWarning`. The [IGEO7 docs](https://dggrid4py.readthedocs.io/en/latest/IGEO7.html) describe the background on orientation and ellipsoid, the explicit configuration and the migration.

### Portable DGGRID binary

> **Warning:** the precompiled portable DGGRID binaries are still **experimental** and are built without GDAL support. For a reliable setup, use the pixi or conda-forge installation above.

If DGGRID can neither be installed from conda-forge nor compiled, dggrid4py can download a portable binary from the [DGGRID_portables](https://github.com/allixender/DGGRID_portables) releases (Linux, macOS and Windows, x86_64 and arm64 each):

```python
import tempfile

from dggrid4py import DGGRIDv8, tool

# downloads and unpacks the binary for the current platform into the given folder (only once)
dggrid_exec = tool.get_portable_executable("dggrid_portable")
dggrid_instance_portable = DGGRIDv8(executable=dggrid_exec, working_dir=tempfile.mkdtemp(), capture_logs=False, silent=True, has_gdal=False)

# other lines: "edge" follows the DGGRID development version, or any release tag of DGGRID_portables
dggrid_exec_edge = tool.get_portable_executable("dggrid_portable", line="edge")
```

The default line `stable` is the DGGRID release that this version of dggrid4py is tested with (currently DGGRID 8.44). The archive is verified against the SHA256 checksums of the release. The binaries are built without GDAL, thus the instance has to be created with `has_gdal=False`, and they should be used with the `DGGRIDv8` class.

## TODO:

- get parent_for_cell_id at coarser resolution
- get children_for_cell_id at finer resolution

Remark: with the IGEO7/Z7 index system, the parent of a cell can be derived from the cell identifier. A function that enumerates the index children is still missing.

## Related work:

Originally inspired by [dggridR](https://github.com/r-barnes/dggridR), Richard Barnes’ R interface to DGGRID. However, dggridR is directly linked via Rcpp to DGGRID and calls native C/C++ functions.

After some unsuccessful trials with ctypes, cython, CFFI, pybind11 or cppyy (rather due to lack of experience) I found [am2222/pydggrid](https://github.com/am2222/pydggrid) ([on PyPI](https://pypi.org/project/pydggrid/)) which made apparently some initial scaffolding for the transform operation with [pybind11](https://pybind11.readthedocs.io/en/master/) including some sophisticated conda packaging for Windows. This might be worth following up. Interestingly, its todos include "Adding GDAL export Geometry Support" and "Support GridGeneration using DGGRID" which this dggrid4py module supports with integration of GeoPandas.


## Bundling for different operating systems

Having to compile DGGRID for Windows can be a bit challenging. We are
working on an updated conda package. DGGRID (currently v8.44) is available on conda-forge, which is what the pixi quickstart above uses:

[![Latest version on conda-forge](https://anaconda.org/conda-forge/dggrid/badges/version.svg)](https://anaconda.org/conda-forge/dggrid)

## greater context DGGS in Earth Sciences and GIS

Some reading to be excited about: [discourse.pangeo.io](https://discourse.pangeo.io/t/discrete-global-grid-systems-dggs-use-with-pangeo/2274)
