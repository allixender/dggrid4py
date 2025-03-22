import geopandas as gpd
import numpy as np
import pandas as pd

from dggrid4py import igeo7

if __name__ == "__main__":
    grid_fname = 'igeo7_res_9.gpkg'
    gdf = gpd.read_file(grid_fname)

    gdf[['z7_string', 'z7_res', 'parent', 'local_pos', 'is_center' ]] = gdf['name'].apply(igeo7.apply_convert_z7hex_to_z7string)

    # do in notebook only
    # gdf.explore()

    gpd_sindex = gdf.sindex.query(gdf.geometry, predicate="intersects")

    z7_hex_str = '0042aad3ffffffff'

    # lower right corner
    n = igeo7.get_neighbours_by_z7(z7_idx=z7_hex_str,
                             gdf=gdf, gpd_sindex=gpd_sindex,
                             z7_col='name')
    
    z7_string = igeo7.z7hex_to_z7string(z7_hex_str)
    print(f"""Z7 string representation of {z7_hex_str}: {z7_string}""")
    print(f"""Resolution of {z7_hex_str}: {igeo7.get_z7hex_resolution(z7_hex_str)}""")
    
    print(f"""Neighbours of {z7_hex_str}: {n}""")

    parent, local_pos, is_center = igeo7.get_z7hex_local_pos(z7_hex_str)
    print(f"""Parent (in z7_string format) cell of {z7_hex_str}: {parent}""")