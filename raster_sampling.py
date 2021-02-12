from pathlib import Path
import os
import io
import sys
import getopt

import math
import numpy as np
import pandas as pd

from shapely.geometry import box, mapping
import rasterio
from rasterio.windows import Window

import logging
import datetime

# create logger
logger = logging.getLogger(__name__)


# define Python user-defined exceptions
class DggsUnknownError(Exception):
    def __init__(self, message="Selected DGGS type unknown or not fully supported!"):
        self.message = message
        super().__init__(self.message)


class NotYetImplemented(Exception):
    def __init__(self, message="Not yet implemented!"):
        self.message = message
        super().__init__(self.message)


try:
    from h3 import h3
except ImportError:
    print("H3 module not available")

try:
    from dggrid4py import DGGRIDv7
except ImportError:
    print("dggrid4py module not avaialble")


try:
    import rhealpixdggs as rhpix
    import rhealpixdggs.dggs as rhpix_dggs
except ImportError:
    print("rhealpixdggs module not avaialble")


def get_dggrid_setup():
    dggrid_exec = None
    try:
        # julia init
        from julia.api import _Julia

        jl = Julia(compiled_modules=False)

        jl.eval("using DGGRID7_jll")
        dirs = jl.eval("DGGRID7_jll.LIBPATH_list")
        # for Windows path separator is ";" and the variable is PATH
        # for linux the path separator is ":" and the variable is LD_LIBRARY_PATH
        path_update = ";".join(dirs)
        os.environ["PATH"] = os.environ["PATH"] + ";" + path_update
        dggrid_exec = jl.eval("DGGRID7_jll.get_dggrid_path()")
        del jl
    except:
        path_update = """C:\\Users\\Alexander\\.julia\\artifacts\\c7ce76fefe6ce8001d5e745e8faaa47ba4228852\\bin;C:\\Users\\Alexander\\.julia\\artifacts\\5ea43503c7796016c3a937d84b6bcf91b1bb608f\\bin;C:\\Users\\Alexander\\.julia\\artifacts\\584ac54e52794a22f516f08a1ca655a7f49e5cae\\bin;C:\\Users\\Alexander\\.julia\\artifacts\\c4c27acf31ad9b69c67a5d1d44910b603936c42a\\bin;C:\\Users\\Alexander\\.julia\\artifacts\\58603237680e9ea3493b336d902e50a332fa1aef\\bin;C:\\Users\\Alexander\\.julia\\artifacts\\d89c49a1ed52a2a6d7c516643590932437d7cbb0\\bin;C:\\Users\\Alexander\\.julia\\artifacts\\a5355b4efc1953f6d15e22618ff5dabbbd0f84dd\\bin;C:\\Users\\Alexander\\.julia\\artifacts\\f70a0c7c2917e5bde1419ffe4fba6e3f5fe78324\\bin;C:\\Users\\Alexander\\.julia\\artifacts\\bdeeb133a67a52079b5206e7f9f3a7cbdf357969\\bin;C:\\Users\\Alexander\\.julia\\artifacts\\00802f9a8a7cef0aa5751ca363507c3fc07de490\\bin;C:\\Users\\Alexander\\.julia\\artifacts\\0b827c9d6fcb7b9f25ede8cf20fcbcab2a497200\\bin;C:\\dev\\Julia_1.5.3\\bin;C:\\dev\\Julia_1.5.3\\bin\\..\\lib\\julia;C:\\dev\\Julia_1.5.3\\bin\\..\\lib;C:\\Users\\Alexander\\.julia\\artifacts\\5c823a847569d0789a68f65aba8aa78b34d4ca3f\\bin;C:\\Users\\Alexander\\.julia\\artifacts\\0407eb942f4c6c53c9bcc1eeccf3def7b46d6cc7\\bin;C:\\Users\\Alexander\\.julia\\artifacts\\0a6628bc2215b9593589daa09bd637c50febfa05\\bin;C:\\Users\\Alexander\\.julia\\artifacts\\42704833bfa671ef47841f64b568a87b9cb25668\\bin;C:\\Users\\Alexander\\.julia\\artifacts\\3f1b12feb59c5ba75346ce47d9ccc46a614c1ca4\\bin;C:\\Users\\Alexander\\.julia\\artifacts\\0f0553e3bdc0535ac962f8b46e21a88536ee2c24\\bin"""
        os.environ["PATH"] = os.environ["PATH"] + ";" + path_update
        dggrid_exec = "C:\\Users\\Alexander\\.julia\\artifacts\\3623dbc61425e9c70e82ed04a4a00af546c1495e\\bin\\dggrid.exe"

    dggrid_instance = DGGRIDv7(
        executable=dggrid_exec, working_dir=os.curdir, capture_logs=False, silent=True
    )

    return dggrid_instance


h3_buf = io.StringIO(
    """H3 Resolution	Average Hexagon Area km2	Average Hexagon Edge Length km	Number of unique indexes
0	4,250,546.8477000	1,107.712591000	122
1	607,220.9782429	418.676005500	842
2	86,745.8540347	158.244655800	5,882
3	12,392.2648621	59.810857940	41,162
4	1,770.3235517	22.606379400	288,122
5	252.9033645	8.544408276	2,016,842
6	36.1290521	3.229482772	14,117,882
7	5.1612932	1.220629759	98,825,162
8	0.7373276	0.461354684	691,776,122
9	0.1053325	0.174375668	4,842,432,842
10	0.0150475	0.065907807	33,897,029,882
11	0.0021496	0.024910561	237,279,209,162
12	0.0003071	0.009415526	1,660,954,464,122
13	0.0000439	0.003559893	11,626,681,248,842
14	0.0000063	0.001348575	81,386,768,741,882
15	0.0000009	0.000509713	569,707,381,193,162"""
)

h3_res = pd.read_csv(h3_buf, sep="\t", thousands=",")
h3_res = h3_res.rename(
    columns={col: col.lower().replace(" ", "_") for col in h3_res.columns}
).set_index("h3_resolution")
h3_res["average_hexagon_area_m2"] = np.float32(
    h3_res["average_hexagon_area_km2"] * 1000000
)
h3_res["average_hexagon_edge_length_m"] = np.float32(
    h3_res["average_hexagon_edge_length_km"] * 1000
)


def get_rhpix_geoid():
    rhpix_geoid = rhpix_dggs.WGS84_003
    return rhpix_geoid


def __haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = (
        math.sin(dlat / 2) ** 2
        + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    )
    c = 2 * math.asin(math.sqrt(a))
    r = 6371  # Radius of earth in kilometers. Use 3956 for miles
    return c * r * 1000


def guess_chunk_size(rs_src, mem_use_mb):
    print("########### guess_chunk_size ##########")
    shapes = []
    for bs in rs_src.block_shapes:
        print("block shapes:" + str(bs))
        shapes.append(bs)
    block_shape = shapes[0]

    # approx chuck size 1000 -> blockshape 256 -> 4 times 256
    # 500 MB mem with float32 byte type, ca 4000x4000 -> find best multiples of blockshape
    # mem_use_mb = 500

    # lat_lon sampling will take extra mem, 2 or 3 times while shifting around
    baselength = math.sqrt(mem_use_mb / 2 * 1024 * 1024 / 32)
    x = int(math.floor(baselength / 2 / block_shape[0]))

    window_blocks_per_chunk_size = x

    print(
        "suggested chunk size -> {} (at approx. window width {} pixels)".format(
            window_blocks_per_chunk_size, str(int(baselength))
        )
    )
    return window_blocks_per_chunk_size


def get_raster_window_defs(rs_src, window_blocks_per_chunk_size):

    block_wins = []

    for ji, window in rs_src.block_windows(1):
        j = ji[0]
        i = ji[1]
        if i == 0:
            block_wins.append([])

        block_wins[j].append((ji, window))

    # print(block_wins[0][0])
    # print(block_wins[-1][-1])
    # print(len(block_wins))

    wj = block_wins
    j_len = len(wj)
    i_len = max([len(wlist) for wlist in block_wins])

    # print(f"w {len(wj)} - j_len ={j_len}, i_len={i_len}")

    # window_blocks_per_chunk_size = 31
    ul_j = 0
    ul_i = 0
    lr_j = 0
    lr_i = window_blocks_per_chunk_size - 1

    j_wise_steps = [
        x
        for x in range(
            window_blocks_per_chunk_size - 1, j_len, window_blocks_per_chunk_size
        )
    ]
    if j_wise_steps[-1] < j_len - 1:
        j_wise_steps.append(j_len - 1)

    i_wise_steps = [
        x
        for x in range(
            window_blocks_per_chunk_size - 1, i_len, window_blocks_per_chunk_size
        )
    ]
    if i_wise_steps[-1] < i_len - 1:
        i_wise_steps.append(i_len - 1)

    for j_counter in j_wise_steps:
        lr_j = j_counter

        for i_counter in i_wise_steps:

            lr_i = i_counter

            # print(f"from [{ul_j}][{ul_i}] to [{lr_j}][{lr_i}]")

            ul = block_wins[ul_j][ul_i][1]
            lr = block_wins[lr_j][lr_i][1]

            # print("-------")
            # print(f"ul: {ul}")
            # print(f"lr: {lr}")
            # print("-------")
            start_c = (ul.col_off, ul.row_off)

            end_c = (
                lr.col_off - ul.col_off + lr.width,
                lr.row_off - ul.row_off + lr.height,
            )

            # print(
            #     f"window(col_off={start_c[0]}, row_off={start_c[1]}, width={end_c[0]}, height={end_c[1]})"
            # )

            cur_win = Window(
                col_off=start_c[0], row_off=start_c[1], width=end_c[0], height=end_c[1]
            )
            yield cur_win

            # print("########")
            ul_i = i_counter + 1

            if ul_i >= max(i_wise_steps):
                ul_i = 0
        # next row
        ul_j = j_counter + 1


def get_extent_for_window(rs_src, window):
    loc_trans = rs_src.window_transform(window)
    extent = rasterio.transform.array_bounds(window.width, window.height, loc_trans)
    return extent


def hexagon_area_from_edgelen(edge_length):
    # Area = (3√3 s2)/ 2
    area = (3 * math.sqrt(3) * edge_length ** 2) / 2
    return area


def hexagon_area_from_radius(radius):
    # apothem is 5√3
    x = radius / math.sqrt(3)
    s = 2 * x
    P = 6 * s
    area = 0.5 * P * radius
    return area


def get_raster_pixel_edge_len(rs_src, adjust_latitudes, pix_size_factor):
    xtransform = rs_src.transform
    xheight = rs_src.height

    a_hor_step = xtransform[0]
    a_vert_step = xtransform[4]

    ax1 = xtransform[2]  # origin
    ay1 = xtransform[5]  # origin
    ax2 = xtransform[2] + a_hor_step  # upper left pixel x2
    ay2 = xtransform[5] + a_vert_step  # upper left pixel y2

    pixel_edge_len = __haversine(ax1, ay1, ax1, ay2)
    if adjust_latitudes:
        widths = []
        lats = []
        ty = ay2
        for i in range(xheight - 1):
            hor_cell_side = __haversine(ax1, ty, ax2, ty)
            widths.append(hor_cell_side)
            lats.append(ty)
            ty = ty + a_vert_step

        dfb = np.array(widths)
        pixel_edge_len = dfb.std() + dfb.min()

    return pixel_edge_len


def suggest_sampling_resolution(
    rs_src, dggstype, adjust_latitudes, pix_size_factor, **kwargs
):
    print("########### suggest_sampling_resolution ##########")
    pixel_edge_len = get_raster_pixel_edge_len(
        rs_src, adjust_latitudes, pix_size_factor
    )
    max_res = None
    if "max_res" in kwargs.keys():
        max_res = kwargs["max_res"]
    resolution = 0
    rhpix_geoid = None
    dggrid_instance = None
    if dggstype.upper() == "H3":

        h3_edge_len = 0
        h3_cell_size = 0

        # Get h3 resolution for raster pixel size
        for idx, row in h3_res.iterrows():
            if row["average_hexagon_edge_length_m"] < pixel_edge_len / pix_size_factor:
                h3_edge_len = row["average_hexagon_edge_length_m"]
                h3_cell_size = row["average_hexagon_area_m2"]
                resolution = idx
                break

        print(
            f"pixel_edge_len: {pixel_edge_len} - resolution: {resolution} - H3 edge len: {h3_edge_len} m - avg cell_size: {h3_cell_size} m2"
        )
    elif dggstype.upper() == "RHEALPIX":
        if "rhpix_geoid" in kwargs.keys():
            rhpix_geoid = kwargs["rhpix_geoid"]
        else:
            rhpix_geoid = get_rhpix_geoid()
        rhpix_res = rhpix_get_res(max_res=max_res) if max_res else rhpix_get_res()

        average_cell_width_m = 0
        # Get resolution for raster pixel size
        for idx, row in rhpix_res.iterrows():
            if row["average_cell_width_m"] < pixel_edge_len / pix_size_factor:
                average_cell_width_m = row["average_cell_width_m"]
                resolution = idx
                break

        print(
            f"pixel_edge_len: {pixel_edge_len} - resolution: {resolution} - RHEALPIX average_cell_width_m: {average_cell_width_m} m"
        )
    elif dggstype.upper() == "DGGRID":
        if "dggrid_instance" in kwargs.keys():
            dggrid_instance = kwargs["dggrid_instance"]
        else:
            dggrid_instance = get_dggrid_setup()

        cls_m = 0
        average_hexagon_area_m2 = 0
        dggrid_res = (
            dggrid_get_res(max_res=max_res, dggrid_instance=dggrid_instance)
            if max_res
            else dggrid_get_res(dggrid_instance=dggrid_instance)
        )
        for idx, row in dggrid_res.iterrows():
            if row["cls_m"] < pixel_edge_len / pix_size_factor:
                cls_m = row["cls_m"]
                average_hexagon_area_m2 = row["average_hexagon_area_m2"]
                resolution = idx
                break

        print(
            f"pixel_edge_len: {pixel_edge_len} - resolution: {resolution} - ISEA7H cls_m: {cls_m} m - avg cell_size: {average_hexagon_area_m2} m2"
        )
    else:
        raise DggsUnknownError()

    return resolution


def h3_gen_ids_points_for_extent(extent, resolution):
    geo_interface = mapping(box(extent[0], extent[1], extent[2], extent[3]))
    h3_idx_gen = h3.polyfill_geojson(geo_interface, res=resolution)
    return h3_idx_gen


# ISEA7H default
def dggrid_gen_ids_points_for_extent(extent, resolution, **kwargs):
    dggrid_instance = None
    if "dggrid_instance" in kwargs.keys():
        dggrid_instance = kwargs["dggrid_instance"]
    else:
        dggrid_instance = get_dggrid_setup()

    geo = box(extent[0], extent[1], extent[2], extent[3])
    df = dggrid_instance.grid_cellids_for_extent("ISEA7H", resolution, clip_geom=geo)
    dggrid_idx_gen = df.values[:, 0]
    _lon = np.asfarray(df.values[:, 1])
    _lat = np.asfarray(df.values[:, 2])
    lat_lon = np.array([_lat, _lon]).T
    return dggrid_idx_gen, lat_lon


def rhpix_suid_to_string(suid):
    return "".join(map(lambda x: str(x), suid))


def rhpix_string_to_suid(s):
    b = s[0]
    rs = map(lambda x: int(x), s[1:])
    return (b, *rs)


def rhpix_gen_ids_points_for_extent(extent, resolution, **kwargs):
    rhpix_geoid = None
    if "rhpix_geoid" in kwargs.keys():
        rhpix_geoid = kwargs["rhpix_geoid"]
    else:
        rhpix_geoid = get_rhpix_geoid()
    se = (extent[1], extent[2])
    nw = (extent[3], extent[0])
    rhpix_cells = list(
        pd.core.common.flatten(
            rhpix_geoid.cells_from_region(resolution, se, nw, plane=False)
        )
    )
    rhpix_idx_gen = np.array([rhpix_suid_to_string(c.suid) for c in rhpix_cells])
    return rhpix_idx_gen, rhpix_cells


def sample_bulk(rs_src, window, lat_lon, full_mask, idx_nodata):
    indices = rs_src.index(lat_lon[:, 1], lat_lon[:, 0])
    xx0_arr = np.array(indices[0]) - window.row_off
    xx1_arr = np.array(indices[1]) - window.col_off

    ins0 = xx0_arr[full_mask]
    ins1 = xx1_arr[full_mask]

    band_array = rs_src.read(1, window=window)
    try:
        return band_array[ins0, ins1]
    except IndexError as ex:
        print(ex)
        return np.full_like(ins0, idx_nodata)


def array_outside_masking(rs_src, window, lat_lon):
    indices = rs_src.index(lat_lon[:, 1], lat_lon[:, 0])
    xx0_arr = np.array(indices[0]) - window.row_off
    xx1_arr = np.array(indices[1]) - window.col_off

    mask0 = xx0_arr < window.height
    mask0_n = xx0_arr >= 0
    mask1 = xx1_arr < window.width
    mask1_n = xx1_arr >= 0
    full_mask = mask0 & mask1 & mask0_n & mask1_n
    return full_mask


def run_pre_check(
    fname, dggstype, adjust_latitudes, pix_size_factor, mem_use_mb, **kwargs
):
    print(
        f"Pre Flight Raster Check: fname {fname}, dggstype {dggstype}, adjust_latitudes {adjust_latitudes}, pix_size_factor {pix_size_factor}, mem_use_mb {mem_use_mb}"
    )
    with rasterio.open(fname, "r") as rs_src:
        print(rs_src.meta)
        target_resolution_base = suggest_sampling_resolution(
            rs_src, dggstype, adjust_latitudes, pix_size_factor, **kwargs
        )
        window_blocks_per_chunk_size = guess_chunk_size(rs_src, mem_use_mb)
        windows = list(get_raster_window_defs(rs_src, window_blocks_per_chunk_size))
        print(f"would create approx. {len(windows)} read/sample/convert windows.")

        # todo
        # if res below x -> consider raster_to_xyz and binning, instead of point query |widht times height points for binning

        return target_resolution_base


def run_binning():
    # todo
    # if res below x -> consider raster_to_xyz and binning, instead of point query |widht times height points for binning
    raise NotYetImplemented()


def run_sampling(
    fname,
    target_resolution,
    dggstype,
    adjust_latitudes,
    pix_size_factor,
    blocks_per_chunk,
    output_format,
    target_dir,
    **kwargs,
):
    p = Path(fname)
    print(f"dir: {p.parent} | file: {p.name}")
    target_dir_path = Path(target_dir)
    print(f"writing all files to {str(target_dir_path.resolve())}")
    target_dcube_data = p.name.split(".")[0]

    # setting up DGGS "engine"
    rhpix_geoid = None
    dggrid_instance = None
    if dggstype.upper() == "H3":
        print("using H3 direct")
    elif dggstype.upper() == "RHEALPIX":
        if "rhpix_geoid" in kwargs.keys():
            rhpix_geoid = kwargs["rhpix_geoid"]
        else:
            rhpix_geoid = get_rhpix_geoid()
        print(f"using RHEALPIX <WGS84_003>")
    elif dggstype.upper() == "DGGRID":
        if "dggrid_instance" in kwargs.keys():
            dggrid_instance = kwargs["dggrid_instance"]
        else:
            dggrid_instance = get_dggrid_setup()
        print(f"using DGGRID <ISEA7H>")
    else:
        raise DggsUnknownError()

    with rasterio.open(fname, "r") as rs_src:
        print(rs_src.meta)
        raster_nodata = rs_src.nodata
        cell_idx_err_nodata = -9999

        target_resolution_base = target_resolution
        if target_resolution is None:
            target_resolution_base = suggest_sampling_resolution(
                rs_src, dggstype, adjust_latitudes, pix_size_factor, **kwargs
            )

        window_iter = get_raster_window_defs(
            rs_src, window_blocks_per_chunk_size=blocks_per_chunk
        )
        print(f"blocks_per_chunk={blocks_per_chunk} fixed")
        # having the windows now  the following loop should be easily
        # parallelizable if enough memory available
        for cur_win in window_iter:
            cur_data_id = f"{target_dcube_data}-{cur_win.col_off}_{cur_win.row_off}-{cur_win.width}_{cur_win.height}-{dggstype}_{target_resolution_base}"

            data_out = str(target_dir_path.resolve() / str(cur_data_id + ".csv"))

            if output_format.upper() == "CSV":
                data_out = str(target_dir_path.resolve() / str(cur_data_id + ".csv"))
            elif output_format.upper() == "PARQUET":
                data_out = str(
                    target_dir_path.resolve() / str(cur_data_id + ".parquet")
                )

            if os.path.isfile(data_out) or os.path.isfile(
                str(target_dir_path.resolve() / str(cur_data_id + ".empty.txt"))
            ):
                print(f"This data slice already sampled ({data_out}/.empty.txt exists)")
                continue
            else:
                print(f"This data slice aimed for {data_out}")

            loc_extent = get_extent_for_window(rs_src, cur_win)

            # nodata short-circuit?
            band_array = rs_src.read(1, window=cur_win)
            if np.any(band_array != raster_nodata):
                del band_array

                print("Generating cell_ids")
                ##### here DGGS specific
                idx_gen = None
                lat_lon = None
                #######
                if dggstype.upper() == "H3":
                    idx_gen = list(
                        h3_gen_ids_points_for_extent(loc_extent, target_resolution_base)
                    )
                    lat_lon = np.array([h3.h3_to_geo(xi) for xi in idx_gen])
                elif dggstype.upper() == "RHEALPIX":
                    idx_gen, rhpix_cells = rhpix_gen_ids_points_for_extent(
                        loc_extent, target_resolution_base, rhpix_geoid=rhpix_geoid
                    )
                    lat_lon = np.array([c.centroid(plane=False) for c in rhpix_cells])
                elif dggstype.upper() == "DGGRID":
                    idx_gen, lat_lon = dggrid_gen_ids_points_for_extent(
                        loc_extent,
                        target_resolution_base,
                        dggrid_instance=dggrid_instance,
                    )
                ###############################

                print(lat_lon.shape)
                print(lat_lon[0])
                # print(rs_src.nodata)

                print("Get raster values for each cell_id")
                positive_mask = array_outside_masking(rs_src, cur_win, lat_lon)
                filtered_cells = np.array(idx_gen)[positive_mask]
                values = sample_bulk(rs_src, cur_win, lat_lon, positive_mask, -9999)

                if filtered_cells.size == values.size:
                    df = pd.DataFrame({"cellids": filtered_cells, "values": values})
                    df["values"] = df["values"].astype(rs_src.meta["dtype"])

                    # Drop nodata
                    df = df[df["values"] != raster_nodata]
                    print(df.head(3))
                    if len(df.index) > 0:
                        idx_err_count = len(
                            df.loc[df.values == cell_idx_err_nodata].index
                        )
                        nan_err_count = np.sum(df.isna().sum())
                        if idx_err_count + nan_err_count > 0:
                            print(
                                f"idx_err_count {idx_err_count} | nan_err_count {nan_err_count}"
                            )

                        df = df[df["values"] != cell_idx_err_nodata].dropna()
                        if output_format.upper() == "CSV":
                            df.to_csv(data_out)
                        elif output_format.upper() == "PARQUET":
                            df.to_parquet(data_out, compression="snappy")
                    else:
                        print(
                            f"This data slice {cur_data_id} is empty (nodata/noindex) after nodata removal"
                        )
                        with open(
                            str(
                                target_dir_path.resolve()
                                / str(cur_data_id + ".empty.txt")
                            ),
                            "w",
                            encoding="utf-8",
                        ) as fh:
                            fh.write(
                                f"This data slice {cur_data_id} is empty (nodata/noindex) after nodata removal"
                            )

                else:
                    print(
                        f"The arrays are not same size for this data slice {cur_data_id} (filter_h3.size {filter_h3.size} != values.size {values.size})"
                    )
            else:
                print(
                    f"This data slice {cur_data_id} is empty (nodata) / all nodata short-ciruit"
                )
                with open(
                    str(target_dir_path.resolve() / str(cur_data_id + ".empty.txt")),
                    "w",
                    encoding="utf-8",
                ) as fh:
                    fh.write(
                        f"This data slice {cur_data_id} is empty (nodata) / all nodata short-ciruit"
                    )


def dggrid_get_res(max_res=16, **kwargs):
    dggrid_instance = None
    if "dggrid_instance" in kwargs.keys():
        dggrid_instance = kwargs["dggrid_instance"]
    else:
        dggrid_instance = get_dggrid_setup()

    isea7h_res = dggrid_instance.grid_stats_table("ISEA7H", max_res)
    isea7h_res = isea7h_res.rename(
        columns={
            "Resolution": "isea7h_resolution",
            "Area (km^2)": "average_hexagon_area_km2",
            "CLS (km)": "cls_km",
        }
    )
    isea7h_res = isea7h_res.rename(
        columns={col: col.lower().replace(" ", "_") for col in isea7h_res.columns}
    ).set_index("isea7h_resolution")
    isea7h_res["average_hexagon_area_m2"] = np.float32(
        isea7h_res["average_hexagon_area_km2"] * 1000000
    )
    isea7h_res["cls_m"] = np.float32(isea7h_res["cls_km"] * 1000)
    return isea7h_res


def rhpix_get_res(max_res=16, **kwargs):
    if "rhpix_geoid" in kwargs.keys():
        rhpix_geoid = kwargs["rhpix_geoid"]
    else:
        rhpix_geoid = get_rhpix_geoid()
    rhpix_resolutions = []
    for i in range(0, max_res, 1):
        rhpix_resolutions.append([i, rhpix_geoid.cell_width(i)])
    rhpix_res = (
        pd.DataFrame(rhpix_resolutions)
        .rename(columns={0: "rhpix_resolution", 1: "average_cell_width_m"})
        .set_index("rhpix_resolution")
    )
    return rhpix_res


def core_dggs_info(dggstype, **kwargs):
    print(f"########### {dggstype} dggs info ##########")
    pd.set_option("display.max_rows", 30)
    pd.set_option("display.max_columns", 6)
    # pd.set_option("precision", 4)
    pd.set_option("display.float_format", lambda n: "{:10.3f}".format(n))

    if dggstype.upper() == "H3":
        print(h3_res)
    elif dggstype.upper() == "RHEALPIX":
        rhpix_res = rhpix_get_res(max_res=16)
        print(rhpix_res)
    elif dggstype.upper() == "DGGRID":
        isea7h_res = dggrid_get_res(max_res=16)
        print(isea7h_res)
    else:
        raise DggsUnknownError()


def info():
    print(
        """
    Supported raster types: technically all GDAL readable which convert into band-arrays, but preferably GeoTIFF/COG (with even blockshapes 128/128 or 256/256)

    - source file: -s <path>

    Supported (in work) DGGS types "-d":
    - Uber H3 <H3>
    - rHealPix DGGS <RHEALPIX>
    - DGGRID ISEA7H <DGGRID>

    run options -o:
    - <info> print this synopsis
    - <precheck> to analyse resolution etc for ingestion, incl. following
      - pixel_size_factor -p 2 (default 2 to capture details)
      - adjust for high latitudes -a y/n (default no)
      - memory usage per window op in MB -m 500 (suggests the number of blocks per chunk)
      - pixel_size_factor -p 2 (default 2 to capture details)
    - <sample> to convert the raster into DGGS type at resolution, list of cell-id and value for your convenience
      - adjust for high latitudes -a y/n (default no)
      - target dggs resolution -r
      - output format -f (<parquet>, <csv>, others maybe too)
      - blocks_per_chunk -b 8 (smaller is usually faster, but too small might cause artifacts)
      - target directory where converted data batches will be deposited -t <target_dir> (must exist, default '.' current)


    sampling workflow (in work):
    ex2: raster_sampling.py -o info
    ex2: raster_sampling.py -o precheck -d H3 -p 3 -a n -s raster.tif
    ex3: raster_sampling.py -o sample -d H3 -p 2 -a y -r 5 -f csv -b 8 -s raster.tif
    """
    )


def main(argv):
    target_resolution = 3
    dggstype = None
    adjust_latitudes = False
    pix_size_factor = 3.0
    output_format = "csv"
    run_option = "info"
    blocks_per_chunk = None
    mem_use_mb = None
    source_file = r"C:\dev\05_geodata\dem\srtm\m30_extracted_eesti\srtm_m30_eesti.tif"
    target_dir = os.curdir

    if len(argv) <= 0:
        info()
        sys.exit(2)

    try:
        opts, args = getopt.getopt(
            argv,
            "o:d:p:a:r:s:b:m:f:t:",
        )
        print(f"opts {opts} args ({args})")
    except getopt.GetoptError:
        print("options error")
        info()
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-o":
            if arg.upper() in ["SAMPLE", "PRECHECK", "INFO"]:
                run_option = arg
            else:
                info()
                sys.exit()
        elif opt == "-d":
            if arg.upper() in ["H3", "DGGRID", "RHEALPIX"]:
                dggstype = arg
            else:
                print(f"opts: dggstype: using default {dggstype}")
        elif opt == "-p":
            pix_size_factor = float(arg)
        elif opt == "-m":
            mem_use_mb = int(arg)
        elif opt == "-b":
            blocks_per_chunk = int(arg)
        elif opt == "-a":
            if arg in ["y", "Y"]:
                adjust_latitudes = True
            else:
                adjust_latitudes = False
        elif opt == "-r":
            target_resolution = int(arg)
        elif opt == "-f":
            if arg.upper() in ["CSV", "PARQUET"]:
                output_format = arg
            else:
                print(f"opts: output_format: using default {output_format}")
        elif opt == "-s":
            if os.path.isfile(arg):
                source_file = arg
            else:
                print(f"opts: source_file: not found")
                sys.exit(2)
        elif opt == "-t":
            if os.path.isdir(arg):
                target_dir = arg
            else:
                print(f"opts: target_dir: not a directory")
                sys.exit(2)

    if run_option.upper() == "INFO":
        if dggstype:
            core_dggs_info(dggstype)
        else:
            info()

    elif run_option.upper() == "SAMPLE":
        if dggstype is None:
            dggstype = "H3"
            print(f"opts: dggstype: using default {dggstype}")

        log_level = logging.INFO
        logger.setLevel(log_level)
        fh = logging.FileHandler(f"raster_sampling-p{my_pid}.log")
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
        )
        fh.setFormatter(formatter)
        # add the handlers to the logger
        logger.addHandler(fh)

        if blocks_per_chunk is None:
            if mem_use_mb is None:
                print("using default approx. memory scaling per window at 500 MB")
                mem_use_mb = 500
            with rasterio.open(source_file, "r") as rs_src:
                blocks_per_chunk = guess_chunk_size(rs_src, mem_use_mb)
        run_sampling(
            source_file,
            target_resolution,
            dggstype,
            adjust_latitudes,
            pix_size_factor,
            blocks_per_chunk,
            output_format,
            target_dir,
        )

    elif run_option.upper() == "PRECHECK":
        if dggstype is None:
            dggstype = "H3"
            print(f"opts: dggstype: using default {dggstype}")
        if mem_use_mb is None:
            print("using default approx. memory scaling per window at 500 MB")
            mem_use_mb = 500
        run_pre_check(
            source_file, dggstype, adjust_latitudes, pix_size_factor, mem_use_mb
        )

    else:
        print("dont know what to do")


if __name__ == "__main__":
    my_pid = os.getpid()

    test_ds = [
        r"C:\dev\05_geodata\dem\srtm\m30_extracted_eesti\srtm_m30_eesti.tif"
        # r"R:\kmoch\datacube_data\ESA_CCI_Landcover\C3S-LC-L4-LCCS-Map-300m-P1Y-2019-v2.1.1.tif",
        # r"R:\kmoch\datacube_data\Chelsea_bioclim\CHELSA_bio10_03_256x256.tif",
        # r"R:\kmoch\datacube_data\HiHydroSoils\top_subsoil\Ksat\Ksat_M_250m_TOPSOIL.tif",
    ]

    main(sys.argv[1:])
