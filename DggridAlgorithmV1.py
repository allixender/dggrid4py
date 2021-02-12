# -*- coding: utf-8 -*-

"""
***************************************************************************
*                                                                         *
*   This program is free software; you can redistribute it and/or modify  *
*   it under the terms of the GNU General Public License as published by  *
*   the Free Software Foundation; either version 2 of the License, or     *
*   (at your option) any later version.                                   *
*                                                                         *
***************************************************************************
"""

from qgis.PyQt.QtCore import QCoreApplication, QVariant

from qgis.core import (
    QgsProcessing,
    Qgis,
    QgsMessageLog,
    QgsRaster,
    QgsRasterLayer,
    QgsProcessingException,
    QgsProcessingAlgorithm,
    QgsProcessingParameterVectorLayer,
    QgsProcessingParameterVectorDestination,
    QgsProcessingParameterExtent,
    QgsProcessingParameterFile,
    QgsProcessingParameterNumber,
    QgsProcessingParameterBoolean,
    QgsProcessingParameterEnum,
    QgsProcessingParameterString,
    QgsVectorLayer,
    QgsField,
)
from qgis import processing

from osgeo import gdal, osr, ogr

import sys
import os

import numpy as np
import numpy.ma as ma

from pathlib import Path
import uuid
import shutil
import os
import sys
import subprocess
import traceback
import tempfile
import platform

# specify a ISEA3H
dggs_types = (
    "CUSTOM",  # parameters will be specified manually
    "SUPERFUND",  # Superfund_500m grid
    "PLANETRISK",
    "ISEA3H",  # ISEA projection with hexagon cells and an aperture of 3
    "ISEA4H",  # ISEA projection with hexagon cells and an aperture of 4
    "ISEA4T",  # ISEA projection with triangle cells and an aperture of 4
    "ISEA4D",  # ISEA projection with diamond cells and an aperture of 4
    "ISEA43H",  # ISEA projection with hexagon cells and a mixed sequence of aperture 4 resolutions followed by aperture 3 resolutions
    "ISEA7H",  # ISEA projection with hexagon cells and an aperture of 7
    "FULLER3H",  # FULLER projection with hexagon cells and an aperture of 3
    "FULLER4H",  # FULLER projection with hexagon cells and an aperture of 4
    "FULLER4T",  # FULLER projection with triangle cells and an aperture of 4
    "FULLER4D",  # FULLER projection with diamond cells and an aperture of 4
    "FULLER43H",  # FULLER projection with hexagon cells and a mixed sequence of aperture 4 resolutions followed by aperture 3 resolutions
    "FULLER7H",  # FULLER projection with hexagon cells and an aperture of 7
)

# control grid generation
clip_subset_types = ("SHAPEFILE", "WHOLE_EARTH", "GDAL", "AIGEN", "SEQNUMS")

# specify the output
cell_output_types = ("AIGEN", "GDAL", "GEOJSON", "SHAPEFILE", "NONE", "TEXT")

gdal_ogr_vector_extension_to_driver = {
    "shp": "ESRI Shapefile",
    "gpkg": "GPKG",
    "geojson": "GeoJSON",
    "json": "GeoJSON",
    "bna": "BNA",
    "dxf": "DXF",
    "csv": "CSV",
    "dbf": "ESRI Shapefile",
    "geojsonl": "GeoJSONSeq",
    "geojsons": "GeoJSONSeq",
    "gml": "GML",
    "xml": "GML",
    "gmt": "OGR_GMT",
    "gdb": "FileGDB",
    "gpx": "GPX",
    "gtm": "GPSTrackMaker",
    "gtz": "GPSTrackMaker",
    "tab": "MapInfo File",
    "mif": "MapInfo File",
    "mid": "MapInfo File",
    "dgn": "DGN",
    "pix": "PCIDSK",
}

output_address_types = (
    "GEO",  # geodetic coordinates -123.36 43.22 20300 Roseburg
    "Q2DI",  # quad number and (i, j) coordinates on that quad
    "SEQNUM",  # DGGS index - linear address (1 to size-of-DGG), not supported for parameter input_address_type if dggs_aperture_type is SEQUENCE
    "INTERLEAVE",  # digit-interleaved form of Q2DI, only supported for parameter output_address_type; only available for hexagonal aperture 3 and 4 grids
    "PLANE",  # (x, y) coordinates on unfolded ISEA plane,  only supported for parameter output_address_type;
    "Q2DD",  # quad number and (x, y) coordinates on that quad
    "PROJTRI",  # PROJTRI - triangle number and (x, y) coordinates within that triangle on the ISEA plane
    "VERTEX2DD",  # vertex number, triangle number, and (x, y) coordinates on ISEA plane
    "AIGEN",  # Arc/Info Generate file format
)

input_address_types = (
    "GEO",  # geodetic coordinates -123.36 43.22 20300 Roseburg
    "Q2DI",  # quad number and (i, j) coordinates on that quad
    "SEQNUM",  # DGGS index - linear address (1 to size-of-DGG), not supported for parameter input_address_type if dggs_aperture_type is SEQUENCE
    "Q2DD",  # quad number and (x, y) coordinates on that quad
    "PROJTRI",  # PROJTRI - triangle number and (x, y) coordinates within that triangle on the ISEA plane
    "VERTEX2DD",  # vertex number, triangle number, and (x, y) coordinates on ISEA plane
    "AIGEN",  # Arc/Info Generate file format
)

### CUSTOM args
dggs_projections = ("ISEA", "FULLER")

dggs_topologies = ("HEXAGON", "TRIANGLE", "DIAMOND")
dggs_aperture_types = ("PURE", "MIXED43", "SEQUENCE")

dggs_res_specify_types = ("SPECIFIED", "CELL_AREA", "INTERCELL_DISTANCE")
dggs_orient_specify_types = ("SPECIFIED", "RANDOM", "REGION_CENTER")

# specify the operation
dggrid_operations = (
    "GENERATE_GRID",
    "TRANSFORM_POINTS",
    "BIN_POINT_VALS",
    "BIN_POINT_PRESENCE",
    "OUTPUT_STATS",
)

"""
class representing a DGGS grid system configuration, projection aperture etc
"""


class Dggs(object):

    """
    dggs_type: str     # = 'CUSTOM'
    projection: str     # = 'ISEA'
    aperture: int      #  = 3
    topology: str      #  = 'HEXAGON'
    res: int           #  = None
    precision: int     #  = 7
    area: float         #   = None
    spacing: float       #  = None
    cls_val: float        #     = None
    resround: str      #  = 'nearest'
    metric: bool        #  = True
    show_info: bool     #  = True
    azimuth_deg: float   #  = 0
    pole_lat_deg: float  #  = 58.28252559
    pole_lon_deg: float  # = 11.25

    mixed_aperture_level:  # e.g. 5 -> dggs_num_aperture_4_res 5  for ISEA_43_H etc
    metafile = []
    """

    def __init__(self, dggs_type, **kwargs):
        self.dggs_type = dggs_type

        self.init_args = kwargs.keys()

        for key, value in kwargs.items():
            self.set_par(key, value)

    # def dgverify(self):
    #     #See page 21 of documentation for further bounds
    #     if not self.projection in ['ISEA','FULLER']:
    #         raise ValueError('Unrecognised dggs projection')
    #
    #     if not self.topology in ['HEXAGON','DIAMOND','TRIANGLE']:
    #         raise ValueError('Unrecognised dggs topology')
    #     if not self.aperture in [ 3, 4 ]:
    #         raise ValueError('Unrecognised dggs aperture')
    #     if self.res < 0:
    #         raise ValueError('dggs resolution must be >=0')
    #     if self.res > 30:
    #         raise ValueError('dggs resolution must be <=30')
    #     if self.azimuth_deg < 0 or self.azimuth_deg > 360:
    #         raise ValueError('dggs azimuth_deg must be in the range [0,360]')
    #     if self.pole_lat_deg < -90  or self.pole_lat_deg > 90:
    #         raise ValueError('dggs pole_lat_deg must be in the range [-90,90]')
    #     if self.pole_lon_deg < -180 or self.pole_lon_deg > 180:
    #         raise ValueError('dggs pole_lon_deg must be in the range [-180,180]')

    def set_par(self, par_key, par_value):
        if par_key == "dggs_type":
            self.dggs_type = par_value
        if par_key == "projection":
            self.projection = par_value
        if par_key == "aperture":
            self.aperture = par_value
        if par_key == "topology":
            self.topology = par_value
        if par_key == "res":
            self.res = par_value
        if par_key == "precision":
            self.precision = par_value
        if par_key == "area":
            self.area = par_value
        if par_key == "spacing":
            self.spacing = par_value
        if par_key == "cls_val":
            self.cls_val = par_value
        if par_key == "resround":
            self.resround = par_value
        if par_key == "metric":
            self.metric = par_value
        if par_key == "show_info":
            self.show_info = par_value
        if par_key == "azimuth_deg":
            self.azimuth_deg = par_value
        if par_key == "pole_lat_deg":
            self.pole_lat_deg = par_value
        if par_key == "pole_lon_deg":
            self.pole_lon_deg = par_value
        if par_key == "mixed_aperture_level":
            self.mixed_aperture_level = par_value

        return self

    def get_par(self, par_key, alternative=None):
        if par_key == "dggs_type":
            try:
                return self.dggs_type
            except AttributeError:
                return alternative
        if par_key == "projection":
            try:
                return self.projection
            except AttributeError:
                return alternative
        if par_key == "aperture":
            try:
                return self.aperture
            except AttributeError:
                return alternative
        if par_key == "topology":
            try:
                return self.topology
            except AttributeError:
                return alternative
        if par_key == "res":
            try:
                return self.res
            except AttributeError:
                return alternative
        if par_key == "precision":
            try:
                return self.precision
            except AttributeError:
                return alternative
        if par_key == "area":
            try:
                return self.area
            except AttributeError:
                return alternative
        if par_key == "spacing":
            try:
                return self.spacing
            except AttributeError:
                return alternative
        if par_key == "cls_val":
            try:
                return self.cls_val
            except AttributeError:
                return alternative
        if par_key == "resround":
            try:
                return self.resround
            except AttributeError:
                return alternative
        if par_key == "metric":
            try:
                return self.metric
            except AttributeError:
                return alternative
        if par_key == "show_info":
            try:
                return self.show_info
            except AttributeError:
                return alternative
        if par_key == "azimuth_deg":
            try:
                return self.azimuth_deg
            except AttributeError:
                return alternative
        if par_key == "pole_lat_deg":
            try:
                return self.pole_lat_deg
            except AttributeError:
                return alternative
        if par_key == "pole_lon_deg":
            try:
                return self.pole_lon_deg
            except AttributeError:
                return alternative
        if par_key == "mixed_aperture_level":
            try:
                return self.mixed_aperture_level
            except AttributeError:
                return alternative
        else:
            return alternative

    def __str__(self):
        txt = []
        for key in self.init_args:
            txt.append(f"{key}:{self.get_par(key)}")
        return "DGGS(" + ", ".join(txt) + ")"

    def dg_closest_res_to_area(self, area, resround, metric, show_info=True):
        raise ValueError("not yet implemented")

    def dg_closest_res_to_spacing(self, spacing, resround, metric, show_info=True):
        raise ValueError("not yet implemented")

    def dg_closest_res_to_cls(self, cls_val, resround, metric, show_info=True):
        raise ValueError("not yet implemented")


class DggridGridGenAlgorithmV1(QgsProcessingAlgorithm):
    """
    runs DGGRID tool and constructs a DGGS as vector layer
    """

    INPUT_DGGS_TYPE = "INPUT_DGGS_TYPE"
    INPUT_DGGS_RESOLUTION = "INPUT_DGGS_RESOLUTION"

    # mixed_aperture_level
    INPUT_DGGS_MIXED_APERTURE_LEVEL = "INPUT_DGGS_MIXED_APERTURE_LEVEL"

    INPUT_CLIP_EXTENT = "INPUT_CLIP_EXTENT"
    INPUT_EXTENT_LAYER = "INPUT_EXTENT_LAYER"

    INPUT_DGGRID_EXECUTABLE = "INPUT_DGGRID_EXECUTABLE"
    INPUT_DGGRID_EXECUTABLE_LIBDIRS = "INPUT_DGGRID_EXECUTABLE_LIBDIRS"

    INPUT_MARK_DATELINE_CROSSED = "INPUT_MARK_DATELINE_CROSSED"

    OUTPUT_GRID = "OUTPUT_GRID"

    def tr(self, string):
        """
        Returns a translatable string with the self.tr() function.
        """
        return QCoreApplication.translate("DGGRID GridGen Algorithm v1", string)

    def createInstance(self):
        return DggridGridGenAlgorithmV1()

    def name(self):
        """
        Returns the algorithm name, used for identifying the algorithm.
        """
        return "dggridgridgenalgorithmv1"

    def displayName(self):
        """
        Returns the translated algorithm name.
        """
        return self.tr("DGGRID Grid Generating Algorithm V1 Script")

    def group(self):
        """
        Returns the name of the group this algorithm belongs to.
        """
        return self.tr("DGGS scripts")

    def groupId(self):
        """
        Returns the unique ID of the group this algorithm belongs to.
        """
        return "dggsscripts"

    def shortHelpString(self):
        """
        Returns a localised short helper string for the algorithm.
        """
        return self.tr(
            "DGGRID GridGen Processing Algorithm short description, constructs DGGS with DGGRID"
        )

    def initAlgorithm(self, config=None):
        """
        Here we define the inputs and output of the algorithm, along
        with some other properties.
        """

        self.INPUT_DGGS_TYPE_OPTIONS_LIST = [
            "SUPERFUND",  # Superfund_500m grid
            "PLANETRISK",
            "ISEA3H",  # ISEA projection with hexagon cells and an aperture of 3
            "ISEA4H",  # ISEA projection with hexagon cells and an aperture of 4
            "ISEA4T",  # ISEA projection with triangle cells and an aperture of 4
            "ISEA4D",  # ISEA projection with diamond cells and an aperture of 4
            "ISEA43H",  # ISEA projection with hexagon cells and a mixed sequence of aperture 4 resolutions followed by aperture 3 resolutions
            "ISEA7H",  # ISEA projection with hexagon cells and an aperture of 7
            "FULLER3H",  # FULLER projection with hexagon cells and an aperture of 3
            "FULLER4H",  # FULLER projection with hexagon cells and an aperture of 4
            "FULLER4T",  # FULLER projection with triangle cells and an aperture of 4
            "FULLER4D",  # FULLER projection with diamond cells and an aperture of 4
            "FULLER43H",  # FULLER projection with hexagon cells and a mixed sequence of aperture 4 resolutions followed by aperture 3 resolutions
            "FULLER7H",  # FULLER projection with hexagon cells and an aperture of 7
        ]

        self.addParameter(
            QgsProcessingParameterEnum(
                self.INPUT_DGGS_TYPE,
                self.tr("Input select DGGRID supported DGGS type"),
                optional=False,
                allowMultiple=False,
                defaultValue=7,
                options=self.INPUT_DGGS_TYPE_OPTIONS_LIST,
            )
        )

        self.addParameter(
            QgsProcessingParameterNumber(
                self.INPUT_DGGS_RESOLUTION,
                self.tr("Input DGGS resolution"),
                defaultValue=4,
                optional=False,
                minValue=1,
                maxValue=22,
            )
        )

        self.addParameter(
            QgsProcessingParameterString(
                self.INPUT_DGGS_MIXED_APERTURE_LEVEL,
                self.tr(
                    "for mixed aperture level (e.g. ISEA43H) provide series of aperture numbers as one string"
                ),
                optional=True,
                multiLine=False,
            )
        )

        self.addParameter(
            QgsProcessingParameterBoolean(
                self.INPUT_CLIP_EXTENT,
                self.tr("Apply layer to clip DGGS grid"),
                defaultValue=False,
                optional=True,
            )
        )

        self.addParameter(
            QgsProcessingParameterVectorLayer(
                self.INPUT_EXTENT_LAYER,
                self.tr(
                    "Input vector layer to constrain DGGS grid (must be in WGS84) and supported OGR/GDAL vector format"
                ),
                optional=True,
            )
        )

        # INPUT_DGGRID_EXECUTABLE
        self.addParameter(
            QgsProcessingParameterFile(
                self.INPUT_DGGRID_EXECUTABLE,
                self.tr("location of DGGRID executable program"),
                optional=False,
                defaultValue="C:\\Users\\Alexander\\.julia\\artifacts\\3623dbc61425e9c70e82ed04a4a00af546c1495e\\bin\\dggrid.exe",
            )
        )

        # INPUT_DGGRID_EXECUTABLE_LIBDIRS
        julia_dirs_alex = [
            "C:\\Users\\Alexander\\.julia\\artifacts\\c7ce76fefe6ce8001d5e745e8faaa47ba4228852\\bin",
            "C:\\Users\\Alexander\\.julia\\artifacts\\5ea43503c7796016c3a937d84b6bcf91b1bb608f\\bin",
            "C:\\Users\\Alexander\\.julia\\artifacts\\584ac54e52794a22f516f08a1ca655a7f49e5cae\\bin",
            "C:\\Users\\Alexander\\.julia\\artifacts\\c4c27acf31ad9b69c67a5d1d44910b603936c42a\\bin",
            "C:\\Users\\Alexander\\.julia\\artifacts\\58603237680e9ea3493b336d902e50a332fa1aef\\bin",
            "C:\\Users\\Alexander\\.julia\\artifacts\\d89c49a1ed52a2a6d7c516643590932437d7cbb0\\bin",
            "C:\\Users\\Alexander\\.julia\\artifacts\\a5355b4efc1953f6d15e22618ff5dabbbd0f84dd\\bin",
            "C:\\Users\\Alexander\\.julia\\artifacts\\f70a0c7c2917e5bde1419ffe4fba6e3f5fe78324\\bin",
            "C:\\Users\\Alexander\\.julia\\artifacts\\bdeeb133a67a52079b5206e7f9f3a7cbdf357969\\bin",
            "C:\\Users\\Alexander\\.julia\\artifacts\\00802f9a8a7cef0aa5751ca363507c3fc07de490\\bin",
            "C:\\Users\\Alexander\\.julia\\artifacts\\0b827c9d6fcb7b9f25ede8cf20fcbcab2a497200\\bin",
            "C:\\dev\\Julia_1.5.3\\bin",
            "C:\\dev\\Julia_1.5.3\\bin\\..\\lib\\julia",
            "C:\\dev\\Julia_1.5.3\\bin\\..\\lib",
            "C:\\Users\\Alexander\\.julia\\artifacts\\5c823a847569d0789a68f65aba8aa78b34d4ca3f\\bin",
            "C:\\Users\\Alexander\\.julia\\artifacts\\0407eb942f4c6c53c9bcc1eeccf3def7b46d6cc7\\bin",
            "C:\\Users\\Alexander\\.julia\\artifacts\\0a6628bc2215b9593589daa09bd637c50febfa05\\bin",
            "C:\\Users\\Alexander\\.julia\\artifacts\\42704833bfa671ef47841f64b568a87b9cb25668\\bin",
            "C:\\Users\\Alexander\\.julia\\artifacts\\3f1b12feb59c5ba75346ce47d9ccc46a614c1ca4\\bin",
            "C:\\Users\\Alexander\\.julia\\artifacts\\0f0553e3bdc0535ac962f8b46e21a88536ee2c24\\bin",
        ]

        self.addParameter(
            QgsProcessingParameterString(
                self.INPUT_DGGRID_EXECUTABLE_LIBDIRS,
                self.tr("additional temporary PATH/LD_LIBRARY_PATH config for DGGRID"),
                optional=True,
                multiLine=True,
                defaultValue=";".join(julia_dirs_alex),
            )
        )

        self.addParameter(
            QgsProcessingParameterBoolean(
                self.INPUT_MARK_DATELINE_CROSSED,
                self.tr("Mark cells that cross the dateline?"),
                defaultValue=False,
                optional=True,
            )
        )

        self.addParameter(
            QgsProcessingParameterVectorDestination(
                self.OUTPUT_GRID,
                self.tr("Output vector layer with DGGS cells"),
            )
        )

    ###########
    # helper functions
    ####################

    def dggrid_init(
        self,
        executable="dggrid",
        working_dir=None,
        capture_logs=True,
        silent=False,
        qgis_feedback=None,
    ):
        """
        initialise stuff to run dggrid smoothly
        """
        self.dggrid_executable = Path(executable).resolve()
        self.dggrid_capture_logs = capture_logs
        self.dggrid_silent = silent
        self.dggrid_last_run_succesful = False
        self.dggrid_last_run_logs = ""
        self.qgis_feedback = qgis_feedback

        if working_dir is None:
            self.dggrid_working_dir = tempfile.mkdtemp(prefix="dggrid_")
        else:
            self.dggrid_working_dir = working_dir

    def dggrid_is_runnable(self, elem):
        """
        testable if DGGRID can be executed
        """
        is_runnable = 0

        swx = Path(elem)
        if swx.exists() and swx.is_file():
            if os.access(elem, os.X_OK):
                is_runnable = 1

        return is_runnable

    def dggrid_run(self, dggs_meta_ops):
        """
        execute DGGRID with prepared metafile info
        """
        curdir = os.getcwd()
        tmp_id = uuid.uuid4()

        # subprocess.call / Popen swat_exec, check if return val is 0 or not
        # yield logs?
        try:
            os.chdir(self.dggrid_working_dir)

            with open("metafile_" + str(tmp_id), "w", encoding="utf-8") as metafile:
                for line in dggs_meta_ops:
                    metafile.write(line + "\n")

            logs = []
            o = subprocess.Popen(
                [
                    os.path.join(self.dggrid_working_dir, self.dggrid_executable),
                    "metafile_" + str(tmp_id),
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )

            while o.poll() is None:
                for b_line in o.stdout:
                    line = b_line.decode().strip()
                    # sys.stdout.write(line)
                    if not self.dggrid_silent:
                        # print(line)
                        self.qgis_feedback.pushInfo(line)
                    if self.dggrid_capture_logs:
                        logs.append(line.strip())

            if o.returncode == 0:
                self.dggrid_last_run_succesful = True
                try:
                    os.remove("metafile_" + str(tmp_id))
                except Exception:
                    pass
            else:
                self.dggrid_last_run_succesful = False

            if self.dggrid_capture_logs:
                self.dggrid_last_run_logs = "\n".join(logs)
            else:
                self.dggrid_last_run_logs = ""

        except Exception as e:
            self.dggrid_last_run_succesful = False
            print(repr(e))
            traceback.print_exc(file=sys.stdout)
            self.dggrid_last_run_logs(repr(e))
        finally:
            os.chdir(curdir)

        return o.returncode

    """
    ##############################################################################################
    # lower level API
    ##############################################################################################
    """

    """
    helper function to generate the metafile from a DGGS config
    """

    def dgselect(self, dggs_type, **kwargs):
        dggs = None

        topo_dict = {"H": "HEXAGON", "T": "TRIANGLE", "D": "DIAMOND"}

        if dggs_type in dggs_types:
            if dggs_type in ["SUPERFUND", "PLANETRISK"]:
                # keep it simple, only that spec
                dggs = Dggs(dggs_type=dggs_type, metric=True, show_info=True)

                for res_opt in [
                    "res",  # dggs_res_spec
                    "precision",  # default 7
                    "area",  # dggs_res_specify_area
                    "spacing",
                    "cls_val",  # dggs_res_specify_intercell_distance
                ]:
                    if res_opt in kwargs.keys():
                        dggs.set_par(res_opt, kwargs[res_opt])

            elif not dggs_type == "CUSTOM":

                # if dggs_type == 'ISEA3H'
                #     projection, aperture, topology = 'ISEA', 3, 'HEXAGON'

                projection, aperture, topology = "ISEA", 3, "HEXAGON"

                if dggs_type.find("ISEA") > -1:
                    projection == "ISEA"
                    sub1 = dggs_type.replace("ISEA", "")
                    topology = topo_dict[sub1[-1]]
                    aperture = int(sub1.replace(sub1[-1], ""))

                elif dggs_type.find("FULLER") > -1:
                    projection == "FULLER"
                    sub1 = dggs_type.replace("FULLER", "")
                    topology = topo_dict[sub1[-1]]
                    aperture = int(sub1.replace(sub1[-1], ""))

                else:
                    raise ValueError("projection not ISEA nor FULLER???")

                dggs = Dggs(
                    dggs_type=dggs_type,
                    projection=projection,  # dggs_projection
                    aperture=aperture,  # dggs_aperture_type / dggs_aperture
                    topology=topology,  # dggs_topology
                    metric=True,
                    show_info=True,
                )

                for res_opt in [
                    "res",  # dggs_res_spec
                    "precision",  # default 7
                    "area",  # dggs_res_specify_area
                    "spacing",
                    "cls_val",  # dggs_res_specify_intercell_distance
                ]:
                    if res_opt in kwargs.keys():
                        dggs.set_par(res_opt, kwargs[res_opt])

                if aperture == 43:
                    if "mixed_aperture_level" in kwargs.keys():
                        dggs.set_par(
                            "mixed_aperture_level", kwargs["mixed_aperture_level"]
                        )

            elif dggs_type == "CUSTOM":
                # load and insert grid definition from dggs obj

                # dggs_projections = ( "ISEA", "FULLER")
                # dggs_res_specify_types = ( "SPECIFIED", "CELL_AREA", "INTERCELL_DISTANCE" )

                # specify_resolution(proj_spec, dggs_res_spec_type)
                """
                proj_spec,
                            dggs_res_spec_type,
                            dggs_res_spec=9,
                            dggs_res_specify_area=120000,
                            dggs_res_specify_intercell_distance=4000,
                            dggs_res_specify_rnd_down=True
                """

                # dggs_topologies = ( 'HEXAGON', 'TRIANGLE', 'DIAMOND')
                # dggs_aperture_types = ( 'PURE', 'MIXED43', 'SEQUENCE')

                # specify_topo_aperture(topology_type, aperture_type, aperture_res)
                """
                specify_topo_aperture(topology_type, aperture_type, aperture_res, dggs_num_aperture_4_res=0, dggs_aperture_sequence="333333333333")
                """

                # dggs_orient_specify_types = ( 'SPECIFIED', 'RANDOM', 'REGION_CENTER' )

                if (
                    "orient_type" in kwargs.keys()
                    and kwargs["orient_type"] in dggs_orient_specify_types
                ):

                    orient_type = kwargs["orient_type"]
                    # specify_orient_type_args(orient_type)
                    """
                                                dggs_vert0_lon=11.25,
                                                dggs_vert0_lat=58.28252559,
                                                dggs_vert0_azimuth=0.0,
                                                dggs_orient_rand_seed=42
                    """

                raise ValueError("custom not yet implemented")

        # dggs.dgverify()

        return dggs

    def dg_grid_meta(self, dggs):

        dggrid_par_lookup = {
            "res": "dggs_res_spec",
            "precision": "precision",
            "area": "dggs_res_specify_area",
            "cls_val": "dggs_res_specify_intercell_distance",
            "mixed_aperture_level": "dggs_num_aperture_4_res",
        }
        metafile = []

        if dggs.dggs_type in ["SUPERFUND", "PLANETRISK"]:
            metafile.append(f"dggs_type {dggs.dggs_type}")

        elif not dggs.dggs_type == "CUSTOM":
            metafile.append(f"dggs_type {dggs.dggs_type}")

        elif dggs_type == "CUSTOM":
            raise ValueError("custom not yet implemented")

        for res_opt in [
            "res",  # dggs_res_spec
            "precision",  # default 7
            "area",  # dggs_res_specify_area
            "cls_val",  # dggs_res_specify_intercell_distance
            "mixed_aperture_level",  # dggs_num_aperture_4_res 5
        ]:
            if not dggs.get_par(res_opt, None) is None:
                opt_val = dggs.get_par(res_opt, None)
                if not opt_val is None:
                    metafile.append(f"{dggrid_par_lookup[res_opt]} {opt_val}")

        return metafile

    def dgapi_grid_gen(self, dggs, subset_conf, output_conf):
        """
        Grid Generation. Generate the cells of a DGG, either covering the complete surface of the earth or covering only a
        specific set of regions on the earth’s surface.
        """
        dggrid_operation = "GENERATE_GRID"
        metafile = []
        metafile.append("dggrid_operation " + dggrid_operation)

        dggs_config_meta = self.dg_grid_meta(dggs)

        for cmd in dggs_config_meta:
            metafile.append(cmd)

        # clip_subset_types
        if subset_conf["clip_subset_type"] == "WHOLE_EARTH":
            metafile.append("clip_subset_type " + subset_conf["clip_subset_type"])
        elif (
            subset_conf["clip_subset_type"] in ["SHAPEFILE", "AIGEN", "GDAL"]
            and not subset_conf["clip_region_files"] is None
        ):
            metafile.append("clip_subset_type " + subset_conf["clip_subset_type"])
            metafile.append("clip_region_files " + subset_conf["clip_region_files"])
        elif (
            subset_conf["clip_subset_type"] in ["SEQNUMS"]
            and not subset_conf["clip_region_files"] is None
        ):
            if not dggs.dggs_type in ["ISEA7H", "FULLER7H", "PLANETRISK"]:
                metafile.append("clip_subset_type " + subset_conf["clip_subset_type"])
                metafile.append("clip_region_files " + subset_conf["clip_region_files"])
            else:
                # if dggs_aperture_type would be SEQUENCE
                # have to reset to WHOLE_EARTH and clip based on
                # output_first_seqnum and output_last_seqnum
                subset_conf["clip_subset_type"] = "WHOLE_EARTH"
                metafile.append("clip_subset_type " + subset_conf["clip_subset_type"])
                # loading seqnums
                files = subset_conf["clip_region_files"].split(" ")
                seqnums = []
                for file in files:
                    seqs = pd.read_csv(file, header=None)[0].values
                    seqnums.append(seqs.min())
                    seqnums.append(seqs.max())

                first_seqnum = np.array(seqnums).min()
                last_seqnum = np.array(seqnums).max()
                subset_conf["output_first_seqnum"] = first_seqnum
                metafile.append(
                    "output_first_seqnum " + str(subset_conf["output_first_seqnum"])
                )
                subset_conf["output_last_seqnum"] = last_seqnum
                metafile.append(
                    "output_last_seqnum " + str(subset_conf["output_last_seqnum"])
                )
        else:
            raise ValueError("something is not correct in subset_conf")

        # join grid gen params add to metafile

        # cell_output_types
        if "cell_output_type" in output_conf.keys():
            if (
                output_conf["cell_output_type"]
                in ["SHAPEFILE", "AIGEN", "GEOJSON", "TEXT"]
                and not output_conf["cell_output_file_name"] is None
            ):
                for elem in filter(
                    lambda x: x.startswith("cell_output_"), output_conf.keys()
                ):
                    metafile.append(f"{elem} " + output_conf[elem])
            elif (
                output_conf["cell_output_type"] in ["GDAL"]
                and not output_conf["cell_output_gdal_format"] is None
                and not output_conf["cell_output_file_name"] is None
            ):
                for elem in filter(
                    lambda x: x.startswith("cell_output_"), output_conf.keys()
                ):
                    metafile.append(f"{elem} " + output_conf[elem])
            elif output_conf["cell_output_type"] in ["NONE"]:
                metafile.append("cell_output_type NONE")

            # check join cell grid params add to metafile

        # point_output_types
        if "point_output_type" in output_conf.keys():
            if (
                output_conf["point_output_type"]
                in ["SHAPEFILE", "AIGEN", "GEOJSON", "TEXT"]
                and not output_conf["point_output_file_name"] is None
            ):
                for elem in filter(
                    lambda x: x.startswith("point_output_"), output_conf.keys()
                ):
                    metafile.append(f"{elem} " + output_conf[elem])
            elif (
                output_conf["point_output_type"] in ["GDAL"]
                and not output_conf["point_output_gdal_format"] is None
                and not output_conf["point_output_file_name"] is None
            ):
                for elem in filter(
                    lambda x: x.startswith("point_output_"), output_conf.keys()
                ):
                    metafile.append(f"{elem} " + output_conf[elem])
            elif output_conf["point_output_type"] in ["NONE"]:
                metafile.append("point_output_type NONE")

            # check join point grid params add to metafile

        result = self.dggrid_run(metafile)

        if not result == 0:
            if self.dggrid_capture_logs == True:
                message = (
                    f"some error happened under the hood of dggrid (exit code {result}): "
                    + str(self.dggrid_capture_logs)
                )
                raise ValueError(message)
            else:
                message = f"some error happened under the hood of dggrid (exit code {result}), try capture_logs=True for dggrid instance"
                raise ValueError(message)

        return {"metafile": metafile, "output_conf": output_conf}

    def dgapi_grid_stats(self, dggs):
        """
        Output Grid Statistics. Output a table of grid characteristics for the specified DGG.
        """
        import copy

        np_table_switch = True
        try:
            import numpy as np
        except ImportError:
            np_table_switch = False

        dggrid_operation = "OUTPUT_STATS"
        metafile = []
        metafile.append("dggrid_operation " + dggrid_operation)

        dggs_config_meta = self.dg_grid_meta(dggs)

        for cmd in dggs_config_meta:
            metafile.append(cmd)

        # we need to capturethe logs for this one:
        save_state = copy.copy(self.dggrid_capture_logs)
        self.dggrid_capture_logs = True

        result = self.run(metafile)

        if not result == 0:
            if self.dggrid_capture_logs == True:
                message = (
                    f"some error happened under the hood of dggrid (exit code {result}): "
                    + self.last_run_logs
                )
                raise ValueError(message)
            else:
                message = f"some error happened under the hood of dggrid (exit code {result}), try capture_logs=True for dggrid instance"
                raise ValueError(message)

        table = []
        earth_line_switch = False
        earth_radius_info = ""
        for line in self.dggrid_last_run_logs.split("\n"):
            if "Earth Radius" in line:
                earth_line_switch = True
                earth_radius_info = line.strip().replace(",", "")

            if earth_line_switch == True:
                table.append(line.strip().replace(",", ""))

        # set capture logs back to original
        self.dggrid_capture_logs == save_state

        if np_table_switch == True:
            np_table = np.genfromtxt(table, skip_header=3)

            return {
                "metafile": metafile,
                "output_conf": {
                    "stats_output": np_table,
                    "earth_radius_info": earth_radius_info,
                },
            }
        else:
            return {
                "metafile": metafile,
                "output_conf": {
                    "stats_output": table,
                    "earth_radius_info": earth_radius_info,
                },
            }

    def driver_from_extension(self, path):
        """
        Attempt to auto-detect driver based on the extension.
        Parameters
        ----------
        path: str or pathlike object
            The path to the dataset to write with.
        Returns
        -------
        str:
            The name of the driver for the extension.
        """
        try:
            # in case the path is a file handle
            # or a partsed path
            path = path.name
        except AttributeError:
            pass

        # basic list for GDAL/OGR
        try:
            return gdal_ogr_vector_extension_to_driver[
                os.path.splitext(path)[-1].lstrip(".").lower()
            ]
        except KeyError:
            raise ValueError("Unable to detect driver.")

    def check_crossing(self, lon1: float, lon2: float, validate: bool = True):
        """
        Assuming a minimum travel distance between two provided longitude coordinates,
        checks if the 180th meridian (antimeridian) is crossed.
        """
        if validate and any(abs(x) > 180.0 for x in [lon1, lon2]):
            raise ValueError("longitudes must be in degrees [-180.0, 180.0]")
        return abs(lon2 - lon1) > 180.0

    def mark_geom_crossed(self, geom):
        crossed = False
        t_geom = geom.asPolygon()

        init_x = t_geom[0][0][0]

        for px in range(1, len(t_geom[0])):
            tx = t_geom[0][px][0]

            if self.check_crossing(init_x, tx):
                crossed = True

        return crossed

    def processAlgorithm(self, parameters, context, feedback):
        """
        Here is where the processing itself takes place.
        """
        """
        Here is where the processing itself takes place.
        """

        input_dggs_resolution = self.parameterAsInt(
            parameters, self.INPUT_DGGS_RESOLUTION, context
        )

        input_dggs_type = self.parameterAsInt(parameters, self.INPUT_DGGS_TYPE, context)

        input_mixed_aperture_level = self.parameterAsString(
            parameters, self.INPUT_DGGS_MIXED_APERTURE_LEVEL, context
        )

        input_clip_extent_switch = self.parameterAsBoolean(
            parameters, self.INPUT_CLIP_EXTENT, context
        )

        input_mark_datecross_switch = self.parameterAsBoolean(
            parameters, self.INPUT_MARK_DATELINE_CROSSED, context
        )

        input_clip_layer = self.parameterAsVectorLayer(
            parameters, self.INPUT_EXTENT_LAYER, context
        )

        # if we want a clip layer but the the layer is invali, we raise error
        if input_clip_extent_switch and not input_clip_layer.isValid():
            raise QgsProcessingException(
                self.invalidVectorError(parameters, self.INPUT_EXTENT_LAYER)
            )

        dggrid_libdirs = self.parameterAsString(
            parameters, self.INPUT_DGGRID_EXECUTABLE_LIBDIRS, context
        )

        path_cache = ""

        # 'Linux', 'Darwin', 'Java', 'Windows'
        if platform.system() == "Windows":
            sep = ";"
            path_cache = os.environ["PATH"]
            os.environ["PATH"] = os.environ["PATH"] + sep + dggrid_libdirs
            feedback.pushInfo(f"adding '{dggrid_libdirs}' to {platform.system()} PATH")
        elif platform.system() == "Linux":
            sep = ":"
            path_cache = os.environ["LD_LIBRARY_PATH"]
            os.environ["LD_LIBRARY_PATH"] = (
                dggrid_libdirs + sep + os.environ["LD_LIBRARY_PATH"]
            )
            feedback.pushInfo(
                f"adding '{dggrid_libdirs}' to {platform.system()} LD_LIBRARY_PATH"
            )
        elif platform.system() == "Darwin":
            sep = ":"
            path_cache = os.environ["LD_LIBRARY_PATH"]
            os.environ["LD_LIBRARY_PATH"] = (
                dggrid_libdirs + sep + os.environ["LD_LIBRARY_PATH"]
            )
            feedback.pushInfo(
                f"adding '{dggrid_libdirs}' to {platform.system()} LD_LIBRARY_PATH"
            )
        else:
            feedback.pushInfo(
                f"{platform.system()} not recognized, cannot adjust additional paths"
            )

        input_dggrid_executable = self.parameterAsString(
            parameters, self.INPUT_DGGRID_EXECUTABLE, context
        )

        dest_fname = self.parameterAsOutputLayer(parameters, self.OUTPUT_GRID, context)

        self.OUTPUT_GRID

        outpath = Path(dest_fname)
        out_drivertype = self.driver_from_extension(outpath)

        if self.dggrid_is_runnable(input_dggrid_executable):
            self.dggrid_init(
                executable=input_dggrid_executable,
                working_dir=outpath.parent,
                capture_logs=True,
                silent=False,
                qgis_feedback=feedback,
            )
        else:
            raise ValueError("can't seem to find/execute dggrid executable")

        # Send some information to the user
        feedback.pushInfo(f"aiming to write grid to {dest_fname} / {out_drivertype}")

        # Update the progress bar
        feedback.setProgress(5)

        ## set everything up and run DGGRID
        mixed_aperture_level = input_mixed_aperture_level

        if (
            not input_mixed_aperture_level is None
            and len(input_mixed_aperture_level) <= 0
        ):
            mixed_aperture_level = None
        elif isinstance(input_mixed_aperture_level, str):
            mixed_aperture_level = input_mixed_aperture_level
        else:
            mixed_aperture_level = None

        dggs_select_option = self.INPUT_DGGS_TYPE_OPTIONS_LIST[int(input_dggs_type)]

        feedback.pushInfo(
            f"selected grid options: dggs_type = {dggs_select_option}, res= {input_dggs_resolution}, mixed_aperture_level={mixed_aperture_level}"
        )

        dggs = self.dgselect(
            dggs_type=dggs_select_option,
            res=input_dggs_resolution,
            mixed_aperture_level=mixed_aperture_level,
        )

        feedback.pushInfo(f"construction info {dggs}")

        subset_conf = {"update_frequency": 100000, "clip_subset_type": "WHOLE_EARTH"}

        if input_clip_extent_switch:

            provider = input_clip_layer.dataProvider()

            subset_conf.update(
                {
                    "clip_subset_type": "GDAL",
                    "clip_region_files": str(Path(input_clip_layer.source()).resolve()),
                }
            )

        # GDAL: "GeoJSON" , "ESRI Shapefile", "GPKG"
        # or just "SHAPEFILE"

        output_conf = {
            "cell_output_type": "GDAL",
            "cell_output_gdal_format": out_drivertype,
            "cell_output_file_name": str(Path(dest_fname).resolve()),
        }

        dggs_ops = self.dgapi_grid_gen(dggs, subset_conf, output_conf)

        # resetting env
        if platform.system() == "Windows":
            os.environ["PATH"] = path_cache
        elif platform.system() == "Linux":
            os.environ["LD_LIBRARY_PATH"] = path_cache
        elif platform.system() == "Darwin":
            os.environ["LD_LIBRARY_PATH"] = path_cache
        else:
            pass

        counter = 0
        if input_mark_datecross_switch:
            feedback.pushInfo(f"calculating dateline crossings")
            dest_layer = QgsVectorLayer(dest_fname, "ogr")
            dataProvider = dest_layer.dataProvider()
            caps_string = dataProvider.capabilitiesString()
            feedback.pushInfo(f"{caps_string}")
            dest_layer.startEditing()
            crossed_field = QgsField("dt_cross", QVariant.Int)
            res = dataProvider.addAttributes([crossed_field])
            # dest_layer.addAttribute(crossed_field)
            dest_layer.updateFields()
            # dest_layer.commitChanges()
            field_id = dest_layer.fields().indexFromName("dt_cross")
            features = dataProvider.getFeatures()
            for feat in features:
                crossed = self.mark_geom_crossed(feat.geometry())
                attVal = 1 if crossed else 0
                id = feat.id()
                attr_value = {field_id: attVal}
                # feat.setAttribute(field_id, attVal)
                dataProvider.changeAttributeValues({id: attr_value})
                counter = counter + 1

            # write results
            dest_layer.commitChanges()

        feedback.pushInfo(f"processed {counter} feature cells")

        # Return the results of the algorithm, all be included in the returned
        # dictionary, with keys matching the feature corresponding parameter
        # or output names.
        return {self.OUTPUT_GRID: dest_fname}
