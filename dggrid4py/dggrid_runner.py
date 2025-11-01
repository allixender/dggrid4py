# -*- coding: utf-8 -*-
#
# Copyright (c) 2022-2025 - Alexander Kmoch
# Licenced under GNU AFFERO GENERAL PUBLIC LICENSE. Please consult the LICENCE
# file for details.
#
# Authors:
# Alexander Kmoch (alexander.kmoch@ut.ee)
# Wai Tik Chan (wai.tik.chan@ut.ee)
#
import abc
import copy
import dataclasses
import decimal
import uuid
import shutil
import os
import sys
import subprocess
import traceback
import tempfile
import warnings
from pathlib import Path
from typing import Iterable, Literal, Sequence, TypedDict, cast, get_args, TYPE_CHECKING

import numpy as np
import pandas as pd
import fiona

import geopandas as gpd
from .interrupt import crosses_interruption, interrupt_cell, get_geom_coords

if TYPE_CHECKING:
    # place here to avoid imposing 'shapely' if not installed, but supported interface
    from geopandas.base import GeometryArray
    from shapely.geometry.base import BaseGeometry  # noqa

    AnyGeometry = BaseGeometry | GeometryArray


fiona_drivers = fiona.supported_drivers

# SHAPEFILE is included more natively in DGGRID than GeoJSON (esp. as clip_subset_type)
def get_geo_out(legacy=True, has_gdal=True):
    if legacy is True and has_gdal is False:
        return { "driver": "ESRI Shapefile", "legacy_driver": "Shapefile", "ext": "shp"}

    # TODO in future would be great to have GeoArrow without GDAL
    if legacy is False and has_gdal is False:
        return { "driver": "ESRI Shapefile", "legacy_driver": "Shapefile", "ext": "shp"}

    if legacy is True and has_gdal is True:
        return { "driver": "ESRI Shapefile", "legacy_driver": "Shapefile", "ext": "shp"}

    if legacy is False and has_gdal is True:
        if "FlatGeobuf" in fiona_drivers.keys() and "w" in fiona_drivers["FlatGeobuf"]:
            return { "driver": "FlatGeobuf", "legacy_driver": "FlatGeobuf", "ext": "fgb"}

        return { "driver": "GPKG", "legacy_driver": "GPKG", "ext": "gpgk"}

    # Failsafe
    return { "driver": "ESRI Shapefile", "legacy_driver": "Shapefile", "ext": "shp"}



# specify a ISEA3H
DggsTypeT = Literal[
    'CUSTOM',  # parameters will be specified manually
    'SUPERFUND', # Superfund_500m grid
    'PLANETRISK',
    'ISEA3H', # ISEA projection with hexagon cells and an aperture of 3
    'ISEA4H', # ISEA projection with hexagon cells and an aperture of 4
    'ISEA4T',  # ISEA projection with triangle cells and an aperture of 4
    'ISEA4D', # ISEA projection with diamond cells and an aperture of 4
    'ISEA43H', # ISEA projection with hexagon cells and a mixed sequence of aperture 4 resolutions followed by aperture 3 resolutions
    'ISEA7H', # ISEA projection with hexagon cells and an aperture of 7
    'IGEO7', # ISEA projection with hexagon cells and an aperture of 7 and Z7 address type
    'FULLER3H', # FULLER projection with hexagon cells and an aperture of 3
    'FULLER4H', # FULLER projection with hexagon cells and an aperture of 4
    'FULLER4T', # FULLER projection with triangle cells and an aperture of 4
    'FULLER4D', # FULLER projection with diamond cells and an aperture of 4
    'FULLER43H', # FULLER projection with hexagon cells and a mixed sequence of aperture 4 resolutions followed by aperture 3 resolutions
    'FULLER7H', # FULLER projection with hexagon cells and an aperture of 7
]
dggs_types = get_args(DggsTypeT)

# control grid generation
DggsClipSubsetTypeT = Literal[
    'SHAPEFILE',
    'WHOLE_EARTH',
    'GDAL',
    'AIGEN',
    'SEQNUMS',
    'COARSE_CELLS',
    'INPUT_ADDRESS_TYPE'
]
clip_subset_types = get_args(DggsClipSubsetTypeT)

# specify the output
DggsFeatureOutputTypeT = Literal[
    'AIGEN',
    'GDAL',
    'GDAL_COLLECTION',
    'GEOJSON',
    'KML',
    'SHAPEFILE',
    'NONE',
    'TEXT'
]
cell_output_types = get_args(DggsFeatureOutputTypeT)
point_output_types = get_args(DggsFeatureOutputTypeT)

DggsZoneRelationOutputTypeT = Literal[
    'GDAL_COLLECTION',
    'NONE',
    'TEXT'
]
neighbor_output_types = get_args(DggsZoneRelationOutputTypeT)
children_output_types = get_args(DggsZoneRelationOutputTypeT)

DggsOutputAddressTypeV7T = Literal[
    'GEO', # geodetic coordinates -123.36 43.22 20300 Roseburg
    'Q2DI', # quad number and (i, j) coordinates on that quad
    'SEQNUM', # DGGS index - linear address (1 to size-of-DGG), not supported for parameter input_address_type if dggs_aperture_type is SEQUENCE
    'INTERLEAVE', # digit-interleaved form of Q2DI, only supported for parameter output_address_type; only available for hexagonal aperture 3 and 4 grids
    'PLANE', # (x, y) coordinates on unfolded ISEA plane,  only supported for parameter output_address_type;
    'Q2DD', # quad number and (x, y) coordinates on that quad
    'PROJTRI', # PROJTRI - triangle number and (x, y) coordinates within that triangle on the ISEA plane
    'VERTEX2DD', # vertex number, triangle number, and (x, y) coordinates on ISEA plane
    'AIGEN',  # Arc/Info Generate file format
    'Z3', # hexadecimal characters index system Z3 especially useful for ISEA3H
    'Z3_STRING',  # numerical digits representation of Z3 (as characters, not an integer)
    'Z7', # hexadecimal characters index system Z7 especially useful for ISEA7H, also in preset IGEO7
    'Z7_STRING', # numerical digits representation of Z7 (as characters, not an integer)
    'ZORDER', # index system ZORDER especially useful for ISEA3H, ISEA4H and mixed aperture
    'ZORDER_STRING'  # numerical digits representation of ZORDER (as characters, not an integer)
]
output_address_types_v7 = get_args(DggsOutputAddressTypeV7T)


DggsOutputAddressTypeV8T = Literal[
    'GEO', # geodetic coordinates -123.36 43.22 20300 Roseburg
    'Q2DI', # quad number and (i, j) coordinates on that quad
    'SEQNUM', # DGGS index - linear address (1 to size-of-DGG), not supported for parameter input_address_type if dggs_aperture_type is SEQUENCE
    'INTERLEAVE', # digit-interleaved form of Q2DI, only supported for parameter output_address_type; only available for hexagonal aperture 3 and 4 grids
    'PLANE', # (x, y) coordinates on unfolded ISEA plane,  only supported for parameter output_address_type;
    'Q2DD', # quad number and (x, y) coordinates on that quad
    'PROJTRI', # PROJTRI - triangle number and (x, y) coordinates within that triangle on the ISEA plane
    'VERTEX2DD', # vertex number, triangle number, and (x, y) coordinates on ISEA plane
    'AIGEN',  # Arc/Info Generate file format
    'HIERNDX',
]
output_address_types_v8 = get_args(DggsOutputAddressTypeV8T)

DggsOutputCellLabelTypeV7T = Literal[
    'GLOBAL_SEQUENCE',
    'ENUMERATION',
    'SUPERFUND'
]
output_cell_label_type_v7 = get_args(DggsOutputCellLabelTypeV7T)

DggsOutputCellLabelTypeV8T = Literal[
    'GLOBAL_SEQUENCE',
    'ENUMERATION',
    'SUPERFUND',
    'OUTPUT_ADDRESS_TYPE',
]
output_cell_label_type_v8 = get_args(DggsOutputCellLabelTypeV8T)

output_extra_fields_v7 = {
    'output_cell_label_type': output_cell_label_type_v7
}

DggsOutputHierNdxSystemT = Literal['Z3', 'Z7', 'ZORDER']
output_hier_ndx_systems = get_args(DggsOutputHierNdxSystemT)

DggsOutputHierNdxFormT = Literal['INT64', 'DIGIT_STRING']
output_hier_ndx_forms = get_args(DggsOutputHierNdxFormT)

output_extra_fields_v8 = {
    'output_cell_label_type': output_cell_label_type_v8,
    'output_hier_ndx_system': output_hier_ndx_systems,
    'output_hier_ndx_form': output_hier_ndx_forms,
}

DggsOutputAddressTypeT = DggsOutputAddressTypeV7T | DggsOutputAddressTypeV8T
DggsOutputCellLabelTypeT = DggsOutputCellLabelTypeV7T | DggsOutputCellLabelTypeV8T

DggsCellOutputControlT = Literal["OUTPUT_ALL", "OUTPUT_OCCUPIED"]
cell_output_controls = get_args(DggsCellOutputControlT)

DggsInputAddressTypeV7T = Literal[
    'GEO', # geodetic coordinates -123.36 43.22 20300 Roseburg
    'Q2DI', # quad number and (i, j) coordinates on that quad
    'SEQNUM', # DGGS index - linear address (1 to size-of-DGG), not supported for parameter input_address_type if dggs_aperture_type is SEQUENCE
    'Q2DD', # quad number and (x, y) coordinates on that quad
    'PROJTRI', # PROJTRI - triangle number and (x, y) coordinates within that triangle on the ISEA plane
    'VERTEX2DD', # vertex number, triangle number, and (x, y) coordinates on ISEA plane
    'AIGEN',  # Arc/Info Generate file format
    'Z3', # hexadecimal characters index system Z3 especially useful for ISEA3H
    'Z3_STRING',  # numerical digits representation of Z3 (as characters, not an integer)
    'Z7', # hexadecimal characters index system Z7 especially useful for ISEA7H, also in preset IGEO7
    'Z7_STRING', # numerical digits representation of Z7 (as characters, not an integer)
    'ZORDER', # index system ZORDER especially useful for ISEA3H, ISEA4H and mixed aperture
    'ZORDER_STRING'  # numerical digits representation of ZORDER (as characters, not an integer)
]
input_address_types_v7 = get_args(DggsInputAddressTypeV7T)

DggsInputAddressTypeV8T = Literal[
    'GEO', # geodetic coordinates -123.36 43.22 20300 Roseburg
    'Q2DI', # quad number and (i, j) coordinates on that quad
    'SEQNUM', # DGGS index - linear address (1 to size-of-DGG), not supported for parameter input_address_type if dggs_aperture_type is SEQUENCE
    'INTERLEAVE', # digit-interleaved form of Q2DI, only supported for parameter output_address_type; only available for hexagonal aperture 3 and 4 grids
    'PLANE', # (x, y) coordinates on unfolded ISEA plane,  only supported for parameter output_address_type;
    'Q2DD', # quad number and (x, y) coordinates on that quad
    'PROJTRI', # PROJTRI - triangle number and (x, y) coordinates within that triangle on the ISEA plane
    'VERTEX2DD', # vertex number, triangle number, and (x, y) coordinates on ISEA plane
    'AIGEN',  # Arc/Info Generate file format
    'HIERNDX',
]
input_address_types_v8 = get_args(DggsInputAddressTypeV8T)

DggsInputAddressTypeT = DggsInputAddressTypeV7T | DggsInputAddressTypeV8T

input_extra_fields_v7 = {}

DggsInputHierNdxSystemT = Literal['Z3', 'Z7', 'ZORDER']
input_hier_ndx_systems = get_args(DggsInputHierNdxSystemT)

DggsInputHierNdxFormT = Literal['INT64', 'DIGIT_STRING']
input_hier_ndx_forms = get_args(DggsInputHierNdxFormT)

input_extra_fields_v8 = {
    'input_hier_ndx_system': input_hier_ndx_systems,
    'input_hier_ndx_form': input_hier_ndx_forms,
}

### CUSTOM args
DggsProjectionT = Literal["ISEA", "FULLER"]
dggs_projections = get_args(DggsProjectionT)

DggsTopologyTypeT = Literal['HEXAGON', 'TRIANGLE', 'DIAMOND']
dggs_topologies = get_args(DggsTopologyTypeT)

DggsApertureT = Literal[3, 4, 7]
dggs_apertures = get_args(DggsApertureT)

DggsApertureTypeT = Literal['PURE', 'MIXED43', 'SEQUENCE']
dggs_aperture_types = get_args(DggsApertureTypeT)

DggsResSpecifyTypeT = Literal["SPECIFIED", "CELL_AREA", "INTERCELL_DISTANCE"]
dggs_res_specify_types = get_args(DggsResSpecifyTypeT)

DggsOrientSpecifyTypesT = Literal["SPECIFIED", "RANDOM", "REGION_CENTER"]
dggs_orient_specify_types = get_args(DggsOrientSpecifyTypesT)

DggsBinCoverageT = Literal["GLOBAL", "PARTIAL"]
bin_coverages = get_args(DggsBinCoverageT)

# specify the operation
DggridOperationsT = Literal[
    'GENERATE_GRID',
    'TRANSFORM_POINTS',
    'BIN_POINT_VALS',
    'BIN_POINT_PRESENCE',
    'OUTPUT_STATS'
]
dggrid_operations = get_args(DggridOperationsT)

DggsDegree = decimal.Decimal | float | str   # allow string to hardcode strictly

DggridMetaConfigParameterT = str | float | int
DggridMetaConfigT = TypedDict(
    "DggridMetaConfigT",
    {
        "dggrid_operation": DggridOperationsT,
        "dggs_type": DggsTypeT,
        "dggs_aperture": DggsApertureT,
        "dggs_aperture_type": DggsApertureTypeT,
        "dggs_aperture_sequence": str,
        "dggs_num_aperture_4_res": int,
        "dggs_proj": DggsProjectionT,
        "dggs_topology": DggsTopologyTypeT,
        "dggs_res_spec": int,
        "dggs_res_specify_area": float,
        "dggs_res_specify_type": DggsResSpecifyTypeT,
        "dggs_res_specify_rnd_down": bool,
        "dggs_res_specify_intercell_distance": float,
        "dggs_orient_specify_type": DggsOrientSpecifyTypesT,
        "dggs_orient_rand_seed": int,
        "dggs_vert0_lon": DggsDegree,
        "dggs_vert0_lat": DggsDegree,
        "dggs_vert0_azimuth": DggsDegree,
        "geodetic_densify": float,
        "densification": int,
        "precision": int,
        "bin_coverage": DggsBinCoverageT,
        "clip_subset_type": DggsClipSubsetTypeT,
        "clip_cell_res": int,
        "clip_cell_addresses": str,
        "clip_region_files": str,
        "clip_cell_densification": int,
        "clipper_scale_factor": float,
        "input_address_type": DggsInputAddressTypeT,
        "input_hier_ndx_form": DggsInputHierNdxFormT,
        "input_hier_ndx_system": DggsInputHierNdxSystemT,
        "input_delimiter": str,
        "input_file_name": str,
        "input_files": str,
        "kml_default_color": str,
        "kml_default_width": str,
        "kml_description": str,
        "kml_name": str,
        "cell_output_type": DggsFeatureOutputTypeT,
        "cell_output_gdal_format": str,
        "cell_output_file_name": str,
        "cell_output_control": DggsCellOutputControlT,
        "collection_output_gdal_format": str,
        "collection_output_file_name": str,
        "children_output_type": DggsZoneRelationOutputTypeT,
        "children_output_file_name": str,
        "neighbor_output_type": DggsZoneRelationOutputTypeT,
        "neighbor_output_file_name": str,
        "point_output_type": DggsFeatureOutputTypeT,
        "point_output_gdal_format": str,
        "point_output_file_name": str,
        "output_address_type": DggsOutputAddressTypeT,
        "output_hier_ndx_system": DggsOutputHierNdxSystemT,
        "output_hier_ndx_form": DggsOutputHierNdxFormT,
        "output_cell_label_type": DggsOutputCellLabelTypeT,
        "output_count": bool,
        "output_count_field_name": str,
        "output_num_classes": bool,
        "output_delimiter": str,
        "output_file_name": str,
        "output_first_seqnum": int,
        "output_last_seqnum": int,
        "max_cells_per_output_file": int,
        "shapefile_id_field_length": int,
        "update_frequency": int,
        "verbosity": int,
        # unofficial, only to document (since commented in metafile)
        "#short_name": str,
        "#dggs_aperture": str,
    },
    total=False,
)

DggridMetafileT = list[str]  # concatenated 'parameter value' entries
DggridApiOutputT = TypedDict(
    "DggridApiOutputT",
    {
        "metafile": DggridMetafileT,
        "output_conf": DggridMetaConfigT | dict[str, DggridMetaConfigParameterT],
    },
    total=True,
)


def __getattr__(name):  # for backward compatibility helper/warning
    if name == 'output_address_types':
        warnings.warn(
            "Definitions 'output_address_types' has been deprecated. "
            "Use the 'DGGRID.output_address_types' of desired version!",
            DeprecationWarning
        )
        return output_address_types_v7
    if name == 'input_address_types':
        warnings.warn(
            "Definitions 'input_address_types' has been deprecated. "
            "Use the 'DGGRID.input_address_types' of desired version!",
            DeprecationWarning
        )
        return output_address_types_v7
    raise AttributeError(f"module {__name__} has no attribute {name}")


def dgselect(dggs_type: DggsTypeT, **kwargs) -> "Dggs":
    """
    helper function to create a DGGS config quickly
    """
    dggs = None

    topo_dict = {
        'H' : 'HEXAGON',
        'T' : 'TRIANGLE',
        'D' : 'DIAMOND'
    }

    if dggs_type in dggs_types:
        if dggs_type in ['SUPERFUND', 'PLANETRISK', 'IGEO7']:
            # keep it simple, only that spec
            dggs = Dggs(dggs_type=dggs_type,
                        metric = True,
                        show_info = True)

            for res_opt in [ 'res', # dggs_res_spec
                             'precision', # default 7
                             'area', # dggs_res_specify_area
                             'spacing' ,
                             'cls_val' # dggs_res_specify_intercell_distance
                        ]:
                if res_opt in kwargs.keys():
                    dggs.set_par(res_opt, kwargs[res_opt])

            if dggs_type == 'IGEO7':
                dggs.set_par('aperture', 7)

        elif not dggs_type == 'CUSTOM':

            # if dggs_type == 'ISEA3H'
            #     projection, aperture, topology = 'ISEA', 3, 'HEXAGON'

            projection, aperture, topology = 'ISEA', 3, 'HEXAGON'

            if dggs_type.find('ISEA') > -1:
                projection = 'ISEA'
                sub1 = dggs_type.replace('ISEA','')
                topology = topo_dict[sub1[-1]]
                aperture = int(sub1.replace(sub1[-1], ''))

            elif dggs_type.find('FULLER') > -1:
                projection = 'FULLER'
                sub1 = dggs_type.replace('FULLER','')
                topology = topo_dict[sub1[-1]]
                aperture = int(sub1.replace(sub1[-1], ''))

            else:
                raise ValueError(f'projection {dggs_type} not ISEA nor FULLER???')

            dggs = Dggs(dggs_type=dggs_type,
                        projection=projection,  # dggs_projection
                        aperture=aperture,  # dggs_aperture_type / dggs_aperture
                        topology=topology, # dggs_topology
                        metric = True,
                        show_info = True)

            for res_opt in [ 'res', # dggs_res_spec
                             'precision', # default 7
                             'area', # dggs_res_specify_area
                             'spacing' ,
                             'cls_val' # dggs_res_specify_intercell_distance
                        ]:
                if res_opt in kwargs.keys():
                    dggs.set_par(res_opt, kwargs[res_opt])

            if aperture == 43:
                if 'mixed_aperture_level' in kwargs.keys():
                    dggs.set_par('mixed_aperture_level', kwargs['mixed_aperture_level'])


        elif dggs_type == 'CUSTOM':
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

            # specify_topo_aperture(topology_type, aperture_type, dggs_aperture_res)
            """
            specify_topo_aperture(topology_type, aperture_type, dggs_aperture_res, dggs_num_aperture_4_res=0, dggs_aperture_sequence="333333333333")
            """

            # dggs_orient_specify_types = ( 'SPECIFIED', 'RANDOM', 'REGION_CENTER' )

            if 'orient_type' in kwargs.keys() and kwargs['orient_type'] in dggs_orient_specify_types:

                orient_type = kwargs['orient_type']
                # specify_orient_type_args(orient_type)
                """
                                            dggs_vert0_lon=11.25,
                                            dggs_vert0_lat=58.28252559,
                                            dggs_vert0_azimuth=0.0,
                                            dggs_orient_rand_seed=42
                """

            raise ValueError('custom not yet implemented')

    # dggs.dgverify()

    return dggs


def dg_grid_meta(dggs: "Dggs") -> DggridMetafileT:
    """
    helper function to generate the metafile from a DGGS config
    """
    if dggs.dggs_type == 'CUSTOM':
        raise ValueError('custom not yet implemented')
    dggs_meta = dggs.to_dict()
    dggs_meta.update(specify_resolution(**dggs_meta))
    dggs_meta.update(specify_topo_aperture(**dggs_meta))
    dggs_meta.update(specify_orient_type_args(**dggs_meta))
    metafile = [
        f"{param} {val}"
        for param, val in dggs_meta.items()
        if val is not None
    ]
    return metafile


@dataclasses.dataclass
class Dggs:
    dggs_type: DggsTypeT = 'CUSTOM'
    projection: DggsProjectionT = 'ISEA'
    aperture: DggsApertureT = 3
    topology: DggsTopologyTypeT = 'HEXAGON'
    resolution: int | None = None
    precision: int = 7
    densification: int | None = None
    geodetic_densify: float | None = None
    area: float | None = None
    spacing: float | None = None
    cls_val: float | None = None
    resround: str | None = 'nearest'
    metric: bool = True
    show_info: bool = True
    azimuth_deg: DggsDegree = 0.0
    pole_lat_deg: DggsDegree = 58.28252559
    pole_lon_deg: DggsDegree = 11.25
    mixed_aperture_level: int | None = None

    _hide_metafile = [
        'metric',
        'show_info',
        'resround',
        'spacing',
        '_aliases',
    ]

    # Aliases mapping: alias -> canonical
    _aliases: dict[str, str] = dataclasses.field(default_factory=lambda: {
        'dggs_type': 'dggs_type',
        'dggs_proj': 'projection',
        'dggs_aperture': 'aperture',
        'dggs_topology': 'topology',
        'dggs_res_spec': 'resolution',
        'dggs_res_specify_area': 'area',
        'dggs_res_specify_intercell_distance': 'cls_val',
        'dggs_vert0_azimuth': 'azimuth_deg',
        'dggs_vert0_lat': 'pole_lat_deg',
        'dggs_vert0_lon': 'pole_lon_deg',
        'dggs_num_aperture_4_res': 'mixed_aperture_level',
        'res': 'resolution',
        # Add more aliases as needed
    }, init=False, repr=False)

    def _resolve_key(self, key: str, strict: bool = False) -> str | None:
        # If key is an alias, return canonical; else if it's a canonical, return as is
        if strict:
            if key in self._aliases:
                return self._aliases[key]
            elif key in self.__dataclass_fields__:
                return key
            return None
        return self._aliases.get(key, key)

    def set_par(self, par_key: str, par_value: DggridMetaConfigParameterT):
        key = self._resolve_key(par_key)
        setattr(self, key, par_value)
        return self

    def get_par(self, par_key: str, alternative=None) -> DggridMetaConfigParameterT | None:
        key = self._resolve_key(par_key)
        return getattr(self, key, alternative)

    def to_dict(self) -> dict[str, DggridMetaConfigParameterT]:
        # use the first alias with 'dggs_' prefix if available, or use canonical name
        # Build reverse mapping: canonical -> [aliases]
        result = {}
        reverse_aliases = {}
        for alias, canonical in self._aliases.items():
            reverse_aliases.setdefault(canonical, []).append(alias)
        for field_name in self.__dataclass_fields__:
            if field_name in self._hide_metafile:
                continue
            value = getattr(self, field_name)
            alias_list = reverse_aliases.get(field_name, [])
            dggs_alias = next((a for a in alias_list if a.startswith('dggs_')), None)
            key = dggs_alias if dggs_alias else field_name
            result[key] = value
        return result

    @property
    def metafile(self) -> DggridMetafileT:
        return copy.copy(dg_grid_meta(self))

    def update(self, strict=False, **kwargs):
        """
        Back propagate keyword parameters if they can be mapped to an attribute handled by this class.

        :param strict:
            If True, only keys that match existing attributes will be set.
            If False, unknown keys will be set with provided name.
        """
        for key, value in kwargs.items():
            found = self._resolve_key(key, strict=strict)
            if found:
                self.set_par(found, value)

    def dg_closest_res_to_area (self, area, resround,metric,show_info=True):
        raise ValueError('not yet implemented')

    def dg_closest_res_to_spacing(self, spacing,resround,metric,show_info=True):
        raise ValueError('not yet implemented')

    def dg_closest_res_to_cls (self, cls_val, resround,metric,show_info=True):
        raise ValueError('not yet implemented')


class DGGRID(abc.ABC):
    """
    Abstract instance that needs to be instantiated once to tell where to use and execute the dggrid cmd tool.

    It must be instantiated by one of its specific version to validate respective fields.
    """
    def __init__(self, executable = 'dggrid',
                 working_dir = None,
                 capture_logs=True,
                 silent=False,
                 tmp_geo_out_legacy=False,
                 has_gdal=True,
                 debug=False):
        self.executable = Path(executable).resolve()
        self.capture_logs=capture_logs
        self.silent=silent
        self.last_run_successful = False
        self.last_run_logs = ''
        self.last_ops_meta = {}
        self.tmp_geo_out = get_geo_out(legacy=tmp_geo_out_legacy)
        self.has_gdal = has_gdal
        self.debug = debug

        if working_dir is None:
            self.working_dir = tempfile.mkdtemp(prefix='dggrid_')
        else:
            self.working_dir = working_dir

    @property
    @abc.abstractmethod
    def version(self) -> int:
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def output_address_types(self) -> Iterable[DggsOutputAddressTypeT]:
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def output_extra_fields(self) -> dict[str, Iterable[DggridMetaConfigParameterT]]:
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def input_address_types(self) -> Iterable[DggsInputAddressTypeT]:
        raise NotImplementedError

    @property
    @abc.abstractmethod
    def input_extra_fields(self) -> dict[str, Iterable[DggridMetaConfigParameterT]]:
        raise NotImplementedError

    def check_output_extra_fields(
        self,
        output_conf: DggridMetaConfigT | dict[str, DggridMetaConfigParameterT],
    ) -> DggridMetaConfigT:
        output_conf_extra: DggridMetaConfigT = {}
        for conf, allowed in self.output_extra_fields.items():
            if conf in output_conf:
                out = output_conf[conf]
                if out not in allowed:
                    raise ValueError(f"output_conf '{conf}' has invalid value '{out}'. Allowed are {allowed}")
                output_conf_extra[conf] = out
        return output_conf_extra

    def check_input_extra_fields(
        self,
        input_conf: DggridMetaConfigT | dict[str, DggridMetaConfigParameterT],
    ) -> DggridMetaConfigT:
        input_conf_extra: DggridMetaConfigT = {}
        for conf, allowed in self.input_extra_fields.items():
            if conf in input_conf:
                _in = input_conf[conf]
                if _in not in allowed:
                    raise ValueError(f"input_conf '{conf}' has invalid value '{_in}'. Allowed are {allowed}")
                input_conf_extra[conf] = _in
        return input_conf_extra

    def check_gdal_support(self):
        if self.has_gdal:
            print(f"GDAL types should be possible: has GDAL={self.has_gdal}")
            print("check your dggrid binary | mac: otool -L | Linux ldd ")
            print(fiona_drivers)

        else:
            print(f"GDAL types should not be used: has GDAL={self.has_gdal}")


    def is_runnable(self):
        is_runnable = 0

        takes = []
        take_1 = shutil.which(self.executable)
        if not take_1 is None:
            takes.append(take_1)

        take_2 = shutil.which(os.path.join(self.working_dir, self.executable))
        if not take_2 is None:
            takes.append(take_2)

        if len(takes) < 1:
            print(f"{self.executable} not in executable paths")
        else:
            for elem in takes:
                swx = Path(elem)
                if swx.exists() and swx.is_file():
                    if os.access(elem, os.X_OK):
                        # print(f"{elem} is executable")
                        self.executable = str(elem)
                        is_runnable = 1

        return is_runnable


    def run(self, dggs_meta_ops):

        curdir = os.getcwd()
        tmp_id = uuid.uuid4()

        # subprocess.call / Popen swat_exec, check if return val is 0 or not
        # yield logs?
        o = None
        try:
            os.chdir(self.working_dir)

            with open('metafile_' + str(tmp_id), 'w', encoding='utf-8') as metafile:
                for line in dggs_meta_ops:
                    metafile.write(line + '\n')

            self.last_ops_meta = dggs_meta_ops
            if self.debug:
                print(dggs_meta_ops)

            logs = []
            o = subprocess.Popen([os.path.join(self.working_dir, self.executable), 'metafile_' + str(tmp_id)], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

            while o.poll() is None:
                for b_line in o.stdout:
                    line = b_line.decode().strip()
                    # sys.stdout.write(line)
                    if not self.silent:
                        print(line)
                    if self.capture_logs:
                        logs.append(line.strip())

            if o.returncode == 0:
                self.last_run_successful = True
                if not self.debug:
                    try:
                        os.remove( 'metafile_' + str(tmp_id) )
                    except Exception:
                        pass
            else:
                self.last_run_successful = False

            if self.capture_logs:
                self.last_run_logs = '\n'.join(logs)
            else:
                self.last_run_logs = ''

        except Exception as e:
            self.last_run_successful = False
            print(repr(e))
            traceback.print_exc(file=sys.stdout)
            self.last_run_logs = repr(e)
        finally:
            os.chdir(curdir)

        if not o:
            return -1
        return o.returncode

    ##############################################################################################
    # lower level API
    ##############################################################################################

    def dgapi_grid_gen(
        self,
        dggs: Dggs,
        subset_conf: DggridMetaConfigT,
        output_conf: DggridMetaConfigT,  # contrary to other methods, it should be pre-resolved to avoid duplicate checks
    ) -> DggridApiOutputT:
        """
        Grid Generation. Generate the cells of a DGG, either covering the complete surface of the earth or covering only a
        specific set of regions on the earth’s surface.
        """
        dggrid_operation = 'GENERATE_GRID'
        metafile = []
        metafile.append("dggrid_operation " + dggrid_operation)

        dggs.update(**subset_conf, **output_conf)
        dggs_config_meta = dggs.metafile

        for cmd in dggs_config_meta:
            metafile.append(cmd)

        # clip_subset_types
        if subset_conf['clip_subset_type'] == 'WHOLE_EARTH':
            metafile.append("clip_subset_type " + subset_conf['clip_subset_type'])
        elif subset_conf['clip_subset_type'] in [ 'SHAPEFILE' , 'AIGEN', 'GDAL'] and not subset_conf['clip_region_files'] is None:
            metafile.append("clip_subset_type " + subset_conf['clip_subset_type'])
            metafile.append("clip_region_files " + subset_conf['clip_region_files'])
        elif subset_conf['clip_subset_type'] in [ 'SEQNUMS', 'COARSE_CELLS', 'INPUT_ADDRESS_TYPE'] and not subset_conf['clip_region_files'] is None:
            if not dggs.dggs_type in ['PLANETRISK']:
                # metafile.append("clip_subset_type " + subset_conf['clip_subset_type'])
                # metafile.append("clip_region_files " + subset_conf['clip_region_files'])
                # COARSE_CELLS needs clip_cell_res
                for elem in filter(lambda x: x.startswith('clip_') , subset_conf.keys()):
                    metafile.append(f"{elem} " + str(subset_conf[elem]))

            else:
                # if dggs_aperture_type would be SEQUENCE
                # have to reset to WHOLE_EARTH and clip based on
                # output_first_seqnum and output_last_seqnum
                # this should now almost never happen
                subset_conf['clip_subset_type'] = 'WHOLE_EARTH'
                metafile.append("clip_subset_type " + subset_conf['clip_subset_type'])
                # loading seqnums
                files = subset_conf['clip_region_files'].split(' ')
                seqnums = []
                for file in files:
                    seqs = pd.read_csv( file, header=None)[0].values
                    seqnums.append(seqs.min())
                    seqnums.append(seqs.max())

                first_seqnum = np.array(seqnums).min()
                last_seqnum = np.array(seqnums).max()
                subset_conf['output_first_seqnum'] = first_seqnum
                metafile.append("output_first_seqnum " + str(subset_conf['output_first_seqnum']))
                subset_conf['output_last_seqnum'] = last_seqnum
                metafile.append("output_last_seqnum " + str(subset_conf['output_last_seqnum']))
        else:
            raise ValueError('something is not correct in subset_conf')

        if 'input_address_type' in subset_conf.keys() and subset_conf.get('input_address_type', 'NOPE') in self.input_address_types:
            metafile.append("input_address_type " + subset_conf['input_address_type'])
            for elem in self.input_extra_fields:
                if elem in subset_conf:
                    metafile.append(f"{elem} {subset_conf[elem]}")

        # join grid gen params add to metafile

        # cell_output_types
        if 'cell_output_type' in output_conf.keys():
            if (
                output_conf['cell_output_type'] in [ 'SHAPEFILE' , 'AIGEN', 'GEOJSON', 'TEXT']
                and output_conf['cell_output_file_name'] is not None
            ):
                for elem in filter(lambda x: x.startswith('cell_output_') , output_conf.keys()):
                    metafile.append(f"{elem} " + str(output_conf[elem]))
            elif (
                output_conf['cell_output_type'] in ['GDAL', 'GDAL_COLLECTION'] and
                output_conf['cell_output_gdal_format'] is not None and
                output_conf['cell_output_file_name'] is not None
            ):
                for elem in filter(lambda x: x.startswith('cell_output_') , output_conf.keys()):
                    metafile.append(f"{elem} " + str(output_conf[elem]))
            elif output_conf['cell_output_type'] in [ 'NONE']:
                metafile.append("cell_output_type NONE")

            # check join cell grid params add to metafile

        # collection output
        if 'collection_output_file_name' in output_conf:
            metafile.append("collection_output_file_name " + str(output_conf['collection_output_file_name']))
        if 'collection_output_gdal_format' in output_conf:
            metafile.append("collection_output_gdal_format " + str(output_conf['collection_output_gdal_format']))

        # point_output_types
        if 'point_output_type' in output_conf.keys():
            if (
                output_conf['point_output_type'] in ['SHAPEFILE' , 'AIGEN', 'GEOJSON', 'TEXT'] and
                output_conf['point_output_file_name'] is not None
            ):
                for elem in filter(lambda x: x.startswith('point_output_') , output_conf.keys()):
                    metafile.append(f"{elem} " + str(output_conf[elem]))
            elif (
                output_conf['point_output_type'] in ['GDAL', 'GDAL_COLLECTION'] and
                output_conf['point_output_gdal_format'] is not None and
                output_conf['point_output_file_name'] is not None
            ):
                for elem in filter(lambda x: x.startswith('point_output_') , output_conf.keys()):
                    metafile.append(f"{elem} " + str(output_conf[elem]))
            elif output_conf['point_output_type'] in [ 'NONE']:
                metafile.append("point_output_type NONE")

        # neighbor_output_types
        if 'neighbor_output_type' in output_conf.keys():
            if (
                output_conf['neighbor_output_type'] in neighbor_output_types and
                output_conf['neighbor_output_file_name'] is not None
            ):
                for elem in filter(lambda x: x.startswith('neighbor_output_'), output_conf.keys()):
                    metafile.append(f"{elem} " + str(output_conf[elem]))
            elif output_conf['neighbor_output_type'] in ['NONE']:
                metafile.append("neighbor_output_type NONE")

        # children_output_types
        if 'children_output_type' in output_conf.keys():
            if (
                output_conf['children_output_type'] in children_output_types and
                output_conf['children_output_file_name'] is not None
            ):
                for elem in filter(lambda x: x.startswith('children_output_'), output_conf.keys()):
                    metafile.append(f"{elem} " + str(output_conf[elem]))
            elif output_conf['children_output_type'] in ['NONE']:
                metafile.append("children_output_type NONE")

        # check join point grid params add to metafile

        if 'output_address_type' in output_conf.keys() and output_conf.get('output_address_type', 'NOPE') in self.output_address_types:
            metafile.append("output_address_type " + output_conf['output_address_type'])
            for elem in self.output_extra_fields:
                if elem in output_conf:
                    metafile.append(f"{elem} {output_conf[elem]}")

        result = self.run(metafile)

        if result != 0:
            if self.capture_logs:
                message = f"some error happened under the hood of dggrid (exit code {result}): " + self.last_run_logs
                raise ValueError(message)
            else:
                message = f"some error happened under the hood of dggrid (exit code {result}), try capture_logs=True for dggrid instance"
                raise ValueError(message)

        return { 'metafile': metafile, 'output_conf': output_conf }


    def dgapi_grid_transform(
        self,
        dggs: Dggs,
        subset_conf: DggridMetaConfigT,
        **conf_extra: DggridMetaConfigParameterT,
    ) -> DggridApiOutputT:
        """
        Address Conversion. Transform a file of locations from one address form (such as longitude/latitude) to another (such as DGG cell indexes).
        """
        dggrid_operation = 'TRANSFORM_POINTS'
        metafile = []
        metafile.append("dggrid_operation " + dggrid_operation)

        dggs_config_meta = dg_grid_meta(dggs)

        for cmd in dggs_config_meta:
            metafile.append(cmd)

        # transform input_types
        if (
            'input_file_name' in subset_conf.keys() and
            'input_address_type' in subset_conf.keys() and
            subset_conf['input_address_type'] in self.input_address_types
        ):
            for elem in filter(lambda x: x.startswith('input_') , subset_conf.keys()):
                metafile.append(f"{elem} " + str(subset_conf[elem]))
        else:
            raise ValueError('no input filename or type given')

        # join grid gen params add to metafile
        subset_conf.update(specify_resolution(**conf_extra))
        subset_conf.update(specify_orient_type_args(**conf_extra))
        subset_conf.update(specify_topo_aperture(**conf_extra))

        input_extras = self.check_input_extra_fields(conf_extra)
        subset_conf.update(input_extras)
        if subset_conf:
            for elem, value in subset_conf.items():
                metafile.append(f"{elem} " + str(value))

        # transform output_types
        output_conf = {}
        if (
            'output_file_name' in conf_extra.keys() and
            'output_address_type' in conf_extra.keys() and
            conf_extra['output_address_type'] in self.output_address_types
        ):
            for elem in filter(lambda x: x.startswith('output_') , conf_extra.keys()):
                output_conf[elem] = conf_extra[elem]
                metafile.append(f"{elem} " + str(conf_extra[elem]))
        else:
            raise ValueError('no output filename or type given')

        output_extras = self.check_output_extra_fields(conf_extra)
        if output_extras:
            for elem, value in output_extras.items():
                output_conf[elem] = value
                metafile.append(f"{elem} " + str(value))

        result = self.run(metafile)

        if result != 0:
            if self.capture_logs:
                message = f"some error happened under the hood of dggrid (exit code {result}): " + self.last_run_logs
                raise ValueError(message)
            else:
                message = f"some error happened under the hood of dggrid (exit code {result}), try capture_logs=True for dggrid instance"
                raise ValueError(message)

        return { 'metafile': metafile, 'output_conf': output_conf }


    def dgapi_point_value_binning(
        self,
        dggs: Dggs,
        subset_conf: DggridMetaConfigT,
        **conf_extra: DggridMetaConfigParameterT,
    ) -> DggridApiOutputT:
        """
        Point Value Binning. Bin a set of floating-point values associated with point locations into the cells of a DGG,
        by assigning to each DGG cell the arithmetic mean of the values which are contained in that cell.

        # specify the operation
        dggrid_operation BIN_POINT_VALS

        # specify the DGG

        dggs_type ISEA3H
        dggs_res_spec 9

        # specify bin controls

        bin_coverage PARTIAL
        input_files inputfiles/20k.txt inputfiles/50k.txt inputfiles/100k.txt inputfiles/200k.txt
        input_delimiter " "

        # specify the output

        output_file_name outputfiles/popval3h9.txt
        output_address_type SEQNUM
        output_delimiter ","
        output_count TRUE
        cell_output_control OUTPUT_OCCUPIED

        """
        dggrid_operation = 'BIN_POINT_VALS'
        metafile = []
        metafile.append("dggrid_operation " + dggrid_operation)

        dggs_config_meta = dg_grid_meta(dggs)

        for cmd in dggs_config_meta:
            metafile.append(cmd)

        # transform input_types
        if (
            'input_file_name' in subset_conf.keys() and
            'input_address_type' in subset_conf.keys() and
            subset_conf['input_address_type'] in self.input_address_types
        ):
            for elem in filter(lambda x: x.startswith('input_') , subset_conf.keys()):
                metafile.append(f"{elem} " + subset_conf[elem])
        else:
            raise ValueError('no input filename or type given')

        # join grid gen params add to metafile
        subset_conf.update(specify_resolution(**conf_extra))
        subset_conf.update(specify_orient_type_args(**conf_extra))
        subset_conf.update(specify_topo_aperture(**conf_extra))

        input_extras = self.check_input_extra_fields(conf_extra)
        subset_conf.update(input_extras)
        if subset_conf:
            for elem, value in subset_conf.items():
                metafile.append(f"{elem} " + str(value))

        # transform output_types
        output_conf = {}
        if (
            'output_file_name' in conf_extra.keys() and
            'output_address_type' in conf_extra.keys() and
            conf_extra['output_address_type'] in self.output_address_types
        ):
            for elem in filter(lambda x: x.startswith('output_') , conf_extra.keys()):
                output_conf[elem] = conf_extra[elem]
                metafile.append(f"{elem} " + conf_extra[elem])
        else:
            raise ValueError('no output filename or type given')

        output_extras = self.check_output_extra_fields(conf_extra)
        if output_extras:
            for elem, value in output_extras.items():
                output_conf[elem] = conf_extra[elem]
                metafile.append(f"{elem} " + str(value))

        result = self.run(metafile)

        if result != 0:
            if self.capture_logs:
                message = f"some error happened under the hood of dggrid (exit code {result}): " + self.last_run_logs
                raise ValueError(message)
            else:
                message = f"some error happened under the hood of dggrid (exit code {result}), try capture_logs=True for dggrid instance"
                raise ValueError(message)

        return { 'metafile': metafile, 'output_conf': output_conf }


    def dgapi_pres_binning(
        self,
        dggs: Dggs,
        subset_conf: DggridMetaConfigT,
        **conf_extra: DggridMetaConfigParameterT,
    ) -> DggridApiOutputT:
        """
        Presence/Absence Binning. Given a set of input files, each containing point locations associated with a particular class, DGGRID outputs,
        for each cell of a DGG, a vector indicating whether or not each class is present in that cell.
        """
        dggrid_operation = 'BIN_POINT_PRESENCE'
        metafile = []
        metafile.append("dggrid_operation " + dggrid_operation)

        dggs_config_meta = dg_grid_meta(dggs)

        for cmd in dggs_config_meta:
            metafile.append(cmd)

        # transform input_types
        if (
            'input_file_name' in subset_conf.keys() and
            'input_address_type' in subset_conf.keys() and
            subset_conf['input_address_type'] in self.input_address_types
        ):
            for elem in filter(lambda x: x.startswith('input_') , subset_conf.keys()):
                metafile.append(f"{elem} " + subset_conf[elem])
        else:
            raise ValueError('no input filename or type given')

        # join grid gen params add to metafile
        subset_conf.update(specify_resolution(**conf_extra))
        subset_conf.update(specify_orient_type_args(**conf_extra))
        subset_conf.update(specify_topo_aperture(**conf_extra))

        input_extras = self.check_input_extra_fields(conf_extra)
        subset_conf.update(input_extras)
        if subset_conf:
            for elem, value in subset_conf.items():
                metafile.append(f"{elem} " + str(value))

        # pre-patch elements before integrating them to the metafile
        if "output_delimiter" in conf_extra:
            delim = conf_extra["output_delimiter"] or " "
            delim = f'"{delim}"' if '"' not in delim else delim
            conf_extra["output_delimiter"] = delim
        if "cell_output_control" in conf_extra and conf_extra["cell_output_control"] not in cell_output_controls:
            raise ValueError(f"cell_output_control must be one of: {cell_output_controls}")
        for field in ["output_count", "output_num_classes"]:
            if field in conf_extra:
                conf_extra[field] = str(conf_extra[field]).upper()  # bool or string-like bool

        # transform output_types
        output_conf = {}
        if (
            'output_file_name' in conf_extra.keys() and
            'output_address_type' in conf_extra.keys() and
            conf_extra['output_address_type'] in self.output_address_types
        ):
            for elem in filter(lambda x: x.startswith('output_') or x.startswith('cell_output_'), conf_extra.keys()):
                output_conf[elem] = conf_extra[elem]
                metafile.append(f"{elem} " + conf_extra[elem])
        else:
            raise ValueError('no output filename or type given')

        output_extras = self.check_output_extra_fields(conf_extra)
        if output_extras:
            for elem, value in output_extras.items():
                output_conf[elem] = conf_extra[elem]
                metafile.append(f"{elem} {value!s}")

        result = self.run(metafile)

        if result != 0:
            if self.capture_logs:
                message = f"some error happened under the hood of dggrid (exit code {result}): " + self.last_run_logs
                raise ValueError(message)
            else:
                message = f"some error happened under the hood of dggrid (exit code {result}), try capture_logs=True for dggrid instance"
                raise ValueError(message)

        return { 'metafile': metafile, 'output_conf': output_conf }


    def dgapi_grid_stats(self, dggs: Dggs) -> DggridApiOutputT:
        """
        Output Grid Statistics. Output a table of grid characteristics for the specified DGG.
        """

        dggrid_operation = 'OUTPUT_STATS'
        metafile = []
        metafile.append("dggrid_operation " + dggrid_operation)

        dggs_config_meta = dg_grid_meta(dggs)

        for cmd in dggs_config_meta:
            metafile.append(cmd)

        # we need to capture the logs for this one:
        save_state = copy.copy(self.capture_logs)
        self.capture_logs = True

        result = self.run(metafile)

        if result != 0:
            if self.capture_logs:
                message = f"some error happened under the hood of dggrid (exit code {result}): " + self.last_run_logs
                raise ValueError(message)
            else:
                message = f"some error happened under the hood of dggrid (exit code {result}), try capture_logs=True for dggrid instance"
                raise ValueError(message)

        table = []
        earth_line_switch = False
        earth_radius_info = ''
        for line in self.last_run_logs.split('\n'):
            if "Earth Radius" in line:
                earth_line_switch = True
                earth_radius_info = line.strip().replace(',','')

            if earth_line_switch:
                table.append(line.strip().replace(',',''))

        # set capture logs back to original
        self.capture_logs = save_state

        np_table = np.genfromtxt(table, skip_header=3)
        return { 'metafile': metafile, 'output_conf': {'stats_output': np_table, 'earth_radius_info': earth_radius_info } }


    def post_process_split_dateline(self, gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        cellsNew = gdf.iloc[:0].copy()
        # if we get eventually binning working we will have dynamic additional column for the binned values
        # cols = dict({ (idx, field) for idx, field in enumerate(gdf.columns)})

        name_col = 'name'
        for potential_name_col in ['name', 'Name', 'global_id']:
            if potential_name_col in gdf.columns:
                name_col = potential_name_col
                break

        for row in gdf.itertuples():
            poly = row.geometry
            coords = get_geom_coords(poly)

            if crosses_interruption(coords):
                geoms = interrupt_cell(coords)

                for geom in geoms:
                    cellsNew.loc[len(cellsNew)] = [getattr(row, name_col), geom]
            else:
                cellsNew.loc[len(cellsNew)] = [getattr(row, name_col), poly]

        return cellsNew

    #################################################################################
    # Higher level API
    #################################################################################

    def grid_stats_table(
        self,
        dggs_type: DggsTypeT,
        resolution: int,
        mixed_aperture_level: int = None,
    ) -> pd.DataFrame:
        """
        generates the area and cell statistics for the given DGGS from resolution 0 to the given resolution of the DGGS
        """
        dggs = dgselect(dggs_type = dggs_type, res= resolution, mixed_aperture_level=mixed_aperture_level)

        dggs_ops = self.dgapi_grid_stats(dggs)
        if self.debug:
            print(dggs_ops)

        df = pd.DataFrame(dggs_ops['output_conf']['stats_output'])

        df.rename(columns={0: 'Resolution', 1: "Cells", 2:"Area (km^2)", 3: "CLS (km)"}, inplace=True)
        df['Resolution'] = df['Resolution'].astype(int)
        df['Cells'] = df['Cells'].astype(np.int64)

        return df


    def grid_cell_polygons_for_extent(
        self,
        dggs_type: DggsTypeT,
        resolution: int,
        mixed_aperture_level: int = None,
        clip_geom: "AnyGeometry" = None,
        split_dateline: bool = False,
        output_address_type: DggsOutputAddressTypeT | None = None,
        **conf_extra: DggridMetaConfigParameterT,
    ) -> gpd.GeoDataFrame:
        """
        generates a DGGS grid and returns all the cells as GeoDataFrame with geometry type Polygon
            a) if clip_geom is empty/None: grid cell ids/seqnums for the WHOLE_EARTH
            b) if clip_geom is a shapely shape geometry, takes this as a clip area
        """
        tmp_id = uuid.uuid4()
        tmp_dir = self.working_dir
        dggs = dgselect(dggs_type = dggs_type, res= resolution, mixed_aperture_level=mixed_aperture_level)
        dggs.update(**conf_extra, strict=True)

        subset_conf: DggridMetaConfigT = { 'update_frequency': 100000, 'clip_subset_type': 'WHOLE_EARTH' }

        if not clip_geom is None and clip_geom.area > 0:

            clip_gdf = gpd.GeoDataFrame(pd.DataFrame({'id' : [1], 'geometry': [clip_geom]}), geometry='geometry', crs=4326)
            clip_path = Path(tmp_dir) / f"temp_clip_{tmp_id}.{self.tmp_geo_out['ext']}"
            clip_gdf.to_file(str(clip_path), driver=self.tmp_geo_out['driver'])

            subset_conf.update({
                'clip_subset_type': 'GDAL',
                'clip_region_files': str( (Path(tmp_dir) / f"temp_clip_{tmp_id}.{self.tmp_geo_out['ext']}").resolve()),
                })

            if not self.has_gdal:
                subset_conf.update({
                    'clip_subset_type': cast(DggsClipSubsetTypeT, self.tmp_geo_out['legacy_driver'].upper())
                })

        subset_conf.update(specify_resolution(**conf_extra))
        subset_conf.update(specify_orient_type_args(**conf_extra))
        subset_conf.update(specify_topo_aperture(**conf_extra))

        input_extras = self.check_input_extra_fields(conf_extra)
        subset_conf.update(input_extras)

        output_extras = self.check_output_extra_fields(conf_extra)
        output_conf: DggridMetaConfigT = {
            'cell_output_type': 'GDAL',
            'point_output_type': 'NONE',
            'cell_output_gdal_format' : self.tmp_geo_out['legacy_driver'],
            'cell_output_file_name': str( (Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.{self.tmp_geo_out['ext']}").resolve()),
            **output_extras
        }

        if not self.has_gdal:
            output_conf.update({
                'cell_output_type': cast(DggsFeatureOutputTypeT, self.tmp_geo_out['legacy_driver'].upper()),
                'cell_output_file_name': str( (Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}").resolve())
            })
            output_conf.pop('cell_output_gdal_format', None)

        if not output_address_type is None and output_address_type in self.output_address_types:
            output_conf.update({'output_address_type': output_address_type})
        else:
            if self.debug:
                print(f"ignoring unknown output_address_type: {output_address_type}")

        dggs_ops = self.dgapi_grid_gen(dggs, subset_conf, output_conf)
        if self.debug:
            print(dggs_ops)

        path = Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.{self.tmp_geo_out['ext']}"
        gdf = gpd.read_file(path.resolve(), driver=self.tmp_geo_out['driver'])

        if not self.debug:
            try:
                os.remove(path)
                if "SHAPEFILE" in self.tmp_geo_out['driver'].upper():
                    for ext in ['dbf', 'prj', 'shx', 'cpg']:
                        os.remove( str( Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.{ext}") )
            except Exception:
                pass
            try:
                os.remove( str( Path(tmp_dir) / f"temp_clip_{tmp_id}.{self.tmp_geo_out['ext']}") )
                if "SHAPEFILE" in self.tmp_geo_out['driver'].upper():
                    for ext in ['dbf', 'prj', 'shx', 'cpg']:
                        os.remove( str( Path(tmp_dir) / f"temp_clip_{tmp_id}.{ext}") )
            except Exception:
                pass
            try:
                os.remove( str( Path(tmp_dir) / f"points.gen") )
            except Exception:
                pass

        if split_dateline:
            return self.post_process_split_dateline(gdf)
        return gdf


    def grid_cell_centroids_for_extent(
        self,
        dggs_type: DggsTypeT,
        resolution: int,
        mixed_aperture_level: int = None,
        clip_geom: "AnyGeometry" = None,
        output_address_type: DggsOutputAddressTypeT | None = None,
        **conf_extra: DggridMetaConfigParameterT,
    ) -> gpd.GeoDataFrame:
        """
        generates a DGGS grid and returns all the cell's centroid as GeoDataFrame with geometry type Point
            a) if clip_geom is empty/None: grid cell ids/seqnums for the WHOLE_EARTH
            b) if clip_geom is a shapely shape geometry, takes this as a clip area
        """
        tmp_id = uuid.uuid4()
        tmp_dir = self.working_dir
        dggs = dgselect(dggs_type = dggs_type, res= resolution, mixed_aperture_level=mixed_aperture_level)
        dggs.update(**conf_extra, strict=True)

        subset_conf: DggridMetaConfigT = { 'update_frequency': 100000, 'clip_subset_type': 'WHOLE_EARTH' }

        if not clip_geom is None and clip_geom.area > 0:

            clip_gdf = gpd.GeoDataFrame(pd.DataFrame({'id' : [1], 'geometry': [clip_geom]}), geometry='geometry', crs=4326)
            clip_path = Path(tmp_dir) / f"temp_clip_{tmp_id}.{self.tmp_geo_out['ext']}"
            clip_gdf.to_file(str(clip_path), driver=self.tmp_geo_out['driver'] )

            subset_conf.update({
                'clip_subset_type': 'GDAL',
                'clip_region_files': str( (Path(tmp_dir) / f"temp_clip_{tmp_id}.{self.tmp_geo_out['ext']}").resolve()),
            })

            if not self.has_gdal:
                subset_conf.update({
                    'clip_subset_type': cast(DggsClipSubsetTypeT, self.tmp_geo_out['legacy_driver'].upper())
                })

        subset_conf.update(specify_resolution(**conf_extra))
        subset_conf.update(specify_orient_type_args(**conf_extra))
        subset_conf.update(specify_topo_aperture(**conf_extra))

        input_extras = self.check_input_extra_fields(conf_extra)
        subset_conf.update(input_extras)

        output_extras = self.check_output_extra_fields(conf_extra)
        output_conf: DggridMetaConfigT = {
            'point_output_type': 'GDAL',
            'cell_output_type': 'NONE',
            'point_output_gdal_format' : self.tmp_geo_out['legacy_driver'],
            'point_output_file_name': str( (Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.{self.tmp_geo_out['ext']}").resolve()),
            **output_extras
        }

        if not self.has_gdal:
            output_conf.update({
                'point_output_type': cast(DggsFeatureOutputTypeT, self.tmp_geo_out['legacy_driver'].upper()),
                'point_output_file_name': str( (Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}").resolve())
            })
            output_conf.pop('point_output_gdal_format', None)

        if not output_address_type is None and output_address_type in self.output_address_types:
            output_conf.update({'output_address_type': output_address_type})
        else:
            if self.debug:
                print(f"ignoring unknown output_address_type: {output_address_type}")

        dggs_ops = self.dgapi_grid_gen(dggs, subset_conf, output_conf)
        if self.debug:
            print(dggs_ops)

        gdf = gpd.read_file( Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.{self.tmp_geo_out['ext']}", driver=self.tmp_geo_out['driver'] )

        if not self.debug:
            try:
                os.remove( str( Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.{self.tmp_geo_out['ext']}") )
                if "SHAPEFILE" in self.tmp_geo_out['driver'].upper():
                    for ext in ['dbf', 'prj', 'shx', 'cpg']:
                        os.remove( str( Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.{ext}") )
            except Exception:
                pass
            try:
                os.remove( str( Path(tmp_dir) / f"temp_clip_{tmp_id}.{self.tmp_geo_out['ext']}") )
                if "SHAPEFILE" in self.tmp_geo_out['driver'].upper():
                    for ext in ['dbf', 'prj', 'shx', 'cpg']:
                        os.remove( str( Path(tmp_dir) / f"temp_clip_{tmp_id}.{ext}") )
            except Exception:
                pass
            try:
                os.remove( str( Path(tmp_dir) / f"cells.gen") )
            except Exception:
                pass

        return gdf


    def grid_cell_polygons_from_cellids(
        self,
        cell_id_list: Sequence[str],
        dggs_type: DggsTypeT,
        resolution: int,
        mixed_aperture_level: int = None,
        split_dateline: bool = False,
        clip_subset_type: DggsClipSubsetTypeT = 'WHOLE_EARTH',
        clip_cell_res: int = 1,
        input_address_type: DggsInputAddressTypeT = 'SEQNUM',
        output_address_type: DggsOutputAddressTypeT = 'SEQNUM',
        **conf_extra: DggridMetaConfigParameterT,
    ) -> gpd.GeoDataFrame:
        """
        generates a DGGS grid and returns all the cells as GeoDataFrame with geometry type Polygon
            a) if cell_id_list is empty/None: grid cells for the WHOLE_EARTH
            b) if cell_id_list is a list/numpy array, takes this list as seqnums ids (potentially also Z3, Z7, ZORDER ..?) for subsetting
        """
        tmp_id = uuid.uuid4()
        tmp_dir = self.working_dir
        dggs = dgselect(dggs_type = dggs_type, res= resolution, mixed_aperture_level=mixed_aperture_level)
        dggs.update(**conf_extra, strict=True)

        subset_conf: DggridMetaConfigT = { 'update_frequency': 100000, 'clip_subset_type': clip_subset_type }
        seq_df = None

        if not cell_id_list is None and len(cell_id_list) > 0:

            seq_df = pd.DataFrame({ input_address_type: cell_id_list})
            seq_df.to_csv( str( (Path(tmp_dir) / f"temp_clip_{tmp_id}.txt").resolve()) , header=False, index=False, columns=[input_address_type], sep=' ')

            subset_conf.update({
                'clip_subset_type': 'SEQNUMS',
                'clip_region_files': str( (Path(tmp_dir) / f"temp_clip_{tmp_id}.txt").resolve()),
            })

            # TODO, for Z3, Z7, ZORDER can potentially also be COARSE_CELLS / aka parent cells?
            # clip_subset_type should INPUT_ADDRESS_TYPE for the equivalent of SEQNUM (tp use input_address_type Z3 ...), or COARSE_CELLS as an actual parent cell type clip (also for Z3 ..)
            if (
                input_address_type in self.input_address_types
                and output_address_type in self.output_address_types
            ):
                subset_conf.update(
                    {
                        'clip_subset_type': 'INPUT_ADDRESS_TYPE',
                        'input_address_type': input_address_type
                    }
                )

                if not clip_subset_type is None and clip_subset_type in ['COARSE_CELLS']:
                    subset_conf.update(
                        {
                            'clip_subset_type': clip_subset_type,
                            'input_address_type': input_address_type,
                            'clip_cell_res' : clip_cell_res,
                            'clip_cell_addresses' : " ".join([str(address) for address in cell_id_list])
                        }
                    )
                    if "clip_cell_densification" in conf_extra:
                        subset_conf.update(
                            {
                                'clip_cell_densification' : conf_extra['clip_cell_densification']
                            }
                        )

        subset_conf.update(specify_resolution(**conf_extra))
        subset_conf.update(specify_orient_type_args(**conf_extra))
        subset_conf.update(specify_topo_aperture(**conf_extra))

        input_extras = self.check_input_extra_fields(conf_extra)
        subset_conf.update(input_extras)

        output_extras = self.check_output_extra_fields(conf_extra)
        output_conf: DggridMetaConfigT = {
            'cell_output_type': 'GDAL',
            'point_output_type': 'NONE',
            'cell_output_gdal_format' : self.tmp_geo_out['legacy_driver'],
            'cell_output_file_name': str( (Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.{self.tmp_geo_out['ext']}").resolve()),
            **output_extras
        }

        if not self.has_gdal:
            output_conf.update({
                'cell_output_type': cast(DggsFeatureOutputTypeT, self.tmp_geo_out['legacy_driver'].upper()),
                'cell_output_file_name': str( (Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}").resolve())
            })
            output_conf.pop('cell_output_gdal_format', None)

        if not output_address_type is None and output_address_type in self.output_address_types:
            output_conf.update({'output_address_type': output_address_type})
        else:
            if self.debug:
                print(f"ignoring unknown output_address_type: {output_address_type}")

        dggs_ops = self.dgapi_grid_gen(dggs, subset_conf, output_conf)
        if self.debug:
            print(dggs_ops)

        gdf = gpd.read_file( Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.{self.tmp_geo_out['ext']}", driver=self.tmp_geo_out['driver'] )

        if not cell_id_list is None and len(cell_id_list) > 0 and not seq_df is None:
            # we have to adjust the columns formats for the IDs/Seqnums/Name field to ensure they are comparable for the join
            # as they are all seqnum cell IDs they should all be long integers, otherwise keep string type
            if not output_address_type == input_address_type:
                if self.debug:
                    print(f"cannot switch address types on the fly here: {input_address_type} !=  {output_address_type}")

            # seq_df['cell_exists'] = True
            if output_address_type in ['SEQNUM']:
                seq_df[input_address_type] = seq_df[input_address_type].astype(np.int64)

            # seq_df.set_index(input_address_type, inplace=True)
            # possible column names: 'name', 'Name', 'global_id'
            name_col = 'name'
            for potential_name_col in ['name', 'Name', 'global_id']:
                if potential_name_col in gdf.columns:
                    name_col = potential_name_col
                    break
            if output_address_type in ['SEQNUM']:
                gdf[name_col] = gdf[name_col].astype(np.int64)
            # gdf = gdf.join( seq_df, how='inner', left_on=name_col, right_on=input_address_type)
            # gdf = gdf.loc[gdf['cell_exists']].drop(columns=['cell_exists'])

        if not self.debug:
            try:
                os.remove( str( Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.{self.tmp_geo_out['ext']}") )
                if "SHAPEFILE" in self.tmp_geo_out['driver'].upper():
                    for ext in ['dbf', 'prj', 'shx', 'cpg']:
                        os.remove( str( Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.{ext}") )
            except Exception:
                pass
            try:
                os.remove( str( Path(tmp_dir) / f"temp_clip_{tmp_id}.txt") )
                if "SHAPEFILE" in self.tmp_geo_out['driver'].upper():
                    for ext in ['dbf', 'prj', 'shx', 'cpg']:
                        os.remove( str( Path(tmp_dir) / f"temp_clip_{tmp_id}.{ext}") )
            except Exception:
                pass
            try:
                os.remove( str( Path(tmp_dir) / f"points.gen") )
            except Exception:
                pass

        if split_dateline:
            return self.post_process_split_dateline(gdf)
        return gdf


    def grid_cell_centroids_from_cellids(
        self,
        cell_id_list: Sequence[str],
        dggs_type: DggsTypeT,
        resolution: int,
        mixed_aperture_level: int = None,
        clip_subset_type: DggsClipSubsetTypeT = 'WHOLE_EARTH',
        clip_cell_res: int = 1,
        input_address_type: DggsInputAddressTypeT = 'SEQNUM',
        output_address_type: DggsOutputAddressTypeT = 'SEQNUM',
        **conf_extra: DggridMetaConfigParameterT,
    ):
        """
        generates a DGGS grid and returns all the cell's centroid as GeoDataFrame with geometry type Point
            a) if cell_id_list is empty/None: grid cells for the WHOLE_EARTH
            b) if cell_id_list is a list/numpy array, takes this list as seqnums ids (potentially also Z3, Z7, or ZORDER) for subsetting
        """
        tmp_id = uuid.uuid4()
        tmp_dir = self.working_dir
        dggs = dgselect(dggs_type = dggs_type, res= resolution, mixed_aperture_level=mixed_aperture_level)
        dggs.update(**conf_extra, strict=True)

        subset_conf: DggridMetaConfigT = { 'update_frequency': 100000, 'clip_subset_type': clip_subset_type }
        seq_df = None

        if not cell_id_list is None and len(cell_id_list) > 0:

            seq_df = pd.DataFrame({ input_address_type: cell_id_list})
            seq_df.to_csv( str( (Path(tmp_dir) / f"temp_clip_{tmp_id}.txt").resolve()) , header=False, index=False, columns=[input_address_type], sep=' ')

            # TODO, for Z3, Z7, ZORDER can potentially also be COARSE_CELLS / aka parent cells?
            subset_conf.update({
                'clip_subset_type': 'SEQNUMS',
                'clip_region_files': str( (Path(tmp_dir) / f"temp_clip_{tmp_id}.txt").resolve()),
            })

            # TODO, for Z3, Z7, ZORDER can potentially also be COARSE_CELLS / aka parent cells?
            # clip_subset_type should INPUT_ADDRESS_TYPE for the equivalent of SEQNUM (tp use input_address_type Z3 ...), or COARSE_CELLS as an actual parent cell type clip (also for Z3 ..)
            if (
                input_address_type in self.input_address_types
                and output_address_type in self.output_address_types
            ):
                subset_conf.update(
                    {
                        'clip_subset_type': 'INPUT_ADDRESS_TYPE',
                        'input_address_type': input_address_type
                    }
                )

                if not clip_subset_type is None and clip_subset_type in ['COARSE_CELLS']:
                    subset_conf.update(
                        {
                            'clip_subset_type': clip_subset_type,
                            'input_address_type': input_address_type,
                            'clip_cell_res' : clip_cell_res,
                            'clip_cell_addresses' : " ".join([str(address) for address in cell_id_list])
                        }
                    )

        subset_conf.update(specify_resolution(**conf_extra))
        subset_conf.update(specify_orient_type_args(**conf_extra))
        subset_conf.update(specify_topo_aperture(**conf_extra))

        input_extras = self.check_input_extra_fields(conf_extra)
        subset_conf.update(input_extras)
        output_extras = self.check_output_extra_fields(conf_extra)
        output_conf: DggridMetaConfigT = {
            'point_output_type': 'GDAL',
            'cell_output_type': 'NONE',
            'point_output_gdal_format' : self.tmp_geo_out['legacy_driver'],
            'point_output_file_name': str( (Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.{self.tmp_geo_out['ext']}").resolve()),
            **output_extras
        }

        if not self.has_gdal:
            output_conf.update({
                'point_output_type': cast(DggsFeatureOutputTypeT, self.tmp_geo_out['legacy_driver'].upper()),
                'point_output_file_name': str( (Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}").resolve())
            })
            output_conf.pop('point_output_gdal_format', None)


        if not output_address_type is None and output_address_type in self.output_address_types:
            output_conf.update({'output_address_type': output_address_type})
        else:
            if self.debug:
                print(f"ignoring unknown output_address_type: {output_address_type}")

        dggs_ops = self.dgapi_grid_gen(dggs, subset_conf, output_conf )
        if self.debug:
            print(dggs_ops)

        gdf = gpd.read_file( Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.{self.tmp_geo_out['ext']}", driver=self.tmp_geo_out['driver'] )

        if not cell_id_list is None and len(cell_id_list) > 0 and not seq_df is None:
            # we have to adjust the columns formats for the IDs/Seqnums/Name field to ensure they are comparable for the join
            # as they are all seqnum cell IDs they should all be long integers, otherwise keep string type
            if not output_address_type == input_address_type:
                if self.debug:
                    print(f"cannot switch address types on the fly here: {input_address_type} !=  {output_address_type}")

            # seq_df['cell_exists'] = True
            if output_address_type in ['SEQNUM']:
                seq_df[input_address_type] = seq_df[input_address_type].astype(np.int64)

            # seq_df.set_index(input_address_type, inplace=True)
            name_col = 'name' if 'name' in gdf.columns else 'Name'
            if output_address_type in ['SEQNUM']:
                gdf[name_col] = gdf[name_col].astype(np.int64)
            # gdf = gdf.join( seq_df, how='inner', on=name_col)
            # gdf = gdf.loc[gdf['cell_exists']].drop(columns=['cell_exists'])

        if not self.debug:
            try:
                os.remove( str( Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.{self.tmp_geo_out['ext']}") )
                if "SHAPEFILE" in self.tmp_geo_out['driver'].upper():
                    for ext in ['dbf', 'prj', 'shx', 'cpg']:
                        os.remove( str( Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.{ext}") )
            except Exception:
                pass
            try:
                os.remove( str( Path(tmp_dir) / f"temp_clip_{tmp_id}.txt") )
                os.remove( str( Path(tmp_dir) / f"cells.gen") )
            except Exception:
                pass

        return gdf


    def grid_cellids_for_extent(
        self,
        dggs_type: DggsTypeT,
        resolution: int,
        mixed_aperture_level: int = None,
        clip_geom: "AnyGeometry" = None,
        output_address_type: DggsOutputAddressTypeT | None = None,
        **conf_extra: DggridMetaConfigParameterT,
    ) -> pd.DataFrame:
        """
        generates a DGGS grid and returns all the cellids as a pandas dataframe
            a) if clip_geom is empty/None: grid cell ids/seqnums for the WHOLE_EARTH
            b) if clip_geom is a shapely shape geometry, takes this as a clip area
            TODO could cellids be generated for COARSE_CELLS? Generate child id from list of parent ids?
        """
        tmp_id = uuid.uuid4()
        tmp_dir = self.working_dir
        dggs = dgselect(dggs_type = dggs_type, res= resolution, mixed_aperture_level=mixed_aperture_level)
        dggs.update(**conf_extra, strict=True)

        subset_conf = { 'update_frequency': 100000, 'clip_subset_type': 'WHOLE_EARTH' }

        if not clip_geom is None and clip_geom.area > 0:

            clip_gdf = gpd.GeoDataFrame(pd.DataFrame({'id' : [1], 'geometry': [clip_geom]}), geometry='geometry', crs=4326)
            clip_path = Path(tmp_dir) / f"temp_clip_{tmp_id}.{self.tmp_geo_out['ext']}"
            clip_gdf.to_file(str(clip_path), driver=self.tmp_geo_out['driver'] )

            subset_conf.update({
                'clip_subset_type': 'GDAL',
                'clip_region_files': str( (Path(tmp_dir) / f"temp_clip_{tmp_id}.{self.tmp_geo_out['ext']}").resolve()),
                })

            if not self.has_gdal:
                subset_conf.update({
                    'clip_subset_type': cast(DggsClipSubsetTypeT, self.tmp_geo_out['legacy_driver'].upper())
                })

        subset_conf.update(specify_resolution(**conf_extra))
        subset_conf.update(specify_orient_type_args(**conf_extra))
        subset_conf.update(specify_topo_aperture(**conf_extra))

        input_extras = self.check_input_extra_fields(conf_extra)
        subset_conf.update(input_extras)

        output_extras = self.check_output_extra_fields(conf_extra)
        output_conf: DggridMetaConfigT = {
            'point_output_type': 'TEXT',
            'cell_output_type': 'NONE',
            'point_output_file_name': str( (Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}").resolve()),
            **output_extras
        }

        if not output_address_type is None and output_address_type in self.output_address_types:
            output_conf.update({'output_address_type': output_address_type})
        else:
            if self.debug:
                print(f"ignoring unknown output_address_type: {output_address_type}")

        dggs_ops = self.dgapi_grid_gen(dggs, subset_conf, output_conf )
        if self.debug:
            print(dggs_ops)

        datatype = {0: str} if (output_address_type is not ['SEQNUM', 'AIGEN']) else {}
        df = pd.read_csv( Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.txt" , header=None, dtype=datatype)
        df = df.dropna()

        if not self.debug:
            try:
                os.remove( str( Path(tmp_dir) / f"temp_{dggs_type}_{resolution}_out_{tmp_id}.txt") )
                os.remove( str( Path(tmp_dir) / f"temp_clip_{tmp_id}.{self.tmp_geo_out['ext']}") )
                if "SHAPEFILE" in self.tmp_geo_out['driver'].upper():
                    for ext in ['dbf', 'prj', 'shx', 'cpg']:
                        os.remove( str( Path(tmp_dir) / f"temp_clip_{tmp_id}.{ext}") )
            except Exception:
                pass
            try:
                os.remove( str( Path(tmp_dir) / f"points.gen") )
            except Exception:
                pass

        return df


    def cells_for_geo_points(
        self,
        geodf_points_wgs84: gpd.GeoDataFrame,
        cell_ids_only: bool,
        dggs_type: DggsTypeT,
        resolution: int,
        mixed_aperture_level: int = None,
        split_dateline: bool = False,
        output_address_type: DggsOutputAddressTypeT | None = None,
        **conf_extra: DggridMetaConfigParameterT,
    ) -> gpd.GeoDataFrame:
        """
        takes a GeoDataFrame with point geometry and optional additional value columns and returns:
            a) if cell_ids_only == True: the same GeoDataFrame with an additional column with the cell ids
            b) if cell_ids_only == False: a new GeoDataFrame with geometry type Polygon, with column of cell ids and the additional columns
        """
        tmp_id = uuid.uuid4()
        tmp_dir = self.working_dir
        dggs = dgselect(dggs_type = dggs_type, res= resolution, mixed_aperture_level=mixed_aperture_level)
        dggs.update(**conf_extra, strict=True)

        cols = set(geodf_points_wgs84.columns.tolist())
        cols = cols - {'geometry'}
        geodf_points_wgs84['lon'] = geodf_points_wgs84['geometry'].x
        geodf_points_wgs84['lat'] = geodf_points_wgs84['geometry'].y
        cols_ordered = ['lon', 'lat']
        for col_name in cols:
            cols_ordered.append(col_name)

        geodf_points_wgs84[cols_ordered].to_csv( str( (Path(tmp_dir) / f"geo_{tmp_id}.txt").resolve()) , header=False, index=False, columns=cols_ordered, sep=' ')

        subset_conf: DggridMetaConfigT = {
            'input_file_name':  str( (Path(tmp_dir) / f"geo_{tmp_id}.txt").resolve()),
            'input_address_type': 'GEO',
            'input_delimiter': "\" \""
        }
        subset_conf.update(specify_resolution(**conf_extra))
        subset_conf.update(specify_orient_type_args(**conf_extra))
        subset_conf.update(specify_topo_aperture(**conf_extra))

        input_extras = self.check_input_extra_fields(conf_extra)
        subset_conf.update(input_extras)

        output_extras = self.check_output_extra_fields(conf_extra)
        output_conf = {
            'output_file_name': str( (Path(tmp_dir) / f"{output_address_type}_{tmp_id}.txt").resolve()),
            'output_address_type': 'SEQNUM',
            'output_delimiter': "\",\"",
            **output_extras
        }

        if not output_address_type is None and output_address_type in self.output_address_types:
            output_conf.update({'output_address_type': output_address_type})
        else:
            if self.debug:
                print(f"ignoring unknown output_address_type: {output_address_type}")

        dggs_ops = self.dgapi_grid_transform(dggs, subset_conf, **output_conf)
        if self.debug:
            print(dggs_ops)
        datatype = {0: str} if (output_address_type is not ['SEQNUM', 'AIGEN']) else {}
        df = pd.read_csv( dggs_ops['output_conf']['output_file_name'] , header=None, dtype=datatype)
        df = df.dropna()
        cell_id_list = df[0].values

        if not self.debug:
            try:
                os.remove( str( Path(tmp_dir) / f"geo_{tmp_id}.txt") )
                os.remove( str( Path(tmp_dir) / f"{output_address_type}_{tmp_id}.txt") )
                # probably not needed anymore
                os.remove( str( Path(tmp_dir) / f"seqnums_{tmp_id}.txt") )
            except Exception:
                pass

        if cell_ids_only:
            geodf_points_wgs84['name'] = cell_id_list
            return geodf_points_wgs84
        else:
            # grid_gen from seqnums
            dggs_conf = dggs.to_dict()
            gdf = self.grid_cell_polygons_from_cellids(
                cell_id_list=cell_id_list,
                # dggs_type=dggs_type,  passed via dggs_conf
                resolution=resolution,
                mixed_aperture_level=mixed_aperture_level,
                input_address_type=output_address_type,
                output_address_type=output_address_type,
                **dggs_conf,  # ensure any extra parameters are passed on
            )
            try:
                # avoid conflict between input 'name' column and generated 'name' Zone ID
                gdf.rename(columns={'name': 'zone'}, inplace=True)
                for col in cols_ordered:
                    gdf[col] = geodf_points_wgs84[col].values
            except Exception:
                pass

            if split_dateline:
                return self.post_process_split_dateline(gdf)

            return gdf


    def address_transform(
        self,
        cell_id_list: Sequence[str],
        dggs_type: DggsTypeT,
        resolution: int,
        mixed_aperture_level: int = None,
        input_address_type: DggsInputAddressTypeT = 'SEQNUM',
        output_address_type: DggsOutputAddressTypeT = 'SEQNUM',
        **conf_extra: DggridMetaConfigParameterT,
    ) -> pd.DataFrame:
        """
            generates the DGGS for the input cell_ids and returns all the transformed cell_ids
            cell_id_list is a list/numpy array, takes this list as seqnums ids (potentially also Z3, Z7, or ZORDER)
        """
        tmp_id = uuid.uuid4()
        tmp_dir = self.working_dir
        dggs = dgselect(dggs_type = dggs_type, res= resolution, mixed_aperture_level=mixed_aperture_level)
        dggs.update(**conf_extra, strict=True)

        if cell_id_list is None or len(cell_id_list) <= 0:
            raise ValueError("Expecting cell_id_list to transform.")

        if not input_address_type in self.input_address_types:
            raise ValueError(f"unknown input_address_type: {input_address_type}")

        if not output_address_type in self.output_address_types:
            raise ValueError(f"unknown output_address_type: {output_address_type}")

        seq_df = pd.DataFrame({ input_address_type: cell_id_list})
        seq_df.to_csv( str( (Path(tmp_dir) / f"temp_in_{input_address_type}_{tmp_id}.txt").resolve()) , header=False, index=False, columns=[input_address_type], sep=' ')

        subset_conf: DggridMetaConfigT = {
            'input_file_name':  str( (Path(tmp_dir) / f"temp_in_{input_address_type}_{tmp_id}.txt").resolve()),
            'input_address_type': input_address_type,
            'input_delimiter': "\" \""
        }

        input_extras = self.check_input_extra_fields(conf_extra)
        subset_conf.update(input_extras)

        output_extras = self.check_output_extra_fields(conf_extra)
        output_conf: DggridMetaConfigT = {
            'output_file_name': str( (Path(tmp_dir) / f"temp_out_{output_address_type}_{tmp_id}.txt").resolve()),
            'output_address_type': output_address_type,
            'output_delimiter': "\" \"",
            **output_extras,
        }

        dggs_ops = self.dgapi_grid_transform(dggs, subset_conf, **output_conf)
        if self.debug:
            print(dggs_ops)

        df = pd.read_csv( dggs_ops['output_conf']['output_file_name'] , header=None, dtype={0:'str', 1:'str'})
        df = df.dropna()
        seq_df[output_address_type] = df.iloc[:,0]

        if not self.debug:
            try:
                os.remove( str( Path(tmp_dir) / f"temp_in_{input_address_type}_{tmp_id}.txt") )
                os.remove( str( Path(tmp_dir) / f"temp_out_{output_address_type}_{tmp_id}.txt") )
            except Exception:
                pass

        return seq_df


class DGGRIDv7(DGGRID):
    version = 7
    output_address_types: DggsOutputAddressTypeV7T = output_address_types_v7
    output_extra_fields = output_extra_fields_v7
    input_address_types: DggsInputAddressTypeV7T = input_address_types_v7
    input_extra_fields = input_extra_fields_v7


class DGGRIDv8(DGGRID):
    version = 8
    output_address_types: DggsOutputAddressTypeV8T = output_address_types_v8
    output_extra_fields = output_extra_fields_v8
    input_address_types: DggsInputAddressTypeV8T = input_address_types_v8
    input_extra_fields = input_extra_fields_v8


class AnyDGGRID(DGGRID):
    version = 0
    output_address_types: DggsOutputAddressTypeT = list(set(output_address_types_v7) | set(output_address_types_v8))
    output_extra_fields = list(set(output_extra_fields_v7) | set(output_extra_fields_v8))
    input_address_types: DggsInputAddressTypeT = list(set(input_address_types_v7) | set(input_address_types_v8))
    input_extra_fields = list(set(input_extra_fields_v7) | set(input_extra_fields_v8))


#############################################################
"""
below is an earlier attempt on custom DGGS configuration, currently we only barely support the predefined DGGS_TYPES
"""


def specify_orient_type_args(
    # use official DGGRID parameter names to map directly by keyword unpacking
    # note: first 3 are placed in the same order as originals to handle positional invocation seamlessly
    dggs_orient_specify_type=None,
    dggs_vert0_lon: DggsDegree = None,
    dggs_vert0_lat: DggsDegree = None,
    dggs_vert0_azimuth: DggsDegree = None,
    dggs_orient_rand_seed=42,
    # backward compatibility parameters
    orient_type=None,
    **_,  # ignore unknowns
) -> DggridMetaConfigT:
    dggs_orient_specify_type = dggs_orient_specify_type or orient_type

    if dggs_orient_specify_type == 'SPECIFIED' or any(
        val is not None for val in [dggs_vert0_lon, dggs_vert0_lat, dggs_vert0_azimuth]
    ):
        if dggs_vert0_lon is not None:
            dggs_vert0_lon_val = decimal.Decimal(dggs_vert0_lon)
            if dggs_vert0_lon_val < -180.0 or dggs_vert0_lon_val > 180.0:
                raise ValueError('dggs_vert0_lon must be in range [-180,180]')
        if dggs_vert0_lat is not None:
            dggs_vert0_lat_val = decimal.Decimal(dggs_vert0_lat)
            if dggs_vert0_lat_val < -90.0 or dggs_vert0_lat_val > 90.0:
                raise ValueError('dggs_vert0_lat must be in range [-90,90]')
        if dggs_vert0_azimuth is not None:
            dggs_vert0_azimuth_val = decimal.Decimal(dggs_vert0_azimuth)
            if dggs_vert0_azimuth_val < 0.0 or dggs_vert0_azimuth_val > 360.0:
                raise ValueError('dggs_vert0_azimuth must be in range [0,360]')
        return {
            'dggs_orient_specify_type' : 'SPECIFIED',
            'dggs_vert0_lon' : dggs_vert0_lon,
            'dggs_vert0_lat' : dggs_vert0_lat,
            'dggs_vert0_azimuth' : dggs_vert0_azimuth
        }
    if dggs_orient_specify_type == 'RANDOM':
        return { 'dggs_orient_rand_seed' : dggs_orient_rand_seed }

    # else default REGION_CENTER
    return {}


def specify_topo_aperture(
    # use official DGGRID parameter names to map directly by keyword unpacking
    # note: first 3 are placed in the same order as originals to handle positional invocation seamlessly
    dggs_topology_type: DggsTopologyTypeT = None,
    dggs_aperture_type: DggsApertureTypeT = None,
    dggs_aperture: DggsApertureT = None,
    dggs_num_aperture_4_res: int = 0,
    dggs_aperture_sequence: str = "333333333333",
    # backward compatibility parameters
    topology_type: DggsTopologyTypeT = None,
    aperture_type: DggsApertureTypeT = None,
    aperture_res: DggsApertureT = None,
    **_,  # ignore unknowns
) -> DggridMetaConfigT:
    dggs_topology_type = dggs_topology_type or topology_type
    dggs_aperture_type = dggs_aperture_type or aperture_type
    dggs_aperture: DggsApertureT = dggs_aperture or aperture_res

    if dggs_topology_type is None and dggs_aperture_type is None:
        return {}

    if not dggs_topology_type in dggs_topologies or not dggs_aperture_type in dggs_aperture_types:
        raise ValueError('topology or aperture type unknown')

    if dggs_aperture_type == 'PURE':
        if dggs_topology_type == 'HEXAGON':
            if not dggs_aperture in [3, 4, 7]:
                print(f"combo not possible / {dggs_topology_type} {dggs_aperture} / setting 3H")
                return {
                    '#short_name': "3H",
                    'dggs_topology': dggs_topology_type,
                    'dggs_aperture_type': dggs_aperture_type,
                    'dggs_aperture': 3
                }
            else:
                return {
                    '#short_name': f"{dggs_aperture}H",
                    'dggs_topology': dggs_topology_type,
                    'dggs_aperture_type': dggs_aperture_type,
                    'dggs_aperture': dggs_aperture
                }

        if dggs_topology_type in ['TRIANGLE', 'DIAMOND']:
            if not dggs_aperture == 4:
                print(f"combo not possible / {dggs_topology_type} {dggs_aperture} / setting 4{dggs_topology_type[0]}")
                return {
                    '#short_name': f"4{dggs_topology_type[0]}",
                    'dggs_topology': dggs_topology_type,
                    'dggs_aperture_type': dggs_aperture_type,
                    'dggs_aperture': 4
                }
            else:
                return {
                    '#short_name': f"4{dggs_topology_type[0]}",
                    'dggs_topology': dggs_topology_type,
                    'dggs_aperture_type': dggs_aperture_type,
                    'dggs_aperture': dggs_aperture
                }

    elif dggs_aperture_type == 'MIXED43':
        # dggs_aperture is ignored, only HEXAGON can have MIXED34
        if dggs_topology_type == 'HEXAGON':
            # dggs_num_aperture_4_res (default 0)
            return {
                '#short_name': f"43H",
                'dggs_topology': dggs_topology_type,
                'dggs_aperture_type': dggs_aperture_type,
                'dggs_num_aperture_4_res': dggs_num_aperture_4_res,
                '#dggs_aperture': dggs_aperture
            }
        else:
            raise ValueError('not yet implemented')

    elif dggs_aperture_type == "SEQUENCE":
        # dggs_aperture_sequence (default “333333333333”).
        return {
            '#short_name': f"SEQ{dggs_topology_type[0]}",
            'dggs_topology': dggs_topology_type,
            'dggs_aperture_type': dggs_aperture_type,
            'dggs_aperture_sequence': str(dggs_aperture_sequence),
            '#dggs_aperture': dggs_aperture
        }

    return {}


def specify_resolution(
    # use official DGGRID parameter names to map directly by keyword unpacking
    # note: first 2 are placed in the same order as original to handle positional invocation seamlessly
    dggs_proj: DggsProjectionT = None,
    dggs_res_specify_type: DggsResSpecifyTypeT =None,
    dggs_res_spec: int = 9,
    dggs_res_specify_area: float = 120000,
    dggs_res_specify_intercell_distance: float = 4000,
    dggs_res_specify_rnd_down: bool = True,
    # backward compatibility
    proj_spec=None,
    dggs_res_spec_type=None,
    **_,  # ignore unknowns
) -> DggridMetaConfigT:
    dggs_proj = dggs_proj or proj_spec
    dggs_res_specify_type = dggs_res_specify_type or dggs_res_spec_type

    if dggs_proj is None or dggs_res_specify_type is None:
        return {}

    if not dggs_proj in list(dggs_projections) or not dggs_res_specify_type in dggs_res_specify_types:
        raise ValueError("base projection (ISEA or FULLER) or resolution spec unknown")

    if dggs_res_specify_type == "SPECIFIED":
        return {
            'dggs_proj': dggs_proj,
            'dggs_res_spec': dggs_res_spec
        }

    elif dggs_res_specify_type == 'CELL_AREA':
        return {
            'dggs_proj': dggs_proj,
            'dggs_res_specify_area': dggs_res_specify_area, # (in square kilometers)
            'dggs_res_specify_rnd_down' : dggs_res_specify_rnd_down
        }
    elif dggs_res_specify_type == 'INTERCELL_DISTANCE':
        return {
            'dggs_proj': dggs_proj,
            'dggs_res_specify_intercell_distance': dggs_res_specify_intercell_distance, # (in kilometers)
            'dggs_res_specify_rnd_down' : dggs_res_specify_rnd_down
        }

    return {}


def dgconstruct(dggs_type: DggsTypeT        = "CUSTOM", # dggs_type
                projection: DggsProjectionT = 'ISEA',  # dggs_projection
                aperture: DggsApertureT     = 3,  # dggs_aperture_type / dggs_aperture
                topology: DggsTopologyTypeT = 'HEXAGON', # dggs_topology
                res: int                    = None, # dggs_res_spec
                precision: int              = 7,
                area: float                 = None, # dggs_res_specify_area
                spacing: float              = None,
                cls_val: float              = None, # dggs_res_specify_intercell_distance
                resround: str               = 'nearest',
                metric: bool                = True,
                show_info: bool             = True,
                azimuth_deg: float          = 0, # dggs_vert0_azimuth
                pole_lat_deg: float         = 58.28252559, # dggs_vert0_lat
                pole_lon_deg: float         = 11.25 # dggs_vert0_lon
                ):

    if not len(list(filter(lambda x: not x is None, [res,area,spacing,cls_val])))  == 1:
        raise ValueError('dgconstruct(): Only one of res, area, length, or cls can have a value!')

    #Use a dummy resolution, we'll fix it in a moment
    dggs = Dggs(
        dggs_type    = dggs_type,
        pole_lon_deg = pole_lon_deg,
        pole_lat_deg = pole_lat_deg,
        azimuth_deg  = azimuth_deg,
        aperture     = aperture,
        res          = 1,
        topology     = topology.upper(),
        projection   = projection.upper(),
        precision    = precision
    )

    if not res is None:
        dggs.res = res
    elif not area is None:
        dggs.res = dggs.dg_closest_res_to_area (area=area, resround=resround,metric=metric,show_info=True)
    elif not spacing is None :
        dggs.res = dggs.dg_closest_res_to_spacing(spacing=spacing,resround=resround,metric=metric,show_info=True)
    elif not cls_val is None:
        dggs.res = dggs.dg_closest_res_to_cls ( cls_val=cls_val, resround=resround,metric=metric,show_info=True)
    else:
        raise ValueError('dgconstruct(): Logic itself has failed us.')

    # dggs.dgverify()

    return dggs
