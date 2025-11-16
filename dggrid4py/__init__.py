# -*- coding: utf-8 -*-
#
# Copyright (c) 2022 - Alexander Kmoch
# Licenced under GNU AFFERO GENERAL PUBLIC LICENSE. Please consult the LICENCE 
# file for details.
#
# Author: Alexander Kmoch (alexander.kmoch@ut.ee)
# Date: 07-12-2022 
#

from .dggrid_runner import DGGRIDv7, DGGRIDv8, Dggs, dgselect, dggs_types
import dggrid4py.igeo7 as igeo7
import dggrid4py.tool as tool
from .auxlat import geoseries_to_authalic, geoseries_to_geodetic, geodetic_to_authalic, authalic_to_geodetic

__version__ = "0.5.3"
