#!/usr/bin/python3
# coding=utf8
#
# Copyright (c) 2022 - Luís Moreira de Sousa & Alexander Kmoch
# Licenced under GNU AFFERO GENERAL PUBLIC LICENSE. Please consult the LICENCE 
# file for details.
#
# Splits every cell from an input DGG crossing the 180º meridian in two. Meant
# to facilitate the display of a DGG in a cartesian GIS such as QGIS. Presently,
# this is achieved by testing every single cell in a DGGS. Future versions might
# use heuristics to restric the number of cells tested.
#
# Both the input and output of this programme are shapefiles. While not optimal,
# this format roots on the outputs of the DGGRID library.
#
# Author: Luís Moreira de Sousa (luis.de.sousa[@]protonmail.ch)
# Date: 07-10-2022 
#

import geopandas as gpd
from shapely.geometry import Polygon
from shapely.ops import split
from shapely.wkt import loads
from shapely.affinity import translate

intersectWest = Polygon([(180.005,-90),(180.005,90),(360,90),(360,-90),(180.005,-90)])
intersectEast = Polygon([(0,-90),(0,90),(179.995,90),(179.995,-90),(0,-90)])


def get_geom_coords(geometry):

    wkt = geometry.wkt
    wkt = wkt.replace("POLYGON ((","")
    wkt = wkt.replace("))","")
    wkt = wkt.split(", ")

    return wkt


def crosses_interruption(coords):

    lat_0 = float(coords[0].split(" ")[0])

    for i in range(1, len(coords)):
        lat_1 = float(coords[i].split(" ")[0])
        if (abs(lat_0 - lat_1) > 180):
            return True

    return False


def interrupt_cell(coords):

    result = []
    
    # Translate latitudes to east of 180ª 
    for i in range(len(coords)):
        lat = float(coords[i].split(" ")[0])
        if lat < 0:
            coords[i] = str(lat + 360) + " " + coords[i].split(" ")[1]

    coords = "POLYGON((%s))" % str(coords).replace("'","").replace("[","").replace("]","")
    poly = loads(coords)
    
    result.append(intersectEast.intersection(poly))
    result.append(translate(intersectWest.intersection(poly), xoff=-360))

    return result

