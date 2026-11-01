#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
dggs_type IGEO7 means IGEO7v2 (geodetic latitudes in and out) since DGGRID 9.1. The DGGRIDv7 and DGGRIDv8 classes
keep their IGEO7 on any DGGRID version, otherwise the latitudes get converted twice and the cells are wrong, without an error.
"""
import warnings

import pytest

from dggrid4py import DGGRIDv7, DGGRIDv8
from dggrid4py.dggrid_runner import DGGRID, AnyDGGRID

OPS = ['dggrid_operation GENERATE_GRID', 'dggs_type IGEO7', 'dggs_res_spec 5']


def _instance(cls, monkeypatch, version):
    monkeypatch.setattr(DGGRID, 'dggrid_version', version)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', DeprecationWarning)
        return cls(executable='dggrid')


@pytest.mark.parametrize('cls', [DGGRIDv7, DGGRIDv8, AnyDGGRID])
@pytest.mark.parametrize('version', [(9, 1), (9, 2), (10, 0)])
def test_igeo7_runs_as_igeo7v1(cls, version, monkeypatch):
    ops = _instance(cls, monkeypatch, version)._metafile_for_executable(OPS)
    assert ops == ['dggrid_operation GENERATE_GRID', 'dggs_type IGEO7v1', 'dggs_res_spec 5']
    assert OPS[1] == 'dggs_type IGEO7'


@pytest.mark.parametrize('version', [(8, 42), (8, 44), (9, 0), None])
def test_igeo7_unchanged_before_dggrid91(version, monkeypatch):
    assert _instance(DGGRIDv8, monkeypatch, version)._metafile_for_executable(OPS) == OPS


@pytest.mark.parametrize('line', ['dggs_type ISEA7H', 'dggs_type IGEO7v1', 'dggs_type IGEO7v2', 'kml_name dggs_type IGEO7'])
def test_other_metafile_lines_unchanged(line, monkeypatch):
    assert _instance(DGGRIDv8, monkeypatch, (9, 1))._metafile_for_executable([line]) == [line]
