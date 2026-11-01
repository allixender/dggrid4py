#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pytest

from dggrid4py.dggrid_runner import Dggs, dg_grid_meta, dgselect, dggs_types_v9


def test_dgselect_minimal_only_given_parameters():
    assert dgselect("IGEO7v2", minimal=True, res=9).metafile == ["dggs_type IGEO7v2", "dggs_res_spec 9"]


def test_dgselect_minimal_mixed_aperture_level_and_precision():
    meta = dgselect("ISEA43HS", minimal=True, res=5, mixed_aperture_level=2, precision=9).metafile
    assert sorted(meta) == ["dggs_num_aperture_4_res 2", "dggs_res_spec 5", "dggs_type ISEA43HS", "precision 9"]


def test_dgselect_minimal_unknown_type():
    with pytest.raises(ValueError, match="unknown dggs_type: 'ISEA5H'"):
        dgselect("ISEA5H", minimal=True, res=5)


def test_dgselect_default_keeps_explicit_parameters():
    meta = dgselect("ISEA7H", res=5).metafile
    assert "dggs_vert0_lon 11.25" in meta
    assert "dggs_proj ISEA" in meta


def test_extra_params_go_to_the_metafile():
    dggs = dgselect("ISEA7HL", minimal=True, res=5)
    dggs.extra_params["polyclip_algo"] = "CLIPPER1"
    assert "polyclip_algo CLIPPER1" in dg_grid_meta(dggs)
    assert "extra_params" not in dggs.to_dict()


def test_extra_params_orientation_is_checked():
    dggs = Dggs(dggs_type="ISEA7HL", projection=None, aperture=None, topology=None, precision=None,
                azimuth_deg=None, pole_lat_deg=None, pole_lon_deg=None)
    dggs.extra_params["dggs_orient_specify_type"] = "NOT_A_TYPE"
    with pytest.raises(ValueError, match="dggs_orient_specify_type"):
        dg_grid_meta(dggs)


def test_all_v9_types_in_list_are_unique():
    assert len(dggs_types_v9) == len(set(dggs_types_v9))
