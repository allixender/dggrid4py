#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Version of the DGGRID executable, from the output of dggrid -v.
"""
import os
import shutil
import stat

import pytest

from dggrid4py import DGGRIDv8


def _fake_dggrid(tmp_path, output, name="dggrid"):
    path = tmp_path / name
    path.write_text(f'#!/bin/sh\necho "{output}"\necho called >> "{tmp_path}/calls"\n')
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


@pytest.mark.parametrize(
    "output, expected",
    [
        ("DGGRID version 9.1 released October 23, 2026", (9, 1)),
        ("DGGRID version 8.44 released December 1, 2025", (8, 44)),
        ("DGGRID version 9.0b", (9, 0)),
        ("** executing DGGRID version 8.42 with GDAL **", (8, 42)),
        ("usage: dggrid metaFileName", None),
    ],
)
def test_version_from_output(tmp_path, output, expected):
    dggrid = DGGRIDv8(executable=_fake_dggrid(tmp_path, output))
    assert dggrid.dggrid_version == expected


def test_version_missing_executable_is_none(tmp_path):
    assert DGGRIDv8(executable=tmp_path / "not_there").dggrid_version is None


def test_version_asked_once_per_executable(tmp_path):
    dggrid = DGGRIDv8(executable=_fake_dggrid(tmp_path, "DGGRID version 9.1 released"))
    assert dggrid.dggrid_version == (9, 1)
    assert dggrid.dggrid_version == (9, 1)
    assert (tmp_path / "calls").read_text().count("called") == 1

    # is_runnable() and users can change the executable of an instance
    dggrid.executable = _fake_dggrid(tmp_path, "DGGRID version 8.44 released", name="dggrid844")
    assert dggrid.dggrid_version == (8, 44)


def test_version_of_real_executable():
    path = os.getenv("DGGRID_PATH") or shutil.which("dggrid")
    if not path or not os.path.isfile(path):
        pytest.skip("DGGRID executable not available")
    version = DGGRIDv8(executable=path).dggrid_version
    assert version is not None and version[0] >= 8
