#!/usr/bin/env python
# -*- coding: utf-8 -*-
import hashlib
import io
import os
import shutil
import tarfile
import zipfile
from pathlib import Path

import pytest

from dggrid4py import tool


@pytest.mark.parametrize(
    "system, machine, asset",
    [
        ("Linux", "x86_64", "dggrid-linux-x86_64.tar.gz"),
        ("Linux", "aarch64", "dggrid-linux-arm64.tar.gz"),
        ("Darwin", "arm64", "dggrid-macos-arm64.tar.gz"),
        ("Darwin", "x86_64", "dggrid-macos-x86_64.tar.gz"),
        ("Windows", "AMD64", "dggrid-windows-x86_64.zip"),
        ("Windows", "ARM64", "dggrid-windows-arm64.zip"),
    ],
)
def test_portable_asset_name(system, machine, asset):
    assert tool.portable_asset_name(system, machine) == asset


def test_portable_asset_name_unsupported_platform():
    with pytest.raises(ValueError, match="No portable executable available for freebsd riscv64"):
        tool.portable_asset_name("FreeBSD", "riscv64")


def _tar_gz(path, members):
    with tarfile.open(path, "w:gz") as archive:
        for name, content in members.items():
            info = tarfile.TarInfo(name)
            info.size = len(content)
            archive.addfile(info, io.BytesIO(content))


class _Release:
    """a local stand-in for a DGGRID_portables release"""

    def __init__(self, tmp_path, monkeypatch, asset="dggrid-linux-x86_64.tar.gz", members=None):
        self.asset = asset
        self.source = tmp_path / "release"
        self.source.mkdir()
        self.fetched = []
        self.offline = False
        if members is None:
            members = {"dggrid-9.0b-linux-x86_64/dggrid": b"#!/bin/sh\n", "dggrid-9.0b-linux-x86_64/LICENSE": b"AGPL"}
        self.build(members)
        monkeypatch.setattr(tool, "portable_asset_name", lambda: self.asset)
        monkeypatch.setattr(tool, "_fetch", self.fetch)

    def build(self, members):
        archive = self.source / self.asset
        if self.asset.endswith(".zip"):
            with zipfile.ZipFile(archive, "w") as zf:
                for name, content in members.items():
                    zf.writestr(name, content)
        else:
            _tar_gz(archive, members)
        checksum = hashlib.sha256(archive.read_bytes()).hexdigest()
        (self.source / "SHA256SUMS").write_text(f"{'0' * 64}  dggrid-other.tar.gz\n{checksum}  {self.asset}\n")

    def fetch(self, url, local_path):
        if self.offline:
            raise OSError("no network")
        assert url.startswith(f"{tool.PORTABLES_URL}/edge/")
        name = url.split("/")[-1]
        self.fetched.append(name)
        shutil.copyfile(self.source / name, local_path)


def test_get_portable_executable_unpacks_and_caches(tmp_path, monkeypatch):
    release = _Release(tmp_path, monkeypatch)
    folder = tmp_path / "bin"

    executable = tool.get_portable_executable(folder)
    assert executable == str((folder / "dggrid-9.0b-linux-x86_64" / "dggrid").resolve())
    assert os.access(executable, os.X_OK)
    assert release.fetched == ["SHA256SUMS", release.asset]
    # no archive or checksum list is left behind
    assert sorted(p.name for p in folder.iterdir() if not p.name.startswith(".")) == ["dggrid-9.0b-linux-x86_64"]

    # same release checksum: the binary is used again, only the checksums are looked up
    assert tool.get_portable_executable(folder) == executable
    assert release.fetched == ["SHA256SUMS", release.asset, "SHA256SUMS"]

    # force downloads again
    assert tool.get_portable_executable(folder, force=True) == executable
    assert release.fetched[-2:] == ["SHA256SUMS", release.asset]

    # the rolling release was rebuilt: downloaded again
    release.build({"dggrid-9.1-linux-x86_64/dggrid": b"#!/bin/sh\necho new\n"})
    rebuilt = tool.get_portable_executable(folder)
    assert rebuilt == str((folder / "dggrid-9.1-linux-x86_64" / "dggrid").resolve())

    # offline: the binary that is already there is returned
    release.offline = True
    assert tool.get_portable_executable(folder) == rebuilt
    with pytest.raises(OSError):
        tool.get_portable_executable(tmp_path / "empty")


def test_get_portable_executable_zip(tmp_path, monkeypatch):
    members = {"dggrid-9.0b-windows-x86_64/dggrid.exe": b"MZ", "dggrid-9.0b-windows-x86_64/LICENSE": b"AGPL"}
    _Release(tmp_path, monkeypatch, asset="dggrid-windows-x86_64.zip", members=members)
    executable = tool.get_portable_executable(tmp_path / "bin")
    assert Path(executable).name == "dggrid.exe"
    assert Path(executable).read_bytes() == b"MZ"


def test_get_portable_executable_checksum_mismatch(tmp_path, monkeypatch):
    release = _Release(tmp_path, monkeypatch)
    (release.source / "SHA256SUMS").write_text(f"{'0' * 64}  {release.asset}\n")
    with pytest.raises(ValueError, match="checksum of dggrid-linux-x86_64.tar.gz does not match"):
        tool.get_portable_executable(tmp_path / "bin")
    assert not any(p.name.startswith("dggrid") for p in (tmp_path / "bin").iterdir())


def test_get_portable_executable_missing_checksum(tmp_path, monkeypatch):
    release = _Release(tmp_path, monkeypatch)
    (release.source / "SHA256SUMS").write_text(f"{'0' * 64}  dggrid-other.tar.gz\n")
    with pytest.raises(ValueError, match="no checksum for"):
        tool.get_portable_executable(tmp_path / "bin")


@pytest.mark.parametrize("name", ["../evil/dggrid", "/tmp/evil/dggrid"])
def test_get_portable_executable_rejects_paths_outside_folder(tmp_path, monkeypatch, name):
    _Release(tmp_path, monkeypatch, members={name: b"#!/bin/sh\n"})
    with pytest.raises(ValueError, match="unsafe path in archive"):
        tool.get_portable_executable(tmp_path / "bin")
    assert not (tmp_path / "evil").exists()


def test_get_portable_executable_without_dggrid_in_archive(tmp_path, monkeypatch):
    _Release(tmp_path, monkeypatch, members={"dggrid-9.0b-linux-x86_64/LICENSE": b"AGPL"})
    with pytest.raises(ValueError, match="no dggrid executable found"):
        tool.get_portable_executable(tmp_path / "bin")
