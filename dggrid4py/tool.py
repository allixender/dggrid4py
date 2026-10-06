"""
Download of the portable DGGRID binaries from https://github.com/allixender/DGGRID_portables

The portable binaries are built without GDAL, use them with ``has_gdal=False``.
"""
import hashlib
import os
import platform
import shutil
import stat
import tarfile
import tempfile
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

PORTABLES_URL = "https://github.com/allixender/DGGRID_portables/releases/download"
# lines of portable binaries and their release tag in DGGRID_portables:
# 'stable' is the DGGRID release that this dggrid4py version is tested with,
# 'edge' is the rolling pre-release, built from the current DGGRID master
PORTABLE_LINES = {
    "stable": "v8.44",
    "edge": "edge",
}
DEFAULT_LINE = "stable"
CHECKSUMS_FILE = "SHA256SUMS"

_portable_assets = {
    ('linux', 'x86_64'): 'dggrid-linux-x86_64.tar.gz',
    ('linux', 'arm64'): 'dggrid-linux-arm64.tar.gz',
    ('darwin', 'x86_64'): 'dggrid-macos-x86_64.tar.gz',
    ('darwin', 'arm64'): 'dggrid-macos-arm64.tar.gz',
    ('windows', 'x86_64'): 'dggrid-windows-x86_64.zip',
    ('windows', 'arm64'): 'dggrid-windows-arm64.zip',
}

_machine_aliases = {'amd64': 'x86_64', 'x64': 'x86_64', 'aarch64': 'arm64'}


def portable_asset_name(system=None, machine=None):
    """
    Name of the portable DGGRID archive for a platform, by default for the current one.
    """
    system = (system or platform.system()).lower()
    machine = (machine or platform.machine()).lower()
    machine = _machine_aliases.get(machine, machine)
    try:
        return _portable_assets[(system, machine)]
    except KeyError:
        raise ValueError(
            f"No portable executable available for {system} {machine}. "
            "Please use dggrid4py with a local DGGRID installation."
        ) from None


def _fetch(url, local_path):
    with urllib.request.urlopen(url, timeout=60) as response, open(local_path, 'wb') as f:
        shutil.copyfileobj(response, f)


def _sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _make_executable(path):
    current_permissions = os.stat(path).st_mode
    os.chmod(path, current_permissions | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _release_checksum(release, asset, folder):
    sums_path = Path(folder) / f".{CHECKSUMS_FILE}"
    _fetch(f"{PORTABLES_URL}/{release}/{CHECKSUMS_FILE}", sums_path)
    try:
        for line in sums_path.read_text().splitlines():
            parts = line.split()
            if len(parts) == 2 and parts[1].lstrip('*') == asset:
                return parts[0].lower()
    finally:
        sums_path.unlink()
    raise ValueError(f"no checksum for {asset} in {CHECKSUMS_FILE} of release '{release}'")


def _safe_member(name, folder):
    # archive members must stay inside the target folder
    root = Path(folder).resolve()
    target = (root / name).resolve()
    if target != root and root not in target.parents:
        raise ValueError(f"unsafe path in archive: {name}")
    return name


def _extract(archive_path, folder):
    """
    unpacks the archive into folder and returns the path of the dggrid executable in it
    """
    if str(archive_path).endswith('.zip'):
        with zipfile.ZipFile(archive_path) as archive:
            names = [_safe_member(name, folder) for name in archive.namelist()]
            archive.extractall(folder)
    else:
        with tarfile.open(archive_path, 'r:gz') as archive:
            names = [_safe_member(name, folder) for name in archive.getnames()]
            if hasattr(tarfile, 'data_filter'):
                archive.extractall(folder, filter='data')
            else:
                archive.extractall(folder)

    for name in names:
        if Path(name).name in ('dggrid', 'dggrid.exe'):
            return Path(folder) / name
    raise ValueError(f"no dggrid executable found in {Path(archive_path).name}")


def download_executable(url, folder="./"):
    """
    Download a single executable file into folder and return its absolute path.
    """
    os.makedirs(folder, exist_ok=True)

    filename = url.split('/')[-1]
    local_path = os.path.join(folder, filename)

    _fetch(url, local_path)
    _make_executable(local_path)

    return os.path.abspath(local_path)


def portable_release(line=DEFAULT_LINE):
    """
    Release tag in DGGRID_portables for a line of portable binaries.

    ``stable`` and ``edge`` are the lines that dggrid4py knows. Any other value is taken as the release tag
    itself, e.g. ``edge-v91b`` for a development line or ``v8.44`` for a specific DGGRID release.
    """
    return PORTABLE_LINES.get(line, line)


def get_portable_executable(folder="./", line=DEFAULT_LINE, force=False):
    """
    Download the portable DGGRID binary for the current platform and return the absolute path of the executable.

    The archive is taken from the DGGRID_portables release of the given line, verified against the
    ``SHA256SUMS`` file of that release, and unpacked into a subfolder of ``folder`` that is named after the
    release. A binary that is already there is used again as long as its checksum is the one of the release,
    so a rolling release is downloaded again only after it was rebuilt. Without a network connection, a binary
    that is already there is returned.

    Args:
        folder (str): where the binaries are kept, created if needed
        line (str): ``stable`` (default, the DGGRID release that this dggrid4py version is tested with),
            ``edge`` (rolling pre-release of the current DGGRID development version), or any release tag
            of DGGRID_portables, e.g. ``edge-v91b``
        force (bool): download again, also if the binary is already there

    Returns:
        str: absolute path of the ``dggrid`` executable

    Raises:
        ValueError: if there is no portable binary for this platform, no such release, or the checksum does not match
    """
    asset = portable_asset_name()
    release = portable_release(line)
    folder = Path(folder) / release
    folder.mkdir(parents=True, exist_ok=True)

    # remembers checksum and location of the unpacked executable
    marker = folder / f".{asset}.sha256"
    cached_checksum, cached_executable = None, None
    if marker.is_file():
        cached_checksum, _, relative_path = marker.read_text().strip().partition(' ')
        if relative_path and (folder / relative_path).is_file():
            cached_executable = folder / relative_path

    try:
        checksum = _release_checksum(release, asset, folder)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            if not any(folder.iterdir()):
                folder.rmdir()
            raise ValueError(
                f"no portable DGGRID release '{release}' (line '{line}') in DGGRID_portables, "
                f"known lines are {list(PORTABLE_LINES)}"
            ) from None
        raise
    except OSError:
        if cached_executable is not None and not force:
            return str(cached_executable.resolve())
        raise

    if cached_executable is not None and cached_checksum == checksum and not force:
        return str(cached_executable.resolve())

    with tempfile.TemporaryDirectory(dir=folder) as tmp_dir:
        archive_path = Path(tmp_dir) / asset
        _fetch(f"{PORTABLES_URL}/{release}/{asset}", archive_path)
        actual = _sha256(archive_path)
        if actual != checksum:
            raise ValueError(f"checksum of {asset} does not match {CHECKSUMS_FILE}: {actual} != {checksum}")
        executable = _extract(archive_path, folder)

    _make_executable(executable)
    marker.write_text(f"{checksum} {executable.relative_to(folder).as_posix()}\n")
    return str(executable.resolve())
