"""Resolve the Cellulose_interface package version from one authoritative value."""

from __future__ import annotations

import tomllib
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


PACKAGE_NAME = "cellulose-ethanol-tools"


def project_version() -> str:
    try:
        return version(PACKAGE_NAME)
    except PackageNotFoundError:
        pyproject = Path(__file__).resolve().with_name("pyproject.toml")
        return str(tomllib.loads(pyproject.read_text())["project"]["version"])
