#!python3
from __future__ import annotations

import argparse
from pathlib import Path

import tomli
from packaging.requirements import Requirement
from packaging.version import Version


def min_dep(req: Requirement) -> str:
    """
    Given a requirement, return the minimum version specifier.

    Example
    -------

    >>> min_dep(Requirement("numpy>=1.0"))
    "numpy==1.0"
    """
    req_name = req.name
    if req.extras:
        req_name = f"{req_name}[{','.join(req.extras)}]"

    # TODO: Should this be allowed?
    if not req.specifier:
        return req_name

    min_version = Version("0.0.0.a1")
    for spec in req.specifier:
        if spec.operator in [">", ">=", "~-"]:
            min_version = max(min_version, Version(spec.version))

    # TODO: should this return `~=` or `==`?
    return f"{req_name}=={min_version}.*"


def extract_min_deps(dependencies: list[str]) -> list[str]:
    return [min_dep(Requirement(dep)) for dep in dependencies]


def main():
    # TODO: Allow optional dependencies
    parser = argparse.ArgumentParser(
        prog="min-deps",
        description="""Parse a pyproject.toml file and output a list of minimum dependencies.

        Output is directly passable to `pip install`.""",
        usage="pip install `python min-deps.py pyproject.toml`",
    )
    parser.add_argument(
        "path", type=Path, help="pyproject.toml to parse minimum dependencies from"
    )
    parser.add_argument("--extras", type=str, nargs="*", help="extras to install")

    args = parser.parse_args()

    pyproject = tomli.loads(args.path.read_text())

    deps = pyproject["project"]["dependencies"]

    for extra in args.extras:
        deps += pyproject["project"]["optional-dependencies"][extra]

    min_deps = extract_min_deps(deps)

    print(" ".join(min_deps))


if __name__ == "__main__":
    main()
