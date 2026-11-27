#!/usr/bin/env python3
"""Find entries in `[tool.ruff.lint.extend-per-file-ignores]` that are no longer needed.

That table lists the missing-annotation (ANN) codes each file still violates.
`ruff check` already fails for violations that aren’t listed,
this script fails for listed codes that aren’t violated anymore,
so the list only ever shrinks.
"""

from __future__ import annotations

import json
import subprocess
import sys
import tomllib
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).parent.parent.parent
TABLE = "tool.ruff.lint.extend-per-file-ignores"


def main() -> int:
    ignores: dict[str, list[str]] = tomllib.loads((ROOT / "pyproject.toml").read_text())
    for key in TABLE.split("."):
        ignores = ignores[key]  # type: ignore[assignment]  # nested tables

    # which of the listed codes each file actually violates, without the table applied.
    # `--config` overrides can only add to the table, so skip the config file entirely
    # (none of its other settings affect ANN rules).
    codes = sorted({c for cs in ignores.values() for c in cs})
    proc = subprocess.run(
        [
            *("ruff", "check", "--isolated", "--exit-zero", "--output-format=json"),
            *(f"--select={','.join(codes)}", "."),
        ],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    )
    violated: defaultdict[Path, set[str]] = defaultdict(set)
    for diagnostic in json.loads(proc.stdout):
        violated[Path(diagnostic["filename"]).relative_to(ROOT)].add(diagnostic["code"])

    problems: list[str] = []
    for pattern, listed in ignores.items():
        needed = set().union(
            *(cs for path, cs in violated.items() if path.full_match(pattern))
        )
        if not needed:
            problems.append(f"{pattern}: fully annotated, remove its entry")
        elif stale := sorted(set(listed) - needed):
            problems.append(f"{pattern}: remove {', '.join(stale)}")

    if problems:
        print(f"Stale entries in [{TABLE}]:", *problems, sep="\n  ", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
