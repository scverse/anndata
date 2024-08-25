from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import TypeAlias

    Version: TypeAlias = tuple[int, int] | tuple[int, int, int]

MARKER_F = "<!-- towncrier release notes start: feature -->"
MARKER_P = "<!-- towncrier release notes start: patch -->"

# A new minor release needs a new header
TEMPLATE_F = f"""\
{MARKER_F}

(v{{version}})=
## Version {{version}}

{MARKER_P}\
"""

# A new patch release needs just an include directive
TEMPLATE_P = f"""\
{MARKER_P}

```{{{{include}}}} /release-notes/{{version}}.md
```\
"""


def checked_replace(rel_notes: str, marker: str, replacement: str) -> str:
    if marker not in rel_notes:
        msg = rf"Missing {marker!r} in release notes index!"
        raise RuntimeError(msg)
    return rel_notes.replace(marker, replacement)


def insert_patch(rel_notes: str, version: Version) -> str:
    r"""Insert a new patch version into `rel_notes`.

    >>> rel_notes = f"START\n\n{MARKER_P}\n\nEND"
    >>> print(insert_patch(rel_notes, version=(1, 0, 1)))
    START
    <BLANKLINE>
    <!-- towncrier release notes start: patch -->
    <BLANKLINE>
    ```{include} /release-notes/1.0.1.md
    ```
    <BLANKLINE>
    END
    """
    assert len(version) == 3
    return checked_replace(
        rel_notes, MARKER_P, TEMPLATE_P.format(version=".".join(map(str, version)))
    )


def insert_feature(rel_notes: str, version: Version) -> str:
    r"""Insert a new feature version into `rel_notes`.

    >>> rel_notes = f"START\n\n{MARKER_F}\n\nOLD:\n\n{MARKER_P}\n\nEND"
    >>> print(insert_feature(rel_notes, version=(1, 1)))
    START
    <BLANKLINE>
    <!-- towncrier release notes start: feature -->
    <BLANKLINE>
    (v1.1)=
    ## Version 1.1
    <BLANKLINE>
    <!-- towncrier release notes start: patch -->
    <BLANKLINE>
    ```{include} /release-notes/1.1.0.md
    ```
    <BLANKLINE>
    OLD:
    <BLANKLINE>
    END
    """
    assert len(version) == 2
    # delete old patch marker
    rel_notes = checked_replace(rel_notes, f"{MARKER_P}\n\n", "")
    # replace feature marker with header and new patch marker
    rel_notes = checked_replace(
        rel_notes, MARKER_F, TEMPLATE_F.format(version=".".join(map(str, version)))
    )
    # insert include for new .0 patch version
    rel_notes = insert_patch(rel_notes, (*version, 0))
    return rel_notes
