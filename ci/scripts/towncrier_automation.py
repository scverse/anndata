#!/usr/bin/env python3
# /// script
# dependencies = [ "towncrier", "packaging" ]
# ///
from __future__ import annotations

import argparse
import re
import subprocess
from functools import cache
from typing import TYPE_CHECKING

from packaging.version import Version

if TYPE_CHECKING:
    from collections.abc import Sequence


class BumpVersion(Version):
    def __init__(self, version: str) -> None:
        super().__init__(version)

        if len(self.release) != 3:
            msg = f"{version} must contain major, minor, and patch version."
            raise argparse.ArgumentTypeError(msg)

        base_branch = get_base_branch()
        patch_branch_pattern = re.compile(r"\d+\.\d+\.x")
        if self.micro != 0 and not patch_branch_pattern.fullmatch(base_branch):
            msg = (
                f"{version} is a patch release, but "
                f"you are trying to release from a non-patch release branch: {base_branch}."
            )
            raise argparse.ArgumentTypeError(msg)

        if self.micro == 0 and base_branch != "main":
            msg = (
                f"{version} is a minor or major release, "
                f"but you are trying to release not from main: {base_branch}."
            )
            raise argparse.ArgumentTypeError(msg)


class Args(argparse.Namespace):
    version: BumpVersion
    dry_run: bool


def parse_args(argv: Sequence[str] | None = None) -> Args:
    parser = argparse.ArgumentParser(
        prog="towncrier-automation",
        description=(
            "This script runs towncrier for a given version, "
            "creates a branch off of the current one, "
            "and then creates a PR into the original branch with the changes. "
            "The PR will be backported to main if the current branch is not main."
        ),
    )
    parser.add_argument(
        "version",
        type=BumpVersion,
        help=(
            "The new version for the release must have at least three parts, like `major.minor.patch` and no `major.minor`. "
            "It can have a suffix like `major.minor.patch.dev0` or `major.minor.0rc1`."
        ),
    )
    parser.add_argument(
        "--dry-run",
        help="Whether or not to dry-run the actual creation of the pull request",
        action="store_true",
    )
    args = parser.parse_args(argv, Args())
    return args


def main(argv: Sequence[str] | None = None) -> None:
    args = parse_args(argv)

    # Run towncrier
    subprocess.run(
        ["towncrier", "build", f"--version={args.version}", "--yes"], check=True
    )

    # Check if we are on the main branch to know if we need to backport
    base_branch = get_base_branch()
    pr_description = "" if base_branch == "main" else "@meeseeksdev backport to main"
    branch_name = f"release_notes_{args.version}"

    # Create a new branch + commit
    subprocess.run(["git", "switch", "-c", branch_name], check=True)
    subprocess.run(["git", "add", "docs/release-notes"], check=True)
    pr_title = f"(chore): generate {args.version} release notes"
    subprocess.run(["git", "commit", "-m", pr_title], check=True)

    # push
    if not args.dry_run:
        subprocess.run(
            ["git", "push", "--set-upstream", "origin", branch_name], check=True
        )
    else:
        print("Dry run, not pushing")

    # Create a PR
    subprocess.run(
        [
            "gh",
            "pr",
            "create",
            f"--base={base_branch}",
            f"--title={pr_title}",
            f"--body={pr_description}",
            "--label=skip-gpu-ci",
            *(["--label=no milestone"] if base_branch == "main" else []),
            *(["--dry-run"] if args.dry_run else []),
        ],
        check=True,
    )

    # Enable auto-merge
    if not args.dry_run:
        subprocess.run(
            ["gh", "pr", "merge", branch_name, "--auto", "--squash"], check=True
        )
    else:
        print("Dry run, not merging")


@cache
def get_base_branch():
    return subprocess.run(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()


if __name__ == "__main__":
    main()
