#!python3
from __future__ import annotations

import argparse
import subprocess

from packaging.version import parse


def main():
    parser = argparse.ArgumentParser(
        prog="towncrier-automation",
        description="This script runs townncrier for a given version, creates a branch off of the current one, and then creates a PR into the original branch with the changes.  The PR will be backported to main if the current branch is not main.",
        usage="python towncrier_automation.py --version <version> [--dry-run]",
    )
    parser.add_argument("--version", type=str, help="What the new version is")
    parser.add_argument(
        "--dry-run",
        help="Whether or not to dry-run the actual creation of the pull request",
        action="store_true",
    )
    args = parser.parse_args()
    version = args.version
    parse(args.version)

    # Run towncrier
    if subprocess.call(["towncrier", "build", f"--version={version}", "--yes"]):
        raise RuntimeError("Failed to build towncrier")

    # Check if we are on the main branch to know if we need to backport
    base_branch = subprocess.run(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    pr_description = ""
    if base_branch != "main":
        pr_description = "on-merge: backport to main"
    branch_name = f"release_notes_{args.version}"

    # Create a new branch + commit
    subprocess.call(["git", "checkout", "-b", branch_name])
    subprocess.call(["git", "add", "docs/release-notes"])
    pr_title = f"(chore): generate {version} release notes"
    subprocess.call(["git", "commit", "-m", pr_title])

    # Create a PR
    subprocess.call(
        [
            "gh",
            "pr",
            "create",
            "--base",
            base_branch,
            "--head",
            branch_name,
            "--title",
            pr_title,
            "--body",
            pr_description,
            "dry-run" if args.dry_run else "",
        ]
    )

    # Enable auto-merge
    subprocess.call(["gh", "pr", "merge", branch_name, "--auto", "--squash"])


if __name__ == "__main__":
    main()
