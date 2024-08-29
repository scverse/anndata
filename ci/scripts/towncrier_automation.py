#!python3
from __future__ import annotations

import argparse
import subprocess

from packaging.version import parse


def main():
    parser = argparse.ArgumentParser(
        prog="towncrier-automation",
        description=(
            "This script runs towncrier for a given version, "
            "creates a branch off of the current one, "
            "and then creates a PR into the original branch with the changes. "
            "The PR will be backported to main if the current branch is not main."
        ),
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
    if subprocess.run(
        ["towncrier", "build", f"--version={version}", "--yes"], check=False
    ).returncode:
        raise RuntimeError("Failed to build towncrier")

    # Check if we are on the main branch to know if we need to backport
    base_branch = subprocess.run(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    pr_description = "" if base_branch == "main" else "on-merge: backport to main"
    branch_name = f"release_notes_{args.version}"

    # Create a new branch + commit
    subprocess.run(["git", "checkout", "-b", branch_name], check=False)
    subprocess.run(["git", "add", "docs/release-notes"], check=False)
    pr_title = f"(chore): generate {version} release notes"
    subprocess.run(["git", "commit", "-m", pr_title], check=False)

    # Create a PR
    subprocess.run(
        [
            "gh",
            "pr",
            "create",
            f"--base={base_branch}",
            f"--head={branch_name}",
            f"--title={pr_title}",
            f"--body={pr_description}",
            "--dry-run" if args.dry_run else "",
        ],
        check=False,
    )

    # Enable auto-merge
    if not args.dry_run:
        subprocess.run(
            ["gh", "pr", "merge", branch_name, "--auto", "--squash"], check=False
        )
    else:
        print("Dry run, not merging")


if __name__ == "__main__":
    main()
