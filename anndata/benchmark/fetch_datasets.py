import re
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).parent))

from benchmarks.datasets import DATASETS


def main():
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Setup datasets for benchmarking.")
    parser.add_argument(
        "--pattern",
        help="Regex to match dataset names against.",
        type=str,
        default=".*",
    )

    args = parser.parse_args()
    pat = re.compile(args.pattern)

    for dataset in DATASETS:
        if pat.match(dataset.name) is None:
            continue
        else:
            dataset.setup()


if __name__ == "__main__":
    main()
