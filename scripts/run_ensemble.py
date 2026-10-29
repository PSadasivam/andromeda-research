"""Run a Monte Carlo ensemble and write results to runs/<name>/.

Example:
    python scripts/run_ensemble.py --config wu2026 --n 10000 --workers 12
"""

import argparse
import json
from pathlib import Path

from andromeda.ensemble import run_ensemble, save, summarize

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", choices=["sawala2025", "wu2026"], required=True)
    parser.add_argument("--galaxies", default="MW,M31,LMC,M33")
    parser.add_argument("--n", type=int, default=10000)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--friction", default="proportional", choices=["proportional", "hierarchical", "none"])
    parser.add_argument("--softening", type=float, default=20.0)
    parser.add_argument("--merger-threshold", type=float, default=20.0)
    parser.add_argument("--name", help="output folder name under runs/")
    args = parser.parse_args()

    galaxies = args.galaxies.split(",")
    name = args.name or f"{args.config}_{'-'.join(galaxies)}_n{args.n}"
    result = run_ensemble(args.config, galaxies, args.n, workers=args.workers, friction=args.friction,
                          softening=args.softening, merger_threshold=args.merger_threshold)
    summary = summarize(result)
    save(result, summary, ROOT / "runs" / name)
    print(json.dumps({"runtime_seconds": result["metadata"]["runtime_seconds"], **summary},
                     indent=2, default=float))


if __name__ == "__main__":
    main()
