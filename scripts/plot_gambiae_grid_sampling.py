#!/usr/bin/env python3

"""Plot gambiae grid coordinates over sampling locations."""

from argparse import ArgumentParser
from pathlib import Path
from warnings import warn

import matplotlib.pyplot as plt
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GRID_PATH = PROJECT_ROOT / "data" / "vo_agam_release" / "gambiae_grid_coords.csv"
DEFAULT_SAMPLING_PATH = PROJECT_ROOT / "data" / "vo_agam_release" / "gambiae_latlon.txt"
DEFAULT_OUTPUT_PATH = PROJECT_ROOT / "gambiae_grid_sampling_locations.png"


def read_coordinates(path: Path, column_names: list[str]) -> pd.DataFrame:
    coordinates = pd.read_csv(
        path,
        header=None,
        names=column_names,
        usecols=[0, 1],
        dtype="string",
    )
    blank = coordinates.isna().all(axis=1)
    if blank.any():
        warn(f"{path}: skipping {int(blank.sum())} blank coordinate row(s)")
        coordinates = coordinates.loc[~blank].copy()
    coordinates = coordinates.apply(pd.to_numeric, errors="coerce")
    valid = coordinates.notna().all(axis=1)
    valid &= coordinates["latitude"].between(-90, 90)
    valid &= coordinates["longitude"].between(-180, 180)
    if not valid.all():
        invalid_count = int((~valid).sum())
        raise ValueError(f"{path} contains {invalid_count} invalid coordinate row(s)")
    if coordinates.empty:
        raise ValueError(f"{path} contains no coordinate rows")
    return coordinates


def parse_arguments() -> object:
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        "--grid",
        type=Path,
        default=DEFAULT_GRID_PATH,
        help="CSV containing grid coordinates in latitude,longitude order",
    )
    parser.add_argument(
        "--sampling",
        type=Path,
        default=DEFAULT_SAMPLING_PATH,
        help="Text file containing sampling coordinates in latitude,longitude order",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_PATH,
        help="Path for the output PNG",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    grid = read_coordinates(arguments.grid, ["latitude", "longitude"])
    sampling = read_coordinates(arguments.sampling, ["latitude", "longitude"])

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(figsize=(10, 8), constrained_layout=True)
    axes.scatter(
        grid["longitude"],
        grid["latitude"],
        s=18,
        color="tab:blue",
        alpha=0.75,
        label=f"Grid coordinates ({len(grid):,})",
        zorder=1,
    )
    axes.scatter(
        sampling["longitude"],
        sampling["latitude"],
        s=8,
        color="tab:orange",
        alpha=0.35,
        edgecolors="none",
        label=f"Sampling locations ({len(sampling):,})",
        zorder=2,
    )
    axes.set_xlabel("Longitude")
    axes.set_ylabel("Latitude")
    axes.set_title("Anopheles gambiae grid and sampling locations")
    axes.grid(True, linestyle=":", linewidth=0.7, alpha=0.7)
    axes.legend(frameon=True)
    figure.savefig(arguments.output, dpi=200)
    plt.close(figure)
    print(f"Saved {arguments.output}")


if __name__ == "__main__":
    main()