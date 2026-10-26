#!/usr/bin/env python3

"""Create an interactive map of fitted, directed FRAME migration rates."""

import argparse
import math
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.colors import sample_colorscale


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GRAPH_PATH = PROJECT_ROOT / "out" / "all_resolution_3_occupied_spdigraph.pkl"
DEFAULT_NODES_PATH = PROJECT_ROOT / "out" / "all_resolution_3_occupied_nodes.csv"
DEFAULT_OUTPUT_PATH = PROJECT_ROOT / "out" / "all_resolution_3_occupied_migration.html"


def read_nodes(path: Path, n_nodes: int) -> pd.DataFrame:
    nodes = pd.read_csv(path)
    required_columns = {"node_id", "longitude", "latitude", "haplotype_count"}
    missing_columns = required_columns.difference(nodes.columns)
    if missing_columns:
        raise ValueError(
            f"{path} is missing required column(s): {', '.join(sorted(missing_columns))}"
        )

    for column in required_columns:
        nodes[column] = pd.to_numeric(nodes[column], errors="raise")
    if nodes[list(required_columns)].isna().any().any():
        raise ValueError(f"{path} contains missing node data")
    if not np.equal(nodes["node_id"], np.floor(nodes["node_id"])).all():
        raise ValueError(f"{path} contains non-integer node IDs")
    nodes["node_id"] = nodes["node_id"].astype(int)
    nodes = nodes.sort_values("node_id").reset_index(drop=True)
    if not np.array_equal(nodes["node_id"].to_numpy(), np.arange(n_nodes)):
        raise ValueError(
            f"{path} must contain exactly one row for each graph node ID "
            f"from 0 to {n_nodes - 1}"
        )
    if not nodes["longitude"].between(-180, 180).all():
        raise ValueError(f"{path} contains invalid longitude values")
    if not nodes["latitude"].between(-90, 90).all():
        raise ValueError(f"{path} contains invalid latitude values")
    if (nodes["haplotype_count"] < 0).any():
        raise ValueError(f"{path} contains negative haplotype counts")
    return nodes


def edge_path(
    source_longitude: float,
    source_latitude: float,
    target_longitude: float,
    target_latitude: float,
) -> tuple[list[float | None], list[float | None]]:
    """Return a slightly curved line and an arrowhead in local map coordinates."""
    mean_latitude = math.radians((source_latitude + target_latitude) / 2)
    longitude_scale = math.cos(mean_latitude)
    start = np.array([source_longitude * longitude_scale, source_latitude])
    end = np.array([target_longitude * longitude_scale, target_latitude])
    delta = end - start
    distance = float(np.linalg.norm(delta))
    if distance == 0:
        return [source_longitude, target_longitude], [source_latitude, target_latitude]

    direction = delta / distance
    perpendicular = np.array([-direction[1], direction[0]])
    control = (start + end) / 2 + perpendicular * min(distance * 0.08, 0.15)
    t = np.linspace(0, 1, 13)
    curve = (
        ((1 - t) ** 2)[:, None] * start
        + (2 * (1 - t) * t)[:, None] * control
        + (t**2)[:, None] * end
    )

    arrow_direction = curve[-1] - curve[-2]
    arrow_direction /= np.linalg.norm(arrow_direction)
    arrow_perpendicular = np.array([-arrow_direction[1], arrow_direction[0]])
    arrow_length = min(distance * 0.2, 0.5)
    arrow_width = arrow_length * 0.55
    arrow_base = end - arrow_direction * arrow_length
    arrow_left = arrow_base + arrow_perpendicular * arrow_width
    arrow_right = arrow_base - arrow_perpendicular * arrow_width

    def to_longitude(local_x: float) -> float:
        return float(local_x / longitude_scale)

    longitudes = [to_longitude(point[0]) for point in curve]
    latitudes = [float(point[1]) for point in curve]
    longitudes.extend(
        [
            None,
            to_longitude(arrow_left[0]),
            target_longitude,
            None,
            to_longitude(arrow_right[0]),
            target_longitude,
        ]
    )
    latitudes.extend(
        [
            None,
            float(arrow_left[1]),
            target_latitude,
            None,
            float(arrow_right[1]),
            target_latitude,
        ]
    )
    return longitudes, latitudes


def make_figure(graph: object, nodes: pd.DataFrame) -> go.Figure:
    if not hasattr(graph, "M") or not hasattr(graph, "node_pos"):
        raise TypeError("The pickle does not contain a fitted FRAME spatial digraph")

    positions = np.asarray(graph.node_pos, dtype=float)
    if positions.shape != (len(nodes), 2):
        raise ValueError("The node coordinate CSV does not match the graph dimensions")
    csv_positions = nodes[["longitude", "latitude"]].to_numpy(dtype=float)
    if not np.allclose(positions, csv_positions, atol=1e-6, rtol=0):
        raise ValueError("The node coordinate CSV does not match the fitted graph")

    matrix = graph.M.tocoo(copy=True)
    if matrix.shape != (len(nodes), len(nodes)):
        raise ValueError("The fitted migration matrix does not match the node data")
    if len(matrix.data) == 0:
        raise ValueError("The fitted graph contains no migration rates")
    if (matrix.row == matrix.col).any():
        raise ValueError("The fitted graph contains self-migration rates")
    if not np.isfinite(matrix.data).all() or (matrix.data <= 0).any():
        raise ValueError("The fitted graph contains non-positive or invalid migration rates")

    log_rates = np.log10(matrix.data)
    log_relative_rates = log_rates - log_rates.mean()
    color_limit = max(float(np.max(np.abs(log_relative_rates))), 0.01)
    color_scale = "RdBu"
    figure = go.Figure()

    order = np.lexsort((matrix.col, matrix.row))
    for index in order:
        source = int(matrix.row[index])
        target = int(matrix.col[index])
        source_longitude, source_latitude = csv_positions[source]
        target_longitude, target_latitude = csv_positions[target]
        relative_log_rate = float(log_relative_rates[index])
        color_position = (relative_log_rate + color_limit) / (2 * color_limit)
        color = sample_colorscale(color_scale, [color_position])[0]
        longitudes, latitudes = edge_path(
            source_longitude,
            source_latitude,
            target_longitude,
            target_latitude,
        )

        if -300 <= relative_log_rate <= 300:
            relative_rate_text = f"{10**relative_log_rate:.4g}x"
        else:
            relative_rate_text = f"10^{relative_log_rate:.2f}x"
        hover_text = (
            f"From node {source} to node {target}"
            f"<br>Migration rate: {matrix.data[index]:.6g}"
            f"<br>Relative rate (geometric mean = 1): {relative_rate_text}"
            f"<br>log10(relative rate): {relative_log_rate:+.3f}"
        )
        figure.add_trace(
            go.Scattergeo(
                lon=longitudes,
                lat=latitudes,
                mode="lines",
                line=dict(
                    color=color,
                    width=1.2 + 1.8 * abs(relative_log_rate) / color_limit,
                ),
                opacity=0.85,
                hovertemplate=f"{hover_text}<extra></extra>",
                showlegend=False,
            )
        )

    figure.add_trace(
        go.Scattergeo(
            lon=[None, None],
            lat=[None, None],
            mode="markers",
            marker=dict(
                size=0,
                color=[-color_limit, color_limit],
                cmin=-color_limit,
                cmax=color_limit,
                colorscale=color_scale,
                showscale=True,
                colorbar=dict(
                    title=dict(text="log10(relative rate)"),
                    tickvals=[-color_limit, 0, color_limit],
                    ticktext=[
                        f"{-color_limit:.2f}",
                        "0",
                        f"+{color_limit:.2f}",
                    ],
                ),
            ),
            hoverinfo="skip",
            showlegend=False,
        )
    )

    max_haplotypes = max(float(nodes["haplotype_count"].max()), 1.0)
    node_sizes = 7 + 13 * np.sqrt(nodes["haplotype_count"] / max_haplotypes)
    figure.add_trace(
        go.Scattergeo(
            lon=nodes["longitude"],
            lat=nodes["latitude"],
            mode="markers",
            marker=dict(
                size=node_sizes,
                color="#292929",
                line=dict(color="white", width=1),
                opacity=0.9,
            ),
            customdata=np.column_stack(
                [
                    nodes["node_id"],
                    nodes["haplotype_count"],
                    nodes["longitude"],
                    nodes["latitude"],
                ]
            ),
            hovertemplate=(
                "Node %{customdata[0]:.0f}"
                "<br>Haplotypes: %{customdata[1]:.0f}"
                "<br>Longitude: %{customdata[2]:.4f}"
                "<br>Latitude: %{customdata[3]:.4f}"
                "<extra></extra>"
            ),
            name="Nodes",
            showlegend=False,
        )
    )
    figure.update_layout(
        title="Relative migration rates between spatial nodes",
        height=850,
        margin=dict(l=20, r=90, t=60, b=20),
        geo=dict(
            projection_type="natural earth",
            fitbounds="locations",
            showland=True,
            landcolor="#f3f4f1",
            showocean=True,
            oceancolor="#eaf2f8",
            showcoastlines=True,
            coastlinecolor="#858585",
            showcountries=True,
            countrycolor="#b8b8b8",
            showframe=False,
        ),
        uirevision="migration-rates",
    )
    return figure


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--graph",
        type=Path,
        default=DEFAULT_GRAPH_PATH,
        help="Trusted fitted FRAME graph pickle (loading a pickle can execute code)",
    )
    parser.add_argument(
        "--nodes",
        type=Path,
        default=DEFAULT_NODES_PATH,
        help="CSV node summary written by a FRAME graph-fitting script",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT_PATH,
        help="Output path for the standalone interactive HTML map",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()
    graph_size_gib = arguments.graph.stat().st_size / (1024**3)
    print(
        f"Loading trusted fitted graph from {arguments.graph} "
        f"({graph_size_gib:.1f} GiB on disk); deserialization may require "
        "substantial memory."
    )
    with arguments.graph.open("rb") as graph_file:
        graph = pickle.load(graph_file)

    nodes = read_nodes(arguments.nodes, graph.M.shape[0])
    figure = make_figure(graph, nodes)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    figure.write_html(arguments.output, include_plotlyjs=True, full_html=True)
    print(f"Saved interactive migration map to {arguments.output}")


if __name__ == "__main__":
    main()
