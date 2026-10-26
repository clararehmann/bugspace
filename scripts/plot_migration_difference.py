#!/usr/bin/env python3

"""Create an interactive map of asymmetric FRAME migration rates."""

import argparse
import math
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.colors import sample_colorscale

from plot_migration_rates import edge_path, read_nodes


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GRAPH_PATH = PROJECT_ROOT / "out" / "all_resolution_3_occupied_spdigraph.pkl"
DEFAULT_NODES_PATH = PROJECT_ROOT / "out" / "all_resolution_3_occupied_nodes.csv"
DEFAULT_OUTPUT_PATH = PROJECT_ROOT / "out" / "all_resolution_3_occupied_migration_difference.html"


def make_figure(graph: object, nodes: pd.DataFrame) -> go.Figure:
    if not hasattr(graph, "M") or not hasattr(graph, "node_pos"):
        raise TypeError("The pickle does not contain a fitted FRAME spatial digraph")

    positions = np.asarray(graph.node_pos, dtype=float)
    if positions.shape != (len(nodes), 2):
        raise ValueError("The node coordinate CSV does not match the graph dimensions")
    csv_positions = nodes[["longitude", "latitude"]].to_numpy(dtype=float)
    if not np.allclose(positions, csv_positions, atol=1e-6, rtol=0):
        raise ValueError("The node coordinate CSV does not match the fitted graph")

    matrix = graph.M.tocsr(copy=True)
    matrix.sum_duplicates()
    if matrix.shape != (len(nodes), len(nodes)):
        raise ValueError("The fitted migration matrix does not match the node data")
    if not np.isfinite(matrix.data).all() or (matrix.data <= 0).any():
        raise ValueError("The fitted graph contains non-positive or invalid migration rates")
    if matrix.diagonal().any():
        raise ValueError("The fitted graph contains self-migration rates")

    asymmetric_edges = []
    row, column = matrix.nonzero()
    node_pairs = sorted({
        (min(int(source), int(target)), max(int(source), int(target)))
        for source, target in zip(row, column)
    })
    for source, target in node_pairs:
        forward_rate = float(matrix[source, target])
        reverse_rate = float(matrix[target, source])
        if forward_rate <= 0 or reverse_rate <= 0:
            raise ValueError(
                f"Missing reciprocal migration rate for node pair {source}, {target}"
            )
        log_rate_difference = math.log10(forward_rate) - math.log10(reverse_rate)
        if log_rate_difference == 0:
            continue
        if log_rate_difference > 0:
            stronger_source, stronger_target = source, target
        else:
            stronger_source, stronger_target = target, source
        asymmetric_edges.append(
            (
                stronger_source,
                stronger_target,
                forward_rate if stronger_source == source else reverse_rate,
                reverse_rate if stronger_source == source else forward_rate,
                abs(log_rate_difference),
            )
        )

    color_limit = max(
        (edge[4] for edge in asymmetric_edges),
        default=0.0,
    )
    figure = go.Figure()
    if asymmetric_edges:
        color_limit = max(color_limit, 0.01)
        for source, target, stronger_rate, weaker_rate, log_ratio in asymmetric_edges:
            source_longitude, source_latitude = csv_positions[source]
            target_longitude, target_latitude = csv_positions[target]
            color = sample_colorscale(
                "YlOrRd",
                [log_ratio / color_limit],
            )[0]
            longitudes, latitudes = edge_path(
                source_longitude,
                source_latitude,
                target_longitude,
                target_latitude,
            )
            ratio_text = (
                f"{10**log_ratio:.4g}x"
                if log_ratio <= 300
                else f"10^{log_ratio:.2f}x"
            )
            hover_text = (
                f"Stronger migration: node {source} to node {target}"
                f"<br>Migration rate: {stronger_rate:.6g}"
                f"<br>Reverse migration rate: {weaker_rate:.6g}"
                f"<br>Stronger/weaker rate: {ratio_text}"
                f"<br>log10(stronger/weaker): {log_ratio:.3f}"
            )
            figure.add_trace(
                go.Scattergeo(
                    lon=longitudes,
                    lat=latitudes,
                    mode="lines",
                    line=dict(
                        color=color,
                        width=1.2 + 1.8 * log_ratio / color_limit,
                    ),
                    opacity=0.9,
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
                    color=[0, color_limit],
                    cmin=0,
                    cmax=color_limit,
                    colorscale="YlOrRd",
                    showscale=True,
                    colorbar=dict(
                        title=dict(text="log10(stronger/weaker rate)"),
                        tickvals=[0, color_limit / 2, color_limit],
                        ticktext=["0", f"{color_limit / 2:.2f}", f"{color_limit:.2f}"],
                    ),
                ),
                hoverinfo="skip",
                showlegend=False,
            )
        )
    else:
        figure.add_annotation(
            text="No asymmetric migration rates found",
            x=0.5,
            y=0.02,
            xref="paper",
            yref="paper",
            showarrow=False,
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
        title="Asymmetric migration between spatial nodes",
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
        uirevision="migration-difference",
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
    print(f"Saved interactive migration difference map to {arguments.output}")


if __name__ == "__main__":
    main()
