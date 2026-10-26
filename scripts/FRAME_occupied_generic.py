import os

import allel
import cartopy.crs as ccrs
import malariagen_data
import numpy as np
import pandas as pd
from frame.spatial_digraph import SpatialDiGraph
from frame.utils import prepare_graph_inputs
from matplotlib import pyplot as plt
from pyproj import Transformer
from scipy.spatial import Delaunay, QhullError
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import minimum_spanning_tree

from ag3_FRAME import fit_digraph, genotype_array_to_FRAME, lambda_cv, plot_digraph
from gambiae_grid import generate_h3_grid, generate_hull, get_sampling_radius
from frame.visualization import Vis
import argparse

SAMPLE_SETS = [f'3.{i}' for i in range(17)]
REGION = '3L:1000000-2000000'


class OccupiedSpatialDiGraph(SpatialDiGraph):
    def fit(self, *args, **kwargs):
        kwargs['lb'] = max(kwargs.get('lb', -10.0), -10.0)
        kwargs['ub'] = min(kwargs.get('ub', 10.0), 10.0)
        return super().fit(*args, **kwargs)


def connect_occupied_nodes(node_coordinates, lat_center, lon_center):
    """Build a connected, undirected Delaunay graph over occupied node positions."""
    node_count = len(node_coordinates)
    if node_count < 2:
        raise ValueError('FRAME needs at least two sample-occupied nodes.')

    local_crs = (
        f'+proj=aeqd +lat_0={lat_center} +lon_0={lon_center} '
        '+datum=WGS84 +units=m +no_defs'
    )
    transformer = Transformer.from_crs('EPSG:4326', local_crs, always_xy=True)
    x, y = transformer.transform(node_coordinates[:, 0], node_coordinates[:, 1])
    projected_coordinates = np.column_stack((x, y))

    edges = set()
    if node_count >= 3:
        try:
            triangles = Delaunay(projected_coordinates).simplices
            for triangle in triangles:
                edges.update(
                    tuple(sorted(edge))
                    for edge in (
                        (triangle[0], triangle[1]),
                        (triangle[1], triangle[2]),
                        (triangle[2], triangle[0]),
                    )
                )
        except QhullError:
            pass

    # The Euclidean minimum spanning tree connects collinear or degenerate sites.
    spanning_tree = minimum_spanning_tree(cdist(projected_coordinates, projected_coordinates))
    spanning_tree = spanning_tree.tocoo()
    edges.update(
        tuple(sorted((int(first), int(second))))
        for first, second in zip(spanning_tree.row, spanning_tree.col)
    )

    return np.asarray(sorted(edges), dtype=int).reshape(-1, 2)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--year', type=int)
    parser.add_argument('--resolution', type=int)
    args = parser.parse_args()
    YEAR = args.year
    RESOLUTION = args.resolution

    ag3 = malariagen_data.Ag3()
    sample_query = f"year=={YEAR} & taxon=='gambiae'"
    sample_metadata = ag3.sample_metadata(
        sample_sets=SAMPLE_SETS,
        sample_query=sample_query,
    )
    genotype_data = ag3.snp_calls(
        region=REGION,
        sample_query=sample_query,
        site_mask='gamb_colu',
        sample_query_options={'engine': 'python'},
    )

    sample_ids = genotype_data['sample_id'].values
    genotypes = allel.GenotypeArray(genotype_data['call_genotype'].data)
    sample_metadata = sample_metadata.set_index('sample_id').reindex(sample_ids)
    valid_coordinates = sample_metadata[['longitude', 'latitude']].notna().all(axis=1).to_numpy()
    sample_metadata = sample_metadata.loc[valid_coordinates]
    genotypes = genotypes[:, valid_coordinates, :]
    coordinates = sample_metadata[['longitude', 'latitude']].to_numpy()
    genotypes = genotype_array_to_FRAME(genotypes)

    lat_center, lon_center, max_distance = get_sampling_radius(sample_metadata)
    projection = ccrs.EquidistantConic(
        central_longitude=lon_center,
        central_latitude=lat_center,
    )

    os.makedirs('data', exist_ok=True)
    os.makedirs('out', exist_ok=True)
    gridpath = f'data/{YEAR}_occupied'
    generate_hull(sample_metadata, gridpath, x_fact=1.5, y_fact=1.5)
    generate_h3_grid(
        lat_center,
        lon_center,
        max_distance,
        resolution=RESOLUTION,
        outpath=gridpath,
    )

    hull = np.genfromtxt(f'{gridpath}_hull.csv', delimiter=',')
    _, _, mesh_nodes, sample_node_ids = prepare_graph_inputs(
        coord=coordinates,
        ggrid=f'{gridpath}_grid_resolution_{RESOLUTION}.shp',
        buffer=0,
        outer=hull,
    )
    occupied_node_ids = np.unique(sample_node_ids)
    occupied_node_coordinates = mesh_nodes[occupied_node_ids]
    occupied_edges = connect_occupied_nodes(
        occupied_node_coordinates,
        lat_center,
        lon_center,
    )

    spdigraph = OccupiedSpatialDiGraph(
        genotypes,
        coordinates,
        occupied_node_coordinates,
        occupied_edges,
    )
    if spdigraph.n_observed_nodes != len(occupied_node_coordinates):
        raise RuntimeError(
            'FRAME reassigned samples differently from the mesh mapping; '
            'not every retained node is occupied.'
        )
    if spdigraph.n_observed_nodes < 3:
        raise ValueError('At least three occupied nodes are needed for cross-validation.')

    print(
        f'Resolution {RESOLUTION}: {len(mesh_nodes)} mesh nodes, '
        f'{len(occupied_node_coordinates)} occupied nodes, '
        f'{len(occupied_edges)} undirected edges, '
        f'{spdigraph.n_observed_nodes} observed nodes'
    )

    outpath = f'out/{YEAR}_resolution_{RESOLUTION}_occupied'
    n_folds = min(10, spdigraph.n_observed_nodes - 1)
    lamb_m_opt = lambda_cv(spdigraph, outpath, N=n_folds)
    print(f'Optimal lambda: {lamb_m_opt}')

    spdigraph_fit = fit_digraph(spdigraph, lamb_m_opt)
    plot_digraph(spdigraph_fit, outpath, projection)

    undirected_degree = np.bincount(
        occupied_edges.ravel(),
        minlength=len(occupied_node_coordinates),
    )
    node_summary = pd.DataFrame({
        'node_id': np.arange(len(occupied_node_coordinates)),
        'longitude': occupied_node_coordinates[:, 0],
        'latitude': occupied_node_coordinates[:, 1],
        'haplotype_count': spdigraph.n_hap,
        'undirected_degree': undirected_degree,
    })
    node_summary.to_csv(f'{outpath}_nodes.csv', index=False)
    pd.DataFrame(occupied_edges, columns=['node_a', 'node_b']).to_csv(
        f'{outpath}_edges.csv',
        index=False,
    )

    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    visualizer = Vis(
        ax,
        spdigraph_fit,
        projection=projection,
        edge_width=1,
        edge_alpha=1,
        edge_zorder=100,
        sample_pt_size=20,
        obs_node_size=5,
        sample_pt_color='black',
        cbar_font_size=5,
        cbar_ticklabelsize=5,
        cbar_bbox_to_anchor=(0.05, 0.2),
        cbar_width='15%',
        cbar_height='5%',
        compass_bbox_to_anchor=(0, 0),
        compass_font_size=5,
        compass_radius=0.2,
        compass_pad=0.2,
        mutation_scale=8,
    )
    visualizer.draw_migration_rates(
        ax,
        mode='Difference',
        draw_map=True,
        draw_nodes=True,
    )
    fig.savefig(f'{outpath}_difference.png')
    plt.close(fig)

    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    visualizer = Vis(ax, spdigraph_fit, projection=projection)
    visualizer.draw_migration_rates(
        ax,
        mode='Full',
        draw_map=True,
        draw_nodes=True,
    )
    fig.savefig(f'{outpath}_full.png')
    plt.close(fig)


if __name__ == '__main__':
    main()