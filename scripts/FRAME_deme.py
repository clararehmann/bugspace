import numpy as np, pandas as pd
import os, sys
import pickle
import argparse
from itertools import combinations

import dask
import dask.array as da
from dask.diagnostics.progress import ProgressBar
# silence some warnings
dask.config.set(**{'array.slicing.split_large_chunks': False})
import matplotlib
matplotlib.use('Agg')  # no display on the cluster
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import allel
import malariagen_data
import geopandas as gpd
from shapely.geometry import Polygon

# Vis, initialize_digraph, lambda_cv, fit_digraph, plot_digraph and
# genotype_array_to_FRAME are expected to come from ag3_FRAME
from ag3_FRAME import *
from plot_migration_rates import make_figure
from gambiae_grid import generate_hull

ag3 = malariagen_data.Ag3()
projection = ccrs.PlateCarree()


def get_query(sample_sets, min_lon, max_lon, min_lat, max_lat):
    sample_query = (f"taxon=='gambiae' and longitude < {max_lon} and longitude > {min_lon} "
                    f"and latitude < {max_lat} and latitude > {min_lat}")
    df_samples = ag3.sample_metadata(sample_sets=sample_sets, sample_query=sample_query)
    return df_samples, sample_query


def build_deme_grid(demefile, grid_path):
    """Write a FRAME/FEEMS triangle shapefile whose nodes are the deme centers.

    Returns the deme table (only demes with at least one edge), in the same
    order as the nodes in the shapefile.
    """
    nodes = pd.DataFrame(gpd.read_file(demefile, layer='demes', ignore_geometry=True))
    edge_df = pd.DataFrame(gpd.read_file(demefile, layer='edges', ignore_geometry=True))
    if nodes['deme_id'].duplicated().any():
        raise SystemExit('Duplicate deme IDs in the demes layer')

    ## map edges to node positions
    pos_of = {d: i for i, d in enumerate(nodes['deme_id'])}
    src = edge_df['source_deme'].map(pos_of)
    dst = edge_df['destination_deme'].map(pos_of)
    bad = src.isna() | dst.isna()
    if bad.any():
        print(f'WARNING: {bad.sum()} edge(s) reference unknown deme IDs; skipped')
    pairs = np.column_stack([src[~bad], dst[~bad]]).astype(int)
    pairs = pairs[pairs[:, 0] != pairs[:, 1]]  # drop self-loops

    ## make graph undirected: A->B and B->A become one edge
    pairs = np.unique(np.sort(pairs, axis=1), axis=0)

    ## keep demes with at least one edge, re-index from 0
    used = np.unique(pairs)
    if len(used) < len(nodes):
        print(f'WARNING: {len(nodes) - len(used)} deme(s) have no edges and are dropped')
    remap = {old: new for new, old in enumerate(used)}
    edges = np.vectorize(remap.get)(pairs)
    kept = nodes.iloc[used].reset_index(drop=True)
    node_pos = kept[['lon', 'lat']].to_numpy(dtype=float)

    ## triangles of mutually connected demes
    nbrs = {i: set() for i in range(len(node_pos))}
    for i, j in edges:
        nbrs[i].add(j)
        nbrs[j].add(i)
    tris = sorted({tuple(sorted((i, j, k)))
                   for i in nbrs for j, k in combinations(nbrs[i], 2) if k in nbrs[j]})
    tri_edges = {e for t in tris for e in combinations(t, 2)}
    missing = sum(tuple(e) not in tri_edges for e in edges)
    if missing:
        print(f'WARNING: {missing} edge(s) are in no triangle and will be missing from the FRAME graph')

    ## write shapefile
    t = np.array(tris)
    ids = kept['deme_id'].to_numpy()
    gpd.GeoDataFrame({'tri_id': np.arange(len(t)),
                      'v1_id': ids[t[:, 0]],
                      'v2_id': ids[t[:, 1]],
                      'v3_id': ids[t[:, 2]]},
                     geometry=[Polygon(node_pos[list(tr)]) for tr in tris],
                     crs='EPSG:4326').to_file(grid_path)
    print(f'Grid: {len({v for tr in tris for v in tr})} nodes, {len(tri_edges)} edges, '
          f'{len(tris)} triangles -> {grid_path}')
    return kept


parser = argparse.ArgumentParser()
parser.add_argument('--demefile', type=str, help='path to deme and edge geopackage')
parser.add_argument('--min_lon', type=float)
parser.add_argument('--max_lon', type=float)
parser.add_argument('--min_lat', type=float)
parser.add_argument('--max_lat', type=float)
parser.add_argument('--outpath', type=str)
parser.add_argument('--hull_from', choices=['samples', 'demes', 'both'], default='samples',
                    help='points the outer hull is built around; demes outside the hull '
                         'are dropped from the FRAME graph')
args = parser.parse_args()

# ---- samples and genotypes
sample_sets = [f'3.{i}' for i in range(17)]
sample_metadata, sample_query = get_query(sample_sets, args.min_lon, args.max_lon,
                                          args.min_lat, args.max_lat)
genotypes = ag3.snp_calls(
    region='3L:1000000-2000000',
    sample_query=sample_query,
    site_mask='gamb_colu',
    sample_query_options={'engine': 'python'}
)
samples = genotypes['sample_id'].values
genotypes = allel.GenotypeArray(genotypes['call_genotype'].data)
sample_metadata.set_index('sample_id', inplace=True)
sample_metadata = sample_metadata.reindex(samples)

# ---- deme grid
grid_path = f'{args.outpath}_grid.shp'
kept_demes = build_deme_grid(args.demefile, grid_path)

# ---- outer hull
# generate_hull is given a table with longitude/latitude columns; for 'demes' or
# 'both' this assumes it only needs those two columns.
deme_pts = kept_demes.rename(columns={'lon': 'longitude', 'lat': 'latitude'})[['longitude', 'latitude']]
sample_pts = sample_metadata[['longitude', 'latitude']]
hull_pts = {'samples': sample_metadata,
            'demes': deme_pts,
            'both': pd.concat([sample_pts, deme_pts], ignore_index=True)}[args.hull_from]
hull = generate_hull(hull_pts, args.outpath, x_fact=1.1, y_fact=1.1)

# ---- FRAME inputs
coordinates = sample_metadata[['longitude', 'latitude']].to_numpy()
genotypes = genotype_array_to_FRAME(genotypes)
spdigraph = initialize_digraph(genotypes, coordinates, f'{args.outpath}_hull.csv',
                               grid_path, digraph_class=BoundedSpatialDiGraph)

n_graph_nodes = len(spdigraph.node_pos)
print(f'FRAME graph has {n_graph_nodes} of {len(kept_demes)} demes '
      f'({spdigraph.n_observed_nodes} with samples)')
if n_graph_nodes < len(kept_demes):
    print('NOTE: some demes fell outside the hull and were dropped; '
          'use --hull_from both or demes to keep them')

N = min(10, spdigraph.n_observed_nodes - 1)

# ---- lambda cross-validation and fit
lamb_m_opt = lambda_cv(spdigraph, args.outpath, N=N)
print(f'Optimal lambda: {lamb_m_opt}')
spdigraph_fit = fit_digraph(spdigraph, lamb_m_opt)

fitted_graph_path = f'{args.outpath}_spdigraph.pkl'
with open(fitted_graph_path, 'wb') as fitted_graph_file:
    pickle.dump(spdigraph_fit, fitted_graph_file, protocol=pickle.HIGHEST_PROTOCOL)
print(f'Saved fitted spatial digraph: {fitted_graph_path}')

node_positions = np.asarray(spdigraph_fit.node_pos, dtype=float)
node_summary = pd.DataFrame({
    'node_id': np.arange(len(node_positions)),
    'longitude': node_positions[:, 0],
    'latitude': node_positions[:, 1],
    'haplotype_count': spdigraph_fit.n_hap,
})
# attach deme IDs by coordinate match (FRAME's node positions come from the grid)
id_at = {(round(x, 6), round(y, 6)): d
         for x, y, d in zip(kept_demes['lon'], kept_demes['lat'], kept_demes['deme_id'])}
node_summary.insert(1, 'deme_id', [id_at.get((round(x, 6), round(y, 6)))
                                   for x, y in node_positions])
nodes_path = f'{args.outpath}_nodes.csv'
node_summary.to_csv(nodes_path, index=False)

figure = make_figure(spdigraph_fit, node_summary)
html_path = f'{args.outpath}_migration.html'
figure.write_html(html_path, include_plotlyjs=True, full_html=True)
print(f'Saved interactive migration map: {html_path}')

# ---- plots
plot_digraph(spdigraph_fit, args.outpath, projection)

vis_kwargs = dict(projection=projection, edge_width=1,
                  edge_alpha=1, edge_zorder=100, sample_pt_size=20,
                  obs_node_size=5, sample_pt_color='black',
                  cbar_font_size=5, cbar_ticklabelsize=5,
                  cbar_bbox_to_anchor=(0.05, 0.2),
                  cbar_width='15%',
                  cbar_height='5%',
                  compass_bbox_to_anchor=(0, 0),
                  compass_font_size=5,
                  compass_radius=0.2,
                  compass_pad=0.2,
                  mutation_scale=8)

for mode, suffix in [('Difference', 'difference'), ('Full', 'full')]:
    fig = plt.figure(dpi=300)
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    v = Vis(ax, spdigraph_fit, **vis_kwargs)
    v.draw_migration_rates(ax, mode=mode, draw_map=True, draw_nodes=True)
    fig.savefig(f'{args.outpath}_{suffix}.png')
    plt.close(fig)
    print(f'Saved {mode.lower()} plot: {args.outpath}_{suffix}.png')