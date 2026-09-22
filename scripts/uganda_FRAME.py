import numpy as np, pandas as pd
import os, sys
import dask
import dask.array as da
from dask.diagnostics.progress import ProgressBar
# silence some warnings
dask.config.set(**{'array.slicing.split_large_chunks': False})
import cartopy.crs as ccrs
import allel
import malariagen_data
ag3=malariagen_data.Ag3()

from ag3_popgen import get_zarr_genotypes
from gambiae_grid import *
from ag3_FRAME import *

def get_month_query(month, sample_sets):
    # generate query string for 2017 samples from Uganda in that month
    df_samples = ag3.sample_metadata(sample_sets=sample_sets, sample_query=f"country=='Uganda' & year==2017 & taxon=='gambiae' & month=={month}")
    samples = df_samples['sample_id'].tolist()
    sample_query = '|'.join([f"sample_id=='{s}'" for s in samples])
    return df_samples, sample_query

sample_sets = [f'3.{i}' for i in range(17)]
months=[3,4,5,6,7,8,9,11]
month=sys.argv[1]
sample_metadata, sample_query = get_month_query(month, sample_sets)
genotypes = ag3.snp_calls(
                        region='3L',
                        sample_query=sample_query,
                        site_mask='gamb_colu',
                        sample_query_options={'engine':'python'}
                        )
samples = genotypes['sample_id'].values
genotypes = allel.GenotypeArray(genotypes['call_genotype'].data)
sample_metadata.set_index('sample_id', inplace=True)
sample_metadata = sample_metadata.reindex(samples)
sample_metadata.to_csv(f"out/Uganda2017/{month}.csv")
# make sampling grid and edges
lat_c, lon_c, max_d = get_sampling_radius(sample_metadata)
projection = ccrs.EquidistantConic(central_longitude = lon_c, central_latitude = lat_c)
resolution=4
gridpath=f'data/Uganda2017/{month}'
grid = generate_h3_grid(lat_c, lon_c, max_d, resolution=resolution, outpath=gridpath)
hull = generate_hull(sample_metadata, gridpath, x_fact=1.5, y_fact=1.5)
# prep coordinates for FRAME
coordinates = sample_metadata[['longitude','latitude']].to_numpy()
# prep genotypes for FRAME
print(sample_query)

genotypes = genotype_array_to_FRAME(genotypes)

# make spatial digraph
spdigraph = initialize_digraph(genotypes, coordinates, f'{gridpath}_hull.csv', f'{gridpath}_grid_resolution_{resolution}.shp')

if len(sample_metadata) < 10:
    N = len(sample_metadata)-1
else:
    N=10

# lambda cross-validation
lamb_m_opt = lambda_cv(spdigraph, f'out/Uganda2017/{month}', N=10)

# fit digraph
spdigraph_fit = fit_digraph(spdigraph, lamb_m_opt)

# plot full results
plot_digraph(spdigraph_fit, f'out/Uganda2017/{month}', projection)

# plot difference graph
fig = plt.figure(dpi=300)
ax= fig.add_subplot(1, 1, 1, projection=projection)
v = Vis(ax, sp_digraph, projection=projection, edge_width=1,
        edge_alpha=1, edge_zorder=100, sample_pt_size=20,
        obs_node_size=5, sample_pt_color="black",
        cbar_font_size=5, cbar_ticklabelsize=5, 
        cbar_bbox_to_anchor=(0.05, 0.2), 
        cbar_width="15%",
        cbar_height="5%",
        compass_bbox_to_anchor=(0, 0),
        compass_font_size=5,
        compass_radius=0.2,
        compass_pad=0.2,
        mutation_scale=8)

v.draw_migration_rates(ax,
                       mode='Difference',
                       draw_map=True,
                       draw_nodes=True,)
plt.savefig(f'out/Uganda2017/{month}_difference.png')

# plot full graph
fig = plt.figure(dpi=300)
ax = fig.add_subplot(1, 1, 1, projection=projection)
v.draw_migration_rates(ax,
                        mode='Full',
                        draw_map=True,
                        draw_nodes=True,)
plt.savefig(f'out/Uganda2017/{month}_full.png')
