import numpy as np, pandas as pd
import os
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

def get_year_query(year, sample_sets):
    # generate query string for 2017 samples from Uganda in that month
    df_samples = ag3.sample_metadata(sample_sets=sample_sets, sample_query=f"country=='Uganda' & year=={year} & taxon=='gambiae'")
    samples = df_samples['sample_id'].tolist()
    sample_query = '|'.join([f"sample_id=='{s}'" for s in samples])
    return df_samples, sample_query

sample_sets = [f'3.{i}' for i in range(17)]
year=2017
sample_metadata, sample_query = get_year_query(year, sample_sets)

# make sampling grid and edges
lat_c, lon_c, max_d = get_sampling_radius(sample_metadata)
projection = ccrs.EquidistantConic(central_longitude = lon_c, central_latitude = lat_c)
resolution=4
gridpath='data/Uganda2017'
grid = generate_h3_grid(lat_c, lon_c, max_d, resolution=resolution, outpath=gridpath)
hull = generate_hull(sample_metadata, gridpath)
# prep coordinates for FRAME
coordinates = sample_metadata[['longitude','latitude']].to_numpy()
# prep genotypes for FRAME
print(sample_query)
genotypes = get_zarr_genotypes(sample_query=f"country=='Uganda' & year=={year} & taxon=='gambiae'")
genotypes = allel.GenotypeArray(genotypes)
genotypes = genotype_array_to_FRAME(genotypes)

# make spatial digraph
spdigraph = initialize_digraph(genotypes, coordinates, f'{gridpath}_hull.csv', f'{gridpath}_grid_resolution_{resolution}.shp')

# lambda cross-validation
lamb_m_opt = lambda_cv(spdigraph, 'out/Uganda2017')
spdigraph_fit = fit_digraph(spdigraph, lamb_m_opt)
plot_digraph(spdigraph_fit, 'out/Uganda2017', projection)