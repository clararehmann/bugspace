import numpy as np, pandas as pd
import os, sys
import dask
import dask.array as da
from dask.diagnostics.progress import ProgressBar
# silence some warnings
dask.config.set(**{'array.slicing.split_large_chunks': False})
import allel
import plotly
import plotly.express as px
from gambiae_grid import *
import h3
from shapely.geometry import Polygon, Point, mapping, LineString
import shapely
from pyproj import Geod
import itertools

import malariagen_data
ag3 = malariagen_data.Ag3()

def get_query(years, sample_sets):
    # generate query string for {year} samples
    years='|'.join([f'year == {year}' for year in years])
    sample_query="taxon=='gambiae' & " + years
    print(sample_query)
    df_samples = ag3.sample_metadata(sample_sets=sample_sets, sample_query=sample_query)
    return df_samples, sample_query

def assign_grid_cells(sample_metadata, resolution=3):
    ## assign each sample to h3 grid cell
    # generate grid
    lat_c, lon_c, max_d = get_sampling_radius(sample_metadata)
    grid = generate_h3_grid(lat_c, lon_c, max_d, resolution=resolution)
    # assign cells to samples
    sample_metadata['h3_index'] = None
    sample_metadata['point'] = [Point(sample_metadata.iloc[i].longitude, sample_metadata.iloc[i].latitude) for i in range(len(sample_metadata))]
    print(len(grid))
    for cell in range(len(grid)):
        mask = [shapely.contains(grid.iloc[cell].geometry, sample_metadata.iloc[i].point) for i in range(len(sample_metadata))]
        sample_metadata.loc[mask, 'h3_index'] = grid.iloc[cell].h3_index
    return sample_metadata

years = [2022]
sample_sets = [f'3.{i}' for i in range(17)]
sample_metadata, sample_query = get_query(years, sample_sets=sample_sets)

#haps = ag3.haplotypes(
#    region='3L',
#    sample_query=sample_query
#)
# reorder sample metadata to match zarr sample indexing
#sample_metadata.set_index('sample_id', inplace=True)
#samples = haps['sample_id'].values
#sample_metadata = sample_metadata.reindex(samples)
#sample_metadata = sample_metadata.reset_index(names=['sample_id'])
#sample_metadata['zarr_index'] = np.arange(len(sample_metadata))
# assign grid cells
sample_metadata = assign_grid_cells(sample_metadata)
sample_metadata.to_csv('data/2022_gambiae_metadata.csv')
# pairwise combinations of samples
indexes = np.arange(len(sample_metadata))
combos = np.array(list(itertools.combinations(indexes, 2)))
print(combos.shape)
print(combos[0])
print(combos[0,1])
