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
import itertools
from pyproj import Geod

import malariagen_data
ag3 = malariagen_data.Ag3()

def get_query(sample_sets):
    # generate query string for 2017 samples from Uganda in that month
    sample_query=f"country=='Uganda' & year==2017 & taxon=='gambiae'"
    df_samples = ag3.sample_metadata(sample_sets=sample_sets, sample_query=sample_query)
    return df_samples, sample_query

def shared_haplotype(gt1, gt2):
    # where genotypes are the same
    gt_match = np.where(gt1 == gt2)[0]
    # position distance between matches
    gt_diffs = np.diff(gt_match)
    # split the array where the difference is not one
    splits = np.split(gt_match, np.where(gt_diffs != 1)[0] + 1)
    return splits

def haplotype_between(ind1, ind2, gt, pos):
    # compare 4 ways between all chromosomes
    maximum = 0
    lengths_all = []
    for chrom_comb in [(0,0), (0,1), (1,0), (1,1)]:
        splits = shared_haplotype(gt[:,ind1,chrom_comb[0]], gt[:,ind2,chrom_comb[1]])
        lengths = [pos[i[-1]] - pos[i[0]] for i in splits]
        max_length = max(lengths)
        lengths_all.append(lengths)
        if max_length > maximum:
            maximum = max_length
    lengths = np.array([l for ls in lengths_all for l in ls])
    med_length = np.median(lengths)
    mean_length = np.mean(lengths)
    return maximum, med_length, mean_length

def assign_grid_cells(sample_metadata, resolution=4):
    ## assign each sample to h3 grid cell
    sample_metadata['h3_index'] = [h3.latlng_to_cell(sample_metadata.iloc[i].latitude,
                                                     sample_metadata.iloc[i].longitude,
                                                     resolution) for i in range(len(sample_metadata))]
    #lat_c, lon_c, max_d = get_sampling_radius(sample_metadata)
    #resolution=4
    #grid = generate_h3_grid(lat_c, lon_c, max_d, resolution=resolution)
    # assign cells to samples
    #sample_metadata['h3_index'] = None
    #sample_metadata['point'] = [Point(sample_metadata.iloc[i].longitude, sample_metadata.iloc[i].latitude) for i in range(len(sample_metadata))]
    #for cell in range(len(grid)):
    #    mask = [shapely.contains(grid.iloc[cell].geometry, sample_metadata.iloc[i].point) for i in range(len(sample_metadata))]
    #    sample_metadata.loc[mask, 'h3_index'] = grid.iloc[cell].h3_index
    return sample_metadata

sample_sets = [f'3.{i}' for i in range(17)]
sample_metadata, sample_query = get_query(sample_sets=sample_sets)

haps = ag3.haplotypes(
    region='3L',
    sample_query=sample_query
)
sample_metadata.set_index('sample_id', inplace=True)
samples = haps['sample_id'].values
sample_metadata = sample_metadata.reindex(samples)
sample_metadata = sample_metadata.reset_index(names='sample_id')
## assign each sample to h3 grid cell
# generate grid
sample_metadata = assign_grid_cells(sample_metadata)

genotypes = haps['call_genotype'].values
positions = haps['variant_position'].values
geod = Geod(ellps="WGS84")

max_lengths = []
med_lengths = []
mean_lengths = []
hex_dists = []
km_dists = []
month_1 = []
month_2 = []
idxs = np.arange(len(sample_metadata))
comb = np.array(list(itertools.combinations(idxs, 2)))
for c in comb:
    i = c[0]
    j = c[1]
    try:
        max_length, med_length, mean_length = haplotype_between(i, j, genotypes, positions)
        max_lengths.append(max_length)
        med_lengths.append(med_length)
        mean_lengths.append(mean_length)
        hex_dists.append(h3.grid_distance(sample_metadata.iloc[i].h3_index, sample_metadata.iloc[j].h3_index))
        km_dists.append(geod.inv(sample_metadata.iloc[i].longitude, 
                                sample_metadata.iloc[i].latitude,
                                sample_metadata.iloc[j].longitude,
                                sample_metadata.iloc[j].latitude)[2]*0.001)
        month_1.append(sample_metadata.iloc[i, sample_metadata.columns.get_loc('month')])
        month_2.append(sample_metadata.iloc[j, sample_metadata.columns.get_loc('month')])
    except:
        pass
df = pd.DataFrame({'max_length':max_lengths,
                   'med_length':med_lengths,
                   'mean_length':mean_lengths,
                   'hex_distance':hex_dists,
                   'km_distance':km_dists,
                   'month_1':month_1,
                   'month_2':month_2})
df.to_csv(f'out/Uganda2017/all_haplotype_distances.csv')

#print(sample_metadata.iloc[i1].longitude, sample_metadata.iloc[i1].latitude)
#print(sample_metadata.iloc[i2].longitude, sample_metadata.iloc[i2].latitude)

