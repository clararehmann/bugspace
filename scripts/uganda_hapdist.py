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

import malariagen_data
ag3 = malariagen_data.Ag3()

def get_month_query(month, sample_sets):
    # generate query string for 2017 samples from Uganda in that month
    sample_query=f"country=='Uganda' & year==2017 & taxon=='gambiae' & month=={month}"
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
    for chrom_comb in [(0,0), (0,1), (1,0), (1,1)]:
        splits = shared_haplotype(gt[:,ind1,chrom_comb[0]], gt[:,ind2,chrom_comb[1]])
        lengths = [pos[i[-1]] - pos[i[0]] for i in splits]
        length = max(lengths)
        if length > maximum:
            maximum = length
    return maximum

sample_sets = [f'3.{i}' for i in range(17)]
month=sys.argv[1]
sample_metadata, sample_query = get_month_query(month, sample_sets=sample_sets)

haps = ag3.haplotypes(
    region='3L',
    sample_query=sample_query
)
sample_metadata.set_index('sample_id', inplace=True)
samples = haps['sample_id'].values
sample_metadata = sample_metadata.reindex(samples)

## assign each sample to h3 grid cell
# generate grid
lat_c, lon_c, max_d = get_sampling_radius(sample_metadata)
resolution=4
gridpath=f'test'
grid = generate_h3_grid(lat_c, lon_c, max_d, resolution=resolution, outpath=gridpath)
# assign cells to samples
sample_metadata['h3_index'] = None
sample_metadata['point'] = [Point(sample_metadata.iloc[i].longitude, sample_metadata.iloc[i].latitude) for i in range(len(sample_metadata))]
for cell in range(len(grid)):
    mask = [shapely.contains(grid.iloc[cell].geometry, sample_metadata.iloc[i].point) for i in range(len(sample_metadata))]
    sample_metadata.loc[mask, 'h3_index'] = grid.iloc[cell].h3_index


genotypes = haps['call_genotype'].values
positions = haps['variant_position'].values
geod = Geod(ellps="WGS84")

lengths = []
hex_dists = []
km_dists = []
for i in range(len(samples)):
    for j in range(len(samples)):
        if i==j:
            pass
        else:
            lengths.append(haplotype_between(i, j, genotypes, positions))
            hex_dists.append(h3.grid_distance(sample_metadata.iloc[i].h3_index, sample_metadata.iloc[j].h3_index))
            km_dists.append(geod.inv(sample_metadata.iloc[i].longitude, 
                                    sample_metadata.iloc[i].latitude,
                                    sample_metadata.iloc[j].longitude,
                                    sample_metadata.iloc[j].latitude)[2]*0.001)
df = pd.DataFrame({'length':lengths,
                   'hex_distance':hex_dists,
                   'km_distance':km_dists})
df.to_csv(f'out/Uganda2017/{month}_haplotype_distances.csv')

#print(sample_metadata.iloc[i1].longitude, sample_metadata.iloc[i1].latitude)
#print(sample_metadata.iloc[i2].longitude, sample_metadata.iloc[i2].latitude)

