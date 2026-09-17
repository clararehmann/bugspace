import numpy as np, pandas as pd
import dask
import dask.array as da
from dask.diagnostics.progress import ProgressBar
# silence some warnings
dask.config.set(**{'array.slicing.split_large_chunks': False})
import allel
import malariagen_data
import plotly.express as px
import itertools

def run_plot_pca(sample_query, region='3L', n_snps=100000, site_mask='gamb_colu'):
    pca_df, evr_bf = ag3.pca(region=region,
                             n_snps=n_snps,
                             site_mask=site_mask,
                             sample_query=sample_query)
    plot = ag3.plot_pca_coords(pca_df, show=False)
    plot.update_layout(
        xaxis={title:{text:f'PC1: {np.round(evr_bf[0]*100, decimals=3)}% variance explained'}},
        yaxis={title:{text:f'PC2: {np.round(evr_bf[1]*100, decimals=3)}% variance explained'}}
    )
    return plot

def run_popgen_stats(sample_query, cohorts='admin2_year', site_mask='gamb_colu'):
    stats_df = ag3.diversity_stats(
        sample_query=sample_query,
        cohorts=cohorts,
        site_mask=site_mask
    )
    return stats_df

def get_zarr_genotypes(sample_query, region='3L', site_mask='gamb_colu'):
    gt = ag3.snp_calls(
        region=region,
        sample_query=sample_query,
        site_mask=site_mask
    )
    gt = allel.GenotypeDaskArray(gt["call_genotype"].data)
    return gt

def label_location_groups(sample_metadata):
    """
    adds another column to the sample metadata with an integer value denoting location group
    """
    unique_locations = df[['longitude', 'latitude']].drop_duplicates()
    unique_locations = list(zip(unique_locations.longitude, unique_locations.latitude))
    sample_groups = range(unique_locations.shape[0])
    sample_groups = dict(zip(unique_locations, sample_groups))
    sample_metadata.assign(sample_group = lambda x: (sample_groups[(x.longitude, x.latitude)]))
    return sample_metadata

def pairwise_fst(sample_metadata, sample_genotypes):
    # sort subpops
    subpops = []
    for i in range(max(sample_metadata.sample_group) + 1):
        subpops.append(list(np.where(sample_metadata.sample_group == i)))
    # allele counts for each subpop
    allele_counts = {}
    for i in range(len(subpops)):
        allele_counts[i] = sample_genotypes.count_alleles(subpop=subpops[i])
    # FST between each location combination
    fst_values = {}
    sample_combinations = list(itertools.product(sample_metadata.sample_group, sample_metadata.sample_group))
    for i in sample_combinations:
        if (i[1], i[0]) in fst_values.keys():
            pass # don't run if that location combination already has been run
        elif i[0] == i[1]:
            pass # don't run between the same sample group
        else:
            num, den = allel.hudson_fst(ac[i[0]], ac[i[1]])
            fst_values[i] = np.sum(num) / np.sum(den)
    return sample_metadata, fst_values
