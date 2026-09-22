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
ag3 = malariagen_data.Ag3()

def run_plot_pca(sample_query, region='3L', n_snps=100000, site_mask='gamb_colu'):
    pca_df, evr_bf = ag3.pca(region=region,
                             n_snps=n_snps,
                             site_mask=site_mask,
                             sample_query=sample_query,
                             sample_query_options = {'engine':'python'})
    plot = ag3.plot_pca_coords(pca_df, show=False)
    plot.update_layout(
        xaxis_title=f'PC1: {np.round(evr_bf[0]*100, decimals=3)}% variance explained',
        yaxis_title=f'PC2: {np.round(evr_bf[1]*100, decimals=3)}% variance explained'
    )
    return plot

def run_popgen_stats(sample_query, cohort_size, cohorts='admin2_year', site_mask='gamb_colu', region='3L'):
    stats_df = ag3.diversity_stats(
        sample_query=sample_query,
        cohorts=cohorts,
        cohort_size=cohort_size,
        region=region,
        site_mask=site_mask,
        sample_query_options = {'engine':'python'}
    )
    return stats_df

def get_zarr_genotypes(sample_query, region='3L', site_mask='gamb_colu'):
    gt = ag3.snp_calls(
        region=region,
        sample_query=sample_query,
        site_mask=site_mask,
        sample_query_options={'engine':'python'}
    )
    gt = allel.GenotypeDaskArray(gt["call_genotype"].data)
    return gt

def between_cohort_fst(sample_metadata, region='3L', site_mask='gamb_colu'):
    country = np.unique(sample_metadata.country)
    year = np.unique(sample_metadata.year)
    month = np.unique(sample_metadata.month)
    cohorts = np.unique(sample_metadata.cohort_admin2_year)
    cohort_combinations = list(itertools.product(cohorts, cohorts))
    fsts = np.empty(len(cohort_combinations))
    stderrs = np.empty(len(cohort_combinations))
    cohort_sizes = [len(sample_metadata[sample_metadata.cohort_admin2_year==c]) for c in np.unique(sample_metadata.cohort_admin2_year)]
    for i in range(len(cohort_combinations)):
        if cohort_combinations[i][0] == cohort_combinations[i][1]:
            fsts[i] = 0
            stderrs[i] = 0
        else:
            ch1 = f"taxon=='gambiae' & country=={country} & year=={year} & month=={month} & cohort_admin2_year=='{cohort_combinations[i][0]}'"
            ch2 = f"taxon=='gambiae' & country=={country} & year=={year} & month=={month} & cohort_admin2_year=='{cohort_combinations[i][1]}'"
            fst, stderr = ag3.average_fst(region=region,
                                      cohort1_query=ch1,
                                      cohort2_query=ch2,
                                      site_mask=site_mask,
                                      min_cohort_size=min(cohort_sizes))
            fsts[i] = fst
            stderrs[i] = stderr
    return cohort_combinations, fsts, stderrs

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
