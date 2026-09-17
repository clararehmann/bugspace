import argparse
import numpy as np, pandas as pd
import allel
from frame.spatial_digraph import SpatialDiGraph
from frame.visualization import Vis
from frame.cross_validation import run_cv
from frame.digraphstats import Digraphstats
from frame.utils import prepare_graph_inputs
from ag3_popgen import get_zarr_genotypes
from matplotlib import pyplot as plt

def genotype_array_to_FRAME(gt):
    ac = gt.to_allele_counts()
    # restrict to only biallelic sites
    bi = ac.is_biallelic()
    gt = gt[np.sum(bi, axis=1) > 0]
    #convert to alternate allele=1
    gt = gt > 0
    # reduce to one dimension for each sample (1=het, 2=hom alt) and reshape for FRAME
    gt = np.sum(gt, axis=2)
    gt = gt.T
    return gt

def initialize_digraph(coordinates, edges, grid):
    """
    inputs
    coordinates: string of path to coordinates file (in same order as samples) or numpy array of coordinates (in same order as samples)
    edges: string of path to edges file or numpy array of edge coordinates
    grid: path to shapefile of grid nodes
    """
    # read spatial files
    if typeof(coordinates) == str:
        coordinates = np.genfromtxt(coordinates, delimiter=',')
    if typeof(edges) == str:
        edges = np.genfromtxt(edges, delimiter=',')
    outer, edges, grid, _ = prepare_graph_inputs(coord=coordinates,
                                                ggrid=args.grid, 
                                                buffer=0,
                                                outer=edges)


    # construct spatial digraph
    sp_digraph = SpatialDiGraph(gt, 
                                coordinates,
                                grid,
                                edges)

def lambda_cv(N=10, lamb_m_warmup=1e3, outpath):
    """
    run lambda cross validation
    inputs
    N: integer number of cross-validation runs (default 10)
    """
    # N fold cross validation to find lambda
    lamb_m_grid = np.geomspace(1e-3,1e3,20)[::-1]
    cv_errs,node_train_idxs=run_cv(sp_digraph,
                                lamb_m_grid=lamb_m_grid,
                                lamb_m_warmup=lamb_m_warmup,
                                n_folds=N,
                                factr=1e10,
                                random_state=200,)

    plt.figure(figsize=(8, 6))
    plt.plot(np.log10(lamb_m_grid), cv_errs, 'bo')   
    plt.xlabel(r"$\mathrm{log}_{10}(\mathrm{\lambda_m})$")
    plt.ylabel('CV Error')
    plt.savefig(f"{outpath}_lambda_cv.png")

    # use lambda to fit spatial digraph
    lamb_m_opt=lamb_m_grid[np.argmin(cv_errs)]
    lamb_m_opt=float("{:.3g}".format(lamb_m_opt))
    return lamb_m_opt

def fit_digraph(sp_digraph, lamb_m_opt, lamb_m_warmup=1e3, outpath)
    sp_digraph.fit(lamb_m=lamb_m_warmup, factr=1e10)
    logm = np.log(sp_digraph.m)
    logc = np.log(sp_digraph.c)
    trans_alpha=-np.log((1/sp_digraph.alpha)-1)

    sp_digraph.fit(lamb_m=lamb_m_opt,
                factr=1e7,
                logm_init=logm,
                logc_init=logc,
                trans_alpha_init=trans_alpha,
                )
    return sp_digraph

def plot_digraph(sp_digraph):
    fig, axs= plt.subplots(2, 4, figsize=(16, 5), dpi=300,
                            subplot_kw={'projection': projection})

    v = Vis(axs[0,0], sp_digraph, projection=projection, edge_width=1,
            edge_alpha=1, edge_zorder=100, sample_pt_size=20,
            obs_node_size=5, sample_pt_color="black",
            cbar_font_size=5, cbar_ticklabelsize=5, 
            cbar_bbox_to_anchor=(0.05, 0.2), 
            cbar_width="15%",
            cbar_height="5%",
            compass_bbox_to_anchor=(0, 0),
            compass_font_size=5,
            compass_radius=0.2,
            mutation_scale=6)

    v.digraph_wrapper(axs, node_scale=[5, 5, 5])
    plt.subplots_adjust(hspace=0)
    plt.savefig(f"{outpath}_results.png")