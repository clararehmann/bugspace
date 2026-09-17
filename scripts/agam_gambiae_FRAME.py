import argparse
import numpy as np, pandas as pd
import allel
from frame.spatial_digraph import SpatialDiGraph
from frame.visualization import Vis
from frame.cross_validation import run_cv
from frame.digraphstats import Digraphstats
from frame.utils import prepare_graph_inputs
from matplotlib import pyplot as plt

parser=argparse.ArgumentParser(description="run FRAME (optionally on a window) of a VCF file")
parser.add_argument("--vcf", type=str, required=True)
parser.add_argument("--chromosome", type=str,  help="chromosome to analyze")
parser.add_argument("--start", type=int, help="start position of the window")
parser.add_argument("--end", type=int, help="end position of the window")
parser.add_argument("--coordinates", type=str, help="path to a file of lat,lon coordinates for each sample in the VCF")
parser.add_argument("--grid", type=str, help="path to a file of lat,lon grid coordinates in the form of a shapefile to run FRAME over")
parser.add_argument("--edges", type=str, help="path to a file of edges to the sampling region")
parser.add_argument("--output", type=str, help="path to output directory")
args=parser.parse_args()

if args.start and args.end and args.chromosome:
    vcf = allel.read_vcf(args.vcf, region=f"{args.chromosome}:{args.start}-{args.end}", 
                         tabix='/hps/software/users/jlees/rehmann/pixi/envs/bugspace-7036013565992295780/envs/default/bin/tabix')
else:
    vcf = allel.read_vcf(args.vcf)
gt = allel.GenotypeArray(vcf['calldata/GT'])
ac = gt.to_allele_counts()
# restrict to only biallelic sites
bi = ac.is_biallelic()
gt = gt[np.sum(bi, axis=1) > 0]
#convert to alternate allele=1
gt = gt > 0
# reduce to one dimension for each sample (1=het, 2=hom alt) and reshape for FRAME
gt = np.sum(gt, axis=2)
gt = gt.T

# read spatial files
coordinates = np.genfromtxt(args.coordinates, delimiter=',')
central_latitude = np.mean(coordinates[:,0])
central_longitude = np.mean(coordinates[:,1])
projection = ccrs.EquidistantConic(central_longitude = central_longitude, central_latitude = central_latitude)
outer = np.genfromtxt(args.edges, delimiter=',')
outer, edges, grid, _ = prepare_graph_inputs(coord=coordinates,
                                            ggrid=args.grid, 
                                            buffer=0,
                                            outer=outer)


# construct spatial digraph
sp_digraph = SpatialDiGraph(gt, 
                            coordinates,
                            grid,
                            edges)

# N fold cross validation to find lambda
N = 10 
lamb_m_warmup = 1e3
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
plt.savefig(f'{args.output}_lambda_cv.png')

# use lambda to fit spatial digraph
lamb_m_opt=lamb_m_grid[np.argmin(cv_errs)]
lamb_m_opt=float("{:.3g}".format(lamb_m_opt))

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
plt.savefig(f'{args.output}_results.png')
