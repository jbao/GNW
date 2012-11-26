#!/usr/bin/env python

from xstats.MINE import analyze_pair
import numpy as np
import networkx as nx

wd = '/Users/bao/geronto/data/DREAM/gnw/scalefree2/gnw/Size1000/'
mic = np.zeros(5000)
ii = 0
for net in range(100):
    g = nx.read_edgelist(wd+'input/scalefree2-'+str(net+1)+'_1000.tsv', create_using=nx.DiGraph())
    for perturb in range(50):
        genes = np.loadtxt(wd+'output/mds/scalefree2-'+str(net+1)+'_perturbation-'+str(perturb+1)+'_1000_normed.dat', delimiter='|', skiprows=1, usecols=[0], dtype=str) 
        pvals = np.loadtxt(wd+'output/mds/pval_scalefree2-'+str(net+1)+'_perturbation-'+str(perturb+1)+'_1000.dat') 

        res = np.zeros(g.order())
        deg = np.zeros(g.order())
        for i in range(g.order()):
            gene = g.nodes()[i]
            idx = genes.tolist().index(gene)
            res[i] = -np.log10(pvals[idx])
            deg[i] = g.degree()[gene]

        mic[ii] = analyze_pair(res, deg)['MIC']
        ii += 1
