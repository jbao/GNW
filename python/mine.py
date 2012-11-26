#!/usr/bin/env python

from xstats.MINE import analyze_pair
import numpy as np
import networkx as nx

wd = '/home/jbao/data/DREAM/gnw/scalefree2/gnw/Size1000/'
g = nx.read_edgelist(wd+'input/scalefree2-1_1000.tsv', create_using=nx.DiGraph())
genes = np.loadtxt(wd+'output/mds/scalefree2-1_perturbation-1_1000_normed.dat', delimiter='|', skiprows=1, usecols=[0], dtype=str) 
pvals = np.loadtxt(wd+'output/mds/pval_scalefree2-1_perturbation-1_1000.dat') 

res = np.zeros(g.order())
deg = np.zeros(g.order())
for i in range(g.order()):
    gene = g.nodes()[i]
    idx = genes.tolist().index(gene)
    res[i] = -np.log10(pvals[idx])
    deg[i] = g.degree()[gene]

mine = analyze_pair(res, deg)
