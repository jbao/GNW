#!/usr/bin/env python

from xstats.MINE import analyze_pair
import numpy as np
import networkx as nx
import itertools
import multiprocessing as mp
from joblib import Parallel,delayed

wd = '/home/jbao/data/DREAM/gnw/scalefree2/gnw/Size1000/'
net = range(100)
perturb = range(50)
all_net,all_perturb = zip(*itertools.product(net, perturb))

def get_mic(net, perturb):
    g = nx.read_edgelist(wd+'input/scalefree2-'+str(net+1)+'_1000.tsv', create_using=nx.DiGraph())
    genes = np.loadtxt(wd+'output/mds/scalefree2-'+str(net+1)+'_perturbation-'+str(perturb+1)+'_1000_normed.dat', delimiter='|', skiprows=1, usecols=[0], dtype=str) 
    pvals = np.loadtxt(wd+'output/mds/pval_scalefree2-'+str(net+1)+'_perturbation-'+str(perturb+1)+'_1000.dat') 

    res = np.zeros(g.order())
    deg = np.zeros(g.order())
    for i in range(g.order()):
        gene = g.nodes()[i]
        idx = genes.tolist().index(gene)
        res[i] = -np.log10(pvals[idx])
        deg[i] = g.degree()[gene]

    return analyze_pair(res, deg)['MIC']

def get_mic_helper(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return get_mic(*a_b)

result_list = []
def log_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    result_list.append(result)

n_cpu = 1
mic = Parallel(n_jobs=n_cpu, verbose=10)(delayed(get_mic)(net,perturb) for net,perturb in itertools.izip(all_net,all_perturb))

#pool = mp.Pool()
#for i in range(len(all_net)):
#    pool.apply_async(get_mic, args = (all_net[i],all_perturb[i]), callback = log_result)
#pool.close()
#pool.join()
#print(result_list)
#p = Pool()
#mic = p.map(get_mic_helper, itertools.izip(all_net, all_perturb))
#np.savetxt(wd+'output/all_mic.dat', mic)
