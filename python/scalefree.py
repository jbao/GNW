#!/usr/bin/env python

import os
import numpy as np
import networkx as nx
from randht import randht
from plfit import plfit

#g = nx.barabasi_albert_graph(size, degree)
# in-degree: Poisson; out-degree: power-law
size = 1000
for n in range(100):
    connected = False
    #while not connected:
    #    print 'is connected',connected
    #    check = 0
    #    while check < 100 or check > 300:
    #        z = nx.utils.create_degree_sequence(size, nx.utils.powerlaw_sequence,\
    #            exponent=2, max_tries=100)
    #        check = max(z)
    #    print 'is valid',nx.is_valid_degree_sequence(z)
    #    g = nx.directed_configuration_model(np.random.poisson(3,size).tolist(), z)
    #    connected = nx.algorithms.components.is_connected(g.to_undirected())
    #in_ds = np.random.poisson(3.5,size).tolist()
    #valid_ds = False
    #while not valid_ds:
    #in_ds = np.random.poisson(3,size).tolist()
    #out_ds = []
    #for i in range(size):
    #    powerlaw = randht(1,'xmin',5,'powerlaw',2.5)[0]
    #    out_sum = sum(out_ds) + powerlaw
    #    if out_sum <= sum(in_ds):
    #        out_ds.append(powerlaw)
    #    elif i == size-1:
    #        out_ds.append(sum(in_ds)-sum(out_ds))
    #    else:
    #        out_ds.append(0)
    #    #ds = np.array(ds)-min(ds)+1
    #    #valid_ds = nx.is_valid_degree_sequence(ds)
    #    #print sum(in_ds),sum(out_ds)
    #    #print plfit(out_ds)
    #    #    print sum(in_ds),sum(out_ds)
    #    #    in_ds = np.random.poisson(3,size).tolist()
    #    #    out_ds = randht(size,'powerlaw',2)
       
    while not connected:
        #out_ds = np.random.permutation(out_ds)
        #print 'is connected?',connected
        #g = nx.directed_configuration_model(in_ds, out_ds)
        #g = nx.configuration_model(ds)
        g = nx.scale_free_graph(size, 0.1, 0.1, 0.8)
        print 'max in-degree',max(g.in_degree().values())
        print 'plfit',plfit(g.degree().values())
        connected = nx.algorithms.components.is_connected(g.to_undirected())
    g = nx.DiGraph(g)
    g.remove_edges_from(g.selfloop_edges())
    #dg = nx.DiGraph()
    #dg.add_nodes_from(g.nodes())
    #for e in g.edges():
    #    if np.random.rand() < 0.5:
    #        dg.add_edge(e[0],e[1])
    #    else:
    #        dg.add_edge(e[1],e[0])
    g = nx.relabel_nodes(g, dict(zip(g.nodes(), ['G'+str(i+1) for i in \
            range(size)])))
    wd = '/home/jbao/data/DREAM/gnw/scalefree2/gnw/Size'+str(size)+'/'
    if not os.path.isdir(wd):
        os.makedirs(wd)
    nx.write_edgelist(g, wd+'scalefree2-'+str(n+1)+'_'+str(size)+'.tsv', \
            data=False, delimiter='\t')
    #pylab.savetxt(wd+'hub-'+str(n+1)+'_'+str(size)+'_degree.tsv', 
    #        pylab.c_[g.in_degree().values(),g.out_degree().values()], '%d', 
    #        delimiter='\t')
    #g.clear()
