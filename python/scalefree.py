#!/usr/bin/env python

import numpy as np
import networkx as nx
from randht import randht
from plplot import plplot
from plfit import plfit

#g = nx.barabasi_albert_graph(size, degree)
# in-degree: Poisson; out-degree: power-law
size = 1000
connected = False
for i in range(1):
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
    while not connected:
        print 'is connected?',connected
        in_ds = np.random.poisson(3,size).tolist()
        out_ds = randht(size,'powerlaw',2)
        print sum(in_ds),sum(out_ds)
        print plfit(out_ds)
        #    print sum(in_ds),sum(out_ds)
        #    in_ds = np.random.poisson(3,size).tolist()
        #    out_ds = randht(size,'powerlaw',2)
       
        #g = nx.directed_configuration_model(in_ds, np.array(out_ds)-1)
        g = nx.directed_configuration_model(in_ds, out_ds)
        connected = nx.algorithms.components.is_connected(g.to_undirected())
    #g = nx.relabel_nodes(g, dict(zip(g.nodes(), ['G'+str(i+1) for i in \
    #        range(size)])))
    #wd = '/home/jbao/data/DREAM/gnw/hub/gnw/Size'+str(size)+'/'
    #if not os.path.isdir(wd):
    #    os.makedirs(wd)
    #nx.write_edgelist(g, wd+'hub-'+str(n+1)+'_'+str(size)+'.tsv', \
    #        data=False, delimiter='\t')
    #pylab.savetxt(wd+'hub-'+str(n+1)+'_'+str(size)+'_degree.tsv', 
    #        pylab.c_[g.in_degree().values(),g.out_degree().values()], '%d', 
    #        delimiter='\t')
    #g.clear()
