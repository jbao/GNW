#!/usr/bin/env python

import networkx as nx

#g = nx.barabasi_albert_graph(size, degree)
# in-degree: Poisson; out-degree: power-law
size = 1000
connected = False
while not connected:
    check = 0
    while check < 100 or check > 300:
        z = nx.create_degree_sequence(size, nx.utils.powerlaw_sequence,\
            exponent=2.5)
        check = max(z)
    nx.is_valid_degree_sequence(z)
    g = nx.directed_configuration_model(np.random.poisson(3,size).tolist(), z)
    connected = nx.algorithms.components.is_connected(g.to_undirected())
