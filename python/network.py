#!/usr/bin/env python
"""
plot graph from a edgelist file.

$Id: network.py 307 2012-04-16 13:04:29Z jbao $
"""

from pdb import set_trace
import sys,os,gc
import networkx as nx
from networkx.algorithms.core import *
#from hcluster import squareform,linkage,dendrogram
#import pygraphviz as pgv
from numpy import loadtxt
import numpy as np
from scipy import polyval, polyfit, optimize
from scipy.stats import pearsonr, spearmanr
from re import sub,finditer
from pylab import *
from matplotlib.font_manager import FontProperties
from mpl_toolkits.mplot3d import Axes3D
from os import system, getcwd
#import community
#import linreg
import data
reload(data)
import imp
#f, filename, desc = imp.find_module('param_loader', [getcwd()+'/gnw/'])
#pl = imp.load_module('param_loader', f, filename, desc)
#from gnw.param_loader import ParamLoader
#reload(param_loader)
#from tree import rtree
import half_life 
reload(half_life)

#rcParams['axes.labelsize'] = 30
#rcParams['axes.titlesize'] = 32
#rcParams['xtick.labelsize'] = 20
#rcParams['ytick.labelsize'] = 20
#rcParams['xtick.major.pad'] = 10
#rcParams['ytick.major.pad'] = 10
#rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

class Graph:

    def __init__(self, filename):
            
        self.g_d = nx.read_edgelist(sub(r'\\','',filename), \
                delimiter='\t',\
                create_using=nx.DiGraph())
        self.g = nx.read_edgelist(sub(r'\\','',filename), \
                delimiter='\t',\
                create_using=nx.Graph())
        self.filename = filename

            #outdeg = g_d.out_degree()
            #indeg = g_d.in_degree()
        """
            c = nx.betweenness_centrality(g)
            figure()
            l = nx.graphviz_layout(g)
            nx.draw(g, pos=l, node_color='w')
            nx.draw_networkx_nodes(g, l, nodelist=c.keys(), \
                    node_color=[i*10 for i in c.values()],\
                    cmap=cm.Reds)
           
            print '%s\t%f\t%f' % (files[i], mean(multiply(deg,deg))-mean(deg)**2,\
                    mean(d[:,i]**2)-mean(d[:,i])**2)
            cc[0,i] = std(deg)/mean(deg)
            cc[1,i] = std(d[:,i])/mean(d[:,i])
            
            #c[i] = nx.average_clustering(g)
            
            get_activation_inhibition()
            x = []
            y = []
            subplot(2,5,i+1)
            #if i < 5:
            #    ax = Axes3D(fig, rect=[i*0.2,0.5,0.15,0.45])
            #else:
            #    ax = Axes3D(fig, rect=[(i-5)*0.2,0,0.15,0.45])

            
            for ii,gene in enumerate(idx):
                #if deg[g.nodes().index('G'+str(gene+1))] is 0:
                #    plot(ii, d[gene,i], '.r')
                #else:
                myx = indeg[g.nodes().index(mds.genes[gene])]
                myy = outdeg[g.nodes().index(mds.genes[gene])]
                myz = d[gene,i]
                #myy = outdeg[g.nodes().index(mds.genes[gene])]/\
                #        (indeg[g.nodes().index(mds.genes[gene])]+1)
                #myy = outdeg[g.nodes().index(mds.genes[gene])]+\
                #        indeg[g.nodes().index(mds.genes[gene])]
                #myy = nx.closeness_centrality(g)[mds.genes[gene]]
                ax.plot([myx], [myy], [myz], '.k')
                hold(True)
                x.append(myx)
                y.append(myy)
            s,p = stats.stats.spearmanr(x,y)
            #yscale('log')
            ylim(ymin=-1)
            title('yeast'+str(i+1)+'\n'+str(s))
                
            
            subplot(1,5,i+1)
            for gene in range(100):
                plot(deg[g.nodes().index('G'+str(gene+1))], d[gene,i], '.k')
            
            
            hist(g.degree(), logspace(0,5.5,base=2), log=True)
            xscale('log')
            
            my_max = max(max(g.in_degree()), max(g.out_degree()))
            all = zeros((my_max+1,my_max+1))
            outlier_degree = zeros((my_max+1,my_max+1))
            for idx in range(100):
                all[g.in_degree()[i],g.out_degree()[i]] += 1
                if g.nodes()[i] in outliers:
                    outlier_degree[g.in_degree()[i],g.out_degree()[i]] += 1

            imshow(all, cm.Greys, origin='lower')
            imshow(outlier_degree, cm.Reds, origin='lower')
            title(files[i])
            xlabel('In-degree')
            ylabel('Out-degree')
            """

    def get_activation_inhibition(self, gene):

        activation_in = 0
        inhibition_in = 0
        activation_out = 0
        inhibition_out = 0

        for incoming in self.g_d.predecessors(gene):
            if self.g_d.edge[incoming][gene]['weight'] > 0:
                activation_in += 1
            elif self.g_d.edge[incoming][gene]['weight'] < 0:
                inhibition_in += 1
        for outgoing in self.g_d.successors(gene):
            if self.g_d.edge[gene][outgoing]['weight'] > 0:
                activation_out += 1
            elif self.g_d.edge[gene][outgoing]['weight'] < 0:
                inhibition_out += 1

        return activation_in,inhibition_in,activation_out,inhibition_out

    def plot_distance(self, dist, gene):
        #for d,g in zip(dist,gene):
        #    plot(d, self.get_activation_inhibition(g), '.k')
        my_max = max(max(self.g_d.in_degree().values()), \
                max(self.g_d.out_degree().values()))
        sum = zeros((my_max+1,my_max+1))
        num = zeros((my_max+1,my_max+1))
        avg = zeros((my_max+1,my_max+1))
        for d,g in zip(dist,gene):
            if g in self.g.nodes():
                a_in,i_in,a_out,i_out = self.get_activation_inhibition(g)
                #sum[a_in,i_in] += d
                #num[a_in,i_in] += 1
                sum[a_out,i_out] += d
                num[a_out,i_out] += 1
                plot(d,a_in-i_in,'.')

        for ii in range(my_max+1):
            for jj in range(my_max+1):
                if num[ii,jj] != 0:
                    avg[ii,jj] = sum[ii,jj]/num[ii,jj]

        #plot([0,my_max], [0,my_max], '--k')

        #pcolor(avg, cmap=cm.Reds)
        #clim(0,0.6)
        #colorbar()
        #xlim(0,15)
        #ylim(0,15)
        #xlabel('activating in-degree')
        #ylabel('inhibiting in-degree')
        show()
        
    def get_dot_graph(self, dist, width=None, pos=None):
        """Return a graph object that can be saved as a dot file.

        Parameters
        ----------
        dist : dictionary
            radial distance of each node according to MDS

        width : dictionary, optional
            edge weight to be plotted

        Returns
        -------
        g : Graph
            graph object
        """
        #l = nx.graphviz_layout(self.g_d)
        g = pgv.AGraph(strict=False, directed=True)
        g.graph_attr['overlap'] = 'compress'
        g.graph_attr['sep'] = '+10'
        # define outliers
        #set_trace()
        layers = {}
        layers['0'] = []
        hierarchy = zeros(len(self.g_d.nodes()))
        for idx,n in enumerate(self.g_d.nodes()):
            if self.g_d.out_degree()[n] == 0:
                layers['0'].append(n)
                hierarchy[idx] = 0

        md = max(dist.values())
        gene_sorted = sorted(dist, key=dist.get, reverse=True)

        for idx,n in enumerate(self.g_d.nodes()):
            """
            if self.g_d.out_degree()[n] != 0:
                # determine the shortest path to the output layer
                spl = 10000
                for i,nn in enumerate(layers['0']):
                    if nx.shortest_path(self.g_d,n,nn):
                        current = nx.shortest_path_length(self.g_d,n,nn)
                        if current < spl:
                            spl = current
                hierarchy[idx] = spl

                if not layers.has_key(str(spl)):
                    layers[str(spl)] = [n]
                else:
                    layers[str(spl)].append(n)
            """ 

            d = dist[n]
            if pos is None:
                g.add_node(n, fillcolor='#%2x%2x%2x'%((1-d/md)*255,(1-d/md)*255,(1-d/md)*255), \
                        size='20,20', fixedsize='True', \
                        height='0.2', width='0.2',
                        shape='circle', fontsize='0', \
                        fontcolor='white', style='filled')
            else:
                g.add_node(n, fillcolor='#%2x%2x%2x'%((1-d/md)*255,(1-d/md)*255,(1-d/md)*255), \
                        size='20,20', fixedsize='True', \
                        height='0.2', width='0.2',
                        shape='circle', fontsize='0', \
                        fontcolor='white', style='filled', \
                        pos=str(pos[n][0])+','+str(pos[n][1]))
            """
            if n in gene_sorted[:5]:   # outliers have different label color
                d = dist[n]
                if pos is None:
                    g.add_node(n, fillcolor='#%2x%2x%2x'%((1-d/md)*255,(1-d/md)*255,(1-d/md)*255), \
                            size='50,50', fixedsize='True', \
                            shape='circle', fontsize='10', \
                            fontcolor='white', style='filled')
                else:
                    g.add_node(n, fillcolor='#%2x%2x%2x'%((1-d/md)*255,(1-d/md)*255,(1-d/md)*255), \
                            size='50,50', fixedsize='True', \
                            shape='circle', fontsize='10', \
                            fontcolor='white', style='filled', \
                            pos=str(pos[n][0])+','+str(pos[n][1]))
                
                elif self.g_d.in_degree()[n] == 0:
                    g.add_node(n, color='red', size='50,50', fixedsize='True', \
                            shape='circle', fontsize='10', \
                            fontcolor='red')
                    top.append(n)
                else:
                    g.add_node(n, color='red', size='50,50', fixedsize='True', \
                            shape='circle', fontsize='10', \
                            fontcolor='red')
                    middle.append(n)
                
            else:                   # non-outliers
                #if self.g_d.out_degree()[n] == 0:
                d = dist[n]
                if pos is None:
                    g.add_node(n, fillcolor='#%2x%2x%2x'%((1-d/md)*255,(1-d/md)*255,(1-d/md)*255), \
                            size='50,50', fixedsize='True', \
                            shape='circle', fontsize='10', style='filled')
                else:
                    g.add_node(n, fillcolor='#%2x%2x%2x'%((1-d/md)*255,(1-d/md)*255,(1-d/md)*255), \
                            size='50,50', fixedsize='True', \
                            shape='circle', fontsize='10', \
                            style='filled', \
                            pos=str(pos[n][0])+','+str(pos[n][1]))
            
                 
                bottom.append(n)
                elif self.g_d.in_degree()[n] == 0:
                    g.add_node(n, color='black', size='50,50', fixedsize='True', \
                            shape='circle', fontsize='10')
                    top.append(n)
                else:
                    g.add_node(n, color='black', size='50,50', fixedsize='True', \
                            shape='circle', fontsize='10')
                    middle.append(n)
            """

        # set rank
        #for k in layers.keys():
        #    g.subgraph(layers[k], rank='same')
        #g.subgraph(middle, rank='same')
        #g.subgraph(bottom, rank='same')

        # differentiate between activating and inhibiting edges
        activation = []
        inhibition = []
        for e in self.g_d.edges():
            if self.g_d.edge[e[0]][e[1]]['weight'] > 0:
                if width is not None:
                    g.add_edge(e, color='firebrick1', \
                            penwidth='%f'%(width[e[0]][e[1]]*1))
                else:
                    g.add_edge(e, color='firebrick1', penwidth='2')
            elif self.g_d.edge[e[0]][e[1]]['weight'] < 0:
                if width is not None:
                    g.add_edge(e, color='gray50', dir='forward', \
                            arrowhead='tee', \
                            penwidth='%f'%(width[e[0]][e[1]]*1))
                else:
                    g.add_edge(e, color='gray50', dir='forward', \
                            arrowhead='tee', penwidth='2')
        #nx.draw_networkx_edges(self.g_d, pos=l, edgelist=activation,\
        #        edge_color='k')
        #nx.draw_networkx_edges(self.g_d, pos=l, edgelist=inhibition,\
        #        style='dashed')
        #nx.draw_networkx_labels(self.g_d, pos=l, font_size=8)

        #axis('off')
        #axis('equal')
        #suptitle(self.filename)
        #show()
        #set_trace()

        return g

    def plot_degree_distribution(self):
        indeg = self.g_d.in_degree().values()
        outdeg = self.g_d.out_degree().values()
        inmax = max(indeg)
        outmax = max(outdeg)
        bins = range(max(inmax,outmax)+2)
        nin,dummy = histogram(indeg, bins, normed=True)
        nout,dummy = histogram(outdeg, bins, normed=True)
        hin  = bar(bins[:-1], nin, width=0.5, color='r')
        hout = bar([i+0.5 for i in bins[:-1]], nout, width=0.5, \
                color='k')
        legend((hin[0],hout[0]),('in-degree','out-degree'))
        xlabel('degree')
        ylabel('Frequency')
        show()

    def plot_ck(self):
        sum = zeros((max(self.g.degree().values())+1))
        num = zeros((max(self.g.degree().values())+1))
        avg = zeros((max(self.g.degree().values())+1))
        for n in self.g_d.nodes():
            sum[self.g.degree(n)] += nx.clustering(self.g,n)
            num[self.g.degree(n)] += 1
        
        for i,n in enumerate(num):
            if n != 0:
                avg[i] = sum[i]/n
            if avg[i] != 0 and i != 0:
                loglog(i, avg[i], '.k')
        plot([1,10],[1,0.1],'--k')
        show()
       
    def plot_degree_correlation(self):
        in_degree = zeros(len(self.g.nodes()))
        out_degree = zeros(len(self.g.nodes()))
        for i,g in enumerate(self.g.nodes()):
            a_in,i_in,a_out,i_out = self.get_activation_inhibition(g)
            in_degree[i] = a_in + i_in
            out_degree[i] = a_out + i_out
        figure()
        plot(in_degree, out_degree, '.k')
        xlim(xmin=-0.1)
        ylim(ymin=-0.1)
        pearson,p = pearsonr(in_degree, out_degree)
        title(r'$r=%.2f$'%pearson+'\n'+r'$p=%.2e$'%p,fontsize=30)
        show()

    def print_degree_distance(self, out_file, mds, pval_file):
        """Use this function to do the weighted least squares fitting

        """
        genes = mds.genes
        pval = loadtxt(pval_file)
        in_degree = zeros(len(genes))
        out_degree = zeros(len(genes))
        d = zeros(len(genes))
        #cores = find_cores(self.g)
        x = []
        y = []
        symbol = []
        for ii in range(len(genes)):
            #d[ii] = -log10(pval[ii])/max(-log10(pval))
            d[ii] = pval[ii]
            if genes[ii] in self.g.nodes():
                a_in,i_in,a_out,i_out = self.get_activation_inhibition(genes[ii])
                in_degree[ii] = a_in + i_in
                out_degree[ii] = a_out + i_out
                x.append(d[ii])
                #y.append(cores[genes[ii]])
                y.append(in_degree[ii]+out_degree[ii])
                symbol.append(genes[ii])
            else:
                in_degree[ii] = 0
                out_degree[ii] = 0
        #set_trace()
        with open(out_file,'w') as f:
            for ss,yy,xx in zip(symbol,y,x):
                f.write(ss+'\t'+str(yy)+'\t'+str(xx)+'\n')

        # regression
        fitfunc = lambda p, x: p[0]*exp(-p[1]*x)+p[2]
        errfunc = lambda p, x, y, err: (fitfunc(p,x) - y) / err          # Distance to the target function
        p0 = [1,1,1]
        #x = d
        #y = (in_degree + out_degree) #/ max(in_degree + out_degree)
        p1,success = optimize.leastsq(errfunc, p0, args=(array(x), array(y), 
            1/(array(x)*max(y)/max(x))**2+array(y)**2), maxfev=10000)
        return p1

    def plot_degree_distance(self, dist, genes):
        """Plot the anti-correlation between degree and radial distance.

        Parameters
        ----------
        dist : list
            radial distance of each node according to MDS

        genes : list
            list of gene names

        """
        #figure(figsize=(10,9))
        
        #all_in = []
        #all_out = []
        #edges = linspace(0,max(dist),num=10)
        #mean_in_degree = zeros((len(edges)-1))
        #mean_out_degree = zeros((len(edges)-1))
        #
        #for i in range(len(edges)-1):
        #    idx = intersect1d(array(nonzero(dist>edges[i])[0]), \
        #            array(nonzero(dist<=edges[i+1])[0]))
        #    #set_trace()
        #    in_degree = []
        #    out_degree = []
        #    for ii in idx:
        #        if genes[ii] in self.g.nodes():
        #            a_in,i_in,a_out,i_out = self.get_activation_inhibition(genes[ii])
        #            in_degree.append(a_in+i_in)
        #            out_degree.append(a_out+i_out)
        #    if len(idx) == 0:
        #        mean_in_degree[i] = 0
        #        mean_out_degree[i] = 0
        #    else:
        #        mean_in_degree[i] = mean(in_degree)
        #        mean_out_degree[i] = mean(out_degree)
        #half = (edges[1]-edges[0])/2
        #plot(edges[:-1]+half, mean_in_degree, '^r')
        #plot(edges[:-1]+half, mean_out_degree, 'or', markersize=10)

        # scatter plot
        all_dist = []
        all_deg = []
        all_core = []
        self.g.remove_edges_from(self.g.selfloop_edges())
        cores = find_cores(self.g)
        for dd,gg in zip(dist, genes):
            if gg in self.g_d.nodes():
                a_in,i_in,a_out,i_out = self.get_activation_inhibition(gg)
                #all_in.append(a_in+i_in)
                #all_out.append(a_out+i_out)
                #hin = plot(dd, a_in+i_in, '^k')
                #hout = plot(dd, a_out+i_out, '.', markersize=15,\
                #        color=str(0.8-dd/max(dist)*0.8))
                all_dist.append(-log10(dd)/max(-log10(dist)))
                all_deg.append(a_out+i_out+a_in+i_in)
                all_core.append(cores[gg])

        #set_trace()
        #fig = figure()
        #ax = Axes3D(fig)
        #ax.plot(all_dist, all_deg, all_core, '.k')
        #ax.set_xlabel('Response strength')
        #ax.set_ylabel('Degree')
        #ax.set_zlabel('k-Core')
        plot(all_dist, all_deg, '.k', alpha=0.2)
        #legend((hin[0],hout[0]), ('in-degree','out-degree'))
        #xlabel('Radial Distance')
        #ylabel('Out-degree')
        #ax = gca()
        #ax.title.set_y(1.05)
        
        # regression
        fitfunc = lambda p, x: p[0]*exp(-p[1]*x)+p[2]
        errfunc = lambda p, x, y, err: (fitfunc(p,x) - y) / err      
        #errfunc = lambda p, x, y, err: sum((fitfunc(p,x) - y)**2 / err**2)       
        p0 = [1,1,1]
        #p1,success = optimize.leastsq(errfunc, p0, args=(edges[:-1]+half/2,\
        #        mean_in_degree))
        #xfit = arange(0.01,max(dist),0.01)
        #a,b,rr = linreg.linreg(dist,all_in)
        #yfit = fitfunc(p1,edges[:-1]+half/2)
        #in_rr = self.get_rr(mean_in_degree,yfit)
        #in_slope = p1[0]
        #plot(edges[:-1]+half,yfit,'o--r')
        p1,success = optimize.leastsq(errfunc, p0, args=(array(all_dist),\
                array(all_deg), 1/((array(all_dist)*max(all_deg)/max(all_dist))\
                **2+array(all_deg)**2)), maxfev=5000)
        #fit = optimize.fmin_bfgs(errfunc, array(p0), args=(array(all_dist),\
        #        array(all_deg), 1/((array(all_dist)*max(all_deg)/max(all_dist))\
        #        **2+array(all_deg)**2)), full_output=True)
        #set_trace()
        xfit = arange(0.01,max(all_dist),0.01)
        yfit = fitfunc(p1, xfit)
        #out_rr = self.get_rr(mean_out_degree,yfit)
        #out_slope = p1[0]
        #pearson,p = pearsonr(edges[:-1]+half,mean_out_degree)
        plot(xfit, yfit, '--r', linewidth=4)
        #text(0.4,max(all_deg)/2,r'$r=%.2f$'%pearson+'\n'+r'$p=%.2e$'%p,fontsize=30)

        #rcParams['text.usetex'] = True
        #title('E. coli (full-scale)', fontweight='bold')
        #ax = gca()
        #ax.title.set_y(1.05)

        # get log-normal distribution
        #set_trace()
        #x = 10**np.random.normal(mean(log10(all_dist)), \
        #        std(log10(all_dist)), 100)
        #y = 10**np.random.normal(mean(log10(array(all_deg)+1)), \
        #            std(log10(array(all_deg)+1)), 100)
        #savetxt('/Users/bao/work/DREAM/DREAM3 in silico challenge/Size100/mds/ecoli1_lognormal', c_[x,y])

        #hAxes = axes()
        #xlim(xmin=-0.01)
        #ylim(-1,max(all_deg)+1)
        #xlabel('Response strength', fontweight='bold')
        #ylabel('Out-degree', fontweight='bold')
        #hAxes.xaxis.set_label_coords(0.5,-0.07) 
        #
        show()
        
        return p1

    def get_rr(self, y, yfit):
        ybar = mean(y)
        sserr = sum([(yfit[i]-y[i])**2 for i in range(len(yfit))])
        sstot = sum([(ybar-yi)**2 for yi in y])
        return 1-sserr/float(sstot)
    
    def plot_degree(self):
        for n in self.g_d.nodes():
            a_in,i_in,a_out,i_out = self.get_activation_inhibition(n)
            plot(a_out+i_out, a_in-i_in, '.k')
        xlim(xmin=-1)
        ylim(-10,10)
        xlabel(r'$k_o$')
        ylabel(r'$k_{ia} - k_{ii}$')
        show()

    def plot_outdeg_zero_hist(self):
        # get degrees for all nodes
        deg = nx.degree(self.g).values()
        node = nx.degree(self.g).keys()
        idx_sorted = sorted(range(len(deg)), key = deg.__getitem__)
        idx_hub = idx_sorted[-len(deg)/4:]
        idx_output = nonzero(array(self.g_d.out_degree().values())==0)[0]

        delta_hub = []
        delta_output = []
        for i in idx_hub:
            a_in,i_in,a_out,i_out = self.get_activation_inhibition(node[i])
            delta_hub.append(a_in-i_in)
        for i in idx_output:
            a_in,i_in,a_out,i_out = self.get_activation_inhibition(node[i])
            delta_output.append(a_in-i_in)
        hist(delta_hub, 20, color='k', histtype='step')
        hist(delta_output, 20, color='r', histtype='step')
        xlabel(r'$k_{ia} - k_{ii}$')
        #show()

    def plot_ck(self):
        c = nx.clustering(self.g)
        k_unique = unique(self.g.degree().values())
        ck = zeros(len(k_unique))
        for ii,k in enumerate(k_unique):
            idx = nonzero(array(self.g.degree().values()).flatten()==k)[0]
            sum = 0
            for i in idx:
                sum += c[self.g.degree().items()[i][0]]
            ck[ii] = sum/len(idx)

        plot(k_unique, ck, '.k', markersize=20)
        xlabel(r'$k$')
        ylabel(r'$C(k)$')

        # regression
        fitfunc = lambda p, x: p[0]*x+p[1]
        errfunc = lambda p, x, y: fitfunc(p,x) -y          # Distance to the target function
        p0 = [-1,1]
        p1,success = optimize.leastsq(errfunc, p0, args=(k_unique,\
                ck))
        yfit = fitfunc(p1, k_unique)
        pearson,p = pearsonr(k_unique, ck)
        plot(k_unique, yfit, '--r', linewidth=4)

        return pearson,p

    def plot_dendrogram(self):
        Y = squareform(self.get_distance_matrix())
        Z = linkage(Y, method='average')
        dendrogram(Z, labels=self.g.nodes(), leaf_font_size=10, \
                orientation='top')

    def get_distance_matrix(self):
        
        L = zeros((self.g_d.order(),self.g_d.order()))
        
        for i,n in enumerate(self.g_d.nodes()):
            for ii,nn in enumerate(self.g_d.nodes()):
                L[i,ii] = nx.shortest_path_length(self.g,n,nn)
                L[ii,i] = L[i,ii]
        
        return L

    def draw(self, dist, genes):
        figure()
        nc = []
        for dd,gg in zip(dist, genes):
            if gg in self.g_d.nodes():
                nc.append(dd)
        nx.draw_graphviz(self.g_d, node_color=nc, cmap=cm.Blues)

        """
        boxplot(d)
        label = []
        for i in range(len(files)):
            label.append(files[i]+'\n%.3f'%(c[i]))
        xticks(range(1,6), label)
        xlabel('Average clustering coefficient')
        ylabel('Radial distance')
        show()

        f.set_size_inches((10,8))
        savefig(wd+'mds/distance_cc.png', dpi=800)


        f = figure()

        f.set_size_inches((16,10))
        savefig(wd+sub('top','layout',outlier_file)+'.pdf', dpi=800)
"""

def plot_param_dist(params, mds):
    sum = {}
    num = {}
    param = params.n
    for pre in param.keys():
        for post in param[pre].keys():
            #set_trace()
            if not sum.has_key(post):
                sum[post] = param[pre][post]
                num[post] = 1
            else:
                sum[post] += param[pre][post]
                num[post] += 1


    dist = mds.get_distance()
    for i in range(len(dist)):
        if sum.has_key(mds.genes[i]):
            p = sum[mds.genes[i]]/num[mds.genes[i]]
            plot(dist[i], p, '.k')
    xlabel('Radial Distance', fontweight='bold')
    ylabel('Hill coefficient', fontweight='bold')
    #rc('text', usetex=True)
    font0 = FontProperties()
    font = font0.copy()
    font.set_style('italic')
    font.set_weight('bold')
    font.set_size('xx-large')
    title(r"E. coli", fontproperties=font)

    ax = gca()
    ax.title.set_y(1.05)

    hAxes = axes()    
    hAxes.xaxis.set_label_coords(0.5,-0.078) 
    hAxes.yaxis.set_label_coords(-0.1,0.5) 
    
    show()


def generate_graph(type):
    size = 1000
    degree = int(size*0.0122)
    n_networks = 100
    #type = 'random'
    for n in range(n_networks):
        print n,
        sys.stdout.flush()
            
        if type == 'uniform':
            connected = False
            while not connected:
                g = nx.directed_configuration_model([int(i) for i in \
                    rand(size)*degree],\
                    [int(i) for i in rand(size)*degree], \
                    create_using=nx.DiGraph())
                connected = nx.algorithms.components.is_connected(g.to_undirected())
        #g = nx.connected_watts_strogatz_graph(size, degree, 0.1)
        elif type == 'scalefree':
            #g = nx.barabasi_albert_graph(size, degree)
            # in-degree: Poisson; out-degree: power-law
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
        elif type == 'random':
            connected = False
            while not connected:
                g = nx.gnp_random_graph(size, 0.003, directed=True)
                connected = nx.algorithms.components.is_connected(g.to_undirected())
        elif type == 'tree':
            g = nx.random_powerlaw_tree(size, tries=5000, seed=42)
        elif type == 'bowtie':
            connected = False
            while not connected:
                # bow-tie graph
                G1 = nx.gnr_graph(size/2,0.1,create_using=nx.DiGraph())
                G2 = G1.reverse()
                G2 = nx.relabel_nodes(G2, lambda n: n+G1.order())
                g = nx.union(G1, G2)
                g.add_edge(0,G1.order())
                connected = nx.algorithms.components.is_connected(g.to_undirected())
        elif type == 'smallworld':
            g = nx.connected_watts_strogatz_graph(size, size/200, 0.1)
            g = g.to_directed()
            for u in g.nodes():
                for v in g.nodes():
                    if g.has_edge(u,v) and g.has_edge(v,u):
                        if rand() > 0.5:
                            g.remove_edge(u, v)
                        else:
                            g.remove_edge(v, u)
        

        g = nx.relabel_nodes(g, dict(zip(g.nodes(), ['G'+str(i+1) for i in \
                range(size)])))
        wd = '/home/jbao/data/DREAM/gnw/'+type+'2/gnw/Size'+str(size)+'/input/'
        if not os.path.isdir(wd):
            os.makedirs(wd)
        #nx.write_dot(g, wd+type+'_'+str(size)+'.dot')
        nx.write_edgelist(g, wd+type+'2-'+str(n+1)+'_'+str(size)+'.tsv', \
                data=False, delimiter='\t')
        #return g
        g.clear()

def remove_node():
    topology_file = '/Users/bao/work/DREAM/DREAM3_in_silico_challenge/Size100/Networks/InSilicoSize100-Ecoli1.tsv'
    topology = loadtxt(topology_file, dtype=str)

    mds = data.Data('/Users/bao/work/DREAM/DREAM3_in_silico_challenge/Size100/mds/pca_ecoli1_100_normed_euc','/Users/bao/work/DREAM/DREAM3_in_silico_challenge/Size100/mds/ecoli1_100_normed.dat')
    genes = mds.get_ranking()

    for i,g in enumerate(genes[:-3]):
        idx = union1d(nonzero(topology[:,0]==g)[0], \
                nonzero(topology[:,1]==g)[0])
        topology = delete(topology, idx, 0)
        savetxt(topology_file.split('.')[0]+'_'+str(i+1)+'.tsv', topology,\
                fmt='%s', delimiter='\t')

def workflow():
    wd = '/Users/bao/work/DREAM/gnw/ecoli/gnw/Size1000/'
    edge_file = '/Users/bao/work/DREAM/gnw/ecoli/gnw/Size1000/ecoli-1_goldstandard_signed.tsv.sub'
    #g = nx.read_edgelist(edge_file, delimiter='\t', create_using=nx.DiGraph())
    #nx.draw_graphviz(g, node_color='0.5', with_labels=False, node_size=100, edge_color='k')

    g = Graph(edge_file)
    gene_file = '/Users/bao/work/DREAM/gnw/ecoli/gnw/Size1000/ecoli-1_1000_normed.dat'
    mds_file = '/Users/bao/work/DREAM/gnw/ecoli/gnw/Size1000/pca_ecoli-1_1000_normed_euc'
    mds = data.Data(mds_file, gene_file)
    dist,dummy = mds.get_distance()
    #pos = nx.graphviz_layout(g.g_d)
    #gg = g.get_dot_graph(dist, pos=pos)
    #gg.write(wd+'ecoli-1_1000.dot')

    #idx = argsort(dummy)
    #x = mds.time_points
    #y = range(1000)
    #X,Y = meshgrid(x, y)
    #pcolor(X, Y, mds.time_series[idx,:], cmap=cm.RdBu_r)
    #xlabel('Time (a.u.)', fontsize='xx-large', fontweight='bold')
    #ylabel('# Genes', fontsize='xx-large', fontweight='bold')
    #cb = colorbar()
    #cb.set_label(r'$\log_2$ fold expression', fontsize='xx-large')

    qval = loadtxt(wd+'ecoli-1_1000_pvals.txt', skiprows=1, usecols=[27], delimiter='|')
    mds.plot_2d(label=find(qval<0.01))

    show()

def ecoli():
    edge_file = '/home/jbao/gnw/src/ch/epfl/lis/gnwgui/rsc/net/ecoli_transcriptional_network_regulonDB_6_2.tsv.sub'
    #edge_file = '/home/jbao/data/DREAM/gnw/yeast/gnw/full/yeast-all_goldstandard_signed.tsv.sub'
    #edge_file = '/home/jbao/data/DREAM/gnw/random/gnw/Size1000/random-1_goldstandard_signed.tsv'
    #edge_file = '/home/jbao/data/DREAM/gnw/ecoli/gnw/Size1000/default/rewiring/ecoli-1_1000_1_goldstandard_signed.tsv'
    #edge_file = '/home/jbao/data/DREAM/gnw/smallworld/gnw/Size1000/smallworld-1_goldstandard_signed.tsv'
    #system('sed \"s/+-/{\'weight\':1}/g;s/\t-/\t{\'weight\':-1}/g;'+\
    #        's/\t1/\t{\'weight\':1}/g;s/\t+/\t{\'weight\':1}/g;'+
    #        's/\t?/\t{\'weight\':1}/g\" '+\
    #        '< %s > %s.sub' % (edge_file, edge_file))
    g = Graph(edge_file)
     
    #set_trace()
    #cond = ['cdc15']
    cond = ['cold', 'heat', 'lactose', 'oxidative']
    #cond = ['alpha', 'cdc15', 'cdc28', 'elu']
    wd = '/home/jbao/data/ecoli/'
    #colormap = cm.Greys
    #colors = [colormap(i) for i in np.linspace(0.3, 0.9, 1000)]
    #colors = ['#348ABD', '#7A68A6', '#A60628', '#467821']
    #x = zeros(1000*len(g.g_d.nodes()))
    #y = zeros(1000*len(g.g_d.nodes()))
    #counter = 0
    for c in cond:
        #print i,
        #sys.stdout.flush()
        gene_file = '/home/jbao/data/ecoli/GSE20305_'+c+'.dat'
        mds_file = '/home/jbao/data/ecoli/pca_GSE20305_'+c+'_4400_normed_euc'
        pval_file = '/home/jbao/data/ecoli/mds/pval_GSE20305_'+c+'_4400.dat'
        deg_dist_file = '/home/jbao/data/ecoli/deg_dist_GSE20305-'+c
        #gene_file = '/home/jbao/data/yeast/yeast_'+c+'.dat'
        #mds_file = '/home/jbao/data/yeast/pca_yeast_'+c+'_6178_normed_euc'
        #pval_file = '/home/jbao/data/yeast/pval_yeast_'+c+'_6178.dat'
        #deg_dist_file = '/home/jbao/data/yeast/kshell_dist_yeast_'+c+'_6178'
        #gene_file = '/home/jbao/data/DREAM/gnw/random/gnw/Size1000/mds/random-1_perturbation-'+str(i+1)+'_1000_normed.dat'
        #mds_file = '/home/jbao/data/DREAM/gnw/random/gnw/Size1000/mds/pca_random-1_perturbation-'+str(i+1)+'_1000_normed_euc'
        #perturbation_file = '/home/jbao/data/DREAM/gnw/ecoli/gnw/full/ecoli-full_perturbation-'+str(i+1)+'_multifactorial_timeseries_perturbations.tsv'
        #perturbation = loadtxt(perturbation_file, dtype=str)
        #set_trace()
        mds = data.Data(mds_file, gene_file)
        #dist,dummy = mds.get_distance()
        #g.plot_degree_distance(dummy, mds.genes)
        p = g.print_degree_distance(deg_dist_file, mds, pval_file)
        #savetxt(wd+'param_kshell_ecoli_'+c+'.dat', p)
    
        # plot degree vs. perturbation
        #for gg in g.g_d.nodes():
        #    idx = perturbation[0,:].tolist().index('"'+gg+'"')
        #    x[counter] = float(perturbation[1,idx])
        #    y[counter] = g.g_d.out_degree()[gg] + g.g_d.in_degree()[gg]
        #    counter += 1
    #scatter(x, y)
    #show()
    
    """
        flag = 0
        x = []
        y = []
        for gg in mds.genes:
            if gg in g.g_d.out_degree():
                flag += 1
                x.append(dist[gg])
                y.append(g.g_d.out_degree()[gg]+g.g_d.in_degree()[gg])
        plot(x, y, 'o', color='k')
        #legend(frameon=False, numpoints=1)
        
        # regression
        fitfunc = lambda p, x: p[0]*exp(-p[1]*x+p[2])+p[3]
        errfunc = lambda p, x, y, err: (fitfunc(p,x) - y) / err          # Distance to the target function
        p0 = [1,1,1,1]
        p1,success = optimize.leastsq(errfunc, p0, args=(array(x),
                #array(y), 1/((array(x)*max(y)/max(x))**2+array(y)**2)), maxfev=5000)
                array(y), 1/(array(y)**2)), maxfev=5000)
        xfit = arange(0.01,max(x),0.01)
        yfit = fitfunc(p1, xfit)
        plot(xfit, yfit, '--', color='k', linewidth=2)
        
    #yscale('log')
    xlim(xmin=-0.01)
    ylim(-10, max(y)+50)
    xlabel('Response strength', fontweight='bold', fontsize='xx-large')
    ylabel('Degree', fontweight='bold', fontsize='xx-large')
    title(r'Random Network', 
    #title(r'Time resolved response of $\mathit{E. coli}$ to stress', 
    #title(r'Time resolved response of $\mathit{S. cerevisiae}$ to synchronization', 
            fontweight='bold', fontsize='large')
    ax = gca()
    ax.title.set_y(1.05)
    show()
    """
    #return x,y

def yeast():
    size = 4441
    wd = '/Users/bao/work/DREAM/gnw/yeast/gnw/full/'
    edge_file = wd + 'yeast_' + str(size) + '_goldstandard_signed.tsv'
    system('sed \"s/+-/{\'weight\':1}/g;s/\t-/\t{\'weight\':-1}/g;'+\
            's/\t1/\t{\'weight\':1}/g;s/\t+/\t{\'weight\':1}/g\" '+\
            '< %s > %s.sub' % (edge_file, edge_file))
    gene_file = wd + 'mds/yeast_' + str(size) + '_normed.dat'
    mds_file = wd + 'mds/pca_yeast_' + str(size) + '_normed_euc'

    mds = data.Data(mds_file, gene_file)
    dist,dummy = mds.get_distance()
    g = Graph('%s.sub' % edge_file)

    essential_orf = np.loadtxt('/Users/bao/work/yeast/Essential_ORFs.txt', delimiter='\t', skiprows=2, usecols=[1], dtype=str)
    
    #from mpl_toolkits.mplot3d import Axes3D
    #fig = figure()
    #ax = Axes3D(fig)
    
    od = g.g_d.out_degree()
    id = g.g_d.in_degree()
    genes = sorted(dist,key=dist.get)
    response = []
    out_degree = []
    in_degree = []
    for k in genes:
        response.append(dist[k])
        out_degree.append(od[k])
        in_degree.append(id[k])
        color_scale = 0.8 - dist[k]/max(dist.values())*0.8
        #ax.plot([od[k]], [id[k]], [dist[k]], color=str(color_scale))

    scatter(out_degree, in_degree, c=response, cmap=cm.jet, edgecolor='none')
    cb = colorbar()
    cb.set_label(r'Response strength')
    xlabel('Out-degree')
    ylabel('In-degree')
    title('Yeast')
    """
    degree = out_degree
    edge = arange(0, max(degree)+25, 25)
    #edge = arange(0, max(degree)+2, 2)
    half = (edge[1]-edge[0])/2
    #window_size = max(degree)/5.0
    #step = max(degree)/50.0
    #n_windows = ceil(max(degree)/step)
    #fraction = zeros((len(out_edge)-1,len(in_edge)-1))
    fraction = zeros(len(edge)-1)
    #for i in range(len(out_edge)-1):
    for i in range(len(edge)-1):
        idx = intersect1d(array(nonzero(array(degree)>=edge[i])[0]), \
                array(nonzero(array(degree)<edge[i+1])[0]))
        #idx = intersect1d(array(nonzero(array(degree)>=i*step)[0]), \
        #        array(nonzero(array(degree)<=i*step+window_size)[0]))
        #idx = nonzero(array(degree)>edge[i])[0]
        #idx = nonzero(array(degree)<=edge[i+1])[0]
        #set_trace()
        #idx = intersect1d(out_idx, in_idx)
        if len(idx) == 0:
            fraction[i] = 0
        else:
            fraction[i] = len([gene[ii] for ii in idx if gene[ii] in essential_orf])/\
                    float(len(idx))

    #genes = mds.genes[:size]
    #dist = mds.get_distance()
    #orf,hl,dummy = half_life.get_half_life()
    #figure()
    #common = intersect1d(genes, orf)
    #x = zeros(len(genes))
    #y = zeros(len(genes))
    #color_scale = zeros(len(common))
    #bc = nx.betweenness_centrality(g.g_d)
    #core = find_cores(g.g)
    #ref = sorted(dist,key=dist.get,reverse=True)[0]
    #set_trace()
    #plot([g.g_d.out_degree()[i] for i in common],[hl[orf.tolist().index(i)] for i in common],'.k')
    #for i,c in enumerate(genes):
        #if genes[i] in orf.tolist():
        #diff = log(max(hl)) - log(min(hl))
        #color_scale[i] = hl[orf.tolist().index(c)]
        #x[i] = dist[c]
        #y[i] = core[c]
        #if nx.shortest_path(g.g,ref,c):
        #    y[i] = nx.shortest_path_length(g.g,ref,c)
        #elif nx.shortest_path(g.g,c,ref):
        #    y[i] = nx.shortest_path_length(g.g,c,ref)
        #else:
        #    y[i] = 0
  
    #set_trace()
    #imshow(fraction, cmap=cm.Greys, origin='lower')
    #xticks(arange(len(in_edge)), in_edge)
    #yticks(arange(len(out_edge)), out_edge)
    #axis('auto')
    #plot(x[color_scale>0],color_scale[color_scale>0],'.k')
    x = edge[:-1]+half
    y = fraction
    plot(x, y, '.k', markersize=10)
    #pearson,p = pearsonr(edge[:-1], fraction)
    #title(r'$r=%.2f$'%pearson+'\n'+r'$p=%.2e$'%p,fontsize=30)
    #hAxes = axes()
    xlim(xmin=-max(x)/30,xmax=max(x)+max(x)/30)
    ylim(ymin=-max(y)/30,ymax=max(y)+max(y)/30)
    
    # regression
    fitfunc = lambda p, x: p[0]*x+p[1]
    errfunc = lambda p, x, y: fitfunc(p,x) -y          # Distance to the target function
    p0 = [-1,1]
    p1,success = optimize.leastsq(errfunc, p0, args=(x,y))
    xfit = linspace(min(x), max(x), 20)
    yfit = fitfunc(p1, xfit)
    pearson,p = pearsonr(x, y)
    plot(xfit, yfit, '--r', linewidth=4)
    text(min(x), 0.8*max(y), r'$r=%.2f$'%pearson+'\n'+r'$p=%.2e$'%p,\
            fontsize=30)
    xlabel('Out-degree', fontweight='bold')
    ylabel('Fraction of essential genes', fontweight='bold')
    title('Yeast network', fontweight='bold')
    ax1 = gca()
    ax1.title.set_y(1.05)
    #hAxes.xaxis.set_label_coords(0.5,-0.078) 
    #colorbar()
    """
    show()
    #return fraction

def dream():
    size = 10
    wd = '/Users/bao/work/DREAM/DREAM3_in_silico_challenge/Size'+str(size)+'/'
    edge_file = wd + 'Networks/InSilicoSize'+str(size)+'-Ecoli1.tsv'
    system('sed \"s/+-/{\'weight\':1}/g;s/\t-/\t{\'weight\':-1}/g;'+\
            's/\t1/\t{\'weight\':1}/g;s/\t+/\t{\'weight\':1}/g\" '+\
            '< %s > %s.sub' % (edge_file, edge_file))
    gene_file = wd + 'mds/ecoli-1_' + str(size) + '_normed.dat'
    mds_file = wd + 'mds/pca_ecoli-1_' + str(size) + '_normed_euc'

    mds = data.Data(mds_file, gene_file)
    dist,dummy = mds.get_distance()
    g = Graph('%s.sub' % edge_file)

    # Random rewiring. For each edge u-v, with probability p, randomly
    # replace with edge u-w.
    import random
    e = g.g_d.edges()
    D = nx.DiGraph()
    for (u, v) in e:
        if random.random() < 0.1:
            newv = random.choice(g.g_d.nodes())
            # avoid self-loops and reject if edge u-newv exists
            # is that the correct WS model?
            while newv == u or g.g_d.has_edge(u, newv):
                newv = random.choice(g.g_d.nodes())
            #D.add_edge(u,v)  # conserve number of edges
            D.add_edge(u,newv)
        else:
            D.add_edge(u,v)

    rewire_file = wd + 'Networks/InSilicoSize'+str(size)+'-Ecoli1_rewiring.tsv'
    system('sed \"s/+-/{\'weight\':1}/g;s/\t-/\t{\'weight\':-1}/g;'+\
        's/\t1/\t{\'weight\':1}/g;s/\t+/\t{\'weight\':1}/g\" '+\
        '< %s > %s.sub' % (rewire_file, rewire_file))
    g_rewire = Graph('%s.sub' % rewire_file)

    #g.print_degree_distance(wd+'mds/ecoli-1_'+str(size)+'_deg_dist',mds)
    #g.plot_degree_correlation()
    #g.plot_degree_distance(dummy, mds.genes)
    pos = nx.graphviz_layout(g.g_d)
    #set_trace()
    gg = g.get_dot_graph(dist, pos=pos)
    gg.write(wd+'mds/ecoli1.dot')
    gg = g_rewire.get_dot_graph(dist, pos=pos)
    gg.write(wd+'mds/ecoli1_rewiring.dot')
    # permutation
    perm = permutation(dist.values())
    for i,k in enumerate(dist.keys()):
        dist[k] = perm[i]
    gg = g.get_dot_graph(dist, pos=pos)
    gg.write(wd+'mds/ecoli1_param.dot')
    return g

def plot_input():
    #wd = '/Users/bao/work/DREAM/DREAM3_in_silico_challenge/Size100/gnw/'
    wd = '/home/jbao/data/DREAM/gnw/ecoli/gnw/full/'
    #mds = data.Data(wd+'mds/pca_ecoli-full_perturbation-1_1502_normed_euc',wd+'mds/ecoli-full_perturbation-1_1502_normed.dat')
    #dist,dummy = mds.get_distance()
    #rank = argsort(dummy).tolist()
    #rank.reverse()
    #ranked = sorted(dist,key=dist.get,reverse=True)
    p = pl.ParamLoader()
    influx = p.load_input(wd+'ecoli-full_perturbation-1_inputs.tsv')
    efflux = p.load_input(wd+'ecoli-full_perturbation-1_outputs.tsv')
    #p.load_param(wd+'ecoli-full.xml')
    #p.load_jacobian(wd+'random-1_jacobian.tsv')
    #p.load_timeseries(wd+'ecoli-10_100_1_multifactorial_timeseries.tsv')
    #p.load_timeseries('/Users/bao/work/DREAM/DREAM3_in_silico_challenge/Size100/DREAM3 data/InSilicoSize100-Ecoli1-trajectories.tsv')
    #p.load_timeseries('/Users/bao/test/InSilicoSize100-Ecoli1_multifactorial_timeseries.tsv')
    #p.load_timeseries(wd+'test_multifactorial_timeseries.tsv')
    #set_trace()
    #order = [p.gene.index(r) for r in ranked]
    
    g = Graph(wd+'ecoli-full_goldstandard_signed.tsv.sub')
    degree = sorted(g.g_d.degree(), key=g.g_d.degree().get, reverse=True)
  
    # plot influx/efflux
    flux = influx + efflux
    for i,gg in enumerate(degree):
        idx = p.gene.index(gg)
        flux[0,:,i] = influx[0,:,idx] + efflux[0,:,idx]
    imshow(sqrt(flux[0,:,:100].T), interpolation='nearest', origin='upper',
            aspect='auto')
    xticks(linspace(0,len(p.time_point)-1,5), 
            p.time_point[[int(i) for i in linspace(0,len(p.time_point)-1,5)]])
    #yticks(linspace(0,80,5), array(degree)[[int(i) for i in linspace(0,80,5)]])
    xlabel('Time (a.u.)')
    #axis('auto')
    colorbar()

    # degree vs. jacobian
    #deg = zeros(len(degree))
    #eigval = zeros(len(degree))
    #response = zeros(len(degree))
    #val,vec = eig(p.jacobian)
    #for i,gg in enumerate(degree):
    #    if gg in p.jacobian_gene:
    #        deg[i] = g.g_d.degree()[gg]
    #        eigval[i] = val[p.jacobian_gene.index(gg)]
    #        #response[i] = dist[gg]
    #        #text(eigval[i], deg[i], gg)
    ##set_trace()
    #plot(eigval, deg, '.k')
    ##scatter(eigval, deg, c=response, cmap=cm.hot)
    #return eigval,deg,response

    # moving average
    #ma = zeros((len(g.g_d.nodes()), len(p.time_point)))
    #from scipy.interpolate.polyint import pchip
    #from scipy.signal import cspline1d, cspline1d_eval
    #for i,gg in enumerate(degree):
    #    if gg in p.delta:
    #        delta = 1/p.delta[gg]/2
    #        idx = p.gene.index(gg)
    #        #print gg,g.g_d.degree()[gg]
    #        #c = pchip(p.time_point, p.input[0,:,idx])
    #        cj = cspline1d(p.input[0,:,idx])
    #        for j,tt in enumerate(p.time_point):
    #            if tt < delta:
    #                left = 0
    #            else:
    #                left = tt - delta

    #            if tt + delta > p.time_point[-1]:
    #                right = p.time_point[-1]
    #            else:
    #                right = tt + delta

    #            newx = linspace(left, right, 10)
    #            #interpol = c.__call__(x)
    #            newy = cspline1d_eval(cj, newx, dx=p.time_point[1]-
    #                    p.time_point[0], x0=p.time_point[0])
    #            integral = np.trapz(newy, newx)
    #            ma[i,j] = integral

    #return ma

    #figure()
    #to_plot = mean(p.input,0)
    #var_input = zeros(len(mds.genes))
    #param = zeros(len(mds.genes))
    #auto_corr = zeros((len(mds.genes),41))
    #for i in range(p.timeseries.shape[2]):
    #    color_scale = rank.index(i)
    #    c = cm.hot(color_scale/float(len(mds.genes)),1)
    #    plot(p.time_point, mean(p.timeseries,0)[:,i].T,color=c)
    """
        #set_trace()
        to_plot[:,i]=to_plot[:,i]/np.max(to_plot,0)[i]
        t = acorr(to_plot[:,i],color=c,usevlines=False,maxlags=20,linestyle='-',marker=None,lw=2)
        #set_trace()
        # regression
        fitfunc = lambda p, x: exp(-x/p[0])
        errfunc = lambda p, x, y: fitfunc(p,x) -y          # Distance to the target function
        p0 = [1]
        param[i],success = optimize.leastsq(errfunc, p0, args=(t[0][20:],\
                t[1][20:]))
        #var_input[i] = var(to_plot[:,i])
        #plot(p.time_point,to_plot[:,i],color=c)
    figure()
    plot(dummy,param,'.k',markersize=15)
    xlabel('Response strength',fontweight='bold')
    ylabel(r'$\tau_{rr}$',fontweight='bold')
    title('100-gene Ecoli network',fontweight='bold')
    #ylim(ymin=-0.01)
    #xticks(range(-20,21,5),range(-200,201,50))
    ax1 = gca()
    ax1.title.set_y(1.05)
    hAxes = axes()
    hAxes.xaxis.set_label_coords(0.5,-0.078) 
    #imshow(to_plot.T,origin='upper')
    #axis('auto')
    #colorbar()
    """
    show()
    return flux


class NetworkLoader:
    def __init__(self, dataset, i, size):
        """
        if dataset == 'yeast':
            self.wd = '/Users/bao/work/DREAM/gnw/yeast/'
            self.n_networks = 10
            self.edge_file = 'Yeast-'+str(i+1)+'_goldstandard_signed.tsv'
            self.gene_file = 'mds/yeast-'+str(i+1)+'_%d_normed.dat'%size
            self.mds_file = 'mds/pca_yeast-'+str(i+1)+'_%d_normed_euc'%size
        """
        if dataset == 'dream':
            files = ['ecoli1','ecoli2','yeast1','yeast2','yeast3']
            self.wd = "/Users/bao/work/DREAM/DREAM3\ in\ silico\ challenge/Size100/"
            self.edge_file = "DREAM3\ gold\ standards/"+files[i]+"_100.txt"
            self.gene_file = "mds/"+files[i]+"_100_normed.dat"
            self.mds_file = "mds/pca_"+files[i]+"_100_normed_euc"
            self.param_file = "Networks/InSilicoSize100-"+files[i].capitalize()+".xml"
            self.n_networks = len(files)
        elif dataset == 'full':
            size = 1502
            self.wd = '/home/jbao/data/DREAM/gnw/ecoli/gnw/full/'
            self.edge_file = 'ecoli-full_goldstandard_signed.tsv'
            starts = [match.start() for match in finditer('_',self.edge_file)]
            self.gene_file = 'mds/'+self.edge_file[:starts[-2]]+'_perturbation-'+str(i+1)+'_'+str(size)+'_normed.dat'
            self.mds_file = 'mds/pca_'+self.edge_file[:starts[-2]]+'_perturbation-'+str(i+1)+'_'+str(size)+'_normed_euc'
            self.pval_file = 'mds/pval_'+self.edge_file[:starts[-2]]+'_perturbation-'+str(i+1)+'_'+str(size)+'.dat'

            self.deg_dist_file = 'mds/deg_dist_ecoli-full_perturbation-'+str(i+1)+'_'+\
                    str(size)
            #self.size = int(self.edge_file.split('_')[1])
            self.size = size
            self.n_networks = 1000
        else:
            self.idx = sys.argv[2]
            #self.idx = int(os.getenv('SGE_TASK_ID')) 
            self.wd = '/home/jbao/data/DREAM/gnw/'+dataset+'/gnw/Size'+\
                    str(size)+'/'
            self.edge_file = dataset+'-'+str(self.idx)+'_goldstandard_signed.tsv'
            self.gene_file = 'mds/'+dataset+'-'+str(self.idx)+'_perturbation-'+\
                    str(i+1)+'_'+str(size)+'_normed.dat'
            self.mds_file = 'mds/pca_'+dataset+'-'+str(self.idx)+'_perturbation-'+\
                    str(i+1)+'_'+str(size)+'_normed_euc'
            self.deg_dist_file = 'mds/deg_dist_'+dataset+'-'+str(self.idx)+\
                    '_perturbation-'+str(i+1)+'_'+str(size)
            self.pval_file = 'mds/pval_'+dataset+'-'+str(self.idx)+\
                    '_perturbation-'+str(i+1)+'_'+str(size)+'.dat'
            self.n_networks = 1000
            self.size = size


if __name__=="__main__":

    #close('all')
    #remove_node()
    #workflow()
    #ecoli()
    #yeast()
    #g=dream()
    #flux = plot_input()

    #"""
    #dataset = 'dream'
    #dataset = 'yeast'
    #dataset = 'bowtie'
    dataset = sys.argv[1]

    #wd = '/Users/bao/work/DREAM/gnw/test/'
    #files = ['ecoli1','ecoli2']
    #files = ['ecoli1','ecoli2','yeast1','yeast2','yeast3']
    #type = ['dream','smallworld','scalefree','random','uniform']

    #all_size = [10,100,1000,2000]
    size = 1000
    data0 = NetworkLoader(dataset,0,size)
    #all_wd = ['/Users/bao/work/DREAM/gnw/yeast/size10/',\
    #        '/Users/bao/work/DREAM/gnw/yeast/',\
    #        '/Users/bao/work/DREAM/gnw/yeast/size1000/',\
    #        '/Users/bao/work/DREAM/gnw/yeast/size2000/']

    #close('all')

    #n_networks = 10
    #indeg = zeros((n_nodes,n_networks))
    #outdeg = zeros((n_nodes,n_networks))
    #modularity = zeros((n_networks))
    #rr = zeros((n_networks))
    #n = zeros((n_networks))
    #c = zeros(5)
    #f = figure()
    #cc = zeros((2,5))
    #fig = figure()
   
    a = zeros(data0.n_networks)
    b = zeros(data0.n_networks)
    tau = zeros(data0.n_networks)
    #wd = all_wd[1]
    iwd = 0
    #for wd,size in zip(all_wd,all_size):
    d = zeros((data0.size,data0.n_networks))
    for i in range(data0.n_networks):
        #i = int(os.getenv('SGE_TASK_ID')) - 1
        print i,
        sys.stdout.flush()
            
        mydata = NetworkLoader(dataset,i,size)

        #edge_file = 'uniform-'+str(i+1)+'_goldstandard_signed.tsv'
        system('sed \"s/+-/{\'weight\':1}/g;s/\t-/\t{\'weight\':-1}/g;'+\
                's/\t1/\t{\'weight\':1}/g;s/\t+/\t{\'weight\':1}/g\" '+\
                '< %s > %s.sub' % (\
                mydata.wd+mydata.edge_file, mydata.wd+mydata.edge_file))
        #system('sed \"s/\t-/\t{\'weight\':-1}/g\" < %s.1 > %s.2' % (\
        #        wd+edge_file, wd+edge_file))
        #gene_file = 'uniform-'+str(i+1)+'_100_normed.dat'
        #outlier_file = 'mds/top_'+files[i]+'_100_ellipse'
        #mds_file = 'mds/pca_uniform-'+str(i+1)+'_100_normed_euc'

        mds = data.Data(sub(r'\\','',mydata.wd+mydata.mds_file), \
                sub(r'\\','',mydata.wd+mydata.gene_file))
        mds.genes = mds.genes[:data0.size]
        #dummy,d[:,i] = mds.get_distance('treated','pca')
        #idx = argsort(d[:,i]).tolist()
        #idx.reverse()

        g = Graph('%s.sub' % (mydata.wd+mydata.edge_file))

        #g.draw(d[:,i], mds.genes[:size])

        #figure()
        pval = loadtxt(mydata.wd+mydata.pval_file)
        #p = g.plot_degree_distance(pval, mds.genes)
        p = g.print_degree_distance(sub(r'\\','',mydata.wd+\
                mydata.deg_dist_file), mds, mydata.wd+mydata.pval_file)
        #savetxt('/home/jbao/data/DREAM/gnw/'+dataset+'/gnw/Size1000/mds/param_'+
        #        dataset+'-'+str(mydata.idx)+'_'+str(i+1)+'.dat', p)
        #savetxt('/home/jbao/data/DREAM/gnw/ecoli/gnw/full/mds/param_ecoli-full_'+
        #        str(i+1)+'.dat', p)
        #set_trace()
        #a[i] = p[0]
        #b[i] = p[2]
        #tau[i] = p[1]

        #connected = nx.algorithms.components.is_connected(g.g)
        #if not connected:
        #    print 'Disconnected!'

        #iwd += 1
    
        #indeg[:len(g.g_d.nodes()),i] = g.g_d.in_degree().values()
        #outdeg[:len(g.g_d.nodes()),i] = g.g_d.out_degree().values()

        #n[i] = nx.average_shortest_path_length(g.g)

        #print i,
        #sys.stdout.flush()
        #q = community.detect_communities(g.g)
        #modularity[i] = q[0]
        #print q[0]
        #set_trace()

        #params = pl.ParamLoader('%s%s'%(sub(r'\\','',mydata.wd),\
        #        sub(r'\\','',mydata.param_file)))
        #plot_param_dist(params, mds)
        #G,layers = g.get_dot_graph([mds.genes[gidx] for gidx in idx[:10]], \
        #        d[:,i], mds.genes[:size], params.n)
        #G.write('%s/mds/%s.dot' % (sub(r'\\','',mydata.wd), \
        #        files[i]))
        #G.write('%s/mds/'%sub(r'\\','',wd)+dataset+'-%d.dot'%(i+1))
   
    #np.save('/home/jbao/data/DREAM/gnw/topo/pearson_'+sys.argv[1]+'-'+
    #        sys.argv[2]+'_1000.npy', slope)
    #np.save('/home/jbao/data/DREAM/gnw/topo/pval_'+sys.argv[1]+'-'+
    #        sys.argv[2]+'_1000.npy', pval)
    #np.save('/home/jbao/data/DREAM/gnw/topo/tau_'+sys.argv[1]+'-'+
    #        sys.argv[2]+'_1000.npy', tau)
    #np.save('/home/jbao/data/DREAM/gnw/topo/tau_ecoli-full_1000.npy', tau)
    #"""
    """
        figure(1, figsize=(10,8))
        clf()
        #subplot(2,5,i+1)
        mds.plot_2d()
        #title(files[i], fontsize=50)
        title('yeast-1', fontsize=40)
        xlim(-0.6,0.6)
        ylim(-0.4,0.4)
        savefig('%s.pdf' % sub(r'\\','',wd+mds_file))
        
        #subplot(2,5,i+1)
        #x = d[:,i]
        #y = array([layers[g.g_d.nodes().index(n)] for n in \
        #        mds.genes[:n_nodes]])
        #plot(x, y, '.k')
        #ylim(-1,5)
        #pearson,p = pearsonr(x,y)
        #title('%.2f;%.2e'%(pearson,p))
        #figure(figsize=(15,12))
        #pearson,p = g.plot_ck()
        #clf()
        #g.plot_degree()
        #subplot(2,5,i+6)
        #g.draw(d[:,i], mds.genes[:n_nodes])
        #g.plot_outdeg_zero_hist()
        #g.plot_distance(d[:,i], mds.genes)
        #g.plot_outlier([mds.genes[gidx] for gidx in idx[:10]])
        #g.plot_outlier([])
        #g.plot_degree_distribution()
        #g.plot_ck()
        #g.plot_dendrogram()
        #ylim(0,1)
        #xlim(0,30)
        #title('Yeast', fontweight='bold')
        #text(0.45,3,r'$r=%.2f$'%pearson+'\n'+r'$p=%.3f$'%p,fontsize=30)
        #title('Ecoli-%d'%(i+1))
        #savefig('%s/mds/uniform-%d_degree_distance.pdf' % (sub(r'\\','',wd), \
        #        i+1))
    
    #subplots_adjust(wspace=0.5)
    dispersion = mean(d,0)/var(d,0)
    plot(modularity, dispersion, '.k')
    xlabel('Modularity')
    ylabel('Dispersion of the distance distribution')
    cc,p = pearsonr(modularity, dispersion)
    title('Pearson correlation\n%.2f; %.2e'%(cc,p))
    ax1 = gca()
    ax1.xaxis.labelpad = 10
    ax1.yaxis.labelpad = 10
    """
    #ax2 = ax1.twinx()
    #ax2.hist(rr, bins=range(-20,100,5), normed=True, cumulative=True, \
    #        histtype='step', color='r')
    """
    plot(range(1,len(all_wd)+1),mean(slope,axis=1),'ro')
    errorbar(range(1,len(all_wd)+1),mean(slope,axis=1),std(slope,axis=1),fmt=None,\
            ecolor='k')
    ylabel('Slope of the regression line')
    xticks(range(1,len(all_wd)+1),('n=10','n=100','n=1000','n=2000'))
    xlim(0,len(all_wd)+1)
    """
    #ax1 = gca()
    #ax1.title.set_y(1.05)
    #ax1.xaxis.labelpad = 10
    #ax1.yaxis.labelpad = 10
    """
    slope = zeros((5,5))
    for dummy,t in enumerate(type):
        for i in range(5):
            if t == 'dream':
                wd,n_networks,edge_file,gene_file,mds_file = load(t,i)
            else:
                wd = '/Users/bao/work/DREAM/gnw/'+t+'/'
                edge_file = t+'-'+str(i+1)+'_100_goldstandard_signed.tsv'
                gene_file = 'mds/'+t+'-'+str(i+1)+'_100_normed.dat'
                mds_file = 'mds/pca_'+t+'-'+str(i+1)+'_100_normed_euc'
            system('sed \"s/+/{\'weight\':1}/g;s/\t-/\t{\'weight\':-1}/g;'+\
                    's/\t1/\t{\'weight\':1}/g\" '+\
                    '< %s > %s.sub' % (\
                    wd+edge_file, wd+edge_file))
            g = Graph('%s.sub' % (wd+edge_file))

            mds = data.Data(sub(r'\\','',wd+mds_file), \
                sub(r'\\','',wd+gene_file))
            dist = mds.get_distance()
            #set_trace()
            slope[i,dummy],pearson,p = g.plot_degree_distance(dist, \
                    mds.genes[:100])
    bar(range(5),mean(slope,axis=0))
    """
    #show()
    #generate_graph()

