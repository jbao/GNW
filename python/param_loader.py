"""
read GeneNetWeaver simulation parameters from the xml file.

$Id: param_loader.py 287 2012-01-31 12:53:33Z jbao $
"""

from pdb import set_trace
from linecache import getline
from re import sub,finditer
from libsbml import *
from pylab import *
import sys, os
from collections import defaultdict
import networkx as nx
import csv

rcParams['xtick.labelsize'] = 20
rcParams['axes.titlesize'] = 32

class ParamLoader:
    def __init__(self):
        pass

    def load_param(self, filename):
        reader = SBMLReader()
        document = reader.readSBML(filename)
        assert document.getNumErrors() == 0, '%s not found'%filename

        model = document.getModel()
        reactions = model.getListOfReactions()

        self.max = {}
        self.deltaProtein = {}
        self.maxTranslation = {}
        self.delta = {}
        self.k = {}
        self.n = {}

        #set_trace()
        for r in reactions:
            p_list = r.getKineticLaw().getListOfParameters()
            starts = [match.start() for match in finditer('_',r.getId())]
            node = r.getId()[:starts[-1]]
            #print node
            #sys.stdout.flush()
            inputs = r.getListOfModifiers()
            input_name = []
            for i in inputs:
                input_name.append(i.getSpecies())
    
            # params
            for p in p_list:
                if p.getId() == 'max':
                    self.max[node] = p.getValue()
                elif p.getId() == 'deltaProtein':
                    self.deltaProtein[node] = p.getValue()
                elif p.getId() == 'maxTranslation':
                    self.maxTranslation[node] = p.getValue()
                elif p.getId() == 'delta':
                    self.delta[node] = p.getValue()
                elif 'k_' in p.getId():
                    idx = int(p.getId().split('_')[-1]) - 1
                    if not self.k.has_key(input_name[idx]):
                        # first encounter of the dict [pre][post][k]
                        self.k[input_name[idx]] = {node:p.getValue()}
                    else:
                        self.k[input_name[idx]][node] = p.getValue()
                elif 'n_' in p.getId():
                    idx = int(p.getId().split('_')[-1]) - 1
                    if not self.n.has_key(input_name[idx]):
                        # first encounter of the dict [pre][post][n]
                        self.n[input_name[idx]] = {node:p.getValue()}
                    else:
                        self.n[input_name[idx]][node] = p.getValue()

        # sanity check
        #ls = model.getListOfSpecies()
        #for i in range(ls.size()):
        #    print ls.get(i).getId(), self.max[ls.get(i).getId()]
        #return model

    def plot_polar(self, param):
        values = []
        for k in sort(param.keys()):
            values.append(param[k])
        fig = figure()
        #ax = fig.add_axes([0.1, 0.1, 0.8, 0.8], polar=True)
        ax = fig.add_subplot(111, polar=True)
        ax.plot(linspace(0,2*pi,len(param)), values, color='r')
        xticks(linspace(0,5.0/3*pi,6), range(1,int(5.0/6*len(param)),int(1.0/6*len(param))))
        yticks([], color='0.5')
        show()

    def load_input(self, filename):
        n_perturb = 1
        gene = getline(filename, 1).rstrip().split('\t')[1:]
        self.gene = [sub('"','',g) for g in gene]
        n_gene = len(gene)
        influx = loadtxt(filename, delimiter='\t', comments='%', skiprows=1, \
                usecols=range(1,n_gene+1))
        self.time_point = loadtxt(filename, delimiter='\t', comments='%', 
                skiprows=1, usecols=[0])[:influx.shape[0]/n_perturb]
        n_time_point = len(self.time_point)
        influx = reshape(influx, (n_perturb,n_time_point,n_gene))
        return influx

    def load_jacobian(self, filename):
        gene = getline(filename, 1).rstrip().split('\t')[1:]
        self.jacobian_gene = [sub('"','',g) for g in gene]
        self.jacobian = loadtxt(filename, delimiter='\t', comments='%', 
                skiprows=1)
        
    def load_timeseries(self, filename, normalized=True):
        n_perturb = 1
        gene = getline(filename, 1).rstrip().split('\t')[1:]
        self.gene = [sub('"','',g) for g in gene]
        n_gene = len(gene)
        #set_trace()
        self.timeseries = loadtxt(filename, delimiter='\t', comments='%', skiprows=1, \
                usecols=range(1,n_gene+1))
        self.time_point = loadtxt(filename, delimiter='\t', comments='%', skiprows=1, \
                usecols=[0])[:self.timeseries.shape[0]/n_perturb]
        n_time_point = len(self.time_point)
        self.timeseries = reshape(self.timeseries, (n_perturb,n_time_point,n_gene))
        # normalization
        if normalized:
            for i in range(n_perturb):
                self.timeseries[i,:,:] = log2(self.timeseries[i,:,:]/tile(self.timeseries[i,-1,:],(n_time_point,
                            1)))

    def load_wildtype(self, filename):
        n_perturb = 1
        gene = getline(filename, 1).rstrip().split('\t')[1:]
        self.gene = [sub('"','',g) for g in gene]
        n_gene = len(gene)
        #set_trace()
        self.wildtype = loadtxt(filename, delimiter='\t', comments='%', skiprows=1, \
                usecols=range(1,n_gene+1))
        self.wildtype_time_point = loadtxt(filename, delimiter='\t', comments='%', skiprows=1, \
                usecols=[0])[:self.wildtype.shape[0]/n_perturb]
        n_time_point = len(self.wildtype_time_point)
        self.wildtype = reshape(self.wildtype, (n_perturb,n_time_point,n_gene))

    def load_perturbation(self, filename):
        n_perturb = 1
        gene = getline(filename, 1).rstrip().split('\t')[1:]
        self.gene = [sub('"','',g) for g in gene]
        n_gene = len(gene)
        #set_trace()
        self.perturbation = loadtxt(filename, delimiter='\t', comments='%', skiprows=1, \
                usecols=range(1,n_gene+1))
        self.perturbation = reshape(self.perturbation, (n_perturb,n_gene))

    def load_mds(self, raw_file, pval_file):
        gene = loadtxt(raw_file, delimiter='|', usecols=[0], skiprows=1,
                dtype=str)
        pval = loadtxt(pval_file)
        self.response = {}
        for g,p in zip(gene,pval):
            self.response[g] = -log10(p) / max(-log10(pval))

    def load_sign(self, edge_file):
        edge = loadtxt(edge_file, delimiter='\t', dtype=str)
        self.sign = {}
        for i in range(edge.shape[0]):
            if edge[i,2] == '+':
                if not self.sign.has_key(edge[i,0]):
                    # first encounter of the dict [pre][post][n]
                    self.sign[edge[i,0]] = {edge[i,1]:1}
                else:
                    self.sign[edge[i,0]][edge[i,1]] = 1
            elif edge[i,2] == '-':
                if not self.sign.has_key(edge[i,0]):
                    # first encounter of the dict [pre][post][n]
                    self.sign[edge[i,0]] = {edge[i,1]:-1}
                else:
                    self.sign[edge[i,0]][edge[i,1]] = -1

    def get_concentration_index(self, data):
        flipped = defaultdict(dict)
        for key, val in data.items():
            for subkey, subval in val.items():
                flipped[subkey][key] = subval

        concentration_index = {}
        for n in flipped.keys():
            numerator = 0
            denominator = 0
            for v in flipped[n].values():
                numerator += v
                denominator += v**2
            #if denominator != 0:
            concentration_index[n] = numerator**2 / denominator

        return concentration_index

    def get_control_index(self, data):
        flipped = defaultdict(dict)
        for key, val in data.items():
            for subkey, subval in val.items():
                flipped[subkey][key] = subval

        control_index = {}
        for n in data.keys():   # for source node i
            h = 0
            for k in data[n].keys():  # for target node j of i
                denominator = 0
                for v in flipped[k].values():   # for all source node of j
                    denominator += v**2
                h += data[n][k]**2 / denominator

            control_index[n] = h

        return control_index

def flatten(d):
    ret = []

    for v in d.values():
        if isinstance(v, dict):
            ret.extend(flatten(v))
        else:
            ret.append(v)

    return ret

if __name__ == '__main__':
    #close('all')
    type = 'uniform2'
    wd = "/home/jbao/data/DREAM/gnw/"+type+"/gnw/Size1000/norm_perturbation/"
    #wd = '/home/jbao/data/DREAM/DREAM3/prune/Size100/'
    for i in range(100):
        print i,
        sys.stdout.flush()

        original = ParamLoader()
        original.load_param(wd+type+'-'+str(i+1)+'.xml')
        w = csv.writer(open(wd+type+"-"+str(i+1)+"_delta.csv", "w"))
        for kk,vv in original.delta.items():
            w.writerow([kk, vv])
        w = csv.writer(open(wd+type+"-"+str(i+1)+"_k.csv", "w"))
        for key in original.k.keys():
            for kk,vv in original.k[key].items():
                w.writerow([key, kk, vv])
        w = csv.writer(open(wd+type+"-"+str(i+1)+"_n.csv", "w"))
        for key in original.n.keys():
            for kk,vv in original.n[key].items():
                w.writerow([key, kk, vv])
    #g = nx.read_edgelist("/home/jbao/data/DREAM/gnw/scalefree2/gnw/Size1000/input/scalefree2-1_1000.tsv", create_using=nx.DiGraph())
    #for n in g.nodes():
    #    plot(original.delta[n], g.degree()[n], '.k')
    #original.load_sign(wd+'ecoli-full_goldstandard_signed.tsv')
    #concentration_index = original.get_concentration_index(original.k)
    #control_index = original.get_control_index(original.k)
    #original.load_mds(wd+'mds/ecoli-full_perturbation-1_1502_normed.dat',
    #        wd+'mds/pval_ecoli-full_perturbation-1_1502.dat')
    #for c in concentration_index:
    #    plot(original.response[c], concentration_index[c], '.k')
    #original.load_timeseries(wd+'full_multifactorial_timeseries.tsv',False)
    #original.load_wildtype(wd+'full_wildtype.tsv')
    #original.load_perturbation(wd+'full__perturbations.tsv')
    #figure()
    #plot(original.wildtype_time_point, \
    #        mean(original.wildtype,0)[:,70], 'k', linewidth=4)
    #plot(original.time_point+max(original.wildtype_time_point), \
    #        mean(original.timeseries,0)[:,70], 'k', linewidth=4)
    #plot([max(original.wildtype_time_point),max(original.wildtype_time_point)],\
    #        [0,1], '--r', linewidth=4)
    #plot([max(original.wildtype_time_point),max(original.wildtype_time_point)+max(original.time_point)],\
    #        [original.perturbation[0,70],original.perturbation[0,70]], '-.b', linewidth=4)
    #xlabel('Time (a.u.)', fontweight='bold')
    #ylabel('Normalized expression', fontweight='bold')
    #title('G68', fontweight='bold')
    #ax = gca()
    #ax.title.set_y(1.05)
    
    #imshow(mean(pl.timeseries,0).T, origin='upper')
    #colorbar()
    #yticks(range(len(pl.gene)), pl.gene)
    #axis('auto')
    show()
    #pl.plot_polar(pl.max)

