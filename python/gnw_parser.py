"""
parse the simulation results from GNW and write to a .dat file for MDS

$Id: gnw_parser.py 287 2012-01-31 12:53:33Z jbao $
"""
from linecache import getline
from re import sub,finditer
import sys,os
from pylab import *
from pdb import set_trace

net = sys.argv[1]
#size = 100
#wd = '/home/jbao/DREAM/gnw/'+net+'/gnw/Size'+str(size)+'/default/rewiring/'
wd = '/home/jbao/data/DREAM/gnw/'+net+'/gnw/Size1000/norm_perturbation/'
#wd = '/Users/bao/work/DREAM/DREAM3_in_silico_challenge/Size'+str(size)+'/gnw/'
n_perturb = 1
delim = '|'
n_networks = sys.argv[3]

#i_network = int(os.getenv('SGE_TASK_ID'))
fileList = os.listdir(wd)
#interval = 1000
for i_network in range(n_networks):  #range((int(sys.argv[1])-1)*interval,\
        #int(sys.argv[1])*interval):
    
    print i_network,
    sys.stdout.flush()

    # Parsing, convert .tsv to .dat
    #fname = 'DREAM3 data/InSilicoSize'+str(size)+'-Ecoli1-trajectories.tsv'
    fname = net+'-'+sys.argv[2]+'_perturbation-'+str(i_network+1)+'_multifactorial_timeseries.tsv'
    #fname = 'full_multifactorial_timeseries.tsv'
    #fname = [f for f in fileList if f.find(net+'-'+str(i_network+1)+'_') > -1\
    #        and f.find('multifactorial_timeseries.tsv') > -1]
    #fname = fname[0]
    
    gene = getline(wd+fname, 1).rstrip().split('\t')[1:]
    gene = [sub('"','',g) for g in gene]
    n_gene = len(gene)
    #set_trace()
    data = loadtxt(wd+fname, delimiter='\t', comments='%', skiprows=1, \
            usecols=range(1,n_gene+1))
    time_point = loadtxt(wd+fname, delimiter='\t', comments='%', skiprows=1, \
            usecols=[0])[:data.shape[0]/n_perturb]
    n_time_point = len(time_point)

    # normalization
    data = reshape(data, (n_perturb,n_time_point,n_gene))
    for i in range(n_perturb):
        data[i,:,:] = log2(2**data[i,:,:]/2**tile(data[i,-1,:],(n_time_point,1)))

    if not os.path.isdir(wd+'/mds'):
        os.makedirs(wd+'/mds')
    #f = open(wd+'mds/'+net+'-'+str(i_network+1)+'_'+str(size)+\
    #        '_normed.dat','w')
    starts = [match.start() for match in finditer('_',fname)]
    f = open(wd+'mds/'+fname[:starts[-2]]+'_'+str(n_gene)+'_normed.dat','w')
    #f = open(wd+'mds/ecoli1_'+str(size)+'_normed.dat','w')
    # header
    f.write('Gene')
    for t in time_point:
        f.write(delim+str(t))
    f.write('\n')
    # body
    for i in range(n_perturb):
        for j,g in enumerate(gene):
            # gene name
            f.write(g)
            for k,t in enumerate(time_point):
                f.write(delim+str(data[i,k,j]))
            f.write('\n')
    f.close()

    # MDS
    #os.system('./dtw %s/mds/%s-%d_%d_normed.dat 10' % (wd, net, \
    #        i_network+1, size))
    #os.system('./hitmds2 100 2 0.04 0 dist_matrix.bin %d > %s/mds/mds_'%\
    #        (size,wd)+'%s-%d_%d_normed_euc'%(net,i_network+1,size))
    #os.system('./pca %s/mds/mds_%s-%d_%d_normed_euc %d 2' % (wd, net, \
    #        i_network+1, size, size))
    #os.system('mv recons.dat %s/mds/pca_%s-%d_%d_normed_euc' % (wd, net, \
    #        i_network+1, size))


