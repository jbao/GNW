#!/bin/sh
#
# $Id: mut_corr.sh 213 2010-12-15 09:43:47Z jbao $
#
#$ -cwd
#$ -S /bin/bash
# -m eas
# -M jie.bao@frias.uni-freiburg.de
# -pe mpi 8
#$ -l h_vmem=500M
#$ -l h_cpu=23:00:00
#$ -R y
#$ -j y

#net="random uniform scalefree smallworld bowtie"
net=$1
export PYTHONPATH=$PYTHONPATH:$HOME/tool/networkx-1.5/lib/python2.7/site-packages
python $HOME/github/GNW/python/gnw_parser.py $net $SGE_TASK_ID
#for n in $net; do
#/usr/local/bin/python $HOME/MicroarrayAnalysis/trunk/python/network.py $net $SGE_TASK_ID
#done
