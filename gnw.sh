#!/bin/sh
#
# $Id: gnw.sh 29 2012-01-04 17:06:55Z jbao $
#
#$ -cwd
#$ -S /bin/bash
# -m eas
# -M jie.bao@frias.uni-freiburg.de
# -pe mpi 8
#$ -l h_vmem=2g
#$ -l h_cpu=23:00:00
#$ -R y

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/tool/libsbml-4.2.0/lib/:$HOME/tool/getopt_pp/
net=ecoli
#type=rewiring
size=1502
outdir=/export/work/jbao/data/DREAM/gnw/$net/gnw/full/
indir=/home/jbao/data/DREAM/gnw/$net/gnw/full/
if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi
#for idx in `seq 1 1`; do
#idx=1
$HOME/GNW/branches/c++/gnw -o $outdir -i $indir -t $net -n $idx -s $size
    #mv $outdir/$net-${idx}_perturbation-${SGE_TASK_ID}_* \
    #    $outdir/$net-${idx}_goldstandard_signed.tsv \
    #    $outdir/$net-${idx}.xml $indir
    #mv    $outdir/$net-${SGE_TASK_ID}_jacobian.tsv $indir
#done
mv $outdir/$net-full_* \
    $outdir/$net-full.xml $indir
