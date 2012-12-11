GNW
===

C++ version of GeneNetWeaver

See also http://sourceforge.net/projects/gnw/

Workflow
--------

1. Use GNW to generate the time series. To submit the job to the SGE queueing 
system, use the `gnw.sh` script, or for example

    for i in `seq 1 50`;do qsub -t 1-100 gnw.sh scalefree2 $i;done

Here the `-t` option of `qsub` specifies the range of the array job, i.e. the same
job is submitted in parallel with `SGE_TASK_ID` ranging from 1 to 100, which is
the index of different network topologies. 
`scalefree2` is the name of the network, and the second argument `$i` runs from
1 to 50, which is the index of different perturbations in this case.

For the simulation of only one network, the equivalent is

    qsub -t 1 gnw.sh scalefree2 1

Use `watch qstat` to monitor the progress of all jobs.

2. Convert the GNW output to proper format. `python/gnw_parser.py` converts the 
time series (`*multifactorial_timeseries.tsv`) to format that is readable by the
`hitmds` tool and save it in the `mds/` subdirectory, to run in parallel

    qsub -t 1-100 python/fnet.sh scalefree2 50

The array job range 1-100 defines the number of unique network topologies, and 
the last argument 50 is the number of 
