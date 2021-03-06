GNW
===

C++ version of GeneNetWeaver

See also http://sourceforge.net/projects/gnw/

To simulate time series with constant input, add

    GnwSettings::Instance()->generateTsConstantInput(true);

in `main.cpp`, otherwise keep this flag `false`.

Workflow
--------

-   Generate network topology by the `generate_graph` function in 
    `python/network.py`, implemented with the `NetworkX` library. Within the
    script, one can specify the size of the network, the number of networks
    and the directory to save the edge list `.tsv` file. For example, start
    a python console by `ipython` within the `python/` directory, then
        
        from network import generate_graph
        generate_graph('scalefree')

-   Use GNW to generate the time series. To submit the job to the SGE queueing 
    system, use the `gnw.sh` script, or for example

        for i in `seq 1 50`;do qsub -t 1-100 gnw.sh scalefree2 $i;done

    Here the `-t` option of `qsub` specifies the range of the array job, i.e. 
    the same
    job is submitted in parallel with `SGE_TASK_ID` ranging from 1 to 100, which is
    the index of different network topologies. 
    `scalefree2` is the name of the network, and the second argument `$i` runs from
    1 to 50, which is the index of different perturbations in this case.

    For the simulation of only one network, the equivalent is

        qsub -t 1 gnw.sh scalefree2 1

    Use `watch qstat` to monitor the progress of all jobs.

-   Convert the GNW output to proper format. `python/gnw_parser.py` converts the 
    time series (`*multifactorial_timeseries.tsv`) to format that is readable by 
    the
    `hitmds` tool and save it in the `mds/` subdirectory, to run in parallel

        qsub -t 1-100 python/fnet.sh scalefree2 50

    The array job range 1-100 defines the number of unique network topologies, and 
    the last argument 50 is the number of different perturbations.

    In order to convert the sbml parameter file to a plain text file readable by
    R, use the `python/param_loader.py` script. One can open up the script in any
    text editor and go to the line with `if __name__ == '__main__':`, which is
    the main function. In general, the python script loops through all sbml files
    generated by different topologies and spits out 3 files (`*_delta.csv, 
    *_k.csv, *_n.csv`) saving the decay rate (delta), dissociation concentration
    (k) and Hill coefficient (n) respectively.

-   Run MDS. Code is located in the `github/mds` repository, see the documentation
    there or in a nutshell, one can submit the job to the queueing system by
        
        for i in `seq 1 100`;do qsub -t 1-50 mds.sh scalefree2 $i;done

    Notice that the syntax is similar to the `gnw.sh` script, but the network
    and perturbation index are interchanged, such that the array job index 1-50
    denote different perturbations and the loop index 1-100 denote different
    topologies.

-   Convert MDS coordinates to p-values with `github/mds/R/response.R`. Currently
    this R script cannot be submitted to the queueing system and therefore only 
    runs on the head node. If one has multiple networks, one can use the bash 
    script `github/mds/R/response_loop.sh`, which is basically a sequential loop.
