#!/bin/bash
module load python/anaconda-2020.02
source activate /project2/nobrega/grace/conda/TAD
export LD_LIBRARY_PATH=$HOME/bin/hdf5-1.8.20-linux-centos7-x86_64-gcc485-shared/lib/:$LD_LIBRARY_PATH

######### Author: Grace Hansen #########
#This script runs the Snakemake TAD fusion pipeline on the cluster (for me, midway on RCC at UChicago, a slurm-based cluster)

snakemake \
    --snakefile /project2/nobrega/grace/TAD_fusion/scripts/snake/Snakefile \
    -kp \
    -j 500 \
    --rerun-incomplete \
    --cluster-config /project2/nobrega/grace/TAD_fusion/scripts/snake/cluster.json \
    -c "sbatch \
        --partition={cluster.partition} \
        --mem={cluster.mem} \
        --nodes={cluster.n} \
        --tasks-per-node={cluster.tasks} \
        --time={cluster.time} \
        --job-name={cluster.name} \
	    --output={cluster.logfile}" \
    $*
