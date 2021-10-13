#!/bin/bash

echo "RUNNING"
echo "slurm procid = " $SLURM_PROCID
echo "slurm ntasks = " $SLURM_NTASKS
echo "slurm cpus per task = " $SLURM_CPUS_PER_TASK

#Compute sims
cmd="python ./LS_run_and_metacalibrate_one_sim.py"
date
echo $cmd
$cmd
date
