#!/bin/bash
#SBATCH -J network_optimization     # short name for my job
#SBATCH -t 00:59:00                 # wall-clock time. Use <1hr for testing and faster scheduling
#SBATCH -p genoa                     # choose the partition/node to run the job. Thin=default
#SBATCH -N 1                       # node count (min 1/8 of "thin" node = 16 cores)
#SBATCH --ntasks 1                  # total number of tasks across all nodes
#SBATCH --cpus-per-task 32           # cpu-cores per task (>1 if multi-threaded tasks)

export GRB_LICENSE_FILE=/gpfs/home5/dgkontoras/.gurobi/gurobi.lic

module purge
module load 2023
module load Julia/1.9.2-linux-x86_64

cp -r $HOME/GridOptiPlan "$TMPDIR"      # copy data to scratch for a single node job

cd "$TMPDIR"/GridOptiPlan

julia --project main.jl

cp -r "$TMPDIR"/GridOptiPlan/output $HOME/GridOptiPlan  # copy output data from scratch back to your home folder