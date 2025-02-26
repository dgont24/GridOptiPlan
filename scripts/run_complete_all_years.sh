#!/bin/bash
#SBATCH -J network_optimization     # short name for my job
#SBATCH -t 05:00:00                 # wall-clock time. Use <1hr for testing and faster scheduling
#SBATCH -p genoa                    # choose the partition/node to run the job. Thin=default
#SBATCH -N 1                        # node count (min 1/8 of "thin" node = 16 cores)
#SBATCH --ntasks 1                  # total number of tasks across all nodes
#SBATCH --cpus-per-task 96         # cpu-cores per task (>1 if multi-threaded tasks)

export GRB_LICENSE_FILE=/gpfs/home5/dgkontoras/.gurobi/gurobi.lic
export JULIA_NUM_THREADS=96
export GA_ITERATIONS=50 # Number of iterations for the genetic algorithm.
export GA_POPULATION=95 # Number of individuals in each generation(iteration). Keep it <= JULIA_NUM_THREADS
export TIME_LIMIT_HR=5.4  # Time limit for the optimization in hours.

cd $HOME/GridOptiPlan

$HOME/julia-1.10.4/bin/julia --project -e 'import Pkg; Pkg.instantiate();
include("scripts/run_complete_all_years.jl")'
