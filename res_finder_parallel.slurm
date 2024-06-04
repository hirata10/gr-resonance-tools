#!/bin/bash

#Submit this script with: sbatch filename

#SBATCH --time=16:00:00   # walltime
#SBATCH --nodes=1-5   # minnodes-maxnodes #sets SLURM_NNODES
#SBATCH --ntasks=1  # number of processor cores (i.e. tasks) #sets SLURM_NTASKS
#SBATCH --ntasks-per-node=1   # number of tasks per node #sets SLURM_NTASKS_PER_NODE
#SBATCH --error=error_res_finder_parallel.txt    #print out error message
#SBATCH --output=output_res_finder_parallel.txt   #print output
#SBATCH --job-name=res_finder_parallel   # job name
#SBATCH --mail-user=makanas@lanl.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue   # does not requeue when preempted

# LOAD MODULEFILES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

echo "number of nodes = $SLURM_NNODES"
echo "SLURM says:"
echo "number of CPUs per node = $SLURM_CPUS_ON_NODE"
echo "number of tasks per node = $SLURM_NTASKS_PER_NODE"
echo "I say:"
ranks_per_node=1
echo "number of CPUs per node = $ranks_per_node"
total_ranks=$(($SLURM_NNODES*$ranks_per_node))
echo "number of ranks total = $total_ranks"

#### cmd="srun -N $SLURM_NNODES --tasks-per-node $SLURM_NTASKS_PER_NODE python3 action_angle_generator.py"

cmd="python3 resfinder_parallel.py 1"

echo $cmd
$cmd

