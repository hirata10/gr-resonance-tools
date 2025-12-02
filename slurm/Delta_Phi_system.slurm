#!/bin/bash

#Submit this script with: sbatch filename

#SBATCH --time=16:00:00   # walltime
#SBATCH --nodes=16   # minnodes-maxnodes #sets SLURM_NNODES
#SBATCH --ntasks-per-node=128   # number of tasks per node #sets SLURM_NTASKS_PER_NODE ## CHECK THE PYTHON CODE REQUESTS THE SAME NUMBER OF TOTAL JOBS
#SBATCH --error=error_Delta_Phi_50_200-0.1.txt    #print out error message
#SBATCH --output=output_Delta_Phi_50_200-0.1.txt   #print output
#SBATCH --job-name=Delta_Phi_50_200-0.1  # job name
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
ranks_per_node=128
echo "number of CPUs per node = $ranks_per_node"
total_ranks=$(($SLURM_NNODES*$ranks_per_node))
echo "number of ranks total = $total_ranks"

#### cmd="srun -N $SLURM_NNODES --tasks-per-node $SLURM_NTASKS_PER_NODE python3 action_angle_generator.py"
srun hostname

cmd="srun -N $SLURM_NNODES --tasks-per-node $SLURM_NTASKS_PER_NODE -n 1 -c $total_ranks python3 Delta_Phi_single_parallel.py"

echo $cmd
$cmd

