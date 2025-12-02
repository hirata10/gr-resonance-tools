#!/bin/bash

#Submit this script with: sbatch J_dot_sf.slurm

#SBATCH --time=16:00:00   # walltime
#SBATCH --nodes=1-5   # minnodes-maxnodes #sets SLURM_NNODES
#SBATCH --ntasks=1  # number of processor cores (i.e. tasks) #sets SLURM_NTASKS
#SBATCH --ntasks-per-node=1   # number of tasks per node #sets SLURM_NTASKS_PER_NODE
#SBATCH --error=error_J_dot_sf_50_200-0.1_28_14.txt    #print out error message
#SBATCH --output=output_J_dot_sf_50_200_0.1_28_14.txt   #print output
#SBATCH --job-name=J_dot_sf_50_200_0.1_28_14   # job name
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

echo "Running perl pars.pl ..."
perl pars.pl
if [ $? -ne 0 ]; then
    echo "perl pars.pl failed!"
    exit 1
fi
echo "perl pars.pl finished successfully."

echo "Compiling J_dot_res_case..."
gcc calling.c -DIS_J_DOT resonance_find.c J_dot.c kerrtraj.c kerrmode.c kerrgwem.c -w -o J_dot_sf -lm
if [ $? -ne 0 ]; then
    echo "Compilation for J_dot_res_case failed!"
    exit 1
fi
echo "Compilation for J_dot_res_case successful."

## Studying the case of a strong resonance ##

cmd="./J_dot_sf 2.537191410184998 2.282777584641027 2.2510988170394506"

echo $cmd
$cmd