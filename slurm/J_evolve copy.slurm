#!/bin/bash
# SLURM Job Configuration
#SBATCH --time=16:00:00                 # Walltime (max runtime)
#SBATCH --nodes=1                     # Number of nodes (min-max nodes)
#SBATCH --ntasks-per-node=30            # Number of tasks per node
#SBATCH --job-name=J_evolve_50_200_0.1  # Job name
#SBATCH --output=output_J_evolve_50_200_0.1.txt   # SLURM output file
#SBATCH --error=error_J_evolve_50_200_0.1.txt    # SLURM error file
#SBATCH --mail-user=makanas@lanl.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# Log the start of the job
echo "Job started at $(date)" >> output_J_evolve_50_200_0.1.txt
echo "Job started at $(date)" >> error_J_evolve_50_200_0.1.txt

# Making the parameter files
echo "Running perl pars.pl ..."
perl pars.pl
if [ $? -ne 0 ]; then
    echo "perl pars.pl failed!"
    exit 1
fi
echo "perl pars.pl finished successfully."

# Compiles the executable for the beginning evolution
echo "Compiling J_evolve_single_norestart..."
gcc -DIS_RK4_J_DOT_norestart calling.c resonance_find.c kerrtraj.c kerrmode.c kerrgwem.c J_dot.c J2J_dot.c -w -o J_evolve_single_norestart -lm
if [ $? -ne 0 ]; then
    echo "Compilation for J_evolve_single_norestart failed!"
    exit 1
fi
echo "Compilation for J_evolve_single_norestart successful."

# Compiles the executable that starts at a restart point
echo "Compiling J_evolve_single_restart..."
gcc -DIS_RK4_J_DOT_restart calling.c resonance_find.c kerrtraj.c kerrmode.c kerrgwem.c J_dot.c J2J_dot.c -w -o J_evolve_single_restart -lm
if [ $? -ne 0 ]; then
    echo "Compilation for J_evolve_single_restart failed!"
    exit 1
fi
echo "Compilation for J_evolve_single_restart successful."


# Run the Python script
python3 J_evolve_single.py >> output_J_evolve_50_200_0.1.txt 2>> error_J_evolve_50_200_0.1.txt

# Capture the exit status of the Python script
EXIT_STATUS=$?

# Check the exit status and act accordingly
if [ $EXIT_STATUS -eq 0 ]; then
    # Exit code 0: Job completed successfully
    echo "Job completed successfully at $(date)" >> output_J_evolve_50_200_0.1.txt

elif [ $EXIT_STATUS -eq 137 ]; then
    # Exit code 137: Walltime exceeded
    echo "Walltime exceeded at $(date). Resubmitting the job..." >> error_J_evolve_50_200_0.1.txt
    sbatch J_evolve.slurm   # Resubmit the job if walltime was exceeded

else
    # Handle any other exit codes if necessary
    echo "Job failed with exit code $EXIT_STATUS at $(date). Please check the logs." >> error_J_evolve_50_200_0.1.txt
fi

# Optionally log when the job ends or is resubmitted
echo "Job ended at $(date)" >> output_J_evolve_50_200_0.1.txt
echo "Job ended at $(date)" >> error_J_evolve_50_200_0.1.txt

# #!/bin/bash

# #Submit this script with: sbatch filename

# #SBATCH --time=16:00:00   # walltime
# #SBATCH --nodes=1-5   # minnodes-maxnodes #sets SLURM_NNODES
# #SBATCH --ntasks=50  # number of processor cores (i.e. tasks) #sets SLURM_NTASKS
# #SBATCH --ntasks-per-node=50   # number of tasks per node #sets SLURM_NTASKS_PER_NODE
# #SBATCH --error=error_J_evolve_50_200_0.1.txt    #print out error message
# #SBATCH --output=output_J_evolve_50_200_0.1.txt   #print output
# #SBATCH --job-name=J_evolve_50_200_0.1   # job name
# #SBATCH --mail-user=makanas@lanl.gov   # email address
# #SBATCH --mail-type=BEGIN
# #SBATCH --mail-type=END
# #SBATCH --mail-type=FAIL
# #SBATCH --no-requeue   # does not requeue when preempted

# # LOAD MODULEFILES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

# echo "number of nodes = $SLURM_NNODES"
# echo "SLURM says:"
# echo "number of CPUs per node = $SLURM_CPUS_ON_NODE"
# echo "number of tasks per node = $SLURM_NTASKS_PER_NODE"
# echo "I say:"
# ranks_per_node=100
# echo "number of CPUs per node = $ranks_per_node"
# total_ranks=$(($SLURM_NNODES*$ranks_per_node))
# echo "number of ranks total = $total_ranks"

# #### cmd="srun -N $SLURM_NNODES --tasks-per-node $SLURM_NTASKS_PER_NODE python3 action_angle_generator.py"

# cmd="python3 J_evolve_single.py"

# echo $cmd
# $cmd

