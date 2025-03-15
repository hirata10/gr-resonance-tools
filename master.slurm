#!/bin/bash

#Submit this script with: sbatch filename

#SBATCH --time=00:30:00   # walltime
#SBATCH --nodes=1   # minnodes-maxnodes #sets SLURM_NNODES
#SBATCH --ntasks=1  # number of processor cores (i.e. tasks) #sets SLURM_NTASKS
#SBATCH --ntasks-per-node=1   # number of tasks per node #sets SLURM_NTASKS_PER_NODE
#SBATCH --error=error_master_50_200-0.1.txt    #print out error message
#SBATCH --output=output_master_50_200-0.1.txt   #print output
#SBATCH --job-name=master_script   # job name
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
ranks_per_node=32
echo "number of CPUs per node = $ranks_per_node"
total_ranks=$(($SLURM_NNODES*$ranks_per_node))
echo "number of ranks total = $total_ranks"

echo "Master job starting at $(date)"

# Run perl script at the very beginning
# MAKE SURE YOU SET THE NECESSARY PARAMETERS IN THE FILE pars.txt

CURRENT_DIR=$SLURM_SUBMIT_DIR

echo "Current directory is $CURRENT_DIR"

echo "Running perl pars.pl ..."
perl pars.pl
if [ $? -ne 0 ]; then
    echo "perl pars.pl failed!"
    exit 1
fi
echo "perl pars.pl finished successfully."

#Make directory for J_evolve outputs

# mkdir -p outputs_data

# Compilation step for necessary executables
echo "Compiling J_evolve_single..."
gcc -DIS_RK4_J_DOT calling.c resonance_find.c kerrtraj.c kerrmode.c kerrgwem.c J_dot.c J2J_dot.c -w -o J_evolve_single -lm
if [ $? -ne 0 ]; then
    echo "Compilation for J_evolve_single failed!"
    exit 1
fi
echo "Compilation for J_evolve_single successful."

echo "Compiling Delta_J_single..."
gcc calling.c -DIS_DELTA_J_SINGLE kerrtraj.c kerrmode.c kerrgwem.c resonance_find.c Gamma.c Delta_J.c J_dot.c -o Delta_J_single -w -lm
if [ $? -ne 0 ]; then
    echo "Compilation for Delta_J_single failed!"
    exit 1
fi
echo "Compilation for Delta_J_single successful."

echo "Compiling Delta_Phi_single..."
gcc kerrphase.c kerrtraj.c kerrgwem.c kerrmode.c resonance_find.c J_dot.c -w -o Delta_Phi_single -lm
if [ $? -ne 0 ]; then
    echo "Compilation for Delta_Phi_single failed!"
    exit 1
fi
echo "Compilation for Delta_Phi_single successful."

# Submit action_generator.slurm and capture the job ID
jobid1=$(sbatch action_generator.slurm | awk '{print $4}')

# Submit J_evolve.slurm with a dependency on action_generator.slurm 
jobid2=$(sbatch --dependency=afterok:$jobid1 J_evolve.slurm | awk '{print $4}')

# Submit res_finder_parallel.slurm with a dependency on J_evolve.slurm
jobid3=$(sbatch --dependency=afterany:$jobid2 --array=1-16 res_finder_parallel.slurm | awk '{print $4}')

# jobid3=$(sbatch --array=1-16 res_finder_parallel.slurm | awk '{print $4}')

# Submit Delta_J.slurm with a dependency on sub_job2 (only run after sub_job2 finishes successfully)
jobid4=$(sbatch --dependency=afterany:$jobid3 --array=1-16 Delta_J.slurm | awk '{print $4}')

# # Submit cmd1 as a new job with dependency on Delta_J
# jobid5=$(sbatch --dependency=afterok:$jobid4 --wrap="python3 organize_Delta_J.py" | awk '{print $4}')

jobid5=$(sbatch --dependency=afterany:$jobid4 --wrap="python3 organize_Delta_J.py" \
    --error=error_organize_Delta_J.txt --output=output_organize_Delta_J.txt | awk '{print $4}')

# # Submit cmd2 as a new job with dependency on cmd1
# jobid6=$(sbatch --dependency=afterok:$jobid5 --wrap="python3 organize_complete_Delta_J.py $CURRENT_DIR" | awk '{print $4}')

jobid6=$(sbatch --dependency=afterany:$jobid5 --wrap="bash -c 'python3 organize_complete_Delta_J.py \"$CURRENT_DIR\"'" \
    --error=error_organize_complete_Delta_J.txt --output=output_organize_complete_Delta_J.txt | awk '{print $4}')

# Submit Delta_Phi.slurm with dependency on cmd2
jobid7=$(sbatch --dependency=afterok:$jobid6 Delta_Phi_system.slurm | awk '{print $4}')

# # Submit cmd2 as a new job with dependency on cmd1
# jobid8=$(sbatch --dependency=afterok:$jobid7 --wrap="python3 organize_Phi.py" | awk '{print $4}')

jobid8=$(sbatch --dependency=afterany:$jobid7 --wrap="python3 organize_Phi.py" \
    --error=error_organize_Phi.txt --output=output_organize_Phi.txt | awk '{print $4}')


# Log the job IDs
echo "Submitted action_generator with job ID $jobid1"
echo "Submitted J_evolve with job ID $jobid2"
echo "Submitted res_finder_parallel with job ID $jobid3"
echo "Submitted Delta_J with job ID $jobid4"
echo "Submitted organize_Delta_J with job ID $jobid5"
echo "Submitted organize_complete_Delta_J with job ID $jobid6"
echo "Submitted Delta_Phi with job ID $jobid7" 
echo "Submitted organize_Phi.py with job ID $jobid8"

# Log when the master script finishes
echo "Master job finished at $(date)"