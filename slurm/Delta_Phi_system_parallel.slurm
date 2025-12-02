#!/bin/bash

#Submit this script with: sbatch filename

#SBATCH --time=16:00:00   # walltime
#SBATCH --nodes=16   # minnodes-maxnodes #sets SLURM_NNODES
#SBATCH --ntasks-per-node=1   # number of tasks per node #sets SLURM_NTASKS_PER_NODE ## CHECK THE PYTHON CODE REQUESTS THE SAME NUMBER OF TOTAL JOBS
#SBATCH --array=0-6699%2048   # Job array: tasks 0-6699, limit 2048 tasks running at once
#SBATCH --error=error_Delta_Phi_50_200-0.1.txt    #print out error message
#SBATCH --output=output_Delta_Phi_50_200-0.1.txt   #print output
#SBATCH --job-name=Delta_Phi_50_200-0.1  # job name
#SBATCH --mail-user=makanas@lanl.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue   # does not requeue when preempted

# Read specific parameters from pars.txt
M=$(grep -w "M" pars.txt | awk '{print $2}')       # Mass of BH
astar=$(grep -w "astar" pars.txt | awk '{print $2}')  # Spin of BH

# Export the values as environment variables so they can be used in the executable
export M
export astar

# Print values for debugging
echo "Mass (M): $M"
echo "Spin (astar): $astar"

# Calculate the start and end index based on the job array index (SLURM_ARRAY_TASK_ID)
start_index=$(( SLURM_ARRAY_TASK_ID * 2048 ))  # Start index for this job
end_index=$(( start_index + 2047 ))            # End index for this job (2048 jobs)

# Make sure the end index does not exceed the total number of rows in tot_Delta_J.txt (6700)
total_rows=6700
if [ $end_index -ge $total_rows ]; then
    end_index=$((total_rows - 1))
fi

# Use `sed` to extract the relevant chunk of rows from tot_Delta_J.txt
sed -n "${start_index},${end_index}p" tot_Delta_J.txt | while read -r row; do
    # Parse the row values into individual variables (assuming space-delimited)
    # Example for extracting some of the columns (adjust for your needs)
    col1=$(echo $row | awk '{print $1}')
    col2=$(echo $row | awk '{print $2}')
    col15=$(echo $row | awk '{print $15}')
    col16=$(echo $row | awk '{print $16}')
    col17=$(echo $row | awk '{print $17}')
    col21=$(echo $row | awk '{print $21}')
    col22=$(echo $row | awk '{print $22}')
    col23=$(echo $row | awk '{print $23}')

    # Now call the executable with the appropriate parameters
    ./Delta_Phi_single $M $astar $col15 $col16 $col17 $col21 $col22 $col23 $col1 $col2
done

