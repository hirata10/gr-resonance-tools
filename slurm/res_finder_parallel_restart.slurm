#!/bin/bash

# Submit this script with: sbatch --array=1-16 res_finder_parallel_restart.slurm

#SBATCH --time=16:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --error=error_res_finder_parallel_restart_50_200-0.1_%A_%a.txt
#SBATCH --output=output_res_finder_parallel_restart_50_200-0.1_%A_%a.txt
#SBATCH --job-name=res_finder_parallel_50_200_0.1_restart   # job name
#SBATCH --mail-user=makanas@lanl.gov   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --no-requeue   # does not requeue when preempted

# === Runtime Setup ===
echo "number of nodes = $SLURM_NNODES"
echo "SLURM says:"
echo "number of CPUs per node = $SLURM_CPUS_ON_NODE"
echo "number of tasks per node = $SLURM_NTASKS_PER_NODE"
ranks_per_node=1
total_ranks=$(($SLURM_NNODES * $ranks_per_node))
echo "number of ranks total = $total_ranks"

# === Run Python Script ===
cmd="python3 resfinder_parallel_restart.py $SLURM_ARRAY_TASK_ID"
echo "Running: $cmd"
$cmd

# === Completion Check & Conditional Resubmit ===
sleep 2
OUTPUT_FILE="output_res_finder_parallel_restart_50_200-0.1_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt"

# Track the number of retries
MAX_RETRIES=3
RETRY_FILE="retry_count_50_200-0.1_${SLURM_ARRAY_TASK_ID}.txt"

if [[ -f $RETRY_FILE ]]; then
    RETRY_COUNT=$(cat $RETRY_FILE)
else
    RETRY_COUNT=0
fi

if grep -q "Completed resonance finder" "$OUTPUT_FILE"; then
    echo "Task $SLURM_ARRAY_TASK_ID completed successfully."
    rm -f $RETRY_FILE  # Remove retry file on success
else
    echo "Task $SLURM_ARRAY_TASK_ID did NOT finish — attempt $((RETRY_COUNT+1))/$MAX_RETRIES"

    if (( RETRY_COUNT < MAX_RETRIES )); then
        (( RETRY_COUNT++ ))
        echo $RETRY_COUNT > $RETRY_FILE  # Update retry count

        echo "Resubmitting Task $SLURM_ARRAY_TASK_ID..."
        sbatch --array=${SLURM_ARRAY_TASK_ID}-${SLURM_ARRAY_TASK_ID} $0
    else
        echo "Task $SLURM_ARRAY_TASK_ID exceeded retry limit. Not resubmitting."
        rm -f $RETRY_FILE  # Clean up retry file if no more resubmissions are needed
    fi
fi


