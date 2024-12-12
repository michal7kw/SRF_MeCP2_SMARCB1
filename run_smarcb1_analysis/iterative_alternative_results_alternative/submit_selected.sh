#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 <ending>"
    echo "Example: $0 neu_cs"
    echo "This will submit all files ending with neu_cs.sh"
    exit 1
fi

ending=$1

# Store files in an array and sort them
mapfile -t sorted_files < <(ls -1 run_smarcb1_analysis_*${ending}.sh 2>/dev/null | sort)

if [ ${#sorted_files[@]} -eq 0 ]; then
    echo "No files found matching pattern: run_smarcb1_analysis_*${ending}.sh"
    exit 1
fi

# Print the sorted files that will be submitted
echo "Files to be submitted:"
printf '%s\n' "${sorted_files[@]}"
echo "---"

# Submit first job separately to get its job ID
first_script="${sorted_files[0]}"
if [ -f "$first_script" ]; then
    echo "Submitting $first_script"
    job_id=$(sbatch "$first_script" | awk '{print $4}')
    echo "Submitted job ID: $job_id"
fi

# Submit remaining jobs with dependency on previous job
for ((i=1; i<${#sorted_files[@]}; i++)); do
    script="${sorted_files[i]}"
    if [ -f "$script" ]; then
        echo "Submitting $script with dependency on job $job_id"
        job_id=$(sbatch --dependency=afterok:$job_id "$script" | awk '{print $4}')
        echo "Submitted job ID: $job_id"
    fi
done

# # To submit all neu_cs.sh files:
# bash submit_selected.sh neu_cs

# # To submit all neu_ww.sh files:
# bash submit_selected.sh neu_ww

# # To submit all nsc_cs.sh files:
# bash submit_selected.sh nsc_cs