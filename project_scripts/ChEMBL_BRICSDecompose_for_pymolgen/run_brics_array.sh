#!/bin/bash

#======================================================
#
# Job script for running a parallel job on a single node
#
#======================================================

# Based on example script: anaconda-slurm-serial.sh

#======================================================
# Propogate environment variables to the compute node
#SBATCH --export=ALL
#
# Run in the standard partition (queue)
#SBATCH --partition=standard
#
# Specify project account
#SBATCH --account=palmer-addnm
#
# No. of tasks required (max. of 40) (1 for a serial job)
#SBATCH --ntasks=1
#
# Specify (hard) runtime (HH:MM:SS)
#SBATCH --time=48:00:00
#
# Job name
#SBATCH --job-name=brics
#
# Output file
#SBATCH --output=slurm-%j.out
#
### Kill job 60 s before walltime limit to print output (Not sure this makes much difference)
###SBATCH --signal=15
#======================================================

module purge

module load anaconda/python-3.6.8/2019.03

#=========================================================
# Prologue script to record job details
# Do not change the line below
#=========================================================
/opt/software/scripts/job_prologue.sh 
#----------------------------------------------------------

# Modify the line below to run your program

source activate deepchem

# Run script as array job using command:
#sbatch --array 0-474 run_brics_array.sh

# Run slurm_stats.sh every 10 minutes and append results to slurm_stats-$SLURM_JOB_ID.out:
#/users/xpb20111/scripts/sys/slurm/watch_slurm.sh &

chembl_database_file="/users/xpb20111/ChEMBL/chembl_33_chemreps.txt.gz"
chembl_size=2372675 # $(zcat ${chembl_database_file} | wc -l)
brics_errors_file='brics_errors.txt'
n_mols_per_cpu=5000

brics_out_file="chembl_33_brics_${SLURM_ARRAY_TASK_ID}.txt.gz"
brics_error_file="brics_errors.txt"

let start_line=${SLURM_ARRAY_TASK_ID}*${n_mols_per_cpu}+1
let end_line=${start_line}+${n_mols_per_cpu}

#let max_line=${start_line}+${n_mols_per_cpu}
if [ ${end_line} -gt ${chembl_size} ]
then
    end_line=${chembl_size}
fi

echo "Start line: ${start_line}"
echo "End line: ${end_line}"

if [ -f ${brics_out_file} ]
then
    lines_so_far=$(zcat ${brics_out_file} | wc -l)

    # Check length of error file:
    if [ -f ${brics_error_file} ]
    then
        error_lines=$(awk -v array_task_id=${SLURM_ARRAY_TASK_ID} 'BEGIN{sum=0} $1==array_task_id {sum+=1}END{print sum}' ${brics_errors_file})
    else
        error_lines=0
    fi

    let current_line=${start_line}+${lines_so_far}-1+${error_lines}

    if [ ! ${current_line} -lt ${end_line} ]
    then
        echo "Run complete"
        exit 0
    fi

    # Write next molecule to errors file:
    zcat ${chembl_database_file} | awk -v line_no=${current_line} -v array_task_id=${SLURM_ARRAY_TASK_ID} 'NR==line_no+1 {print array_task_id"\t"line_no"\t"$1"\t"$2}' >> ${brics_errors_file}

    let start_line=${current_line}+1
fi

# Submit new job to continue where this job finishes in case of out of memory error:
dependency_job_id=$(sbatch \
       --parsable \
       --array ${SLURM_ARRAY_TASK_ID} \
       --dependency=afterany:${SLURM_JOB_ID} \
       run_brics_array.sh)

python bricsdecompose_chembl.py ${start_line} ${end_line} &> OUT-$SLURM_JOB_ID


#=========================================================
# Epilogue script to record job endtime and runtime
# Do not change the line below
#=========================================================
/opt/software/scripts/job_epilogue.sh 
#----------------------------------------------------------
