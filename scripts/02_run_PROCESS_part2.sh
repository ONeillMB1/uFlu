#!/bin/bash

#SBATCH --mem=100000
#SBATCH -J uFlu
#SBATCH -o logfile_customAnalyses_%A_%a.txt
#SBATCH -e logfile_customAnalyses_%A_%a.txt
#SBATCH --array=1-12
#SBATCH --qos=fast
#SBATCH --partition=dedicated

PROJECT="NGS11"
TASKDIR="/pasteur/zeus/projets/p01/uFlu/reassortment_project"
SAMPLE_NUM=${SLURM_ARRAY_TASK_ID}
line=`head -n ${SAMPLE_NUM} ${TASKDIR}/02_STARsolo_outputs/${PROJECT}/${PROJECT}_samples.txt | tail -n 1`
ID=`echo $line | cut -d ' ' -f 1`

# load modules you need: 
module purge # remove modules that may have been loaded by mistake
module load R/4.1.0

start=`date +%s`
hostname
${ID}
${SAMPLE_NUM}

Rscript ${TASKDIR}/resources/uFlu/scripts/process_part2_supplement.R ${ID} ${PROJECT}

chmod -R 775 ${TASKDIR}/03_custom_analyses/${PROJECT}/${ID}
