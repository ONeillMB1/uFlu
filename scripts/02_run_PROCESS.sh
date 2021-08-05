#!/bin/bash
#sbatch --array=1-12 --qos=geh -p geh --mem 100000 -J uFlu -o "logfile_custom_analyses_%a.log" 02_run_PROCESS.sh

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

mkdir ${TASKDIR}/03_custom_analyses/${PROJECT}/${ID}
mkdir ${TASKDIR}/03_custom_analyses/${PROJECT}/${ID}/n50
mkdir ${TASKDIR}/03_custom_analyses/${PROJECT}/${ID}/n75

Rscript ${TASKDIR}/resources/uFlu/scripts/process.R ${ID}
