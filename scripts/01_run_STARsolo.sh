#!/bin/bash

#SBATCH --mem=60000
#SBATCH --ntasks=12
#SBATCH -J solo
#SBATCH -o logfile_STARsolo_%A_%a.txt
#SBATCH -e logfile_STARsolo_%A_%a.txt
#SBATCH --array=1-<UPDATE TO NUMBER OF SAMPLES> 

##########################################################################################################
##########################################################################################################
##### This script has been written by Mary O'Neill in August 2021 for the uFlu project.
##### It is designed to work on Pasteur's Maestro server.
##### The script will execute the STARsolo pipeline on the custom in-drop-like design
##### of the project consisting of 4 barcodes with linkers. 
##### Update paths where necessary to execute on new sequencing data.
##### To launch the script on Maestro use the sbatch command commented out on line 2.
##### Modify the '--array' parameter to reflect the number of experiments.
##### The script requires a text file with the root names of the fastq files - 
##### update the name and path if needed.
##########################################################################################################
##########################################################################################################

SAMPLE_NUM=${SLURM_ARRAY_TASK_ID} 

PROJECT="<CHANGE ME>"
FASTQPATH="/pasteur/zeus/projets/p01/uFlu/reassortment_project/01_raw_seq_data/NGS11/Novaseq_150621/fastq" #modify for new sequencing run
REFPATH="/pasteur/zeus/projets/p01/uFlu/reassortment_project/resources/Ref_data"
TASKDIR="/pasteur/zeus/projets/p01/uFlu/reassortment_project/02_STARsolo_outputs"

# Extract the sample root name and expected number of cells from the input file
line=`head -n ${SAMPLE_NUM} ${TASKDIR}/${PROJECT}/${PROJECT}_samples.txt | tail -n 1`
ID=`echo $line | cut -d ' ' -f 1`
cells=`echo $line | cut -d ' ' -f 2`

#load STAR
module purge 
module load STAR/2.7.8a

start=`date +%s`
hostname

STAR --genomeDir ${REFPATH}/Set2/ \
	--readFilesCommand zcat \
	--readFilesIn ${FASTQPATH}/${ID}_R2_001.fastq.gz ${FASTQPATH}/${ID}_R1_001.fastq.gz \
	--outFileNamePrefix ${TASKDIR}/${PROJECT}/${ID}/${ID}. \
	--runThreadN 12 \
	--soloType CB_UMI_Complex \
	--soloCBposition 0_0_0_15 0_20_0_35 0_40_0_55 0_60_0_75 \
	--soloUMIposition 0_94_0_103 \
	--soloCBwhitelist ${REFPATH}/A_whitelist.txt ${REFPATH}/B_whitelist.txt ${REFPATH}/C_whitelist.txt ${REFPATH}/D_whitelist.txt \
	--soloCBmatchWLtype Exact \
	--soloCellFilter None \
	--outFilterScoreMinOverLread 0 \
	--outFilterMatchNminOverLread 0 \
	--outFilterMatchNmin 30 \
	--soloStrand Reverse \
	--outSAMtype BAM SortedByCoordinate \
	--limitBAMsortRAM 28811888672 \
	--outSAMattributes NH HI nM AS CR UR CB UB GX GN sM

end=`date +%s`
runtime=$((end-start))
echo $runtime

start=`date +%s`

STAR --runMode soloCellFiltering ${TASKDIR}/${PROJECT}/${ID}/${ID}.Solo.out/Gene/raw \
	${TASKDIR}/${PROJECT}/${ID}/${ID}.Solo.out/Gene/emptyDrops_50counts/ \
	--soloCellFilter EmptyDrops_CR ${cells} 0.99 10 45000 90000 50 0.01 20000 0.01 10000 \
	--outFileNamePrefix ${TASKDIR}/${PROJECT}/${ID}/${ID}_soloCellFiltering_n50counts. \

STAR --runMode soloCellFiltering ${TASKDIR}/${PROJECT}/${ID}/${ID}.Solo.out/Gene/raw \
	${TASKDIR}/${PROJECT}/${ID}/${ID}.Solo.out/Gene/emptyDrops_75counts/ \
	--soloCellFilter EmptyDrops_CR ${cells} 0.99 10 45000 90000 75 0.01 20000 0.01 10000 \
	--outFileNamePrefix ${TASKDIR}/${PROJECT}/${ID}/${ID}_soloCellFiltering_n75counts. \

mv ${TASKDIR}/${PROJECT}/${ID}/${ID}.Solo.out/Gene/raw/features.tsv ${TASKDIR}/${PROJECT}/${ID}/${ID}.Solo.out/Gene/raw/genes.tsv
mv ${TASKDIR}/${PROJECT}/${ID}/${ID}.Solo.out/Gene/emptyDrops_50counts/features.tsv ${TASKDIR}/${PROJECT}/${ID}/${ID}.Solo.out/Gene/emptyDrops_50counts/genes.tsv
mv ${TASKDIR}/${PROJECT}/${ID}/${ID}.Solo.out/Gene/emptyDrops_75counts/features.tsv ${TASKDIR}/${PROJECT}/${ID}/${ID}.Solo.out/Gene/emptyDrops_75counts/genes.tsv

chmod -R 775 ${TASKDIR}/${PROJECT}/${ID} 

end=`date +%s`
runtime=$((end-start))
echo $runtime
