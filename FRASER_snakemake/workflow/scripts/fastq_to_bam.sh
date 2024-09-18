#!/bin/bash
#SBATCH --job-name=sort
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --partition=batch
#SBATCH --mem=256G
#SBATCH --time=24:00:00
#SBATCH --account=sjaiswal

echo "entering fastq_to_bam.sh"

BWA_GREF=$1
bam_path=$2

line_number=$SLURM_ARRAY_TASK_ID #get index of which file to process from $SLURM_ARRAY_TASK_ID provided by SLURM
fastq_file="${bam_path}/fastq_files" #provide path to file containing list of fastq files
fastq_prefix="$(sed "${line_number}q; d" "${fastq_file}")" #extract only the line number corresponding to $SLURM_ARRAY_TASK_ID

#Collect the sample name from the BAMs file
#Ex: If $bam_prefix is "/atac_seq/data/ATAC_tet2_KO_LDL", then the PREFIX is "ATAC_tet2_KO_LDL"
SAMPLE_NAME=$(basename "${fastq_prefix}")
echo "SAMPLE_NAME: $SAMPLE_NAME"

R1="${SAMPLE_NAME}_R1.fq.gz"
R2="${SAMPLE_NAME}_R2.fq.gz"
READGROUP="@RG\tID:${SAMPLE_NAME}\tLB:${SAMPLE_NAME}\tPL:illumina\tSM:${SAMPLE_NAME}"

echo "fastq_to_bam command used the following parameters: 
$0 $1"

echo "variables used:
SAMPLE_NAME: $SAMPLE_NAME
R1: $R1
R2: $R2
READGROUP: $READGROUP"

if [ ! -f "${SAMPLE_NAME}.bam" ]; then
    echo "Aligning to reference genome..."
    module load bwa
    module load samtools
    bwa mem -R "${READGROUP}" "${BWA_GREF}" "${R1}" "${R2}" | samtools view -b - | samtools sort - -o "${bam_path}/${SAMPLE_NAME}.bam"
    echo "...alignment complete"
else
    echo "Alignment already completed"
fi

if [ ! -f "${SAMPLE_NAME}.bam.bai" ]; then
    echo "Indexing BAM file..."
    module load samtools
    samtools index "${SAMPLE_NAME}.bam" "${SAMPLE_NAM}.bam.bai"
    echo "...indexing complete"
else
    echo "Indexing of BAM already completed"
fi

echo "fastq_to_bam.sh complete"