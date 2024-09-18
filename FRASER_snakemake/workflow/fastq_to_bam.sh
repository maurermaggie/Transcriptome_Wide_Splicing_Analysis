#!/bin/bash
set -euo pipefail

#########################################---Step 1: Define Functions---###########################################################
check_for_star_index() {
    argument_name="${1}"
    directory_path="${2}"
    if [[ ! -d ${directory_path} ]]; then
        echo "Error: filepath ${directory_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

check_for_fastq() {
    argument_name="${1}"
    file_path="${2}"
    if [[ ! -f ${file_path} ]]; then
        echo "Error: filepath ${file_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

#############################################---Step 2: Set Up Parameters---######################################################
options_array=(
    star_index
    fastq1
    fastq2
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g'):

arguments=$(getopt --options a --longoptions "${longoptions}" --name 'MPRA pipeline' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --star_index )
            star_index="${2}"; check_for_star_index "${1}" "${2}"; shift 2 ;;
        --star_out )
            output_directory="$2"; shift 2 ;;
        --fastq1 )
            fastq1="${2}"; check_for_fastq "${1}" "${2}"; shift 2 ;;
        --fastq2 )
            fastq2="${2}"; check_for_fastq "${1}" "${2}"; shift 2 ;;
        -- )
            shift; break ;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

#################################################---Step 2: Run Star---###########################################################
FASTQ_1_ABS=~{fastq1}
FASTQ_2_ABS=~{fastq2}

echo "FASTQs:"
echo "$FASTQ_1_ABS"
echo "$FASTQ_2_ABS"

# extract index
echo "--- $(date "+[%b %d %H:%M:%S]") Extracting STAR index ---"
mkdir star_index
tar -xvf ~{star_index} -C star_index --strip-components=1
echo "--- $(date "+[%b %d %H:%M:%S]")" Done extracting index

mkdir star_out
echo "--- $(date "+[%b %d %H:%M:%S]") Running STAR ---"
STAR --genomeDir star_index \
     --readFilesIn "$FASTQ_1_ABS" "$FASTQ_2_ABS" \
     --outFileNamePrefix star_out/~{SID}. \
     --readFilesCommand zcat \
     --outSAMattributes ~{outSAMattributes} \
     --outFilterType ~{outFilterType} \
     --runThreadN ~{ncpu} \
     --outSAMtype ~{outSAMtype} \
     --quantMode ~{quantMode} \
     --outSAMunmapped Within \
     --outSAMstrandField intronMotif \
     --outSAMattrRGline ID:rg1 PL:Illumina LB:~{SID} SM:~{SID}

cd star_out
ls

echo "--- $(date "+[%b %d %H:%M:%S]") Running samtools index ---"
samtools index ~{SID}.Aligned.sortedByCoord.out.bam

echo "--- $(date "+[%b %d %H:%M:%S]") Finished running samtools index ---"
ls
echo "--- $(date "+[%b %d %H:%M:%S]") Completed task ---"