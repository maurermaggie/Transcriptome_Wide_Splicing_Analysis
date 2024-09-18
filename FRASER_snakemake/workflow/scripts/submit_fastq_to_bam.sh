##################################################################################################################################
#############################################---STEP 1: SET UP PARAMETERS---###################################################### 
##################################################################################################################################
if [ -z $1 ] || [ -z $2 ]; then
    echo "Format: ./submit_fastq_to_bam.sh [data_directory] [output_path]"
    echo "user can use argument --genome_build to specify the genome build used when converting fastqs to bams"
    exit 1
else
#pipeline in shell for ATACseq footprinting
    TEMP=`getopt -o vdm: --long genome_build: \
        -n './submit_fastq_to_bam' -- "$@"`

       if [ $? != 0 ]; then
           echo "Unrecognized argument. Possible arguments: --genome_build path_to_genome_build." >&2 ; exit 1 ; 
       fi
           eval set -- "$TEMP"

            genome_build="/oak/stanford/groups/sjaiswal/Herra/CHIP_Panel_AmpliSeq/GRCh38.p12.genome.u2af1l5_mask.fa"
            
        while true; do
            case "$1" in
                --genome_build ) genome_build="$2"; shift 2 ;;
                -- ) shift; break ;;
                * ) break ;;
            esac
        done

        echo "genome build selected $genome_build"

    data_directory=$1
    output_path=$2
    bwa_gref=$genome_build

    code_directory=$( realpath . )

    module load samtools/1.9

##################################################################################################################################
#######################################---STEP 2: CREATE NECESSARY FOLDERS---#####################################################
##################################################################################################################################
    if [ ! -d "$output_path/Logs" ]; then 
            mkdir -p "$output_path/Logs"
    fi
            
    if [ ! -d "$output_path/Parameters" ]; then
            mkdir -p "$output_path/Parameters"
    fi

    Logs="${output_path}/Logs"
    Parameters="${output_path}/Parameters"

##################################################################################################################################
#######################################---STEP 3: CREATE PARAMETER LOG---#########################################################
##################################################################################################################################
    now=$(date +%m_%d_%H_%M)
    
    parameter_file="$Parameters/${now}_parameters.txt"                                 #add date stamp to parameter files and, if provided, the log name
    
    touch $parameter_file
        
    if [ $genome_build != "/oak/stanford/groups/sjaiswal/Herra/CHIP_Panel_AmpliSeq/GRCh38.p12.genome.u2af1l5_mask.fa" ]; then     
        set_genome_build="--genome_build ${genome_build}"                                                         
    fi

    echo "location of scripts used to run code : $code_directory
        " > $parameter_file 
        
    echo "call made to execute code: $0 $1 $2 $set_genome_build  
    " >> $parameter_file
            
 ##################################################################################################################################
############################################---STEP 4: INDEX GENOME BUILD---###################################################### 
##################################################################################################################################
    if [ ! -f "${genome_build}.bwt" ]; then   
            module load bwa
            cd $genome_build

            echo "indexing genome build"
            bwa index -a bwtsw $genome_build
    fi

##################################################################################################################################
##########################################---STEP 5: CONVERT FASTQS TO BAMS---#################################################### 
##################################################################################################################################
    cd $data_directory 
    
    fastq_list="$data_directory/fastq_files" #give a path to a file to store the paths to the fastq files in $fastq_directory

   find "${data_directory}/" -type f `#list all files in ${fastq_directory}` | \
        grep ".*\.fq.gz$" `#only keep files with FASTQ in name (case insensitive)` | \
        grep -v "Undetermined" `#remove Undetermined FASTQs` | \
        sed -e 's/_R1.*$//g' | sed -e 's/_R2.*$//g' `#remove _R1/2_fastq.gz file extension`| \
        sort -u  `#sort and remove duplicate names` > ${fastq_list} 
    fastq_array_length=$(wc -l < ${fastq_list}) #get the number of FASTQs 
    echo "fastq array length: $fastq_array_length"
     
    sbatch -o "$Logs/%A_%a.log" `#put into log` \
            -a "1-${fastq_array_length}" `#initiate job array equal to the number of fastq files` \
            -W `#indicates to the script not to move on until the sbatch operation is complete` \
            "${code_directory}/fastq_to_bam.sh" \
            "$bwa_gref" \
            "$data_directory"

    echo "all fastqs have been converted to bams"

fi