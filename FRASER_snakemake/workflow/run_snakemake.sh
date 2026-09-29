#!/bin/bash

#########################################---Step 1: Define Functions---###########################################################
check_for_config_file() {
    argument_name="${1}"
    file_path="${2}"
    if [[ ! -f ${file_path} ]]; then
        echo "Error: filepath ${file_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

check_for_profile() {
    argument_name="${1}"
    directory_path="${2}"
    if [[ ! -d ${directory_path} ]]; then
        echo "Error: filepath ${directory_path} passed with ${argument_name} does not exist."
        exit 1
    fi
}

function parse_yaml {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\):|\1|" \
        -e "s|^\($s\)\($w\)$s:$s[\"']\(.*\)[\"']$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/2;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}
#############################################---Step 2: Set Up Parameters---######################################################
options_array=(
    config_file
    profile
)

longoptions=$(echo "${options_array[@]}" | sed -e 's/ /:,/g'):

arguments=$(getopt --options a --longoptions "${longoptions}" --name 'MPRA pipeline' -- "$@")
eval set -- "${arguments}"

while true; do
    case "${1}" in
        --config_file )
            config_file="${2}"; check_for_config_file "${1}" "${2}"; shift 2 ;;
        --profile )
            profile="${2}"; check_for_profile "${1}" "${2}"; shift 2 ;;
        -- )
            shift; break ;;
        * )
            echo "Invalid argument ${1} ${2}" >&2
            exit 1
    esac
done

##########################################---Step 3: Clean yaml file---############################################################
#The script will produce an error if there are double slashes in filepaths. Thus, for all directories, they must NOT end in a /
#This code replaces all instances of /" with " in the yaml file.
sed -i 's/\/"/"/g' ${config_file}

#########################################---Step 3: Make Directories---###########################################################
eval $(parse_yaml ${config_file})
mkdir -p $output_directory
mkdir -p $input_directory

###########################################---Step 4: Run Snakemake---############################################################
if [ -z $config_file ]; then
    echo "no config_file selected; please select config_file to run snakemake on"
    exit 1
elif [ -z $profile ]; then
    echo "no profile selected"
    exit 1
else 
    micromamba run -n snakemake7 snakemake --configfile ${config_file} -c ${no_cores} --profile ${profile} --use-conda
fi
