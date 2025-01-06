# Transcriptome_Wide_Splicing_Analysis

## FRASER_Snakemake
This folder contains my Snakemake pipeline for running FRASER (and not FRASER2 at the moment).
More details coming soon.

### Installation
Please install the conda environments with:
```
micromamba create -f /path/to/fraser1.yml
micromamba create -f /path/to/fraser2.yml
```

The yml files can be found at:
```
FRASER_snakemake/conda_envs
```
### Making file list
-Make a csv with a list of all of the full filepaths of the files you want to be included in your FRASER run. 

-For an example, see symlink_blood.csv in the FRASER_snakemake/config directory

### Setting up configurations
-Make a config.yaml file in the config folder with the following information
```
output_directory: "/file/path/to/output/directory"
metadata_file: "file/path/to/metadata/file"
data_type: "Tissue Type"
mig_file: "/file/path/to/Homo_sapiens_gene.csv"
input_directory: "/file/path/where/you/want/to/symlink/data/to"
file_list: "csv/file/where/you/have/list/of/files/to/run.csv"
```
### Running FRASER
-Go to the FRASER_snakemake/workflow directory

-run the following command:
```
./run_snakemake --config_file "/path/to/config/yaml" --profile "path/to/config/slurm_scg"
```

The files in the slurm_scg directory will allow the snakemake to be run with resources on Stanford's scg/ oak

## run_results_paper
This folder contains the code and images for the main figures of my manuscript. 

## Gene Information
Contains the genesets used in the manuscript. They were originally from Cormier, et al., 2022 (PMID: 36376793) and can also be found at https://github.com/macarthur-lab/gene_lists


