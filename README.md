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
This folder contains the code for the main and supplemental figures of my manuscript (except Figure 1). 

## misc_scripts
This folder contains the code for Figure 1, the creation of the metadata file, and a comparison of novel samples vs samples in Ungar et al., 2024 (PMID: PMC10802764).

## Gene Information
Contains the gene sets used in the manuscript. They were originally from Cormier, et al., 2022 (PMID: 36376793) and can also be found at https://github.com/macarthur-lab/gene_lists

Sources of the gene sets:

| Gene set | Description | Source |
| :---: | :---: | :---: |
| Haploinsufficient | ClinGen dataset (Cormier et al., 2021) | 
| Autosomal recessive | Blekhman et al., 2008; Berg et al., 2013 | 
| Autosomal dominant  Blekhman et al., 2008; Berg et al., 2013 | 
| Olfactory receptor | Mainland, et al., 2015 |
| CRISPR non-essential | Hart et al., 2017 | 
| Developmental delay | https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz > Firth et al., 2011;  Fitzgerald, et al., 2015; Wright et al., 2015; McRae, et al., 2017; Wright, et al., 2018 | 
| OMIM | https://omim.org/downloads > Amberger et al., 2019 |

