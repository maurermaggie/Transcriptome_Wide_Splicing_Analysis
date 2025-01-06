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

| Gene set | Description | URL Source | Source |
| :---: | :---: | :---: | :---: 
| Haploinsufficient | https://github.com/maurermaggie/Transcriptome_Wide_Splicing_Analysis/blob/main/Gene_Information/haploinsufficient.tsv | ClinGen dataset (Cormier et al., 20212) | 
| Autosomal recessive | https://github.com/maurermaggie/Transcriptome_Wide_Splicing_Analysis/blob/main/Gene_Information/autosomal_recessive.tsv | Blekhman et al., 20083; Berg et al., 20134 | 
| Autosomal dominant | https://github.com/maurermaggie/Transcriptome_Wide_Splicing_Analysis/blob/main/Gene_Information/autosomal_dominant.tsv | Blekhman et al., 20083; Berg et al., 20134 | 
| Olfactory receptor | https://github.com/maurermaggie/Transcriptome_Wide_Splicing_Analysis/blob/main/Gene_Information/olfactory_receptors.tsv | Mainland, et al., 20155 |
| CRISPR non-essential | https://github.com/maurermaggie/Transcriptome_Wide_Splicing_Analysis/blob/main/Gene_Information/CRISPR_nonessential_genes.tsv | Hart et al., 20176 | 
| Developmental delay | https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz | Firth et al., 20117;  Fitzgerald, et al., 20158; Wright et al., 20159; McRae, et al., 201710; Wright, et al., 201811 | 
| OMIM | https://omim.org/downloads/ | OMIM dataset (Amberger et al., 201912) |

