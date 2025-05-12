# Transcriptome_Wide_Splicing_Analysis

## FRASER_snakemake
This folder contains my Snakemake pipeline for running FRASER (and not FRASER2 at the moment).

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
input_directory: "/file/path/where/you/want/to/symlink/data/to"
file_list: "/csv/file/where/you/have/list/of/files/to/run.csv"
FRASER_type: "Both"
```
FRASER_type can be "FRASER" (to indicate you want to run FRASER with the outputs of theta, psi3, and psi5), "FRASER2" (to indicate you want to run FRASER with the Jaccard index output), or "Both" (to indicate you want to receive all four outputs)

### Running FRASER
-Go to the FRASER_snakemake/workflow directory

-run the following command:
```
./run_snakemake --config_file "/path/to/config/yaml" --profile "path/to/config/slurm_scg"
```

The files in the slurm_scg directory will allow the snakemake to be run with resources on Stanford's scg/ oak

## run_results_paper
This folder contains the code for the main and supplemental figures of my manuscript (except Figure 1). 

### Installation
Please install the conda environments with:
```
micromamba create -f /path/to/fraser1.yml
micromamba create -f /path/to/fraser2.yml
```

The yml files can be found at:
```
run_results_paper/conda_envs
```

### Setting up configurations
-Make a config.yaml file in the config folder with the following information
```
FRASER1_results_uncompiled: "/file/path/to/raw/FRASER/output.csv"
input_file_FRASER: "/file/path/to/csv/with/filepaths/and/ids/of/all/samples/run/in/FRASER.csv"
FRASER2_results_uncompiled: "/file/path/to/raw/FRASER2/output.csv"
input_file_FRASER2: "/file/path/to/csv/with/filepaths/and/ids/of/all/samples/run/in/FRASER2.csv"
metadata_file: "/file/path/to/metadata/file.csv"
mig_file: "/file/path/to/Homo_sapiens_gene.csv"
low_RIN: "/file/path/to/csv/with/samples/to/exclude/due/to/low/RIN.csv"
output_directory: "/file/path/to/desired/output/directory"
genesets: "/file/path/to/directory/with/genesets"
genes: "/file/path/to/FRASER/output/rds/file.rds"
genes_FRASER2: "/file/path/to/FRASER2/output/rds/file.rds"
missing: "/file/path/to/csv/with/samples/to/exclude/due/to/missing/metadata.csv"
size_run_dir: "/file/path/to/directory/with/different/iterations/of/FRASER/across/different/run/sizes"
```

### Running run_results_paper
-Go to the run_results_paper/workflow directory

-run the following command:
```
./run_snakemake --config_file "/path/to/config/yaml" --profile "path/to/config/slurm_scg"
```

The files in the slurm_scg directory will allow the snakemake to be run with resources on Stanford's scg/ oak

## misc_scripts
This folder contains the code for Figure 1, the creation of the metadata file, and a comparison of novel samples vs samples in Ungar et al., 2024 (PMID: PMC10802764).

## Gene Information
Contains the gene sets used in the manuscript. They were originally from Cormier, et al., 2022 (PMID: 36376793) and can also be found at https://github.com/macarthur-lab/gene_lists

Sources of the gene sets:

| Gene set | Source | Filename | 
| :---: | :---: | :---: |
| Haploinsufficient | ClinGen dataset > Cormier et al., 2021 | haploinsufficient.tsv | 
| Autosomal recessive | Blekhman et al., 2008; Berg et al., 2013 | autosomal_recessive.tsv | 
| Autosomal dominant | Blekhman et al., 2008; Berg et al., 2013 | autosomal_dominant.tsv | 
| Olfactory receptor | Mainland, et al., 2015 | olfactory_receptors.tsv | 
| CRISPR non-essential | Hart et al., 2017 | CRISPR_nonessential_genes.tsv | 
| Developmental delay | https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz > Firth et al., 2011;  Fitzgerald, et al., 2015; Wright et al., 2015; McRae, et al., 2017; Wright, et al., 2018 | developmental_delay_genes.csv | 
| OMIM | https://omim.org/downloads > Amberger et al., 2019 | OMIM_genes.tsv |

