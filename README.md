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

## Paper
This folder contains the code and images for the figures of my upcoming manuscript. 
### Figure 1
Contains the code for:
-Figure 1A: Getting the number fo cases vs controls

-Figure 1B: Boxplot of the splicing metrics

-Figure 1C: Geneset enrichment analysis

Contains the data for the geneset enrichment analysis in the "Gene Information" folder

## ASHG_Poster
This folder contains the code for the poster I am making at ASHG.

## ASHG_Presentation
This folder contains the code for the plenary talk I am giving at ASHG. It contains four folders: 
-MIG_count: code for the MIG boxplot and geom points shown (labeled and unlabeled) 

-Metadata Analysis: code for the age and sex distributions (now removed from the talk)

-Theta vs Psi: code for the excess theta vs psi comparison (now removed from the talk)

-Excess Junctions: code (by gene and by junction) showing the excess outliers in geom_point form

## Projects
This folder contains the code for various images I have created for presentations and posters. 
