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
-Figure 1A: Getting the number of cases vs controls

-Figure 1B: Boxplot of the splicing metrics

-Figure 1C: Geneset enrichment analysis

Contains the data for the geneset enrichment analysis in the "Gene Information" folder

### Figure 2
Contains the code for:
-Figure 2A: By both gene and junction, geom_point plots of the excess outliers from the different FRASER and FRASER2 metrics

-Figure 2B: PCA plot (removed from analysis)

-Figure 2C (now B): Metadata QC plots (Age, Batch, RIN, Sex) as well as Variance Explained. 

### Figure 3
Contains the code for:
-Figure 3A: Geom_point showing the separation of the RNU4ATAC samples from the rest of the cohort 

-Figure 3B: Boxplot of the number of genes impacted by the theta events in the RNU4ATAC samples compared to the rest of the cohort

### Figure 4
Contains the code for:
-Figure 4A: Geom_point showing the separation of the RNU4ATAC and RNU6ATAC samples from the rest of the cohort 

-Figure 4B: Boxplot of the number of genes impacted by the theta events in the RNU4ATAC and RNU6ATAC samples compared to the rest of the cohort

### Figure 5
Contains the code for:
-Plots showing the impact of sample size on the number of excess MIGs splicing efficiency outliers. 

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
