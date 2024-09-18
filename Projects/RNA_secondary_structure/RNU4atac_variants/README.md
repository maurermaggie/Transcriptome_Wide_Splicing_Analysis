# Prediction of the pathogenicity of RNU4ATAC GenomAD variants

These scripts allow to reproduce the analysis on RNU4ATAC variants from paper: [Benoit-Pilven et al, (2020)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0235655).


## Prerequisites

This code run with python 3. The following python libraries are needed:
- Biopython
- re
- argparse
- collections
- pprint

To compute the bimolecule secondary structure, we use the RNAstructure software.
This software can be downloaded [here](https://rna.urmc.rochester.edu/RNAstructure.html).
(Warning: The version of RNAstructure needed is the "text interface" not the "graphical interface".)


## Input files

To run this analysis, 6 files are needed :
- **variants list :** a csv file listing the variants to study (this file comes from GnomAD website)
- **new variants list :** a tab-delimited file containing variants that are not present in GnomAD but are known in patients
- **fasta 1 :** a fasta file containing the sequence of the RNA to mutate (RNU4atac)
- **fasta 2 :** a fasta file containing the sequence of the partner RNA in the bimolecular structure (RNU6atac)
- **structure annotations :** a tab-delimited file containing annotations about the bimolecular structure
- **pathologies annotations :** a tab-delimited file containing annotations about the known pathologies


##### Variants lists :
The 2 variants lists: csv file downloaded from GnomAD website and the list of new variants are concatenated in 1 file 
for the rest of the analysis.

The result file will have the following format:
- the first column must be the chromosome ("Chrom"), 
- the second column the position of the variants on the chomosome 
("Position")
- and the seventh column the denomination of the mutation (eg: n.2A>C).
- the following columns contain the allele frequency information for each population from GnomAD. 

The list of new variants must contain the columns 1, 2, 4, 5 and 7. All others columns can be filled with "0".

##### Fasta 1 :
The fasta file corresponding to the RNA for which we selected the mutation in GnomAD: RNU4atac.

##### Fasta 2 :
The fasta file of the partner RNA with which the first RNA forms a bimolecular structure: RNU6atac.

##### Structure annotations :
This file must have 4 columns. The first column is the position in the RNA ("position"), it must be intergers
starting at 1 (for the first nucleotide) until the length of the RNA. The second column is the annotation of region
in the bimolecular structure. The third column contain the known importance of the position from the litterature.
The possible values for this column is: "Limited", "Variable" and "Major". Here we took the importance of the
base/region for splicing. Finally, the fourth and last column contain information about the conservation of this
position. The possible values are: "Diverged", "Moderate/Low" and "High".

##### Pathologies annotations :
This file must have 3 columns. The first column is the variant (the denomination must be the same as in the 7th
column of the variants list. The second column contain the pathologies that has been associated with this variant.
And the third column contain the known zygosity ("Heterozygous", "Homozygous" or "Homozygous, Heterozygous").



## Run the analysis

To run the analysis, you should be in the "data" directory.

##### First step : Process the file downloaded from GnomAD
```
python ../ProcessGnomAD.py --input "gnomAD_v2.1_ENSG00000264229_2018_12_11_16_56_58.csv" --output "Variants_RNU4atac_GnomAD.tab" --newVariants "new_variants_RNU4atac.tab" --dir "input/"
```

##### Second step : Create the mutated sequence of RNU4atac
```
python ../AllMutationsSeq.py --mutation "input/Variants_RNU4atac_GnomAD.tab" --fasta "sequences_WT/RNU4atac_wild_type.fasta" --outputDir "sequences_Mut/" --outputPrefix "RNU4atac_mut_"
```

##### Third step : Predict the structure of the bimolecule for each mutation
```
python ../bimoleculeStructure.py --RNAstructureDir "../RNAstructure_Text_Interface/" --inputDir "sequences_Mut/" --partnerFasta "sequences_WT/RNU6atac_wild_type.fasta" --outputDir "RNAstructure/"
mkdir -p RNAstructure_WT/
export DATAPATH=../RNAstructure_Text_Interface/data_tables/ && ../RNAstructure_Text_Interface/exe/./bifold -p 2 'sequences_WT/RNU4atac_wild_type.fasta' 'sequences_WT/RNU6atac_wild_type.fasta' 'RNAstructure_WT/RNU4atac_wild_type_RNU6atac_wild_type.ct'
export DATAPATH=../RNAstructure_Text_Interface/data_tables/ && ../RNAstructure_Text_Interface/exe/./draw 'RNAstructure_WT/RNU4atac_wild_type_RNU6atac_wild_type.ct' 'RNAstructure_WT/RNU4atac_wild_type_RNU6atac_wild_type.ps'
```
Warning : don't forget to generate the WT structure...

##### Fourth step : Compare the mutated bimolecule structure with the wild type one
```
python ../CompareRNAStructure.py --inputDir "RNAstructure/" --structureWT "RNAstructure_WT/RNU4atac_wild_type_RNU6atac_wild_type.ct" --outputDir "differences/" --partnerName "RNU6atac"
```

##### Fifth step : Concatenate all annotations and information in one file
```
python ../AnnotateVariant.py --mutation "input/Variants_RNU4atac_GnomAD.tab" --annotation "annotations/annotations_RNU4atac.tab" --patient "annotations/annotations_RNU4atac_patients_mutations.tab" --caddScore "input/GRCh37-v1.4_bc228813829ae6a04c62f53709356e00.tsv" --inputDir "differences/" --output "results/Variants_RNU4atac_GnomAD_annotated.tab"
```
