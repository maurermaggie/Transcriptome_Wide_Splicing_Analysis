#!/usr/bin/env python

import MutateSeq
import argparse
import os


# Function that read the file containing all mutations
# and save them in a list
def readAnnotationFile(annotationFile):
    file = open(annotationFile, "r")
    annotationData = []
    for line in file:
        if "Chromosome" in line:
            continue
        tmp = line.split("\t")
        mut = tmp[5]
        annotationData.append(mut)
    return annotationData


# Function that apply a list of mutation to a fasta file
# all mutated fasta files are created in pathOut folder
def mutateAll(idSeq, seqWT, annotationData, pathOut):
    for mut in annotationData:
        fastaFileOut = pathOut+mut+".fasta"
        MutateSeq.mutateSequence(idSeq, seqWT, fastaFileOut, mut)
    return 1


# ------------- Main functions ------------- #
def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Apply each mutation present in the mutation file in the given fasta "
                                                 "file. All output fasta files are written in the output directory and "
                                                 "start with prefix, followed by the name of the mutation.")
    parser.add_argument("--mutation", action="store", dest="annotationFile", type=str,
                        help="Path and name of the mutation file.",
                        required=True)
    parser.add_argument("--fasta", action="store", dest="fastaFileIn", type=str,
                        help="Path and name of the fasta file that should be mutated.",
                        required=True)
    parser.add_argument("--outputDir", action="store", dest="pathOutput", type=str,
                        help="Path to the output directory.",
                        required=True)
    parser.add_argument("--outputPrefix", action="store", dest="prefixOutput", type=str,
                        help="Prefix of the output file. Each output file name will start by this prefix followed by "
                             "the mutation name.",
                        required=True)
    options = parser.parse_args()
    print(">>> begin mutating all sequences")
    annotationData = readAnnotationFile(options.annotationFile)
    # print annotationData
    (idSeq, seqWT) = MutateSeq.readSeq(options.fastaFileIn)
    # create the output directory if it does not exist
    if not os.path.exists(options.pathOutput):
        os.makedirs(options.pathOutput)
    # mutate all seq
    output = options.pathOutput+options.prefixOutput
    mutateAll(idSeq, seqWT, annotationData, output)
    print(">>> end mutating all sequences")


if __name__ == "__main__":
    main()
