#!/usr/bin/env python

import os
import argparse


# Function that compute the bimolecular structure with the bifold function of RNAstructure
def bifold(file1, file2, output, rnaStructureDir):
    #bashCommand = "export DATAPATH=" + rnaStructureDir + "data_tables/ " \
    #              "&& " + rnaStructureDir + "exe/./bifold -p 2 '" + file1 + "' '" + file2 + "' '" + output + "'"
    bashCommand = "export DATAPATH=" + rnaStructureDir + "data_tables/ " \
                  "&& " + "bifold -p 2 '" + file1 + "' '" + file2 + "' '" + output + "'"
    os.system(bashCommand)


# Function that use RNAstructure draw function to draw all bimolecular function from a CT file
def drawAll(ctFile, imageFile, rnaStructureDir):
    bashCommand = "export DATAPATH=" + rnaStructureDir + "data_tables/ " \
                  "&& " + "draw '" + ctFile + "' '" + imageFile + "'"
    os.system(bashCommand)


# Function that list all fasta file in a directory
def listAllFastaFile(directory):
    fileList = []
    for file in os.listdir(directory):
        if file.endswith(".fasta"):
            fileList.append(file)
    return fileList


# ------------- Main functions ------------- #
def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Compute the bimolecule structure between each fasta file present in "
                                                 "the input directory and the partner fasta file. The ct files and png "
                                                 "files are written in the output directory.")
    parser.add_argument("--RNAstructureDir", action="store", dest="rnaStructureDir", type=str,
                        help="Directory of the RNAstructure tool. This directory should contain an exe directory and a "
                             "data_tables directory",
                        required=True)
    parser.add_argument("--inputDir", action="store", dest="directory", type=str,
                        help="Input directory. This directory should contain all mutated fasta file.",
                        required=True)
    parser.add_argument("--partnerFasta", action="store", dest="filePartner", type=str,
                        help="Path and name of the fasta file of the partner RNA. Example : U6atac",
                        required=True)
    parser.add_argument("--outputDir", action="store", dest="directoryOut", type=str,
                        help="Path to the output directory.",
                        required=True)
    parser.add_argument("--RNU4ATAC", action="store", dest="RNU4ATAC_file", type=str,
                        help="Path and name of the fasta file of RNU4ATAC fasta.",
                        required=True)
    options = parser.parse_args()

    RNU4ATAC_file = options.RNU4ATAC_file
    print(">>> begin computing bimolecule structures")
    files = listAllFastaFile(options.directory)
    namePartner = os.path.basename(options.filePartner).split(".fasta")[0]
    RNU4ATAC = os.path.basename(RNU4ATAC_file).split(".fasta")[0]
    # create the output directory if it does not exist
    if not os.path.exists(options.directoryOut):
        os.makedirs(options.directoryOut)
    #for file in files:
    #    name = os.path.basename(file).split(".fasta")[0]
    #    ctFile = options.directoryOut + name + "_" + namePartner + ".ct"
    #    imageFile = options.directoryOut + name + "_" + namePartner + ".ps"
    #    print("file1 = " + RNU4ATAC_file)
    #    print("file2 = " + options.filePartner)
    #    print("ct file = " + ctFile)
    #    print("image = " + imageFile)
    #    bifold(options.directory+file, options.filePartner, ctFile, options.rnaStructureDir)
        #drawAll(ctFile, imageFile, options.rnaStructureDir)
    
    ctFile = options.directoryOut + RNU4ATAC + "_" + namePartner + ".ct"
    print("ctFile output =" + ctFile)
    bifold(RNU4ATAC_file, options.filePartner, ctFile, options.rnaStructureDir)
    print(">>> end computing bimolecule structures")


if __name__ == "__main__":
    main()
