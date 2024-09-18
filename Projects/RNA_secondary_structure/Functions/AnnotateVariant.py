#!/usr/bin/env python

import collections
import os
import argparse
import glob


# Function that take into input a file containing all variants
# and save this data in a dictionary
def readVariantFile(variantFile):
    file = open(variantFile, "r")
    variantData = collections.defaultdict(dict)
    header = ""
    for line in file:
        if "Chrom" in line:
            line = line.rstrip()
            header = line
            continue
        line = line.rstrip()
        tmp = line.split("\t")
        mut = tmp[5]
        pos = tmp[1]
        variantData[int(pos)][mut] = tmp
    file.close()
    return variantData, header


# Function that take into input a file containing the annotation of the gene of interest
# (columns of this file : position; Region; Importance_for_splicing; Conservation)
# and save this info in a dictionary
def readAnnotationFile(annotationFile):
    file = open(annotationFile, 'r')
    annotationData = {}
    header = ""
    for line in file:
        if "position" in line:
            line = line.rstrip()
            header = line
            continue
        line = line.rstrip()
        tmp = line.split("\t")
        pos = int(tmp[0])
        annotationData[pos] = tmp
    file.close()
    return annotationData, header


# Function that take into input a file containing the annotation of the mutations found in patients
# and save this data in a dictionary
def readVariantPatientFile(variantPatientFile):
    file = open(variantPatientFile, 'r')
    variantPatientData = collections.defaultdict(dict)
    header = ""
    for line in file:
        if "mutation" in line:
            line = line.rstrip()
            tmp = line.split("\t")
            header = "\t".join(tmp[1:])
            continue
        line = line.rstrip()
        tmp = line.split("\t")
        variantPatientData[tmp[0]] = tmp[1:]
    file.close()
    return variantPatientData, header


# Function that take into input the file containing the CADD scores
# and save the cadd score in a dictionary 
def readCADDscoreFile(caddScoreFile, posInit):
    file = open(caddScoreFile, 'r')
    caddScoreData = collections.defaultdict(dict)
    for line in file:
        if "##" in line:
            continue
        elif "#CHROM" in line:
            continue
        line = line.rstrip()
        tmp = line.split("\t")
        pos = int(tmp[1])
        ref = tmp[2]
        alt = tmp[3]
        caddScoreData[pos][ref+":"+alt] = tmp
    file.close()
    return caddScoreData


# Functions that take into input a comparison structure file and an annotation data
# and compute a score : high score => important differences in structure
#                       low score => small or not important differences in structure
# First score method : this is the method currently used
def computeScore(structureData, annotationData):
    valueMajor = 3
    valueVariable = 1
    valueLimited = 0.5
    major = 0
    variable = 0
    limited = 0
    for pos in structureData:
        if annotationData[int(pos)][2] == "Major":
            major += 1
        elif annotationData[int(pos)][2] == "Variable":
            variable += 1
        elif annotationData[int(pos)][2] == "Limited":
            limited += 1
    score = major*valueMajor + variable*valueVariable + limited*valueLimited
    return score


# Second score method : this method take into account the conservation of the base to compute the score
def computeScoreConservation(structureData, annotationData):
    valueMajor = 3
    valueVariable = 1
    valueLimited = 0
    valueConsHigh = 3
    valueConsModerate = 2
    valueConsDiverged = 1
    score = 0
    for pos in structureData:
        importanceSplicing = annotationData[int(pos)][2]
        conservation = annotationData[int(pos)][3]
        if importanceSplicing == "Major":
            if conservation == "High":
                score += valueMajor * valueConsHigh
            elif conservation == "Moderate/Low":
                score += valueMajor * valueConsModerate
            elif conservation == "Diverged":
                score += valueMajor * valueConsDiverged
            else:
                print("Warning : Conservation = " + conservation + " not defined !")
        elif importanceSplicing == "Variable":
            if conservation == "High":
                score += valueVariable * valueConsHigh
            elif conservation == "Moderate/Low":
                score += valueVariable * valueConsModerate
            elif conservation == "Diverged":
                score += valueVariable * valueConsDiverged
            else:
                print("Warning : Conservation = " + conservation + " not defined !")
        elif importanceSplicing == "Limited":
            if conservation == "High":
                score += valueLimited * valueConsHigh
            elif conservation == "Moderate/Low":
                score += valueLimited * valueConsModerate
            elif conservation == "Diverged":
                score += valueLimited * valueConsDiverged
            else:
                print("Warning : Conservation = " + conservation + " not defined !")
        else:
            print("Warning : Importance for splicing = " + importanceSplicing + " not defined !")
    return score


# Function that take into input a comparison structure file
# output "Yes" if the file contain data except the header and "No" otherwise
# and a score : high score => important differences in structure
#               low score => small or not important differences in structure
def analyzeDifferencesFile(diffFile, annotationData, energyCtrl):
    file = open(diffFile, 'r')
    comparisonStructureData = []
    structureDiff = "No"
    energy = 0
    ctrlEnergyDiff = 999
    for line in file:
        if "position" in line:
            element = line.rsplit()
            energy = (element[2].split("="))[1]
            ctrlEnergyDiff = float(energyCtrl) - float(energy)
            continue
        structureDiff = "Yes"
        line = line.rstrip()
        tmp = line.split("\t")
        comparisonStructureData.append(tmp[0])
    file.close()
    score = computeScore(comparisonStructureData, annotationData)
    return structureDiff, score, energy, ctrlEnergyDiff


# Function that take into input a comparison structure file
# and get the "best" structure (closest structure to the control)
# output "Yes" if this structure is exactly the same as the control and "No" otherwise
# and a score : high score => important differences in structure
#               low score => small or not important differences in structure
def analyzeAllDifferencesFile(diffFile, annotationData, energyCtrl):
    file = open(diffFile, 'r')
    comparisonStructureData = {}
    comparisonEnergyData = {}
    structureNb = 0
    passLines = 0
    # read each line of diffFile and save the structure
    for line in file:
        if "position" in line:
            structureNb += 1
            comparisonStructureData[structureNb] = []  # initialize the structure in the dictionary
            passLines = 0
            element = line.rsplit()
            energyTmp = element[2].split("=")[1]
            ctrlEnergyDiffTmp = float(energyCtrl) - float(energyTmp)
            comparisonEnergyData[structureNb] = [energyTmp, ctrlEnergyDiffTmp]
            continue
        if passLines == 1:
            continue
        line = line.rstrip()
        tmp = line.split("\t")
        comparisonStructureData[structureNb].append(tmp[0])
    file.close()
    # for each structure saved, compute a score
    bestScore = 999
    bestStructure = 0
    structureDiff = "Yes"
    firstEnergy = comparisonEnergyData[1][0]
    firstEnergyDiff = "NA"
    ctrlEnergyDiff = "NA"
    for structure in comparisonStructureData.keys():
        score = computeScore(comparisonStructureData[structure], annotationData)
        energy = comparisonEnergyData[structure][0]
        if score < bestScore:
            bestScore = score
            bestStructure = structure
            firstEnergyDiff = float(firstEnergy) - float(energy)
            ctrlEnergyDiff = comparisonEnergyData[structure][1]
            if bestScore == 0:
                structureDiff = "No"
    return structureDiff, bestScore, bestStructure, firstEnergyDiff, ctrlEnergyDiff


# Function that extract the mutation name from the diff file name
def getMutation(fileName, partner="", type="bimolecule"):
    mut = ""
    if type == "bimolecule":
        toSplit = "_"+partner
        mut = fileName.split("mut_")[1].split(toSplit)[0]
    elif type == "monomer":
        mut = fileName.split("mut_")[1].split(".tab")[0]
    else:
        print("WARNING: invalid type " + type)
    return mut


# Function that extract the "best structure" from all the sub-optimal structure analyzed
# The "best structure" is the one with the smallest score (the structure more similar to the WT one)
def extractBestStructureData(structureDirectory, annotationData, energyCtrl, partner="", type="bimolecule"):
    structureData = {}
    # analyze all comparison structure files from the structure directory
    for file in glob.glob(structureDirectory+"*.tab"):
        mut = getMutation(file, partner=partner, type=type)
        filePath = file
        structureDiff, score, bestStructure, firstEnergyDiff, ctrlEnergyDiff = analyzeAllDifferencesFile(filePath,
                                                                                                         annotationData,
                                                                                                         energyCtrl)
        structureData[mut] = [structureDiff, str(score), str(bestStructure)]#, str(firstEnergyDiff), str(ctrlEnergyDiff)]
    return structureData


# Function that merge all information in one dict
# It takes in input all available data (annotations, structure, pathologies, ...)
def mergeData(variantData, annotationData, caddScoreData, structureDataBimolecule, posInit,
              structureDataMonomer=None, variantPatientData=None):
    # assign the default value {} in case these arguments are equal to None
    if structureDataMonomer is None:
        structureDataMonomer = {}
    if variantPatientData is None:
        variantPatientData = {}
    mergedData = {}
    for position in variantData.keys():
        #print(position)
        if position in caddScoreData.keys():
            if position-posInit in annotationData.keys():
                for mutation in variantData[position].keys():
                    #print(mutation)
                    if mutation in structureDataBimolecule.keys():
                        phenotype = ["-", ""]
                        if mutation in variantPatientData.keys():
                            phenotype = variantPatientData[mutation]
                        refalt = variantData[position][mutation][3] + ":" + variantData[position][mutation][4]
                        if refalt in caddScoreData[position].keys():
                            if mutation in structureDataMonomer.keys():
                                annotatedMutation = variantData[position][mutation] + annotationData[position - posInit] + \
                                                    phenotype + structureDataBimolecule[mutation] + \
                                                    structureDataMonomer[mutation] + caddScoreData[position][refalt][4:]
                                mergedData[mutation] = annotatedMutation
                            else:
                                annotatedMutation = variantData[position][mutation] + annotationData[position-posInit] + \
                                                    phenotype + structureDataBimolecule[mutation] + caddScoreData[position][refalt][4:]
                                mergedData[mutation] = annotatedMutation
                        else:
                            print("Mutation ("+refalt+") not present in CADD score file.")
                    else:
                        print("Mutation ("+mutation+") not present in structure data of the bimolecule.")
            else:
                print("position ("+str(position)+") not found in the annotation file.")
        else:
            print("position ("+str(position)+") not found in the CADD score file.")
    return mergedData


# Function that write the merged data in the final result file
def writeMergeData(mergedData, fileName, header):
    file = open(fileName, "w")
    file.write(header)
    for mutation in mergedData.keys():
        lineToWrite = "\t".join(mergedData[mutation])+"\n"
        file.write(str(lineToWrite))
    file.close()


# ------------- Main functions ------------- #
def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Annotate the variants with information about:"
                                                 "pathogenicity, bimolecule variation score, region...")
    parser.add_argument("--mutation", action="store", dest="variantFile", type=str,
                        help="Path and name of the mutation file.",
                        required=True)
    parser.add_argument("--annotation", action="store", dest="annotationFile", type=str,
                        help="Path and name of the annotation file.",
                        required=True)
    parser.add_argument("--patient", action="store", dest="variantPatientFile", type=str,
                        help="Path and name of the patient file.",
                        required=True)
    parser.add_argument("--caddScore", action="store", dest="caddScoreFile", type=str,
                        help="Path and name of the file containing the CADD score.",
                        required=True)
    parser.add_argument("--inputDir", action="store", dest="structureDirectoryBimolecule", type=str,
                        help="Input directory. This directory should contain all comparison files for the bimolecule "
                             "structure.",
                        required=True)
    parser.add_argument("--energyWT", action="store", dest="energyWT", type=str,
                        help="Free energy value for the wild type structure.",
                        default='-88.1')
    parser.add_argument("--posInit", action="store", dest="posInit", type=int,
                        help="Start position of the gene on chromosome By default for U4atac: 122288455.",
                        default=122288455)
    parser.add_argument("--output", action="store", dest="outputFile", type=str,
                        help="Path and name of the output file.",
                        required=True)
    parser.add_argument("--partnerName", action="store", dest="partnerName", type=str,
                        help="Name of the partner RNA. Example : RNU6atac",
                        default="RNU6atac")

    options = parser.parse_args()
    print(">>> begin RNU4atac annotate variant")
    (variantData, headerVariantFile) = readVariantFile(options.variantFile)
    (annotationData, headerAnnotationFile) = readAnnotationFile(options.annotationFile)
    (variantPatientData, headerPatientFile) = readVariantPatientFile(options.variantPatientFile)
    caddScoreData = readCADDscoreFile(options.caddScoreFile, options.posInit)
    structureDataBimolecule = extractBestStructureData(options.structureDirectoryBimolecule, annotationData,
                                                       options.energyWT, partner=options.partnerName)
    header = headerVariantFile + "\t" + headerAnnotationFile + "\t" + headerPatientFile + "\t" \
             + "Modif_structure_bimolecule" + "\t" + "Score_modif_structure_bimolecule" + "\t" \
             + "Best_structure_bimolecule" + "\t" + "CADD rawScore" + "\t" + "CADD phredScore" + "\n"
             #+ "Energy_diff_structure_bimolecule" + "\t" + "Energy_diff_WT_structure_bimolecule" 
    mergedData = mergeData(variantData, annotationData, caddScoreData, structureDataBimolecule, options.posInit,
                           variantPatientData=variantPatientData)
    writeMergeData(mergedData, options.outputFile, header)
    print(">>> end RNU4atac annotate variant")


if __name__ == "__main__":
    main()
