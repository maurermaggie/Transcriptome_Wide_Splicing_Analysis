#!/usr/bin/env python

import re
import collections
import os
import argparse


# Function that parse the first structure in a CT file
# Outputs : Position and match of each base in a dict.
#           Energy value of the structure.
def parseOneStructureCTFile(CTFile):
    file = open(CTFile, "r")
    dataCT = {}
    energy = 0
    # get position and match of each base of RNU4atac (1 to 130)
    for line in file:
        if "ENERGY" in line:
            element = line.split()
            energy = element[3]
            continue
        element = line.split()
        if int(element[5]) == 0:
            break
        dataCT[int(element[5])] = int(element[4])
    file.close()
    return dataCT, energy


# Function that parse all structures in a CT file
# Outputs : Position and match of each base of each structure in a dict of dict.
#           Energy value of all structure in a dict.
def parseAllStructuresCTFile(CTFile):
    file = open(CTFile, "r")
    dataCT = collections.defaultdict(dict)
    energy = {}
    structure = 0
    passLines = 0
    # get position and match of each base of RNU4atac (1 to 130)
    for line in file:
        if "ENERGY" in line:
            structure += 1
            element = line.split()
            energy[structure] = element[3]
            passLines = 0
            continue
        if passLines == 1:
            continue
        element = line.split()
        if int(element[5]) == 0:
            passLines = 1
        dataCT[structure][int(element[5])] = int(element[4])
    file.close()
    return dataCT, energy


# Function that compute the offset needed to match each base of the WT sequence and the mutant sequence
# Outputs : value of the offset (nbr of base of the difference)
#           start of the offset (at which base should the offset be applied)
def computeOffset(mutation):
    offset = 0
    startOffset = 0
    mut = mutation[2:]
    if '>' in mut:
        pass
    else:
        match = re.search('^\d+', mut)
        if 'del' in mut:
            if match:
                startOffset = int(match.group(0))
            else:
                print("WARNING : no position for this deletion (" + mut + ")")
            matchSeq = re.search('(?<=del)[ATGC]+$', mut)
            sequence = ""
            if matchSeq:
                sequence = matchSeq.group(0)
            else:
                print("WARNING : no sequence for this deletion (" + mut + ")")
            offset = -len(sequence)
        elif 'ins' in mut:
            if match:
                startOffset = int(match.group(0))
            else:
                print("WARNING : no position for this insertion (" + mut + ")")
            matchSeq = re.search('(?<=ins)[ATGC]+$', mut)
            sequence = ""
            if matchSeq:
                sequence = matchSeq.group(0)
            else:
                print("WARNING : no sequence for this insertion (" + mut + ")")
            offset = len(sequence)
        elif 'dup' in mut:
            if match:
                startOffset = int(match.group(0))
            else:
                print("WARNING : no position for this duplication (" + mut + ")")
            matchPosEnd = re.search('(?<=_)\d+', mut)
            posEnd = startOffset
            if matchPosEnd:
                posEnd = int(matchPosEnd.group(0))
            offset = posEnd - startOffset + 1
            startOffset = startOffset + offset - 1
        else:
            raise Exception("Problem with mutation "+mutation)
    return offset, startOffset


# Function that compare one CT structure to another one
# and outputs 2 lists of the differences found in the WT structure and in the mutant structure
def compareOneVsOneCTdata(dataCtrl, dataTest, mutation):
    diffCtrl = []
    diffTest = []
    offset, startOffset = computeOffset(mutation)
    for nucleotide in dataCtrl.keys():
        if nucleotide >= startOffset:
            if dataCtrl[nucleotide] < startOffset:
                if dataTest[nucleotide+offset] != dataCtrl[nucleotide]:
                    diffCtrl.append(nucleotide)
                    diffTest.append(nucleotide+offset)
            elif dataCtrl[nucleotide] >= startOffset:
                if dataTest[nucleotide+offset] != dataCtrl[nucleotide]+offset:
                    diffCtrl.append(nucleotide)
                    diffTest.append(nucleotide+offset)
            else :
                raise Exception("Should never happen !")
        else:
            if dataCtrl[nucleotide] < startOffset:
                if dataTest[nucleotide] != dataCtrl[nucleotide]:
                    diffCtrl.append(nucleotide)
                    diffTest.append(nucleotide)
            elif dataCtrl[nucleotide] >= startOffset:
                if dataTest[nucleotide] != dataCtrl[nucleotide]+offset:
                    diffCtrl.append(nucleotide)
                    diffTest.append(nucleotide)
            else:
                raise Exception("Should never happen !")
    return diffCtrl, diffTest


# Function that compare all mutant structure to a single WT structure
# and outputs 2 dict of lists of the differences found in the WT structure and in the mutant structures
def compareAllvsOneCTdata(dataCtrl, allDataTest, mutation):
    allDiffCtrl = {}
    allDiffTest = {}
    for structure in allDataTest.keys():
        dataTest = allDataTest[structure]
        diffCtrl, diffTest = compareOneVsOneCTdata(dataCtrl, dataTest, mutation)
        allDiffCtrl[structure] = diffCtrl
        allDiffTest[structure] = diffTest
    return allDiffCtrl, allDiffTest


# Function that write the differences found for one comparison in a tab-delimited file
def writeDifferences(nameFileOutput, diffCtrl, diffTest, energyTest):
    fileOutput = open(nameFileOutput, 'w')
    header = "position_Ctrl" + "\t" + "position_Test" + "\t" + "ENERGY=" + energyTest + "\n"
    fileOutput.write(header)
    i = 0
    for position in diffCtrl:
        lineTowrite = str(position)+"\t"+str(diffTest[i])+"\n"
        fileOutput.write(lineTowrite)
        i += 1
    fileOutput.close()


# Function that write the differences found for all comparisons in a tab-delimited file
def writeAllDifferences(nameFileOutput, diffCtrl, diffTest, energyTest):
    fileOutput = open(nameFileOutput, 'w')
    for structure in diffCtrl.keys():
        header = "position_Ctrl_" + str(structure) + "\t" + "position_Test_" + str(structure) + "\t" + "ENERGY=" +\
                 energyTest[structure] + "\n"
        fileOutput.write(header)
        i = 0
        for position in diffCtrl[structure]:
            lineTowrite = str(position) + "\t" + str(diffTest[structure][i]) + "\n"
            fileOutput.write(lineTowrite)
            i += 1
    fileOutput.close()


# Function that extract the mutation name from the CT file name
def getMutation(fileName, partner="", type="bimolecule"):
    mut = ""
    if type == "bimolecule":
        toSplit = "_"+partner
        mut = fileName.split("mut_")[1].split(toSplit)[0]
    elif type == "monomer":
        mut = fileName.split("mut_")[1].split(".ct")[0]
    else:
        raise Exception("WARNING: invalid type " + type)
    return mut


# Function that launch the whole comparisons between all the CT files from a directory and a control CT file
# and output the number of mutant for which no structure had 0 differences (nbrDiff)
# and the number of mutant for which a structure had 0 differences (nbrNoDiff)
def compareAll(searchDirectory, outputDirectory, fileCtrl, partner="", mode="AllVsOne", type="bimolecule"):
    dataCtrl, energyCtrl = parseOneStructureCTFile(fileCtrl)
    nbrDiff = 0
    nbrNoDiff = 0
    for file in os.listdir(searchDirectory):
        if file.endswith(".ct") and "mut" in file:
            mutation = getMutation(file, partner=partner, type=type)
            fileTest = os.path.join(searchDirectory, file)
            outputFile = str(file.split('.ct')[0]) + ".tab"
            nameFileOutput = os.path.join(outputDirectory, outputFile)
            if mode == "OneVsOne":
                dataTest, energyTest = parseOneStructureCTFile(fileTest)
                diffCtrl, diffTest = compareOneVsOneCTdata(dataCtrl, dataTest, mutation)
                writeDifferences(nameFileOutput, diffCtrl, diffTest, energyTest)
                # increment the number of diff if diffCtrl is not empty
                if diffCtrl:
                    nbrDiff += 1
                else:
                    nbrNoDiff += 1
            elif mode == "AllVsOne":
                dataTest, energyTest = parseAllStructuresCTFile(fileTest)
                diffCtrl, diffTest = compareAllvsOneCTdata(dataCtrl, dataTest, mutation)
                writeAllDifferences(nameFileOutput, diffCtrl, diffTest, energyTest)
            else:
                raise Exception("WARNING: invalid mode. Mode should be equal to OneVsOne or AllVsOne.")
    return nbrDiff, nbrNoDiff


# ------------- Main functions ------------- #
def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Compare all mutated bimolecule structure to the wild type one.")
    parser.add_argument("--inputDir", action="store", dest="inputDirectory", type=str,
                        help="Input directory. This directory should contain all bimolecule prediction files.",
                        required=True)
    parser.add_argument("--structureWT", action="store", dest="CTFileWT", type=str,
                        help="Path and name of the wild type structure file.",
                        required=True)
    parser.add_argument("--outputDir", action="store", dest="outputDirectoryDiff", type=str,
                        help="Path to the output directory (where all diff files will be written).",
                        required=True)
    parser.add_argument("--partnerName", action="store", dest="partnerName", type=str,
                        help="Name of the partner RNA. Example : RNU6atac",
                        required=True)
    parser.add_argument("--chosenMode", action="store", dest="chosenMode", type=str,
                        help="Comparison mode: OneVsOne or AllVsOne.",
                        default="AllVsOne")
    options = parser.parse_args()
    print(">>> begin RNU4atac compare structure")
    # create the output directory if it does not exist
    if not os.path.exists(options.outputDirectoryDiff):
        os.makedirs(options.outputDirectoryDiff)
    nbrDiff, nbrNoDiff = compareAll(options.inputDirectory, options.outputDirectoryDiff, options.CTFileWT,
                                    partner=options.partnerName, mode=options.chosenMode)
    #print(str(nbrDiff)+" variants with differences.")
    #print(str(nbrNoDiff)+" variants with no differences")
    print(">>> end RNU4atac compare structure")


if __name__ == "__main__":
    main()
