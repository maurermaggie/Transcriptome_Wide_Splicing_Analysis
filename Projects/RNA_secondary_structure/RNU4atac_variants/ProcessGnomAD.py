#!/usr/bin/env python

import argparse
import csv
import re


# Function that read the input file
def readGnomAD(gnomADFile):
    gnomADdata = {}
    header = []
    # open file
    with open(gnomADFile, "r") as file:
        csv_reader = csv.reader(file, delimiter=',')
        header = next(csv_reader)
        for element in csv_reader:
            # pass header
            if element[0] == "Chromosome":
                header = element[0:5]+[element[8]]+element[13:]
                continue
            if isInsideRNU4atac(element):
                if element[8] == "":
                    element[8] = "n."+element[1]+element[3]+">"+element[4]
                gnomADdata[element[8]] = element[0:5]+[element[8]]+element[13:]
    return(gnomADdata, header)


# Function that select the variants inside RNU4ATAC gene
def isInsideRNU4atac(variant):
    isInside = False
    # if variant is annotated inside the gene with Ensembl annotation
    if variant[13] == "non_coding_transcript_exon_variant":
        isInside = True
    # if variant is between position 122288456 & 122288586, the variant is inside the gene with RefSeq annotation
    elif int(variant[2]) >= 122288456 and int(variant[2]) <= 122288585:
        # if deletion and it starts after coord 122288585, don't select it
        if len(variant[4]) > 1 and int(variant[2]) < 122288585:
            isInside = True
        elif len(variant[4]) == 1:
            isInside = True
    return(isInside)


# Function that read the file containing the new variants (= variants not present in GnomAD file)
#def readNewVar(newVariantsFile):
#    newVariantsData = {}
#    # open file
#    with open(newVariantsFile, "r") as file:
#        for line in file:
#            line = line.rstrip()
#            element = line.split('\t')
#            newVariantsData[element[5]] = element
#    return(newVariantsData)


# Function that process the GnomADdata
#def processGnomADdata(gnomADdata, newVariantsData):

def processGnomADdata(gnomADdata): 
    processedData = {}
    # first, modify the data that comes from GnomAD
    # for each variant
    for variant in gnomADdata:
        # keep only the columns of interest
        data = gnomADdata[variant]
        # correct the consequence on the consequence
        # get the old coordinates on RNU4ATAC and add 1 to them
        coords = re.findall('\\d+', data[5])
        newConseq = data[5]
        for coord in sorted(coords, reverse=True):
            if int(coord) > 130 and len(data[3]) > 1 and len(data[4]) == 1:
                # deletion outside gnomAD annotation
                newCoord1 = str(int(data[1]) - 122288456 + 2)
                newCoord2 = str(int(newCoord1) + len(data[3]) - 1)
            if int(newCoord2) > 130:
                newCoord2 = str(130)
                lenDel = int(newCoord2) - int(newCoord1) + 1
                newConseq = "n."+newCoord1+"_"+newCoord2+"del"+data[3][1:lenDel]
            elif int(coord) > 130 and len(data[3]) == 1 and len(data[4]) > 1:
                # insertion outside gnomAD annotation
                newCoord1 = str(int(data[1]) - 122288456 + 1)
                newCoord2 = str(int(newCoord1) + 1)
                newConseq = "n."+newCoord1+"_"+newCoord2+"ins"+data[4][1:]
            elif int(coord) > 130:
                # substitution outside gnomAD annotation
                newCoord = str(int(data[1]) - 122288456 + 1)
                newConseq = newConseq.replace(coord, newCoord, 1)     
            else:
                # inside gnomAD annotation
                newCoord = str(int(coord)+1)
                newConseq = newConseq.replace(coord, newCoord, 1)
        data[5] = newConseq
        processedData[data[5]] = data
    return(processedData)
    # then, merge the data with the new variants
    #for variant in newVariantsData:
    #    data = newVariantsData[variant] + ([0] * 37)
    #    processedData[variant] = data

# Function that write the output file
def writeOutput(data, header, outputFile):
    # open the file
    with open(outputFile, "w") as file:
        # write header
        toWrite = "\t".join(str(i) for i in header)
        toWrite += "\n"
        file.write(toWrite)
        # for each variant, write data
        for variant in data:
            toWrite = "\t".join(str(i) for i in data[variant])
            toWrite += "\n"
            file.write(toWrite)
    return(1)


# ------------- Main functions ------------- #
def main():
    # Parse arguments
    parser = argparse.ArgumentParser(description="Process GnomAD file to convert it to the format needed by"
                                                 " the rest of the pipeline")
    parser.add_argument("--input", action="store", dest="input", type=str,
                        help="Name of the input csv file from GnomAD.",
                        required=True)
    parser.add_argument("--output", action="store", dest="output", type=str,
                        help="Name of the output tab file.",
                        required=True)
    #parser.add_argument("--newVariants", action="store", dest="newVariants", type=str,
                        #help="Name of file containing the new variants.",
                        #required=True)
    parser.add_argument("--dir", action="store", dest="dir", type=str,
                        help="Name of the directory where the input is and where the output will be created.",
                        required=True)
    options = parser.parse_args()
    print(">>> begin processing GnomAD file")
    print("> read GnomAD file")
    inputFile = options.dir+options.input
    (gnomADdata, header) = readGnomAD(inputFile)
    print(str(len(gnomADdata)) + " variants")
    print("> read new variants file")
    #newVariantsFile = options.dir + options.newVariants
    #newVariantsData = readNewVar(newVariantsFile)
    print("> process GnomAD data and merge with new variants")
    #processedData = processGnomADdata(gnomADdata, newVariantsData)
    
    processedData = processGnomADdata(gnomADdata)
    print ("> write output file")
    outputFile = options.dir + options.output
    writeOutput(processedData, header, outputFile)
    print(">>> end processing GnomAD file")


if __name__ == "__main__":
    main()
