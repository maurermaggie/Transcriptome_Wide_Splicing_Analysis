#!/usr/bin/env python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re


# Function that read the fasta file and give the sequence and its id
def readSeq(file):
    fastaRecord = SeqIO.read(file, "fasta")
    return fastaRecord.id, fastaRecord.seq


# Function that write a fasta file for a given sequence
def writeSeq(fileName, sequence, idSeq, mut):
    sequenceFasta = SeqRecord(sequence, id=idSeq, description=mut)
    SeqIO.write(sequenceFasta, fileName, "fasta")
    return 1


# Function that get into input the annotation of the variant
# and output the action to apply on the sequence
def processMutation(mut):
    mut = mut[2:]
    action = ""
    pos = 0
    sequence = ""
    # if mut contain ">" : action = substitution
    if '>' in mut:
        action = "substitution"
        match = re.search('^\d+', mut)
        if match:
            pos = match.group(0)
        else:
            print("WARNING : no position for this substitution (" + mut + ")")
        matchSeq = re.search('(?<=>)[ATGC]$', mut)
        if matchSeq:
            sequence = matchSeq.group(0)
        else:
            print("WARNING : no alternative sequence for this substitution (" + mut + ")")
    # elsif mut contain "del" : action = deletion
    elif 'del' in mut:
        action = "deletion"
        match = re.search('^\d+', mut)
        posStart = 0
        if match:
            posStart = match.group(0)
        else:
            print("WARNING : no position for this deletion (" + mut + ")")
        matchPosEnd = re.search('(?<=_)\d+', mut)
        if matchPosEnd:
            posEnd = matchPosEnd.group(0)
            pos = [posStart, posEnd]
        else:
            pos = [posStart]
        matchSeq = re.search('(?<=del)[ATGC]+$', mut)
        if matchSeq:
            sequence = matchSeq.group(0)
        else:
            print("WARNING : no sequence for this deletion (" + mut + ")")
    # elsif mut contain "ins" : action = insertion
    elif 'ins' in mut:
        action = "insertion"
        match = re.search('^\d+', mut)
        posStart = 0
        posEnd = 0
        if match:
            posStart = match.group(0)
        else:
            print("WARNING : no position for this insertion (" + mut + ")")
        matchPosEnd = re.search('(?<=_)\d+', mut)
        if matchPosEnd:
            posEnd = matchPosEnd.group(0)
        else:
            print("WARNING : no end position for this insertion (" + mut + ")")
        pos = (posStart, posEnd)
        matchSeq = re.search('(?<=ins)[ATGC]+$', mut)
        if matchSeq:
            sequence = matchSeq.group(0)
        else:
            print("WARNING : no sequence for this insertion (" + mut + ")")
    # elsif mut contain "dup" : action = duplication
    elif 'dup' in mut:
        action = "duplication"
        match = re.search('^\d+', mut)
        posStart = 0
        if match:
            posStart = match.group(0)
        else:
            print("WARNING : no position for this duplication (" + mut + ")")
        matchPosEnd = re.search('(?<=_)\d+', mut)
        if matchPosEnd:
            posEnd = matchPosEnd.group(0)
            pos = [posStart, posEnd]
        else:
            pos =  [posStart]
        matchSeq = re.search('(?<=dup)[ATGC]+$', mut)
        if matchSeq:
            sequence = matchSeq.group(0)
        else:
            print("WARNING : no sequence for this duplication (" + mut + ")")
    return pos, action, sequence


# Function that takes a sequence and a variation (mutation, indel, ins)
# and apply this variation to the sequence
def applyMutation(seqInit, mut):
    (pos, mutType, sequenceAlt) = processMutation(mut)
    seqMut = ""
    if mutType == "substitution":
        seqMut = seqInit[:(int(pos)-1)]+sequenceAlt+seqInit[int(pos):]
    elif mutType == "deletion":
        if len(pos) == 2:
            seqMut = seqInit[:int(pos[0])]+seqInit[int(pos[1])+1:]
        elif len(pos) == 1:
            lenSeqDel = len(sequenceAlt)
            if lenSeqDel == 1:
                seqMut = seqInit[:int(pos[0])]+seqInit[int(pos[0])+1:]
            elif lenSeqDel > 1 :
                seqMut = seqInit[:int(pos[0])]
        else:
            print("WARNING : problem with deletion, more than 2 positions...")
    elif mutType == "insertion":
        if len(pos) == 2:
            seqMut = seqInit[:int(pos[0])]+sequenceAlt+seqInit[int(pos[1])-1:]
        else:
            print("WARNING : missing start or end for the insertion OR too many positions for this insertion...")
    elif mutType == "duplication":
        if len(pos) == 2:
            seqMut = seqInit[:int(pos[1])]+sequenceAlt+seqInit[int(pos[1]):]
        elif len(pos) == 1:
            seqMut = seqInit[:int(pos[0])]+sequenceAlt+seqInit[int(pos[0]):]
        else:
            print("WARNING : problem with duplication, more than 2 positions...")
    return seqMut


# Function that write the mutated sequence
def mutateSequence(idSeq, seqWT, fastaOutput, mut):
    seqMut = applyMutation(seqWT, mut)
    writeSeq(fastaOutput, seqMut, idSeq, mut)
    return 1
