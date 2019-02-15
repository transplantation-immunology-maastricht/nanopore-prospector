from os.path import splitext, isdir, split, join
from os import makedirs, system, listdir

from Bio.Seq import Seq
from Bio.SeqIO import write, parse

import csv

def consensusSequenceFromSNPFiles(referenceFasta, snpTable, outputDirectory):
    # The reference sequence is a fasta file with just the unmodified reference.
    referenceFastaSequence = list(parse(referenceFasta, 'fasta'))[0]

    # I want to be sure I am changing the correct base, so I include the reference base in the input.
    # snpTable has this format:
    # [ X ,refpos1 ,refpos2 ...]
    # [ X ,refbase1,refbase2...]
    # [nm1,snpbase1,snpbase2...]
    # [nm2,snpbase1,snpbase2...]

    datafile = open(snpTable, 'r')
    datareader = csv.reader(datafile, delimiter=',')
    allSnpData = []
    for row in datareader:
        allSnpData.append(row)

    #print (allSnpData[1:4])

    # Assign Reference SNP Positions. Ignore the first column, it's just blank.
    # I am reversing the lists. Because I want to iterate them in reverse order.
    referenceSNPPositions = allSnpData[0][1:]
    referenceSNPPositions.reverse()
    #print(referenceSNPPositions)
    
    # Assign Reference SNP Values
    # I am reversing the lists. Because I want to iterate them in reverse order.
    referenceSNPValues = allSnpData[1][1:]
    referenceSNPValues.reverse()
    #print(referenceSNPValues)

    # Assign novel SNP data to new list, get rid of reference snp data.
    novelSNPData = allSnpData[2:]

    outputFastaFile = createOutputFile(join(outputDirectory,'GeneratedSequences.fasta'))

    # Loop Sequences. Don't include the reference sequence.
    for sequenceIndex in range(0,len(novelSNPData)):

        sequenceName = novelSNPData[sequenceIndex][0]
        # SNP list for novel sequence. Reversed to iterate in reverse order.
        novelSequenceSNPList = novelSNPData[sequenceIndex][1:]
        novelSequenceSNPList.reverse()

        print('Checking Novel Sequence ' +  str(sequenceName))

        # Create a new sequence, start with the reference sequence.
        newSequence = referenceFastaSequence.seq

        # Loop SNP Positions, In Reverse. Reverse is important because forward loops might mess up indices.
        for index, snpPosition in enumerate(referenceSNPPositions):
            # Get Reference Sequence (and length) from SNP Table
            referenceSNPSequence = referenceSNPValues[index]
            refSequenceLength = len(referenceSNPSequence)

            # Compare with sequence extracted from Reference sequence, this is a sanity check.
            # Subtract 1 because lists are 0-based, and snps are 1-based
            extractedReferenceSNPSequence = newSequence[int(snpPosition)-1:int(snpPosition) + refSequenceLength -1]

            # Get SNP sequence
            novelSNPSequence = novelSequenceSNPList[index]

            #print('checking SNP position:' + str(snpPosition))

            # Final check to make sure we're changing the right sequence.
            if(referenceSNPSequence == extractedReferenceSNPSequence):

                # Replace with the SNP sequence.
                if(referenceSNPSequence != novelSNPSequence):
                    # = s[:index] + newstring + s[index + 1:]
                    print('Position ' +  str(snpPosition) + ' Replacing ' + referenceSNPSequence + ' with ' + novelSNPSequence)
                    newSequence = newSequence[:int(snpPosition)-1] + novelSNPSequence + newSequence[int(snpPosition) + refSequenceLength -1:]
                else:
                    #print('Snp Matches reference, not replacing.')
                    pass

            else:
                print('PROBLEM DETECTED:')
                print ('Position: ' + str(snpPosition)
                    + ', RefSNPSequence: ' + str(referenceSNPSequence)
                    + ', ExtractedReferenceSNPSequence: ' + str(extractedReferenceSNPSequence)
                    + ', NovelSNPSequence: ' + str(novelSNPSequence))

                raise ValueError('Reference sequence does not match extracted reference sequence')

        # Removing the deletion character for the final sequence....:
        newSequence = str(newSequence).replace('-','')
        #print('SEQUENCE:\n' + newSequence)

        # Add novel sequence to Fasta output.
        # This is a super cheap way to make a fasta, i could use biopython but it's harder.

        outputFastaFile.write('>' + str(sequenceName) + '\n')
        outputFastaFile.write(str(newSequence) + '\n')

    outputFastaFile.close()
    
    
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not isdir(tempDir):
        print('Making Directory:' + tempDir)
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput     

