# This file is part of HLA-allele-analysis.
#
# HLA-allele-analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HLA-allele-analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with HLA-allele-analysis. If not, see <http://www.gnu.org/licenses/>.

#import xml.etree.ElementTree
import sys
import os
import getopt
#from os import listdir
from os.path import split
#from HLA_Allele import HLA_Allele
#from HLA_Allele_Group import HLA_Allele_Group
#from Bio.Align.Applications import ClustalwCommandline
#from Bio.Align.Applications import MuscleCommandline
#from Bio import AlignIO
from Bio import SeqIO
#from Bio.Align import AlignInfo
#import traceback
#import subprocess

SoftwareVersion = "Polymorphism Analyzer Version 1.0"

def usage():
    print("You are using this program incorrectly.  Program better Noob.")    

# This method is a directory-safe way to open up a write file.
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not os.path.isdir(tempDir):
        os.mkdir(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput

# Read Commandline Arguments.  Return true if everything looks okay for read extraction.
def readArgs():
    # Default to None.  So I can easily check if they were not passed in.
    global inputFileName
    global outputFileName

    
    inputFileName  = None
    outputFileName = None

    if(len(sys.argv) < 3):
        print ('I don\'t think you have enough arguments.\n')
        usage()
        raise()
        return False    

    # getopt.getopt(..) is a function for parsing args the way smart people do it
    # For More info, use google or 
    # https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    try:
        opts, args = getopt.getopt(sys.argv[1:]
            ,"hvi:o:"
            ,["help", "version", "input=","output="])

        for opt, arg in opts:

            if opt in ('-h', '--help'):
                print (SoftwareVersion)
                usage()
                return False

            elif opt in ('-v', '--version'):
                print (SoftwareVersion)
                return False

            elif opt in ("-i", "--input"):
                inputFileName = arg
            elif opt in ("-o", "--output"):
                outputFileName = arg
                
            else:
                print('Unknown commandline option: ' + opt)
                raise()

    except getopt.GetoptError, errorMessage:
        print ('Something seems wrong with your commandline parameters.')
        print (errorMessage)
        usage()
        return False

    print('Input File:' + str(inputFileName))
    print('Output File:' + str(outputFileName))

    # Quick sanity check.
    if(len(inputFileName) < 4):
        print('Input directory is too short:' + str(inputFileName))
        return False
    if(len(outputFileName) < 4):
        print('Output directory is too short:' + str(outputFileName))
        return False

    return True

# Look at the variation in the alleles.  
def analyzePolymorphism():
    print ('Parsing input file:' + inputFileName)
    
    outputFile = createOutputFile(outputFileName)    
    outputFile.write('Input file,' + inputFileName + '\n')
    
    numberSequences = len(list(enumerate(SeqIO.parse(inputFileName, 'fasta'))))
    print('I found this many sequences to analyze:' + str(numberSequences))
    print('I\'ll get started.')
    
    outputFile.write('Input Sequence Count,' + str(numberSequences) + '\n')
    
    
    matchingAlleles = {}
    
    enumeratedRecords = enumerate(SeqIO.parse(inputFileName, 'fasta'))
    # Each record represents an HLA element in the input fasta file.
    for index, record in enumeratedRecords:
        sequenceID = str(record.id)
        #print ('Read ID:' + currentReadID)
        sequenceSequence = str(record.seq)
        
        #print('ID:' + sequenceID)
        #print('SQ:' + sequenceSequence)
        
        
        #Is this sequence familiar to us already?
        if(sequenceSequence in matchingAlleles):
            #Yes.  Add this sequence to the list.
            matchingAlleles[sequenceSequence].append(sequenceID)
        else:
            # No.  Create a dictionary entry for this sequence, add the sequence to the list.
            matchingAlleles[sequenceSequence] = [sequenceID]

    outputFile.write('Distinct Sequence Count,' + str(len(matchingAlleles)) + '\n')
    
    # Loop through my matching allele dictionary.
    outputFile.write('\n')
    for currentKey in matchingAlleles:
        print ('\n\nSequence:\n' + currentKey)
        currentSequenceIDs = sorted(matchingAlleles[currentKey])
        outputFile.write('\n' + str(len(currentSequenceIDs)) + ',,' + currentKey + '\n')

        print (str(len(currentSequenceIDs)) + ' sequences found.')
        
        for currentSequenceID in currentSequenceIDs:
            print(currentSequenceID)
            outputFile.write(currentSequenceID + '\n')
            
    outputFile.close()
        


if __name__=='__main__':
    try:
        readArgs()
        
        analyzePolymorphism()


        print('Done.  Ben did a great job.')

    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print 'Unexpected problem during execution:'
        print sys.exc_info()[1]
        raise
