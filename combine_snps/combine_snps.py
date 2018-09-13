from os.path import splitext, isdir, split, join
from os import makedirs, system, listdir

#from csv import reader

#from Bio.SeqIO import parse, write

#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord

#from pysam import AlignmentFile

#from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA



class SnpArrayStructure:

    
    def __init__(self):
        self.sequenceName = None
        
        
        self.snpPositions = []
        self.refBases = []
        self.snpBases = []

def combineSNPLists(inputDirectory):  
    # This method extracts a specific sequence from a larger amplicon.
    # For example, we are extracting the HLA-DRA coding sequence from our DRA amplicon.

    mySNPArrays=[]

    print('I am inside extract sequences method.')
    
    for fileNameShort in sorted(listdir(inputDirectory)):
        if fileNameShort.endswith(".csv"):
            print('Found CSV File:' + str(fileNameShort))
            fullPath = join(inputDirectory, fileNameShort)
            
            currentSNPArray = SnpArrayStructure()
            currentSNPArray.sequenceName = fileNameShort.replace('.csv','')

            csvLines = open(fullPath,'rb').read().splitlines()
            
            row1Tokens = str(csvLines[0]).split(',')
            row2Tokens = str(csvLines[1]).split(',')
            row3Tokens = str(csvLines[2]).split(',')

            for index, token in enumerate(row1Tokens):
                # Exclude the ffirst and last, there are no snp data there.
                if(index==0 or index==(len(row1Tokens)-1)):
                    pass
                else:
                    currentSNPArray.snpPositions.append(int(row1Tokens[index]))
                    currentSNPArray.refBases.append(str(row2Tokens[index]))
                    currentSNPArray.snpBases.append(str(row3Tokens[index]))
                    
                #print(token)

            #print('printing the stored snps:')
            #for index, snpPos in enumerate(currentSNPArray.snpPositions):
            #    print(currentSNPArray.snpPositions[index] + '\n')
            #    print(currentSNPArray.refBases[index] + '\n')
            #    print(currentSNPArray.snpBases[index] + '\n')
                
            # Store our new snp array
            mySNPArrays.append(currentSNPArray)


    return mySNPArrays


def printCombinedSnpLists(mySNPArrays, outputDirectory):
    print ('Printing combined snp lists')
    print ('I have ' + str(len(mySNPArrays)) + ' snp arrays.')
    
    outputFileName = join(outputDirectory, 'CombinedSnps.csv')
    outputFile = createOutputFile(outputFileName)
    
    # Get list of all snp positions
    # The key is the reference position
    # The value is the reference base.
    snpDict = {}
    for snpArray in mySNPArrays:
        for index, snpPos in enumerate(snpArray.snpPositions):
            #print ('appending this snp position:' + str(snpPos))
            snpDict[snpPos] = snpArray.refBases[index]
    # convert to set to make distinct
    snpList = sorted(snpDict.keys())
    
    # Burn a column
    outputFile.write(',')
    
    # Write position and reference
    for snpPos in snpList:
        outputFile.write(str(snpPos) + ',')
    outputFile.write('\n')
    
    outputFile.write(',')
    for snpPos in snpList:
        outputFile.write(str(snpDict[snpPos]) + ',')
    outputFile.write('\n')
    
    # for each Sequence
    for snpArray in mySNPArrays:
        
        outputFile.write(str(snpArray.sequenceName) + ',')

        # For each snp position
        for snpPos in snpList:            
            # if sequence has this snp
            if(snpPos in snpArray.snpPositions):
                # print snp
                outputFile.write(str(snpArray.snpBases[snpArray.snpPositions.index(snpPos)]) + ',')
            else:
                outputFile.write(',')
        
        # newline
        outputFile.write('\n')
    
    outputFile.close()
    
    
    

    
    
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not isdir(tempDir):
        print('Making Directory:' + tempDir)
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput     

