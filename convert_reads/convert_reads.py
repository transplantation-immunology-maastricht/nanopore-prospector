from sys import exc_info
#from Bio.SeqIO import parse as parseReads, write
from Bio.SeqIO import parse as parseReads#, FastaIO
from os.path import split, isdir
from os import makedirs


def fastqToFasta(inputFile, outputFile, truncateIDLength):

    try:
        parsedReads = list(parseReads(inputFile, 'fastq'))

        # Create output file for fasta files, with a wrap=None.
        tempDir, tempFilename = split(outputFile)
        if not isdir(tempDir):
            makedirs(tempDir)
        fastaOut = open(outputFile, 'w')
        #fastaOut = FastaIO.FastaWriter(resultsOutput, wrap=None)
        
        for read in parsedReads:
            #read.id = read.id[0:truncateIDLength]
            #read.description = ''
            #fastaOut.write_record(read)
            fastaOut.write('>' + str(read.id) + ' ' + str(read.description) + '\r\n')
            fastaOut.write(str(read.seq) + '\r\n')



    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print ('Fatal problem during read extraction:')
        print (exc_info())
        raise
   
