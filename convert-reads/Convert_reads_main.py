#import sys
#from os import makedirs
##from os.path import isfile, isdir, exists
from sys import argv, exc_info
#from getopt import getopt, GetoptError


#from collect_reads.collect_reads import collectReads

from Bio.SeqIO import parse as parseReads, write

#SoftwareVersion = "nit-picker Version 1.0"



if __name__=='__main__':
    
    try:    
        # convert fastq to fasta.
        
        inputFile = argv[1]
        outputFile = argv[2]

        #inputFile = '/home/ben/ben_share/nanopolish_test/2018.Jan10.ValidationData/A1/A1Reads.fastq'
        #outputFile = '/home/ben/ben_share/nanopolish_test/2018.Jan10.ValidationData/A1/A1Reads.fasta'
        
        
                    
        
        parsedReads = list(parseReads(inputFile, 'fastq'))
        
        for read in parsedReads:
            #read.id = 'Hello'
            read.id = read.id[0:36]
            read.description = ''
        
        write(parsedReads, outputFile, 'fasta')
    #minionReadRecords = enumerate(parsedReads)
    #readRecordList = list(minionReadRecords)
    #readRecordList = list(parsedReads)

    #readCount = len(readRecordList)
        
        
    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print ('Fatal problem during read extraction:')
        print (exc_info())
        raise
   
