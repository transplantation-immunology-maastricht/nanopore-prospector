#import sys
from os import makedirs
from os.path import isfile, isdir, exists
from sys import argv, exc_info
from getopt import getopt, GetoptError


from collect_reads.collect_reads import collectReads

#SoftwareVersion = "nit-picker Version 1.0"



# Read Commandline Arguments.  Return true if everything looks okay for read extraction.
def readArgs():
    
    # Default to None.  So I can easily check if they were not passed in.
   
    global fastqReadFile
    global fast5ReadDirectory
    global outputDirectory
            
    fastqReadFile      = None
    fast5ReadDirectory = None
    outputDirectory    = None

    # https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    try:
        opts, args = getopt(argv[1:]
            ,"f:F:o:"
            ,["fastq=", "fast5=", "output="])

        for opt, arg in opts:

            if opt in ("-o", "--output"):
                outputDirectory = arg

            elif opt in ("-f", "--fastq"):
                fastqReadFile = arg
            elif opt in ("-F", "--fast5"):
                fast5ReadDirectory = arg
   
            else:
                print('Unknown Commandline Option:' + str(opt) + ':' + str(arg))
                raise Exception('Unknown Commandline Option:' + str(opt) + ':' + str(arg))
            
        if(len(argv) < 3):
            print ('I don\'t think you have enough arguments.\n')
            #usage()
            return False     

    except GetoptError:        
        print ('Something seems wrong with your commandline parameters.')
        print(exc_info())
        #usage()
        return False

    # Sanity Checks. The only required parameters are input dir/file and output.    
    if (isfile(fastqReadFile)):
        print ('Fastq input is a file that exists.')
    elif (isdir(fast5ReadDirectory)):
        print ('Fast5 Read input is a directory that exists.')
    else :
        print ('I don\'t understand the read input specified, it is not a file or directory:' + readInput)
        return False
    
    # This output directory should exist
    if not exists(outputDirectory):
        makedirs(outputDirectory)


    return True

if __name__=='__main__':
    
    try:    
        if(readArgs()):
            print('Commandline arguments look fine.\n I will collect the reads that are in the fastq file.')
            
            collectReads(fastqReadFile, fast5ReadDirectory, outputDirectory)
            #prepareReads(readInput, outputResultDirectory, sampleID, barcodeFileLocation, referenceInput, minimumReadLength, maximumReadLength, minimumQuality, maximumQuality )
            
            print ('Done collecting the reads. Have a nice day.')    
        else:
            print('\nI\'m giving up because I was not satisfied with your commandline arguments.')  
            
    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print ('Fatal problem during read extraction:')
        print (exc_info())
        raise
   