from sys import argv, exc_info
from getopt import getopt, GetoptError

from os.path import isdir, isfile

#from combine_snps.combine_snps import combineSNPLists, printCombinedSnpLists
from combine_snps.find_files_with_pattern import combineFastaFiles, collectFilesWithPattern

# TODO: I think this whole file is broken, probably needs some maintenance in order to use it.
def usage():
    print("usage:\n" + 
    "\tThis script is written for python 3.6\n" + 
    "\tI haven't written the usage tutorial yet.  Oops.  Do this now please."
    )      

# Read Commandline Arguments.  Return true if everything looks okay for read extraction.
def readArgs():
    # Default to None.  So I can easily check if they were not passed in.
  
    global input
    global output
    
    input        = None
    output        = None

    # TODO: splitHeterozygotes should be commandline parameter, not just a True
    #splitHeterozygotes      = True

    if(len(argv) < 2):
        print ('I don\'t think you have enough arguments.\n')
        usage()
        return False    

    # https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    try:
        opts, args = getopt(argv[1:]
            ,"i:o:"
            ,['input=','output='])

        for opt, arg in opts:

            if opt in ('-h', '--help'):
                #print (SoftwareVersion)
                usage()
                return False

            elif opt in ('-v', '--version'):
                #print (SoftwareVersion)
                return False

            elif opt in ("-i", "--input"):
                input = arg
            elif opt in ("-o", "--output"):
                output = arg

            else:
                print('Unknown Commandline Option:' + str(opt) + ':' + str(arg))
                raise Exception('Unknown Commandline Option:' + str(opt) + ':' + str(arg))
            

    except GetoptError:
        print ('Something seems wrong with your commandline parameters.')
        print (exc_info())
        usage()
        return False

    # Consensus,threads is optional, the rest are necessary.
    # Sanity Checks
   
    
    #if (isfile(sequenceInput)):
    #    print ('Read input is a file that exists.')
    #elif (isdir(readInput)):
    #    print ('Read input is a directory that exists.')
    #else :
    #    print ('I don\'t understand the read input specified, it is not a file or directory:' + readInput)
    #    return False

    return True

if __name__=='__main__':

    try:    
        if(readArgs()):
            print('Commandline arguments look fine.\nThe hour is at hand.')
 
            #extractGeneSequencesFromLongerSequences(sequenceInput, geneReference)
            #snpLists = combineSNPLists(input)
            #printCombinedSnpLists(snpLists,output)

            #TODO: Pass in the file description parameter as a string, and move this to nanopore prospector.
            #combineFastaFiles(input, output, 'NewConsensus.fasta')

            collectFilesWithPattern(input, output, 'AlignmentSummary.SNPs.csv')

            
            print ('I am done for now, have a nice day.')    
        else:
            print('\nI\'m giving up because I was not satisfied with your commandline arguments.')  
            
    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print ('Fatal problem during sequence extraction:')
        print (exc_info())
        raise


