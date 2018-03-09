# This file is part of nanopore_prospector.
#
# nanopore_prospector is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# nanopore_prospector is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with nanopore_prospector. If not, see <http://www.gnu.org/licenses/>.


#import sys
from os import makedirs
#from os.path import isfile, isdir, exists
from os.path import exists
from sys import argv, exc_info
from getopt import getopt, GetoptError


from allele_wrangler.allele_wrangler import AlleleWrangler





# Read Commandline Arguments.  Return true if everything looks okay for read extraction.
def readArgs():
    
    # Default to None.  So I can easily check if they were not passed in.
   
    global inputFastaSequences
    global predictedClusterCount
    global outputDirectory
    
    inputFastaSequences   = None
    predictedClusterCount = None
    outputDirectory       = None
            
    try:
        opts, args = getopt(argv[1:]
            ,"c:i:o:"
            ,["clusters=", "inputfile=", "outputdir="])

        for opt, arg in opts:

            if opt in ("-o", "--outputdir"):
                outputDirectory = arg

            elif opt in ("-i", "--inputfile"):
                inputFastaSequences = arg
                
            elif opt in ("-c", "--clusters"):
                predictedClusterCount = arg
   
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
    #if (isfile(fastqReadFile)):
    #    print ('Fastq input is a file that exists.')
    #elif (isdir(fast5ReadDirectory)):
    #    print ('Fast5 Read input is a directory that exists.')
    #else :
    #    print ('I don\'t understand the read input specified, it is not a file or directory:' + readInput)
    #    return False
    
    # are my input or output parameters equal to None?
    
    # This output directory should exist
    if not exists(outputDirectory):
        makedirs(outputDirectory)


    return True

if __name__=='__main__':
    
    try:    
        if(readArgs()):
            print('Commandline arguments look fine.\n Time to cluster the sequences.')
            # create an allele_wrangler object.
            
            # Do I make a consensus sequence? 
            # Or should I align my reads straight against the first sequence?
            # Probably just align the consensuses to the first sequnce
            # Name this alignment like it's a "final" alignment.
            
            # Cluster reads, be able to pass in a parameter with how many clusters.

            # TODO Parameters
            numberIterations = 6
            numberThreads = 4
            splitHeterozygotes = True
            
            myAlleleWrangler = AlleleWrangler(inputFastaSequences, outputDirectory, None, numberIterations, numberThreads, splitHeterozygotes)
            currentAssemblyResults = myAlleleWrangler.analyzeReads()
           
            print ('Done. Have a nice day.')    
        else:
            print('\nI\'m giving up because I was not satisfied with your commandline arguments.')  
            
    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print ('Fatal problem during sequence clustering:')
        print (exc_info())
        raise
   