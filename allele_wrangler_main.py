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



#SoftwareVersion = "Allele-Wrangler Version 1.0"

#import os

#import sys
#import getopt


# TODO: I think this whole file is broken, probably needs some maintenance in order to use it.
def usage():
    print("usage:\n" + 
    "\tThis script is written for python 3.6\n" + 
    "\tI haven't written the usage tutorial yet.  Oops.  Do this now please."
    
    )      
    
    
# Read Commandline Arguments.  Return true if everything looks okay for read extraction.
def readArgs():
    # Default to None.  So I can easily check if they were not passed in.
   
    global readInput
    global outputResultDirectory
    global numberIterations
    global consensusFileName
    global numberThreads
    global splitHeterozygotes
    
    readInput        = None
    outputResultDirectory    = None
    consensusFileName        = None
    numberIterations         = 1
    numberThreads            = 1
    # TODO: splitHeterozygotes should be commandline parameter, not just a True
    splitHeterozygotes      = True

    if(len(sys.argv) < 3):
        print ('I don\'t think you have enough arguments.\n')
        usage()
        return False    

    # https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    try:
        opts, args = getopt.getopt(sys.argv[1:]
            ,"hvi:o:r:c:t:"
            ,[ "help", "version", "iterations=","outputdir=","reads=",'consensus=','threads='])

        for opt, arg in opts:

            if opt in ('-h', '--help'):
                print (SoftwareVersion)
                usage()
                return False

            elif opt in ('-v', '--version'):
                print (SoftwareVersion)
                return False
            
            # TODO Commandline args:
            # param:MultipleSequenceAligner Readcount=3
            # Tune it higher for better initial consensus. Maybe longer final consensus?
            # param: Aligment sample readcount=75
            # param: splitheterozygtes

            elif opt in ("-i", "--iterations"):
                numberIterations = arg
            elif opt in ("-o", "--outputdir"):
                outputResultDirectory = arg
            elif opt in ("-r", "--reads"):
                readInput = arg
            elif opt in ("-c", "--consensus"):
                consensusFileName = arg
            elif opt in ("-c", "--threads"):
                numberThreads = arg
            else:
                print('Unknown Commandline Option:' + str(opt) + ':' + str(arg))
                raise Exception('Unknown Commandline Option:' + str(opt) + ':' + str(arg))
            

    except getopt.GetoptError, errorMessage:
        print ('Something seems wrong with your commandline parameters.')
        print (errorMessage)
        usage()
        return False

    # Consensus,threads is optional, the rest are necessary.
    # Sanity Checks
    if (int(numberIterations) < 1):
        print('You must specify an Iteration count >= 1.')
    
    if (os.path.isfile(readInput)):
        print ('Read input is a file that exists.')
    elif (os.path.isdir(readInput)):
        print ('Read input is a directory that exists.')
    else :
        print ('I don\'t understand the read input specified, it is not a file or directory:' + readInput)
        return False

    return True

if __name__=='__main__':

    try:    
        if(readArgs()):
            print('Commandline arguments look fine.\nThe hour is at hand. Let us wrangle the Alleles.')
            
            myAlleleWrangler = AlleleWrangler(readInput, outputResultDirectory, consensusFileName, numberIterations, numberThreads, splitHeterozygotes)
            myAlleleWrangler.analyzeReads()
            
            print ('I am done wrangling alleles for now, have a nice day.')    
        else:
            print('\nI\'m giving up because I was not satisfied with your commandline arguments.')  
            
    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print 'Fatal problem during read extraction:'
        print sys.exc_info()
        raise


