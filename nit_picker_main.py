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

# I'm deprecating this file. This logic is just to accept variables, it is repetition of the Nanopore_Prospector_main.py file.
# So from now on don't use this main method, start with the NanoporeProspectorMain main method.

# TODO Nitpicker: Lots of things:
#nitpicker: put a progress bar on the read extractor and demultiplexer.
#Nitpicker: Feedback from Albacore Q scores: 
#Description
#If you take the qscore for each base, convert it back to an error probability, take the mean of those, and then convert that mean error back into a qscore, you should get something very close to what albacore provides.
#Prospector - do some better statistics on read alignments.  I can do a few things:
#Count statistics for homopolymers.
#TBH this fits in with the nitpicker code. I think nitpicker should have a parameter for just "read analysis".  If that flag is set, we must have 1) reference sequence. 1) Read sequences. Reads are aligned. How many could not align?  For each read: for each position, match or mismatch. Calculate perentage.  There must be papers about analyzing read quality, what are the metrics?
#For each read: for each homopolymer: how many are in read vs sequence? Statistics.
#- Polish consensus sequences using known informatin about the most likely homopolymers.
#-Paul suggested that during nit-picker i can calculate read statistics.  homopolyemers.  nit picker calculates some stats on reads. Get a good idea about reference length vs homopolymer length. I can implement these parameters into allele-wrangler.
#Trim off the barcodes. in the nit picker
# TODO: trim out spaces and special characters from the sample id.
# Check on all the values.

import sys
from os import makedirs
from os.path import isfile, isdir, exists
from getopt import getopt, GetoptError

from nit_picker.nit_picker import prepareReads

SoftwareVersion = "nit-picker Version 1.0"

# TODO Fix usage
def usage():
    print("usage:\n" + 
    "\tThis script is written for python 3.6\n" + 
    "\tInput = fastq of reads, or a directory. Directory might contain multiple fastq files,\n" + 
    "\t\tOR\n" + 
    "\t\tSubdirectories containing .fast5 reads (Barcoded Reads, ex. /BCO1)\n\n" + 
    "\tThe output directory will be filled with .fasta and .fastq reads, sorted by subfolders.\n\n" + 

    "\tOptions:\n" +  
    
    "\t-r\t--reads    \tInput Reads. May be a fastq file, or directory (required)\n" +  
    "\t-o\t--outputdir\tOutput Directory (required)\n" +  
    
    "\t-R\t--reference\tFasta file containing HLA reference sequence(s). Used for calculating actual read quality statistics.\n" +
    # If I use a reference file, this program will assume that my reads all will map to the reference file.
    # All reads should be on-target, more or less.
    # Heterozygous positions will confuse the analysis. 
    
    # TODO: accept "calculateReadStats" as a commandline option. Just passing a reference is not good enough.
    
    "\t-m\t--minlen   \tMinimum Read Length Filter\n" +  
    "\t-M\t--maxlen   \tMaximum Read Length Filter\n" + 
    
    "\t-q\t--minqual  \tMinimum Read Quality (Phred)\n" +  
    "\t-Q\t--maxqual  \tMaximum Read Quality (Phred)\n" +  
           
    "\t-h\t--help     \tPrint this message\n" +   
    "\t-r\t--rundate  \tOptional identifier, provided by you, to be included in filename.\n" +   
    
    "\n\tSee README.MD for instructions on how to set up an anaconda environment for this script\n"
    )   
    # Usage: you must provide a read input and output directory. 
    # reads are fastq or a directory with fastq.
    # The rest are optional.

# Read Commandline Arguments.  Return true if everything looks okay for read extraction.
def readArgs():
    
    # Default to None.  So I can easily check if they were not passed in.
   
    global readInput
    global referenceInput
    global outputResultDirectory
    global barcodeFileLocation    
    global minimumReadLength
    global maximumReadLength
    global minimumQuality
    global maximumQuality    
    global barcodeSampleMapFilename
    global sampleID
            
    readInput                = None
    referenceInput           = None
    outputResultDirectory    = None
    barcodeFileLocation      = None    
    minimumReadLength        = None
    maximumReadLength        = None
    minimumQuality           = None
    maximumQuality           = None    
    barcodeSampleMapFilename = None
    sampleID                 = None

    # https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    try:
        opts, args = getopt(sys.argv[1:]
            ,"m:M:q:Q:hvbo:r:R:s:"
            ,["minlen=", "maxlen=", "minqual=", "maxqual=", "help", "version","barcode=","outputdir=","reads=", "reference=", "sampleid="])

        for opt, arg in opts:

            if opt in ('-h', '--help'):
                print (SoftwareVersion)
                usage()
                return False

            elif opt in ('-v', '--version'):
                print (SoftwareVersion)
                return False

            elif opt in ("-o", "--outputdir"):
                outputResultDirectory = arg
            elif opt in ("-r", "--reads"):
                readInput = arg
                
            # if a .fasta reference file is provided, I can calculate read statistics against it.
            elif opt in ("-R", "--reference"):
                referenceInput = arg
                
            elif opt in ("-b", "--barcode"):
                barcodeFileLocation = arg
                
                
            elif opt in ("-m", "--minlen"):
                minimumReadLength = int(arg)
            elif opt in ("-M", "--maxlen"):
                maximumReadLength = int(arg)
            
            elif opt in ("-q", "--minqual"):
                minimumQuality = int(arg)   
                
            elif opt in ("-Q", "--maxqual"):
                maximumQuality = int(arg)    
            
            elif opt in ("-s", "--sampleid"):
                sampleID = arg
                
            else:
                print('Unknown Commandline Option:' + str(opt) + ':' + str(arg))
                raise Exception('Unknown Commandline Option:' + str(opt) + ':' + str(arg))
            
        if(len(sys.argv) < 3):
            print ('I don\'t think you have enough arguments.\n')
            usage()
            return False     

    except GetoptError:        
        print ('Something seems wrong with your commandline parameters.')
        print(sys.exc_info())
        usage()
        return False

    # Sanity Checks. The only required parameters are input dir/file and output.    
    if (isfile(readInput)):
        print ('Read input is a file that exists.')
    elif (isdir(readInput)):
        print ('Read input is a directory that exists.')
    else :
        print ('I don\'t understand the read input specified, it is not a file or directory:' + readInput)
        return False
    
    # This output directory should exist
    if not exists(outputResultDirectory):
        makedirs(outputResultDirectory)
        
    # TODO: trim out spaces and special characters from the sample id.
    # Check on all the values.

    return True

if __name__=='__main__':
    
    try:    
        if(readArgs()):
            print('Commandline arguments look fine.\n Now I will prepare the reads and calculate quality statistics.')

            prepareReads(readInput, outputResultDirectory, sampleID, barcodeFileLocation, referenceInput,
                         minimumReadLength, maximumReadLength, minimumQuality, maximumQuality)

            print ('Done with nit-picker for now. Have a nice day.')


        else:
            print('\nI\'m giving up because I was not satisfied with your commandline arguments.')  
            
    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print ('Fatal problem during read extraction:')
        print (sys.exc_info())
        raise
   