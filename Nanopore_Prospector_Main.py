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

SoftwareVersion = "nanopore_prospector Version 1.0"
from tkinter import Tk

from sys import argv
from getopt import getopt, GetoptError

from os.path import isfile, isdir, join

from nanopore_prospector.Nanopore_Prospector_Master_Frame import NanoporeProspectorMasterFrame
from nanopore_prospector.common import getDirectoryAlignmentStats
from nit_picker.nit_picker import prepareReads
from allele_wrangler.allele_wrangler import AlleleWrangler
from Parse_IMGT_HLA.HLA_XML_to_fasta import parseImgtHla
from find_homopolymers.find_homopolymers import findHomopolymers
from convert_reads.convert_reads import fastqToFasta
from combine_snps.consensus_from_snp_table import consensusSequenceFromSNPFiles
from combine_snps.name_novels import nameNovels
from combine_DP_region.combine_DP_region import combineDPSequences
from combine_snps.find_files_with_pattern import combineFastaFiles, collectFilesWithPattern, removeDuplicates

def readArgs():
    # Trying to be consistent and global with my parameter inputs.
   
    global readInput
    global referenceInput
    global inputFile
    global inputDirectory
    global outputDirectory
    global outputFile
    global analysisAction
    global minimumReadLength
    global maximumReadLength
    global minimumQuality
    global maximumQuality
    global sampleID
    global barcodeFileLocation
    global snps
    global minimumSnpPctCutoff
    global filePattern

    readInput                = None
    referenceInput           = None
    inputFile                = None
    inputDirectory           = None
    outputDirectory          = None
    outputFile               = None
    analysisAction           = None
    minimumReadLength        = None
    maximumReadLength        = None
    minimumQuality           = None
    maximumQuality           = None  
    sampleID                 = None
    barcodeFileLocation      = None
    snps                     = None
    minimumSnpPctCutoff      = None
    filePattern              = None

    if(len(argv) < 3):
        print ('I don\'t think you have enough arguments.\n')
        #usage()
        return False

    print('Attempting to load commandline arguments')
    # https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    try:
        opts, args = getopt(argv[1:]
            ,"m:M:q:Q:hvb:I:i:o:O:r:R:t:a:s:S:p:P:"
            ,["minlen=", "maxlen=", "minqual=", "maxqual=", "help", "version","barcode=", "inputfile=", "inputdirectory="
                ,"outputdirectory=","outputfile=","reads=","reference=",'threads=','action=', 'sampleid=', 'snps=', 'snpcutoff=', 'pattern='])

        print (str(len(opts)) + ' arguments found.')

        for opt, arg in opts:

            if opt in ('-h', '--help'):
                #print (SoftwareVersion)
                #usage()
                return False

            elif opt in ('-v', '--version'):
                #print (SoftwareVersion)
                return False
            
            elif opt in ('-r', '--reads'):
                #print (SoftwareVersion)
                readInput = arg
                
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
                                        
            elif opt in ('-R', '--reference'):
                #print (SoftwareVersion)
                referenceInput = arg

            # Potential bug: I just changed the -i and -I options to match the input/output file.
            # This removes the --iterations option. Might need to addt hat again later.
            # Might also cause a confusion between inputFile and inputDirectory. watch for that bug.
            elif opt in ("-i", "--inputdirectory"):
                inputDirectory = arg

            elif opt in ("-I", "--inputfile"):
                inputFile = arg

            elif opt in ("-o", "--outputdirectory"):
                outputDirectory = arg

            elif opt in ("-O", "--outputfile"):
                outputFile = arg

            elif opt in ("-a", "--action"):
                analysisAction = arg

            elif opt in ("-s", "--sampleid"):
                sampleID = arg

            elif opt in ("-p", "--snpcutoff"):
                # Why is it -p? polymorphism.
                # This is a tuning parameter, I consider a SNP if the percent of mismatched/deleted is greater than the minimumSnpPctCutoff.
                # So i would use something like .0000001 for consensus, or maybe .20 for minion reads.
                minimumSnpPctCutoff = float(arg)
                if (minimumSnpPctCutoff > 1 or minimumSnpPctCutoff < 0):
                    raise Exception('minimumSnpPctCutoff is the minimum proportion of unmatching reads for a position to be considered polymorphic. It should be between 0 and 1.')

            elif opt in ("-P", "--pattern"):
                filePattern = arg

            elif opt in ("-S", "--snps"):
                # A list of polymorphic positions for analysis. Used mainly (only?) in the heterozygous split.
                # It should be a list of 1-based SNP positions matching the reference sequence, separated by commas, like this:
                # 134, 664, 1233, 1445
                # I choose 1-based because IGV shows 1-based positions.
                # 134,664,1233,1445
                # convert from 1 to 0-based and store as a list.
                snps = str(arg).split(',')
                for index, snp in enumerate(snps):
                    snps[index] = int(snps[index]) - 1


                #snps = arg

            else:
                print('Unknown Commandline Option:' + str(opt) + ':' + str(arg))
                raise Exception('Unknown Commandline Option:' + str(opt) + ':' + str(arg))
            

    except GetoptError as err:
        print ('Something seems wrong with your commandline parameters.')
        print (str(err.msg))
        print ('Unknown option: ' + str(err.opt))
        #print (errorMessage)
        #usage()
        return False

    # TODO: Move the sanity checks inside the use cases.
    #if (isfile(readInput)):
#        print ('Read input is a file that exists.')
#    elif (isdir(readInput)):
#        print ('Read input is a directory that exists.')
#    else :
#        print ('I don\'t understand the read input specified, it is not a file or directory:' + readInput)
#        return False

    return True


if __name__=='__main__':
    # TODO: Read args? What args do i even need?
    readArgs()

    # TODO: how to check parameters in general? Obviously each action has different inputs, can I check that they exist?
    # In a general sense?

    # Quick and Dirty way to do a full analysis. Sort, Split, Assemble all reads.
    if(analysisAction=="fullanalysis"):
        # Quick sanity check.
        if(readInput is not None and outputDirectory is not None):
            # TODO: I don't need the gui to do this, i am cheating here. Separate the logic from the GUI. Maybe even delete the GUI, i hate it.
            root = Tk()
            app = NanoporeProspectorMasterFrame(root)

            # Load config file? 
            # Actually maybe this already happened.
            
            # Set input and output directory of the prospector.
            # Just put the text into the text fields, that's probably easier?
            
            app.setInputDir(readInput, False)
            app.setOutputDir(outputDirectory)
            #app.inputDirectoryText.set(readInput)
            #app.outputDirectoryText.set(outputResultDirectory)
            #assignConfigurationValue('output_directory', self.outputDirectoryText.get())
            
            # Start full analysis
            app.runFullAnalysis()
            
        else:
            raise Exception("Can not do a full analysis without a read input and output directory.")


    elif(analysisAction=="preparereads"):
        # Prepare reads and do simple statistics. Works for Reads or, potentially, HLA alleles from IPD, consensus sequences.
        # If a reference file is provided, I can also calculate information on how polymorphic a region is.
        print('Preparing Reads......')

        prepareReads(readInput, outputDirectory, sampleID, barcodeFileLocation, referenceInput, minimumReadLength, maximumReadLength, minimumQuality, maximumQuality, True, minimumSnpPctCutoff)

        print('Done.')


    elif(analysisAction=="heterosplit"):

        print('Doing the hetero split.')

        # Calculate the prepared reads folder name.
        preparedReadsFolder = join(outputDirectory, '1_prepared_reads')

        # Quality control, reject big reads etc.
        # I will pass "None" as the reference sequence here because
        # I don't need complicated read statistics in this case.
        # TODO: add an "analysis action" for just preparing reads. Programming is fun.
        # TODO: I just switched the reference input from None to referenceInput. This changes...that we will actually calculate read stats during htis step. Probably not an issue but it will be slower.
        # TODO: to fix this, i should refactor this prepareReads method. Somehow pass in yes or no on calculate read stats.
        prepareReads(readInput, preparedReadsFolder, sampleID, barcodeFileLocation, referenceInput, minimumReadLength, maximumReadLength, minimumQuality, maximumQuality, False )

        # If we used length/quality filters, the data is in the "Pass" reads.
        if isfile(join(preparedReadsFolder, 'minion_reads_Pass.fastq')):
            preparedReadsFileName = join(preparedReadsFolder, 'minion_reads_Pass.fastq')
        else:
            preparedReadsFileName = join(preparedReadsFolder, 'minion_reads_All.fastq')
        splitReadsFolder = join(outputDirectory, '2_phased_reads')

        # I don't need to sort the reads in this case....Just hetero split
        # I can hetero split with or without a reference sequence.

        # TODO Parameters don't come from anywhere, i can pass this in.
        # T for threads? Im already passing iteration count.
        #numberIterations = 6
        numberThreads = 8
        splitHeterozygotes = True

        # Call the hetero split method.
        # I guess i need an allele wrangler object
        # TODO: Can I supply a list of polymorphic positions? Yeah pass it in.

        myAlleleWrangler = AlleleWrangler(preparedReadsFileName, splitReadsFolder, referenceInput, 0,
                                          numberThreads, splitHeterozygotes, snps)

        currentAssemblyResults = myAlleleWrangler.analyzeReads()

    elif (analysisAction == "snpanalysis"):
        print('Doing some SNP analysis now.')

        prepareReads(readInput, outputDirectory, sampleID, None, referenceInput, minimumReadLength, maximumReadLength, minimumQuality, maximumQuality, True, minimumSnpPctCutoff)


    elif (analysisAction == "alignmentstatistics"):
        # Search directory for .bam files and output a file with read coverage statistics.
        pass
        getDirectoryAlignmentStats(inputDirectory, outputFile, True)



    # TODO: Analysis Actions to add:
    # Splicing analysis.
    # Demultiplex
    # Sort by Gene.


    elif (analysisAction == 'splicinganalysis'):
        print('I would love to do the MinION Read Splicing Analysis, but Ben has not yet implemented it. Please ask him to.')

    elif (analysisAction == 'collectreads'):
        # Nanopolish will require access to .fast5 read files in order to polish a consensus sequence.
        # After sorting and preparing the reads, a .fastq file should contain the reads that support a consensus.
        # Supply the fastq read file, the directory of the fast5 reads.

        print('Ben has not implemented this functionality yet.')

    elif (analysisAction == 'combinesnps'):
        print('Ben has not implemented this functionality yet.')

    elif (analysisAction == 'fastqtofasta'):
        # Convert a fastq file to a fasta file, effectively deleting the quality scores.
        # Required:
        # -i / --inputfile
        # -O / --outputfile
        fastqToFasta(inputFile, outputFile, 36)

    elif (analysisAction == 'extractsequences'):
        print('The main method is in extract_sequences_main.py, move it here.')

    elif (analysisAction == 'findfiles'):
        collectFilesWithPattern(inputDirectory, outputDirectory, filePattern)

    elif(analysisAction == 'combinefastafiles'):
        #print('This method works similary to collectFilesWithPattern, but specific for making a fasta file.')
        combineFastaFiles(inputDirectory, outputDirectory, filePattern)

    elif (analysisAction == 'removeduplicates'):
        # print('This method works similary to collectFilesWithPattern, but specific for making a fasta file.')
        removeDuplicates(inputFile, outputFile)


    elif (analysisAction == 'namenovels'):
        nameNovels(referenceInput, inputFile, outputDirectory)

    elif (analysisAction == 'consensusfromsnps'):
        consensusSequenceFromSNPFiles(referenceInput, inputFile, outputDirectory)

    elif (analysisAction == 'fillsnptable'):
        print('Ben has not implemented this functionality yet.')

    elif (analysisAction == 'findhomopolymers'):
        # Input a .fasta full of sequences, output some statistics on the homopolymers in these sequences.
        # Required:
        # -i / --inputfile
        # -o / --outputdirectory
        print('Searching for Homopolymers.')

        # TODO: Sanity Checks on required inputs
        findHomopolymers(inputFile, outputDirectory)

    elif (analysisAction == 'hlaalleleanalysis'):
        # To input a IMGT XML file and generate HLA .fasta reference sequences, for a variety of purposes.
        # Required:
        # -i / --inputfile
        # -o / --outputdirectory
        print('Extracting allele sequences from IPD-IMGT/HLA XML file.')

        # TODO: Sanity Checks on required inputs
        parseImgtHla(inputFile, outputDirectory)

    elif (analysisAction == 'combinedpsequences'):
        # This is a specific analysis to combine DPA1, Promoter, DPB1 sequences from 3 amplicons.

        # -i / --inputfile
        #   A fasta file containing, in this exact order:
        #   DP region reference.
        #   DPA1 sequence 5'-3' (I will reverse complement the sequence.)
        #   DPB1 sequence 5'-3'
        #   Promoter sequence 5'-3'
        # -o / --outputdirectory
        #   I write a few intermediate and result files in there.
        print('Combining DP sequences from the fasta file:\n' + str(inputFile) + '\nand write results to:\n' + str(outputDirectory))

        # TODO: Sanity Checks on required inputs
        combineDPSequences(inputFile, outputDirectory)

    elif (analysisAction is None or analysisAction == ''):
        # Start the main GUI, so we can do some analysis steps.
        # Honestly this GUI is a PITA, and doesn't provide much. TODO: Delete the GUI.
        root = Tk()
        app = NanoporeProspectorMasterFrame(root)
        root.mainloop()

    else:
        raise Exception('Unknown analysis action:' + str(analysisAction))

    print('Done.  Yay.')
    
    


