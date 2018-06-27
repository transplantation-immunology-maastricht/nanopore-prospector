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
from nit_picker.nit_picker import prepareReads
from allele_wrangler.allele_wrangler import AlleleWrangler


# TODO Move this code to allele_wrangler. That's where allele calling happens.
#TODO: Talk to halagan about how i should do this now. 
"""
def alleleCall():
    print ('I will attempt an allele call.')
    
    hlaSequence = 'CAGGAGCAGAGGGGTCAGGGCGAAGTCCCAGGGCCCCAGGCGTGGCTCTCAGGGTCTCAGGCCCCGAAGGCGGTGTATGGATTGGGGAGTCCCAGCCTTGGGGATTCCCCAACTCCGCAGTTTCTTTTCTCCCTCTCCCAACCTACGTAGGGTCCTTCATCCTGGATACTCACGACGCGGACCCAGTTCTCACTCCCATTGGGTGTCGGGTTTCCAGAGAAGCCAATCAGTGTCGTCGCGGTCGCTGTTCTAAAGTCCGCACGCACCCACCGGGACTCAGATTCTCCCCAGACGCCGAGGATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTACTCTCGGGGGCCCTGGCCCTGACCCAGACCTGGGCGGGTGAGTGCGGGGTCGGGAGGGAAACCGCCTCTGCGGGGAGAAGCAAGGGGCCCTCCTGGCGGGGGCGCAGGACCGGGGGAGCCGCGCCGGGAGGAGGGTCGGGCAGGTCTCAGCCACTGCTCGCCCCCAGGCTCCCACTCCATGAGGTATTTCTTCACATCCGTGTCCCGGCCCGGCCGCGGGGAGCCCCGCTTCATCGCCGTGGGCTACGTGGACGACACGCAGTTCGTGCGGTTCGACAGCGACGCCGCGAGCCAGAAGATGGAGCCGCGGGCGCCGTGGATAGAGCAGGAGGGGCCGGAGTATTGGGACCAGGAGACACGGAATATGAAGGCCCACTCACAGACTGACCGAGCGAACCTGGGGACCCTGCGCGGCTACTACAACCAGAGCGAGGACGGTGAGTGACCCCGGCCCGGGGCGCAGGTCACGACCCCTCATCCCCCACGGACGGGCCAGGTCGCCCACAGTCTCCGGGTCCGAGATCCACCCCGAAGCCGCGGGACTCCGAGACCCTTGTCCCGGGAGAGGCCCAGGCGCCTTTACCCGGTTTCATTTTCAGTTTAGGCCAAAAATCCCCCCGGGTTGGTCGGGGCGGGGCGGGGCTCGGGGGACTGGGCTGACCGCGGGGTCGGGGCCAGGTTCTCACACCATCCAGATAATGTATGGCTGCGACGTGGGGCCGGACGGGCGCTTCCTCCGCGGGTACCGGCAGGACGCCTACGACGGCAAGGATTACATCGCCCTGAACGAGGACCTGCGCTCTTGGACCGCGGCGGACATGGCAGCTCAGATCACCAAGCGCAAGTGGGAGGCGGTCCATGCGGCGGAGCAGCGGAGAGTCTACCTGGAGGGCCGGTGCGTGGACGGGCTCCGCAGATACCTGGAGAACGGGAAGGAGACGCTGCAGCGCACGGGTACCAGGGGCCACGGGGCGCCTCCCTGATCGCCTATAGATCTCCCGGGCTGGCCTCCCACAAGGAGGGGAGACAATTGGGACCAACACTAGAATATCACCCTCCCTCTGGTCCTGAGGGAGAGGAATCCTCCTGGGTTTCCAGATCCTGTACCAGAGAGTGACTCTGAGGTTCCGCCCTGCTCTCTGACACAATTAAGGGATAAAATCTCTGAAGGAGTGACGGGAAGACGATCCCTCGAATACTGATGAGTGGTTCCCTTTGACACCGGCAGCAGCCTTGGGCCCGTGACTTTTCCTCTCAGGCCTTGTTCTCTGCTTCACACTCAATGTGTGTGGGGGTCTGAGTCCAGCACTTCTGAGTCTCTCAGCCTCCACTCAGGTCAGGACCAGAAGTCGCTGTTCCCTTCTCAGGGAATAGAAGATTATCCCAGGTGCCTGTGTCCAGGCTGGTGTCTGGGTTCTGTGCTCTCTTCCCCATCCCGGGTGTCCTGTCCATTCTCAAGATGGCCACATGCGTGCTGGTGGAGTGTCCCATGACAGATGCAAAATGCCTGAATTTTCTGACTCTTCCCGTCAGACCCCCCCAAGACACATATGACCCACCACCCCATCTCTGACCATGAGGCCACCCTGAGGTGCTGGGCCCTGGGCTTCTACCCTGCGGAGATCACACTGACCTGGCAGCGGGATGGGGAGGACCAGACCCAGGACACGGAGCTCGTGGAGACCAGGCCTGCAGGGGATGGAACCTTCCAGAAGTGGGCGGCTGTGGTGGTGCCTTCTGGAGAGGAGCAGAGATACACCTGCCATGTGCAGCATGAGGGTCTGCCCAAGCCCCTCACCCTGAGATGGGGTAAGGAGGGAGATGGGGGTGTCATGTCTCTTAGGGAAAGCAGGAGCCTCTCTGGAGACCTTTAGCAGGGTCAGGGCCCCTCACCTTCCCCTCTTTTCCCAGAGCTGTCTTCCCAGCCCACCATCCCCATCGTGGGCATCATTGCTGGCCTGGTTCTCCTTGGAGCTGTGATCACTGGAGCTGTGGTCGCTGCCGTGATGTGGAGGAGGAAGAGCTCAGGTGGAGAAGGGGTGAAGGGTGGGGTCTGAGATTTCTTGTCTCACTGAGGGTTCCAAGCCCCAGCTAGAAATGTGCCCTGTCTCATTACTGGGAAGCACCTTCCACAATCATGGGCCGACCCAGCCTGGGCCCTGTGTGCCAGCACTTACTCTTTTGTAAAGCACCTGTTAAAATGAAGGACAGATTTATCACCTTGATTACGGCGGTGATGGGACCTGATCCCAGCAGTCACAAGTCACAGGGGAAGGTCCCTGAGGACAGACCTCAGGAGGGCTATTGGTCCAGGACCCACACCTGCTTTCTTCATGTTTCCTGATCCCGCCCTGGGTCTGCAGTCACACATTTCTGGAAACTTCTCTGGGGTCCAAGACTAGGAGGTTCCTCTAGGACCTTAAGGCCCTGGCTCCTTTCTGGTATCTCACAGGACATTTTCTTCCCACAGATAGAAAAGGAGGGAGTTACACTCAGGCTGCAAGTAAGTATGAAGGAGGCTGATGCCTGAGGTCCTTGGGATATTGTGTTTGGGAGCCCATGGGGGAGCTCACCCACCCCACAATTCCTCCTCTAGCCACATCTTCTGTGGGATCTGACCAGGTTCTGTTTTTGTTCTACCCCAGGCAGTGACAGTGCCCAGGGCTCTGATGTGTCTCTCACAGCTTGTAAAGGTGAGAGCTTGGAGGGCCTGATGTGTGTTGGGTGTTGGGTGGAACAGTGGACACAGCTGTGCTATGGGGTTTCTTTGCGTTGGATGTATTGAGCATGCGATGGGCTGTTTAAGGTGTGACCCCTCACTGTGATGGATATGAATTTGTTCATGAATATTTTTTTCTATAGTGTGAGACAGCTGCCTTGTGTGGGACTGAGAGGCAAGAGTTGTTCCTGCCCTTCCCTTTGTGACTTGAAGAACCCTGACTTTGTTTCTGCAAAGGCACCTGCATGTGTCTGTGTTCGTGTAGGCATAATGTGAGGAGGTGGGGAGAGCACCCCACCCCCATGTCCACCATGACCCTCTTCCCACGCTGACCTGTGCTCCCTCTCCAATCATCTTTCCTGTTCCAGAGAGGTGGGGCTGAGGTGTCTCCATCTCTGTCTCAACTTCATGGTGCACTGAGCTGTAACTTCTTCCTTCCCTATTAAAA'
    gene = 'HLA-A'
    
    print ('this sequence:\n')
    print(hlaSequence)
    
    from Bio import SeqIO
    from BioSQL import BioSeqDatabase
    from seqann.sequence_annotation import BioSeqAnn
    
    import pygfe
    
    seq_file = 'test_dq.fasta'
    gfe = pygfe.pyGFE()
    server = BioSeqDatabase.open_database(driver="pymysql", user="root",
        passwd="", host="localhost",
                                                                                                                     db="bioseqdb")
    seqann = BioSeqAnn(server=server)
    seq_rec = list(SeqIO.parse(seq_file, 'fasta'))[0]
    annotation = seqann.annotate(seq_rec, "HLA-DQB1")
    gfe = gfe.get_gfe(annotation, "HLA-DQB1")
    print(gfe)
    
"""   

def readArgs():
    # Default to None.  So I can easily check if they were not passed in.
   
    global readInput
    global referenceInput
    global outputResultDirectory
    global analysisAction
    global minimumReadLength
    global maximumReadLength
    global minimumQuality
    global maximumQuality
    global sampleID
    global barcodeFileLocation
    
    readInput                = None
    referenceInput           = None
    outputResultDirectory    = None
    analysisAction           = None
    minimumReadLength        = None
    maximumReadLength        = None
    minimumQuality           = None
    maximumQuality           = None  
    sampleID                 = None
    barcodeFileLocation      = None

    if(len(argv) < 3):
        print ('I don\'t think you have enough arguments.\n')
        #usage()
        return False    

    # https://www.tutorialspoint.com/python/python_command_line_arguments.htm
    try:
        opts, args = getopt(argv[1:]
            ,"m:M:q:Q:hvb:i:o:r:R:t:a:s:"
            ,["minlen=", "maxlen=", "minqual=", "maxqual=", "help", "version","barcode=", "iterations=","outputdir=","reads=","reference=",'threads=','action=', 'sampleid='])

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

            elif opt in ("-i", "--iterations"):
                numberIterations = arg
                
            elif opt in ("-o", "--outputdir"):
                outputResultDirectory = arg
                
            elif opt in ("-a", "--action"):
                analysisAction = arg

            else:
                print('Unknown Commandline Option:' + str(opt) + ':' + str(arg))
                raise Exception('Unknown Commandline Option:' + str(opt) + ':' + str(arg))
            

    except GetoptError:
        print ('Something seems wrong with your commandline parameters.')
        #print (errorMessage)
        #usage()
        return False

    

    # Consensus,threads is optional, the rest are necessary.
    # Sanity Checks
    #if (int(numberIterations) < 1):
    ##    print('You must specify an Iteration count >= 1.')
    
    if (isfile(readInput)):
        print ('Read input is a file that exists.')
    elif (isdir(readInput)):
        print ('Read input is a directory that exists.')
    else :
        print ('I don\'t understand the read input specified, it is not a file or directory:' + readInput)
        return False

    return True


if __name__=='__main__':

    # TODO: Read args? What args do i even need?
    readArgs()
    
    # Quick and Dirty way to do a full analysis. Sort, Split, Assemble all reads.
    if(analysisAction=="fullanalysis"):
        if(readInput is not None
            and outputResultDirectory is not None):
            #print ('I am ready to run this from the command line.')
            root = Tk()
            app = NanoporeProspectorMasterFrame(root)
            #root.mainloop()
            
            # Load config file? 
            # Actually maybe this already happened.
            
            # Set input and output directory of the prospector.
            # Just put the text into the text fields, that's probably easier?
            
            app.setInputDir(readInput, False)
            app.setOutputDir(outputResultDirectory)
            #app.inputDirectoryText.set(readInput)
            #app.outputDirectoryText.set(outputResultDirectory)
            #assignConfigurationValue('output_directory', self.outputDirectoryText.get())
            
            # Start full analysis
            app.runFullAnalysis()
            
        else:
            raise Exception("Can not do a full analysis without a read input and output directory.")
        
    elif(analysisAction=="heterosplit"):

        print('Doing the hetero split.')

        # Calculate the prepared reads folder name.
        preparedReadsFolder = join(outputResultDirectory, '1_prepared_reads')

        # Quality control, reject big reads etc.
        # I will pass "None" as the reference sequence here because
        # I don't need complicated read statistics in this case.
        # TODO: add an "analysis action" for just preparing reads. Programming is fun.
        prepareReads(readInput, preparedReadsFolder, sampleID, barcodeFileLocation, None, minimumReadLength, maximumReadLength, minimumQuality, maximumQuality )

        # If we used length/quality filters, the data is in the "Pass" reads.
        if isfile(join(preparedReadsFolder, 'minion_reads_Pass.fastq')):
            preparedReadsFileName = join(preparedReadsFolder, 'minion_reads_Pass.fastq')
        else:
            preparedReadsFileName = join(preparedReadsFolder, 'minion_reads_All.fastq')
        splitReadsFolder = join(outputResultDirectory, '2_phased_reads')

        # I don't need to sort the reads in this case....Just hetero split
        # I can hetero split with or without a reference sequence.

        # TODO Parameters don't come from anywhere, i can pass this in.
        #numberIterations = 6
        numberThreads = 8
        splitHeterozygotes = True

        # Call the hetero split method.
        # I guess i need an allele wrangler object
        myAlleleWrangler = AlleleWrangler(preparedReadsFileName, splitReadsFolder, referenceInput, 0,
                                          numberThreads, splitHeterozygotes)

        currentAssemblyResults = myAlleleWrangler.analyzeReads()

    elif (analysisAction == "snpanalysis"):
        print('Doing some SNP analysis now.')

        prepareReads(readInput, outputResultDirectory, sampleID, None, referenceInput,
            minimumReadLength, maximumReadLength, minimumQuality, maximumQuality)

    # TODO: Analysis Actions to add:
    # HLA Allele Analysis. I can merge that code into prospector.
    # Splicing analysis.
    # Demultiplex
    # Sort by Gene.






    else:
        root = Tk()
        app = NanoporeProspectorMasterFrame(root)
        root.mainloop()

    print('Done.  Yay.')
    
    


