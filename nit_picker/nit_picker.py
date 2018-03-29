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

from os import listdir
from os.path import join, isdir, isfile

from numpy import mean, amin, amax

from .minion_read_collection import minionReadCollection, createCollectionFromReadFile
#, #createScatterPlot, createOutputFile, printSimpleReadStats
from .barcode_demultiplexer import splitByBarcode
from nanopore_prospector.common import getReadFileType

from Bio.SeqIO import write
from allele_wrangler.allele_wrangler import createOutputFile

#nit_picker.prepareReads(readInput, outputResultDirectory, sampleID, barcodeFileLocation, minimumReadLength, maximumReadLength, minimumQuality, maximumQuality )
#inputReads and outputDirectory are necessary.  The rest can be "None" if you would like.            
def prepareReads(inputReads, outputDirectory, sampleID, barcodeFileLocation, referenceFileLocation, minimumReadLength, maximumReadLength, minimumReadQuality, maximumReadQuality ):
    print ('Preparing Reads')
    
    #TODO: Create an output file with read stats, like the read extractor did.
    #Create Read Stat File Here.
    # Just a brief summary. Maybe it's a method that can be returned, so i can include it in the nanopore prospector results.
 
    # Default sample id
    if (sampleID is None):
        sampleID = 'minion_reads'

    # Determine if input is a file, or a directory.
    # Both should work equally well.
    if (isfile(inputReads)):
        print ('Read input is a file that exists.')
        allReads = createCollectionFromReadFile(inputReads)
    elif (isdir(inputReads)):
        print ('Read input is a directory that exists.')
        allReads = minionReadCollection([])
        for currentInputReadFile in listdir(inputReads):
            
            allReads.readInputFormat = getReadFileType(currentInputReadFile)

            print ('loading Reads from:' + str(join(inputReads,currentInputReadFile)))
            newReads = createCollectionFromReadFile(join(inputReads,currentInputReadFile))
            allReads.concatenate(newReads)
            
    else :
        print('I expect a .fasta or .fastq format for read input. Alternatively, specify a directory containing read inputs. Please check your input.')
        raise Exception('Bad Read Input Format')
    
    print ('Total # of allReads found = ' + str(len(allReads.readCollection)))

    passReads=[]
    lengthRejectReads=[]
    qualityRejectReads=[]

    # Iterate all reads
    print ('Rejecting reads for Length and Quality...')
    for currentRead in allReads.readCollection:
        
        currentSeqLength = len(currentRead)
        
        try:
            phredQualities = currentRead.letter_annotations["phred_quality"]                    
            currentAvgPhredQuality = mean(phredQualities)
        except Exception:
            #print ('Cannot find Phred qualities for a read. I think we are looking at a fasta file.')
            currentAvgPhredQuality = None
        
        # Reject allReads that are wrong length
        if( (minimumReadLength is not None and currentSeqLength < minimumReadLength)
            or (maximumReadLength is not None and currentSeqLength > maximumReadLength)
            ):
            lengthRejectReads.append(currentRead)

        # Reject allReads that have the wrong quality
        elif( currentAvgPhredQuality is not None and
            ((minimumReadQuality is not None and currentAvgPhredQuality < minimumReadQuality)
            or (maximumReadQuality is not None and currentAvgPhredQuality > maximumReadQuality))
            ):
            qualityRejectReads.append(currentRead)
      
        # This read is okay.
        else:
            passReads.append(currentRead)
           
    passReadCollection=minionReadCollection(passReads)
    lengthRejectReadCollection=minionReadCollection(lengthRejectReads)
    qualityRejectReadCollection=minionReadCollection(qualityRejectReads)
    
    # TODO I shouldn't have to assign this. But an array of reads does not know what format they used to be in
    # So, i can't detect read type in minion_read_collection, that's dumb, im sure i can think of a better solution.
    passReadCollection.readInputFormat = allReads.readInputFormat
    lengthRejectReadCollection.readInputFormat = allReads.readInputFormat
    qualityRejectReadCollection.readInputFormat = allReads.readInputFormat

    # readStats is a dictionary.
    # a read stats is a tuple, with 2 arrays.
    # readstats is a 2d array, with lengths and qualities.
    readStats = {}

    # Output Reads
    if(len(allReads.readCollection)>0):
        readStats['all_reads'] = allReads.outputReadPlotsAndSimpleStats(outputDirectory, 'All', sampleID, referenceFileLocation)
      
    
    # If Pass reads are not the same size as All reads, 
    # We must print the rejected reads, and the pass reads as well.   
    if(len(allReads.readCollection) != len(passReadCollection.readCollection)):

        if(len(lengthRejectReadCollection.readCollection)>0):
            readStats['length_rejected'] = lengthRejectReadCollection.outputReadPlotsAndSimpleStats(outputDirectory, 'Length Rejected', sampleID, referenceFileLocation)
            
        if(len(qualityRejectReadCollection.readCollection)>0):
            readStats['quality_rejected'] = qualityRejectReadCollection.outputReadPlotsAndSimpleStats(outputDirectory, 'Quality Rejected', sampleID, referenceFileLocation)
    
        if(barcodeFileLocation is None):
            print('No barcode file was provided. I will not attempt to de-multiplex the allReads.')
            readStats['pass_reads'] = passReadCollection.outputReadPlotsAndSimpleStats(outputDirectory, 'Pass', sampleID, referenceFileLocation)
    
        else:
            print('Provided Barcode File:' + str(barcodeFileLocation) + '\nI will de-multiplex the allReads.')
            # Debarcode pass reads
            
            # TODO: This loop is untested. Did I do these stats right? 
            # They should appear in the output file.
            barcodeReadStats = splitByBarcode(passReadCollection, outputDirectory, barcodeFileLocation, sampleID, referenceFileLocation)
            for key in barcodeReadStats.keys():
                readStats[key] = barcodeReadStats[key]

    writeReadStats(readStats, outputDirectory)
    return readStats

def writeReadStats(readStats, outputDirectory):
    outputFileName = join(outputDirectory,'Read_Summary.txt')
    outputFile = createOutputFile(outputFileName)
      
    for key in readStats.keys():
        outputFile.write(writeReadStatsSub(key, readStats[key]))        
    outputFile.close()
     
    
def writeReadStatsSub(readType, readStats):    
    statsSubText=''        
    
    if(len(readStats) > 0):
        readLengths = readStats[0]
        readQualities = readStats[1]
        
        statsSubText += (readType + '\n' 
            + 'Sequences Analyzed :' + str(len(readLengths)) + '\n'
            + 'Minimum Length     :' + str(int(amin(readLengths))) + '\n'
            + 'Maximum Length     :' + str(int(amax(readLengths))) + '\n'
            + 'Mean Length        :' + str(mean(readLengths)) + '\n'
            + 'Mean Quality       :' + str(mean(readQualities)) + '\n\n'
        )
        
    else:
        statsSubText += (readType + '\n' 
            + 'Sequences Analyzed :0\n\n'

        )

    return statsSubText

