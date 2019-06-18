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

from numpy import mean

from .minion_read_collection import MinionReadCollection, createCollectionFromReadFile, writeReadStats
from .barcode_demultiplexer import splitByBarcode
from nanopore_prospector.common import getReadFileType

from Bio.SeqIO import parse as parseRecords

from Bio.SeqIO import write
from allele_wrangler.allele_wrangler import createOutputFile

#nit_picker.prepareReads(readInput, outputResultDirectory, sampleID, barcodeFileLocation, minimumReadLength, maximumReadLength, minimumQuality, maximumQuality )
#inputReads and outputDirectory are necessary.  The rest can be "None" if you would like.            
def prepareReads(inputReads, outputDirectory, sampleID, barcodeFileLocation, referenceFileLocation, minimumReadLength, maximumReadLength, minimumReadQuality, maximumReadQuality, calculateReadStatistics, minimumSnpPctCutoff):
    print ('Preparing Reads:' + inputReads)
    
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

        if(referenceFileLocation is not None):
            allReads.setReferenceSequence(list(parseRecords(referenceFileLocation, 'fasta'))[0])
        else:
            allReads.setReferenceSequence(None)

    elif (isdir(inputReads)):
        # Loop through input files and concatenate them.
        print ('Read input is a directory that exists.')
        allReads = MinionReadCollection([])
        for currentInputReadFile in listdir(inputReads):

            allReads.readInputFormat = getReadFileType(currentInputReadFile)

            print ('loading Reads from:' + str(join(inputReads,currentInputReadFile)))
            newReads = createCollectionFromReadFile(join(inputReads,currentInputReadFile))

            # temp debug stuff...for some reason phred qualities are missing.
            print('Read Format: ' + allReads.readInputFormat   )
            print ('Allreads has ' + str(len(allReads.readCollection)))
            print ('newreads has ' + str(len(newReads.readCollection)))
            print('allreads avgphred= ' + str(len(allReads.readAvgReportedPhredQualities)))
            print('newreads avgphred= ' + str(len(newReads.readAvgReportedPhredQualities)))

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
           
    passReadCollection=MinionReadCollection(passReads)
    lengthRejectReadCollection=MinionReadCollection(lengthRejectReads)
    qualityRejectReadCollection=MinionReadCollection(qualityRejectReads)
    
    # TODO I shouldn't have to assign this. But an array of reads does not know what format they used to be in
    # So, i can't detect read type in minion_read_collection, that's dumb, im sure i can think of a better solution.
    passReadCollection.readInputFormat = allReads.readInputFormat
    lengthRejectReadCollection.readInputFormat = allReads.readInputFormat
    qualityRejectReadCollection.readInputFormat = allReads.readInputFormat

    # readStats is a dictionary.
    # a read stats is a tuple, with 2 arrays.
    # readstats is a 2d array, with lengths and qualities.
    # TODO: I think this information is included in MinIONReadCollection object, I  think i can remove this and just
    # TODO: use the information inside the object.
    # Well....not so sure because this is keeping track of all our sets of reads, even rejected ones.
    # A good idea: TODO: make MinionReadCollection keep track of reads that are rejected for Length, and Quality.
    readStats = {}


    # Output Reads
    # TODO: A potential problem. I'm passing a "minimum read length" as a "minimum alignment length"
    # TODO: I divided it by 2, to help ensure the data will map. Relatively untested.

    # Here's a bandaid, in case the minimum Read Length is None, I can't divide it by 2. Set it as something arbitrarily low.
    if( minimumReadLength is None):
        minimumReadLength = 2

    if(len(allReads.readCollection)>0):
        readStats['all_reads'] = allReads.outputReadPlotsAndSimpleStats(outputDirectory, 'All', sampleID, referenceFileLocation, calculateReadStatistics, minimumReadLength / 2, minimumSnpPctCutoff)
      
    
    # If Pass reads are not the same size as All reads, 
    # We must print the rejected reads, and the pass reads as well.   
    if(len(allReads.readCollection) != len(passReadCollection.readCollection)):

        if(len(lengthRejectReadCollection.readCollection)>0):
            readStats['length_rejected'] = lengthRejectReadCollection.outputReadPlotsAndSimpleStats(outputDirectory, 'Length Rejected', sampleID, referenceFileLocation, False, minimumReadLength / 2, minimumSnpPctCutoff)
            
        if(len(qualityRejectReadCollection.readCollection)>0):
            readStats['quality_rejected'] = qualityRejectReadCollection.outputReadPlotsAndSimpleStats(outputDirectory, 'Quality Rejected', sampleID, referenceFileLocation, False, minimumReadLength / 2, minimumSnpPctCutoff)
    
        if(barcodeFileLocation is None):
            print('No barcode file was provided. I will not attempt to de-multiplex the allReads.')
            readStats['pass_reads'] = passReadCollection.outputReadPlotsAndSimpleStats(outputDirectory, 'Pass', sampleID, referenceFileLocation, calculateReadStatistics, minimumReadLength / 2, minimumSnpPctCutoff)
    
        else:
            print('Provided Barcode File:' + str(barcodeFileLocation) + '\nI will de-multiplex the allReads.')
            # Debarcode pass reads
            
            # TODO: This loop is untested. Did I do these stats right? 
            # They should appear in the output file.
            barcodeReadStats = splitByBarcode(passReadCollection, outputDirectory, barcodeFileLocation, sampleID, referenceFileLocation)
            for key in barcodeReadStats.keys():
                readStats[key] = barcodeReadStats[key]


    # TODO: I don't really want to pass in these stats like this. MinionReadCollection should keep track of pass/fail reads, so it's not necessary.
    # These variables are missing from the read collection until I call setReference.
    #alignedReadCount = allReads.totalAlignedReadCount
    #allReads.totalAlignedReadCount is not assigned yet. I can cheat and just count the # of reads in the collection. This is a bug.
    alignedReadCount = len(allReads.readCollection)
    meanAlignedReadLength = mean(allReads.readLengths)
    meanCalculatedQuality = mean(allReads.readAvgReportedPhredQualities)
    writeReadStats(readStats, outputDirectory, alignedReadCount, meanAlignedReadLength, meanCalculatedQuality)


    return readStats







