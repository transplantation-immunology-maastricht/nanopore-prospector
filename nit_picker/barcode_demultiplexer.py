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


import sys
import os

from .minion_read_collection import minionReadCollection
from Bio.Seq import Seq
#from Bio.SeqIO import write
#from os.path import join

# Read Barcodes from an input file
def readBarcodes(barcodeFileNameWithPath):
    try:
        print ('Reading Barcodes:' + barcodeFileNameWithPath)

        barcodes = {}
        barcodeFile = open(barcodeFileNameWithPath, 'r')
        for i, line in enumerate(barcodeFile):
            if (len(line.strip()) > 0):
                # Filter comments.
                if not (line.strip()[0:1]=='#'):
                    tokens = line.split()
                    barcodes[tokens[0]] = tokens[1]

        barcodeFile.close()
        return barcodes

    except Exception:
        print ('Problem reading barcode file: ' + str(sys.exc_info()[0]))
        raise
        
def splitByBarcode(currentReadCollection, outputDirectory, barcodeFileNameWithPath, sampleID, referenceFileLocation):
    #print('Searching for barcodes in ' + inputDirectory)
    barcodeList = readBarcodes(barcodeFileNameWithPath)

    print('Preparing output files in ' + outputDirectory)

    #allReadFileOutputs = {}
    allBarcodeCollections = {}

    for key in barcodeList.keys():
        #allReadFileName = join(outputDirectory, key + '_pass_reads.fastq')
        #allReadFileOutputs[key] = open(allReadFileName, 'w')
        allBarcodeCollections[key] = minionReadCollection([])
        
        # TODO: This is completely hardcoded. I need to assign the read format to the minion_read_collection.
        # I should look this up based on the input data, i think i have a common method for finding file type.
        allBarcodeCollections[key].readInputFormat = 'fastq'
        
    #unbarcodedOutput = open(join(outputDirectory, 'unbarcoded_pass_reads.fastq'), 'w')
    unbarcodedReadCollection = minionReadCollection([])
    
    # TODO: This is completely hardcoded. I need to assign the read format to the minion_read_collection.
    # I should look this up based on the input data, i think i have a common method for finding file type.
    unbarcodedReadCollection.readInputFormat = 'fastq'

    print('Splitting by barcodes, Just a second...')
     
    # TODO: Maybe I can parse the reads in one step, and write in another
    # Writing a sequence to the hard drive for each iteration is quite slow.
    # But if i keep all the data in memory, it also might be slow.
    for rec in currentReadCollection.readCollection:
        
        #print('Searching Barcodes,id=' + rec.id)
        barcodeFound = False
    
        # Loop through each barcode.
        # Note: It is possible for a read to have multiple barcodes. 
        for barcodeKey in barcodeList.keys():
            barcodeText = barcodeList[barcodeKey]
            
            # TODO: Maybe I can do the revcoms before I look for barcode, instead of EVERY time.
            revcomBarcode = Seq(barcodeText).reverse_complement()
            #print ('barcode:' + barcodeText)
            #print ('reverse complement barcode:' + revcomBarcode)
            # Here is where I check if the barcode exists in the reads.  
       
            # TODO: Trim the barcode out of the read. We don't care about that sequence.
            if barcodeText in rec:
                #print('I found ' + barcodeKey + ' in the record with id:' + rec.id)
                allBarcodeCollections[barcodeKey].readCollection.append(rec)
                barcodeFound = True
                 
            elif revcomBarcode in rec:
                barcodeFound = True
                allBarcodeCollections[barcodeKey].readCollection.append(rec)
                     
            else:
                #Don't write the file here. You're too deep in the loop.
                #SeqIO.write([rec], unbarcodedOutput, 'fastq')
                
                pass
   
        if(barcodeFound):
            # Great.  Do nothing.
            pass
        else:
            unbarcodedReadCollection.readCollection.append(rec)
            #SeqIO.write([rec], unbarcodedOutput, 'fastq')
       
    # this is for a return value
    barcodedReadStats = {}
         
    #close output files. Delete empties and generate scatterplots
    for barcodeKey in barcodeList.keys():
        #print('How many reads for barcode:' + str(barcodeKey) + ':' + str(len(allBarcodeCollections[barcodeKey].readCollection)))
        
        # if barcode file has entries
        if (len(allBarcodeCollections[barcodeKey].readCollection) > 0 ):
            # calculate stats
            allBarcodeCollections[barcodeKey].summarizeSimpleReadStats()            
            
            
            #Bug here, testing...
            #print('Output a barcode file:')
            #print('barcodeKey:' + str(barcodeKey))
            #print('sampleid:' + str(sampleID))
            #print('referencefilelocation:' + str(referenceFileLocation))
            
            barcodedReadStats[barcodeKey] = allBarcodeCollections[barcodeKey].outputReadPlotsAndSimpleStats(outputDirectory, 'Barcode ' + str(barcodeKey), sampleID, referenceFileLocation)
        
        else:
            #print('Barcode:' + str(barcodeKey) + ' has no reads. Nothing to do.')
            pass
        
    # write the unbarcoded reads
    unbarcodedReadCollection.summarizeSimpleReadStats()
    barcodedReadStats['unbarcoded'] = unbarcodedReadCollection.outputReadPlotsAndSimpleStats(outputDirectory, 'Unbarcoded Reads', sampleID, referenceFileLocation)
    
    # TODO: Return a dictionary of read stats. This method will probably crash until i implement this.    
    # Every time i call "outputReadPlotsAndSimpleStats" I can put an entry in the dictionary.    
    #return barcodeList.keys()
    return barcodedReadStats


