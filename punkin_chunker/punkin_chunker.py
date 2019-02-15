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


from datetime import datetime
from os import makedirs, listdir
from os.path import split, isdir, join, splitext

from subprocess import Popen, PIPE, STDOUT, call
from operator import itemgetter

from Bio.SeqIO import parse as parseReads, write
from Bio.Blast.NCBIXML import parse as parseNCBIXML

from io import StringIO

from punkin_chunker.blast_result import blastResult
from nanopore_prospector.common import logMessageToFile, getConfigurationValue

from nit_picker.minion_read_collection import MinionReadCollection

from threading import Thread
from time import sleep

def logBlastProgress(blastLogText):
    # If we have a log file location from Global Variables, we can use that.
    # Otherwise print output to the console.
    if getConfigurationValue('log_file_location') is not None:
        logMessageToFile(blastLogText)
    else:
        print(str(blastLogText))

def createBlastDatabase(HLAReferenceFilename):
    #print ('Creating a blast database...')
    makeBlastDB_cline = ('makeblastdb' 
        + ' -in ' + HLAReferenceFilename 
#        + ' -parse_seqids -dbtype nucl')
        + ' -dbtype nucl')
    #print ('MakeDB Commandline:\n' + makeBlastDB_cline)
    makeDBOutput = call(makeBlastDB_cline, shell=True)
    
    if(makeDBOutput != 0):
        logBlastProgress('An error occured when creating the BLAST database.')
        raise Exception('Cannot create BLAST database')
    

# This method is a directory-safe way to open up a write file.
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not isdir(tempDir):
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput


def sortDirectory(readDirectory, outputDirectory, sortReference, threadCount):
    print ('Sorting this directory of reads:\n' + str(readDirectory))
    
    filesToAnalyze = []
    
    # Which files to process?
    for currentInputReadFile in listdir(readDirectory):
        if(getReadFormat(currentInputReadFile) == 'fastq'):
            # TODO: I don't think I'm naming a file "Unfiltered". I am using "All" reads...
            # Maybe Unfiltered makes more sense...
            if('_Rejected' in currentInputReadFile):
                print ('Skipping Rejected Read File:\n' + currentInputReadFile)
            elif('_Unfiltered' in currentInputReadFile):
                print('Skipping Unfiltered Read File:\n' + currentInputReadFile)
            else:
                # This is a "Pass" or "All" reads
                filesToAnalyze.append(currentInputReadFile)
        elif (getReadFormat(currentInputReadFile) == 'fasta'):
            # TODO Fix this, support fasta for BLAST sorting
            raise Exception('Bad Read Input Format, I dont support fasta in the blast sorter yet')
        else:
            print ('Skipping this file:' + currentInputReadFile)
            
    # Questionable logic:
    # If we have more than one file to analyze
    # That means that I definitely either demultiplexed
    # Or filtered based on length or quality.
    # So I will SKIP processing ALL reads. That will take forever. Still, keep the reads.
    if(len(filesToAnalyze) > 1):
        for currentInputReadFile in filesToAnalyze:
            if ('_All' in currentInputReadFile and '_Pass' not in currentInputReadFile):
                print ('Skipping analysis of ' + str(currentInputReadFile) + '')
                filesToAnalyze.remove(currentInputReadFile)
                
    # A simple dictionary value to return some results from this method.
    # the key is the name of the file analyzed.
    # The value is a list of minion_read_collections, pertaining to what gene they sorted to.    
    sortResults = {}

    for currentInputReadFile in filesToAnalyze:
        fullReadFileName=str(join(readDirectory,currentInputReadFile))
        #print ('loading Reads from:\n' + fullReadFileName)
        #splitext[0] will get the filename without the .fastq extension.
        outputSubDirectory = join(outputDirectory,splitext(currentInputReadFile)[0])
        print('Sorted Reads will go here:\n' + outputSubDirectory)
        
        sortResults[splitext(currentInputReadFile)[0]] = sortMinIONReadsByGene(fullReadFileName, outputSubDirectory, sortReference, threadCount)

    return sortResults



def getReadFormat(fullFilePath):
    # This method only works for fasta and fastq files. 
    # Ok really only fastq files.
    filename, file_extension = splitext(fullFilePath)

    # Trimming off the period here.
    if (".fasta" == file_extension or ".fa" == file_extension):
        return 'fasta'
    elif(".fastq"== file_extension or ".fq" == file_extension):
        return 'fastq'
    else :
        return None

 


def sortMinIONReadsByGene(inputFilename, outputDirectory, HLAReferenceFilename, threadCount):
    #print ('Time to sort our MinION Reads.  Compare against the APD Allele Reference.')

    shortFilename = split(inputFilename)[1]

    print ('HLA APD Reference:\n' + HLAReferenceFilename)
    print ('MinION Read Input file:\n' + inputFilename)
    print ('Output directory:\n' + outputDirectory)
    
    global FileOutputFormat
    if (".fasta" == inputFilename[-6:] or ".fa" == inputFilename[-3:]):
        FileOutputFormat = "fasta"
    elif (".fastq"== inputFilename[-6:] or ".fq" == inputFilename[-3:]):
        FileOutputFormat = "fastq"

    #sortResultsOutput = createOutputFile(join(outputDirectory ,
    #    (shortFilename + '.SortResults.txt'))) 
    
    #shortBlastOutput = createOutputFile(join(outputDirectory ,
    #    (shortFilename + '.ShortBlastOutput.txt'))) 
    
    finalBlastSummaryOutput = createOutputFile(join(outputDirectory ,
        (shortFilename + '.BlastSummary.txt'))) 
     
 
    #print ('Parsing input file.  It\'s format is ' + FileOutputFormat)
    parsedReads = parseReads(inputFilename, FileOutputFormat)
    #minionReadRecords = enumerate(parsedReads)
    #readRecordList = list(minionReadRecords)
    readRecordList = list(parsedReads)

    readCount = len(readRecordList)
    print (str(readCount) + ' reads ready to sort.')
    
    #for i in range(0,10):
    #    print ('index ' + str(i) + 'readrecord at this position is:' + str(readRecordList[i]))
  
    # pause statement for debugging...
    #wait = raw_input("\n\n\n\n**************************PRESS ENTER TO CONTINUE.**************************\n\n\n\n")

    # Print some info to the sorting summary
    finalBlastSummaryOutput.write('HLA APD Reference:' + HLAReferenceFilename + '\n')
    finalBlastSummaryOutput.write('MinION Read Input file:' + inputFilename + '\n')
    finalBlastSummaryOutput.write('Output directory:' + outputDirectory + '\n')
    finalBlastSummaryOutput.write('Read Count: ' + str(readCount) + '\n')

    # Create the Blast Reference before we start.
    createBlastDatabase(HLAReferenceFilename)

    # Split the reads into batches.  How many batches? Same as thread count.
    # Since I'm rounding, I'll add one to the read count so my batch size is right. I'm worried about forgetting a read.
    print ('thread count is ' + str(threadCount))
    batchSize = round(readCount / threadCount) + 1

    
    readBatches = splitReadsIntoBatches(readRecordList, batchSize)
    # After splitting I'll clear the old value, because memory management. This is probably unnecessary, compilers are smart.
    readRecordList = None
    
    print ('I split the ' + str(readCount) + ' reads into ' 
        + str(len(readBatches)) + ' batches of <=' + str(batchSize) + ' reads.')

    threadHandles = []
    threadResults = {}
    combinedResultSets = []
    
    
        # Start a thread to process each batch of reads.
    # I am using 1-based indexing to identify my threads. Why? I just needed some consistency...
    # So far I haven't run into simultaneous access problems. 
    # Either with the BLAST database, or with the RESULTS dictionary.
    for batchIndex, readBatch in enumerate(readBatches):
        print ('Sorting Batch # (' + str(batchIndex + 1) + '/' + str(len(readBatches)) + ') ' + str(datetime.now()))
        
        currentThread = blastThread(batchIndex + 1, readBatch, HLAReferenceFilename, threadResults)
        currentThread.start()
        threadHandles.append(currentThread)

    # Wait for all the threads to finish.
    for currentThread in threadHandles:
        print ('Waiting for thread to finish:' + currentThread.name)
        currentThread.join()
    print ('All the blast threads are done.')
    
    # Combine the results for each thread.
    print ('Time to check the results, this many keys in the results dict:' + str(len(threadResults.keys())))
    for key in threadResults.keys():
        print ('Thread #' + str(key) + ' had ' + str(len(threadResults[key])) + ' BLAST results.')
        combinedResultSets = combinedResultSets + threadResults[key]



    #blastResultSets = sortReadArrayByGene(readRecordList, HLAReferenceFilename, sortResultsOutput, shortBlastOutput)
    # The result is an array of minion_read_collection objects.
    blastResults = writeSortedReadsToFile(combinedResultSets, outputDirectory, finalBlastSummaryOutput)
    
        
        
    #printSortedFastaFiles(blastResultSets, outputDirectory, finalBlastSummaryOutput)
    #sortResultsOutput.close()
    #shortBlastOutput.close()
    
    finalBlastSummaryOutput.close()
    
    # Return results for intepretation downstream.
    return blastResults

    


def splitReadsIntoBatches(recordList, batchSize):
    #print ('I have ' + str(len(recordList)) + ' reads to split into batches of size ' + str(batchSize))
    # remainderRecords keeps track of the "unbatched" reads.
    remainderRecords = recordList
    # readBatches is an array full of read-arrays.  Sort of a 2d array but not quite.
    readBatches = []
        
    while len(remainderRecords) > 0:
        currentBatch = remainderRecords[0:batchSize]
        readBatches.append(currentBatch)
        remainderRecords = remainderRecords[batchSize:]

    return readBatches

def sortReadArrayByGene(minionReadRecordList, HLAReferenceFilename, threadCount):
    blastResultSet = []

    # I group the reads into a single fasta string to pass to BLAST  
    fastaString = ''
        
    for minionReadIndex, minionReadRecord in enumerate(minionReadRecordList):
        currentReadID = str(minionReadRecord.id)
        currentSequence = str(minionReadRecord.seq) 
        
        fastaString += '>' + currentReadID + '\n' + currentSequence + '\n'

    # TODO: What if blastn is not installed?
    # TODO: this thread count is not really working.  
    # TODO: I should do real threading by myself.  
    commandLineArray = ['blastn'
        , '-db' , HLAReferenceFilename               
        , '-num_threads', str(threadCount)         
        , '-outfmt', '5'
        , '-evalue', '0.001'             
        ]
    
    #print('COMMANDLINEARRAY')
    #print(str(commandLineArray))
    
    blastProcessHandle = Popen(commandLineArray, stdout=PIPE, stdin=PIPE, stderr=STDOUT)    
    #grep_stdout = p.communicate(input=fastaString)[0]

    blastXMLText = blastProcessHandle.communicate(input=fastaString.encode())[0].decode()
    
    #print ('resulting blast xml text:')
    #print(blastXMLText)
    
    blastResults = parseXMLForBlastResults(minionReadRecordList,blastXMLText)
    
    # I don't use append, because, blastResults might have more than one minionReadRecord.  Might.
    blastResultSet = blastResultSet + (blastResults)
        
        
    #print('I am done, and I need to do a final check.')
    #print('Length MinION Record List: ' + str(len(minionReadRecordList)))
    #print('Length MblastResultSet: ' + str(len(blastResultSet)))



    return blastResultSet 



           

def parseXMLForBlastResults(readRecords, blastXMLText):
    #Parse XML for blast results. Need a File handle first.
    xmlStringHandle = StringIO(blastXMLText)    
    blastRecords = parseNCBIXML(xmlStringHandle)
    
    blastResults = []
    
    # parse returns multiple blast results.
    for blastRecordIndex, blastRecord in enumerate(blastRecords):  
        
        currentBlastResult = blastResult()   
        currentBlastResult.readRecord = readRecords[blastRecordIndex]
        
        if(len(blastRecord.alignments) < 1):
            #print ('No alignments detected for this read. This is not a problem.')
            pass

        else:
            currentBlastResult.processBlastRecord(blastRecord)
                    
        blastResults.append(currentBlastResult)  

    #print ('Returning Blast Results. Out of ' + str(len(readRecords)) + ' reads i found ' + str(len(blastResults)) +  ' blast result records.')

    return blastResults



#I think this method is deprecated.
def printBlastRecordInformation(blastRecord):
    
    print ('BlastRecord:' + str(blastRecord))
    #print ('this record aligned with: ' + blast)
    
    currentAlignments = blastRecord.alignments
    print ('Number of alignments:' + str(len(currentAlignments)))
    for alignment in currentAlignments:
        print ('****Alignment****')
        print ('sequence:' + str(alignment.title))
        print ('length:' + str(alignment.length))
        # Usually there is only one HSP.
        # hsp stores info about the alignment
        currentHsps = alignment.hsps
        if (len(currentHsps) != 1):
            print ('I only expect exactly one hsp for this alignment. There are ' + str())
            raise Exception('Multiple HSPs. I don\'t know what this means but you have to fix it.')
        #print ('Number of hsps:' + str(len(currentHsps)))
        # I don't need to loop these. Just grab the first.
        for hsp in currentHsps: #print('e value:' + str(hsp.expect))
            print ('score:' + str(hsp.score))
            print ('strand:' + str(hsp.strand))
            print ('frame:' + str(hsp.frame))
            print (hsp.query[0:75] + '...')
            print (hsp.match[0:75] + '...')
            print (hsp.sbjct[0:75] + '...')
            
# readRecordList contains a list of the reads.
# blastResultSets is a list of the blastResults.
# TODO: I really should sync those together.  Maybe put the read record right into the blast result object.  Yeah do that.
# Because It's a pain to keep the two arrays in sync.   With threading it's much harder.
# I Think i did this, but i'll have to watch to see what i broke.
def writeSortedReadsToFile(blastResults, outputDirectory, finalBlastSummaryOutput):    
    print('I will now write the sorted reads to file.')

    # An array of readCollections is enough.
    geneLevelReadCollections = []

    
    global FileOutputFormat

    unsortedReadCollection = MinionReadCollection([])
    unsortedReadCollection.readInputFormat = FileOutputFormat
    unsortedReadCollection.gene = None

    for resultIndex, currentBlastResult in enumerate(blastResults):
        
        #print ('checking this read:' + str(currentBlastResult.readRecord.id))
        
        # if nomenclature fields is empty, this blast result was not assigned to a gene.
        if(len(currentBlastResult.NomenclatureFields) < 1):
            # This corresponds to if there was no blast results for this read. 
            
            #unsortedReadCount += 1
            #write([currentBlastResult.readRecord], unsortedReadOutputFile, FileOutputFormat)
            
            
            #print ('One unsorted read found:' + str(currentBlastResult.readRecord.id))
            
            unsortedReadCollection.readCollection.append(currentBlastResult.readRecord)
            #unsortedReadCollection.readInputFormat = FileOutputFormat
               
        else:
            currentGene = currentBlastResult.Gene
            #currentGroup = currentBlastResult.AssignedAlleleGroup
            
            #print ('This read should belong to this gene:' + str(currentBlastResult.readRecord.id)
            #    + '\n' + currentBlastResult.Gene)
    
            # If the gene is assigned   
            if(currentGene is not None and len(currentGene) > 0 and currentGene != '-1'):
         
                # Print the read to a group level             
                #currentGeneLevelOutputFile = None
                
                # Search for existing Gene level outputffiles.
                foundGeneLevelOutput = False

                
                # If we already have a gene level output file in our list, use that one.
                # I keep track of the outputFileIndex on my own, because enumerate() is smarter than I am.
                #outputFileIndex = 0
                #for outputFileGene, outputFileObject, readCount in geneLevelOutputFiles:
                for readCollection in geneLevelReadCollections:
                    #print ('outputFileObject is this:' + str(outputFileObject))
                    #print ('outputFileObject[0] is this:' + str(outputFileObject[0]))
                    if readCollection.gene == currentGene:
                        #foundGeneLevelOutput = True
                        #currentGeneLevelOutputFile = outputFileObject
                        #readCollection.readCollection.append(currentBlastResult.readRecord)
                        foundGeneLevelOutput = True
                        # Increment Read Count
                        #geneLevelOutputFiles[outputFileIndex] = (
                        #    outputFileGene,
                        #    outputFileObject,
                        #    readCount + 1)
                        if (currentBlastResult.ForwardMatch):
                            readCollection.readCollection.append(currentBlastResult.readRecord)
                            #write([currentBlastResult.readRecord], currentGeneLevelOutputFile, FileOutputFormat)
                        else:
                            #print ('about to check the reverse complement of this read:\n')
                            #print (str(currentBlastResult.readRecord))
                            reverseRecord = currentBlastResult.readRecord
                            forwardRecord = reverseRecord.reverse_complement(id=reverseRecord.id+"_reverse_complement", name=True, description=True)
                            #SeqIO.write([rec], allReadFileOutputs[barcodeKey], 'fastq')
                    
                            #print ('This record is in the reverse direction.')
                            #write([forwardRecord], currentGeneLevelOutputFile, FileOutputFormat)
                            readCollection.readCollection.append(forwardRecord)
                        
                    
                    #outputFileIndex += 1
                # None found, make a new output file.
                #if currentGeneLevelOutputFile is None:
                #    currentGeneLevelOutputFile = createOutputFile(join(outputDirectory,'HLA-' + currentGene + '.' + FileOutputFormat))
                #    geneLevelOutputFiles.append((currentGene,currentGeneLevelOutputFile,1))
                if not foundGeneLevelOutput:
                    newGeneLevelReadCollection = MinionReadCollection([currentBlastResult.readRecord])
                    newGeneLevelReadCollection.gene = currentGene
                    newGeneLevelReadCollection.readInputFormat = FileOutputFormat
                    geneLevelReadCollections.append(newGeneLevelReadCollection)
                
                        
                # Print the sequence to the Gene level output.
                #if currentGeneLevelOutputFile != None:
                    # Is the sequence mapped to the reference in the reverse direction?
                    
                #    if (currentBlastResult.ForwardMatch):
                #        write([currentBlastResult.readRecord], currentGeneLevelOutputFile, FileOutputFormat)
                #    else:
                        #print ('about to check the reverse complement of this read:\n')
                        #print (str(currentBlastResult.readRecord))
                 #       reverseRecord = currentBlastResult.readRecord
                 #       forwardRecord = reverseRecord.reverse_complement(id=reverseRecord.id+"_reverse_complement", name=True, description=True)
                        #SeqIO.write([rec], allReadFileOutputs[barcodeKey], 'fastq')
                
                        #print ('This record is in the reverse direction.')
                  #      write([forwardRecord], currentGeneLevelOutputFile, FileOutputFormat)
             
            #else:
                # In this case we cant find a gene. It's unsorted.
             
             
             
                    
    # New loop for writing objects.
    # Write all reads. 
    for readCollection in geneLevelReadCollections:
        readCollection.summarizeSimpleReadStats()
        readCollection.outputReadPlotsAndSimpleStats(
            outputDirectory
            , 'Reads'
            , 'HLA-' + str(readCollection.gene)
            , None)
        
    # Unsorted reads
    unsortedReadCollection.summarizeSimpleReadStats()
    unsortedReadCollection.outputReadPlotsAndSimpleStats(
        outputDirectory
        , 'Reads' 
        , 'unsorted'
        , None)
    
     
        
    # Write sort results to the output file and close the Gene and Group specific fasta files.    
    finalBlastSummaryOutput.write('\n\nSorting Read Results:\n')

    #for outputFileObject in geneLevelOutputFiles:
    #for outputFileObject in sorted(geneLevelOutputFiles, key=itemgetter(2), reverse=True):
    for readCollection in geneLevelReadCollections: 
        finalBlastSummaryOutput.write('HLA-' + str(readCollection.gene) + ': ' + str(len(readCollection.readCollection)) + ' Reads\n')
        #outputFileObject[1].close()

    finalBlastSummaryOutput.write('Unsorted Reads: ' + str(len(unsortedReadCollection.readCollection)) + ' Reads\n')
    #finalBlastSummaryOutput.close()
    #unsortedReadOutputFile.close()
    
    # Return an array with results. I'll add the unsorted reads to the collection before returning. 
    geneLevelReadCollections.append(unsortedReadCollection)
    return geneLevelReadCollections

class blastThread (Thread):
    def __init__(self, id, readBatch, HLAReferenceFilename, threadResults):
        Thread.__init__(self)
        #self.threadID = threadID
        self.id = id
        self.readBatch = readBatch
        self.HLAReferenceFilename = HLAReferenceFilename
        self.threadResults = threadResults
        
    def run(self):
        print ('Starting ' + self.name)
        # Get lock to synchronize threads
        #acquire()
        #print_time(self.name, self.counter, 3)
        
        sleep(5)
        
        currentResults = sortReadArrayByGene(self.readBatch, self.HLAReferenceFilename, 1)
        print ('Thread ' + self.name + ' found ' + str(len(currentResults)) + ' blast results.')
        print ('Storing them in the results dictionary...')
        self.threadResults[self.id] = currentResults
        #
        
        # Store the results in our global dictionary
        #threadResults[self.id] = currentResults
        #currentResults = sortReadArrayByGene(readBatch, HLAReferenceFilename, threadCount)
       
        # Free lock to release next thread
        #release()
        
        print ('Finishing ' + self.name)

#def print_time(threadName, delay, counter):
#    while counter:
#        sleep(delay)
#        print ("%s: %s" % (threadName, ctime(time())))
#        counter -= 1
    

