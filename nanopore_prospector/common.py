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

from pylab import clf, scatter, xlabel, ylabel, title, savefig
from random import shuffle

from os.path import split, isdir, isfile, join, expanduser, abspath, exists, splitext
from os import makedirs, system, walk

from xml.etree import ElementTree as ET
from xml.dom import minidom as MD

from Bio.SeqIO import parse as parseRecords, write as writeRecords

from datetime import datetime

from pysam import AlignmentFile


def alignReads(referenceLocation, readFileLocation, alignmentOutputDirectory, useReadSubset):
    # TODO: I just removed the last two parameters from this method: , numberThreads, minimumAlignmentLength. They weren't being used.
    # That will probably break things.

    # TODO: Consume the number of threads, pretty sure I can pass this into minimap2.
    # TODO: Maybe use #threads+1 because minimap2 always uses one thread for IO during alignment. Double check that.

    # TODO: Pass in a boolean for whether or not to exclude short alignments. Short alignments are problematic.
    # TODO: They align because (i think) of false homology within the gene. It looks like a hump of coverage with lots of false SNPs.
    # TODO: The solution is in the -m parameter in minimap2. This is a minimum length for smaller aligned chains.
    # TODO: -m 3000 worked for a reference length of about 5600, so maybe 1/2 reference length is a good guess.

    # TODO: this is not a good solution. The major problem is With DP, many of my reads are too short compared \
    # to reference. New solution is to pass a minimum alignment length. i try that, if I have a problem with false, unaligned reads, this is why.
    # If my parameters are wrong, check that I'm passing a minimumReadLength, not a boolean to eliminate the shortReads


    # perform minimap alignment to align reads against a reference.
    # This method returns the # of reads that aligned to this reference.
    print('\nAligning reads against the reference:' + str(referenceLocation) + ' & ' + str(readFileLocation))

    if not isdir(alignmentOutputDirectory):
        makedirs(alignmentOutputDirectory)

    # Part 1 Index the Reference
    try:
        # Copy the reference sequence to the alignment directory. This is a complicated way to do it.
        newReferenceLocation = join(alignmentOutputDirectory, 'AlignmentReference.fasta')
        refSequence = list(parseRecords(referenceLocation, 'fasta'))[0]
        refSequence.id = 'AlignmentReference'
        sequenceWriter = createOutputFile(newReferenceLocation)
        writeRecords([refSequence], sequenceWriter, 'fasta')
        sequenceWriter.close()

        # Index The Reference
        #referenceIndexName = newReferenceLocation.replace('.fasta', '.mmi')
        #cmd = ('minimap2 -d ' + referenceIndexName + ' ' + newReferenceLocation)

        referenceIndexName = newReferenceLocation + '.lastindex'
        cmd = 'lastdb -Q 0 ' + referenceIndexName + ' ' + newReferenceLocation
        system(cmd)

    except Exception:
        print ('Exception indexing alignment reference. Is bwa installed? folder writing permission issue?')
        raise

    # TODO: Make this a commandline parameter.  Lower = faster. Higher = more accurate consensus correction
    # Change the parameter from "usesubset" to "max read count"...something like that.
    alignmentReadSubsetCount = 5000
    try:
        if useReadSubset:
            print ('Using a read Subset')
            # load Reads
            parsedReads = list(parseRecords(readFileLocation, getReadFileType(readFileLocation)))

            # If there aren't enough reads for this
            if (len(parsedReads) < alignmentReadSubsetCount):
                alignmentReadSubsetCount = len(parsedReads)

            # choose random subset
            randomIndexes = list(range(0, len(parsedReads)))
            shuffle(randomIndexes)
            sampledReads = []
            for i in range(0, alignmentReadSubsetCount):
                sampledReads.append(parsedReads[randomIndexes[i]])

            # write random reads to alignment directory
            # Reassign the reads we'll use downstream
            readFileLocation = join(alignmentOutputDirectory, 'ReadSample.fasta')
            readSampleWriter = createOutputFile(readFileLocation)

            writeRecords(sampledReads, readSampleWriter, 'fasta')
            readSampleWriter.close()

        else:
            # We'll use the whole read file.
            print ('Using All reads, not a subset.')
            pass

    except Exception:
        print ('Exception selecting a read subset.')
        raise

    # Part 2 Align
    try:
        alignmentOutputName = join(alignmentOutputDirectory, 'alignment.bam')


        # Parameters depend on "read type"
        # if i have a fastq, i assume these are minion reads.
        # if i have a fasta, i assume these are consensus sequences.

        if (getReadFileType(readFileLocation) == 'fasta'):
            #minimapParams = '-ax asm20'
            lastargs = "-s 2 -T 0 -Q 0 -a 1" # These parameters work for Last.

            #LAST manual says that -D100 might work for shorter alignments. This experiment isn't working.
            #lastargs = "-s 2 -T 1 -Q 0 -a 1 -D 10 -e 5 -d 5"

            #lastargs = "-s 2 -T 0 -Q 0 -a 1 -b 1 -r 1" # Experimental Args for DQA, because DQA1 is full of structural differences.

            # Penalizing Insertions heavily. Extension is still a low penalty. I want to minimize insertions.for consensus sequences. 3 seems to be the sweet spot.
            #lastargs = "-s 2 -T 0 -Q 0 -a 1 -A 3"
        elif (getReadFileType(readFileLocation) == 'fastq'):
            #print ('attempting the consensus parameters instead of ONT settings.')
            #minimapParams = '-ax map-ont'
            lastargs = "-s 2 -T 0 -Q 1 -a 1"
        else:
            raise Exception('Unknown read file type....')

        #if(excludeShortAlignments):
        #if(minimumAlignmentLength is not None):
            # TODO: I'm using -m as a parameter in minimap for minimum bases for aligned reads. Maybe 1/3 length of the reference is a good default value
            # TODO: I should set the file type based on a parameter.
            # TODO: Scrapped that Idea, just going with minimum alignmentLength.
            # TODO: Ok Actually -m is a parameter for "minimum chaining score" in minimap2. This surprises me, i think the parameter was different before.
            # Good to know.
            # TODO: Figure out another minimum alignment length somehow. Short alignments are still getting lost. The minimumAlignmentLength is not actually being used.
            #minimumReadValue = int(round(len(refSequence.seq) / 2))
            #minimapParams += ' -m ' + str(minimumReadValue)
            #minimapParams += ' -m ' + str(minimumAlignmentLength)
        #    pass

        #else:
        #    pass

        # minimap2 is actually pretty bad in my opinion. I'm not using it because I have so many problems with misaligned reads.
        #cmd = ("minimap2 " + minimapParams + " " +
        #       referenceIndexName + " " +
        #       readFileLocation +
        #       " | samtools view -b | samtools sort -o "
        #       + alignmentOutputName)

        # TODO: im calling python2 here, think about a better way to do it.
        lastTxtAlignmentLocation = join(alignmentOutputDirectory, 'LastAlignment.last.txt')
        cmd = ("lastal " + lastargs + " " +
            referenceIndexName + " " +
            readFileLocation + " | python2 /usr/bin/last-map-probs > " +
            lastTxtAlignmentLocation)
        system(cmd)

        # Convert to sam
        samAlignmentName = join(alignmentOutputDirectory, 'LastAlignment.sam')
        cmd = "python2 /usr/bin/maf-convert sam " + lastTxtAlignmentLocation + " > " + samAlignmentName
        system(cmd)


        cmd = "samtools view -T " + newReferenceLocation + " -bS " + samAlignmentName + " | samtools sort -o " + alignmentOutputName
        system(cmd)


        # Convert to bam


        #print ('RunningThisCommand:\n' + cmd + '\n')




              # " | samtools view -b | samtools sort -o "
              # + alignmentOutputName)



        # print ('alignment command:\n' + cmd)


    except Exception:
        print ('Exception aligning reads against reference. Are minimap2 and samtools installed?')
        raise







    # Part 3 Index Alignment
    try:
        cmd = ("samtools index " + alignmentOutputName)
        # print ('alignment index command:\n' + cmd)
        system(cmd)
        # print ('index command:\n' + cmd)
    except Exception:
        print ('Exception indexing alignment reference. Is bwa installed?')
        raise

    alignedReadCount = calculateTotalCoverage(alignmentOutputName)
    return alignedReadCount

def getDirectoryAlignmentStats(analysisDirectory, outputFileName, recursive):
    # Search a directory for bam files. for each bam file, calculate some simple statistics.
    # Most important is total coverage.

    print('Searching ' + str(analysisDirectory) + ' for .bam files.')

    # Create an output file.
    statOutputFile = createOutputFile(outputFileName)
    statOutputFile.write('AlignmentFile,ReadCoverage\n')

    # Loop Bam files
    for root, directories, filenames in walk(analysisDirectory):
        for directory in directories:
            print('Searching ' + str(join(root, directory)) + ' for .bam files.')
            pass
        for filename in filenames:
            fullFileNamePath = join(root, filename)

            if fullFileNamePath.endswith('.bam'):
                readCoverage = calculateTotalCoverage(fullFileNamePath)
                statOutputFile.write(fullFileNamePath + ',' + str(readCoverage) + '\n')


    print('Done. Results were written to ' + outputFileName)

    # Close output file if necessary.
    statOutputFile.close()



def calculateTotalCoverage(alignmentOutputName):
    # The number of reads that are involved in this alignment
    bamfile = AlignmentFile(alignmentOutputName, 'rb')
    alignedReads = list(bamfile.fetch())
    return len(alignedReads)

def getReadFileType(fileNameFull):
    # This method is for deciding if a file is fasta or fastq.
    # It can be extended for other file types.
    fileName, fileExtension = splitext(fileNameFull)
    
    if(fileExtension == '.fasta' or fileExtension == '.fa'):
        return 'fasta'
    elif(fileExtension == '.fastq' or fileExtension == '.fq'):
        return 'fastq'
    else:
        raise Exception('I do not understand this file type: ' + str(fileNameFull))


def createScatterPlot(graphTitleText, xValues, yValues, xAxisName, yAxisName, outputFileName):
    print('Creating a Scatter Plot: ' + outputFileName)
      
    #Clear the figure and start anew.
    clf()
    
    # K is black, we need to repeat N times.
    colors = ['K'] * len(xValues)

    scatter(xValues, yValues, s=1, c=colors, marker='.')
    
    xlabel(xAxisName)
    ylabel(yAxisName)
    title(graphTitleText)
    
    savefig(outputFileName)

# This method is a directory-safe way to open up a write file.

def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not isdir(tempDir):
        print('Making Directory:' + tempDir)
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput     


# I'm storing global variables in a dictionary for now. 
def initializeGlobalVariables():    
    global globalVariables 
    
    if not ("globalVariables" in globals()):
        globalVariables={}
        
def assignConfigurationValue(configurationKey, configurationValue):
    # Overwrite config value without question.
    initializeGlobalVariables()
    globalVariables[configurationKey] = configurationValue
    
def assignIfNotExists(configurationKey, configurationValue):
    # Use this assigner if we want to declare important, new configuration values.
    # Using this method, we will not overwrite custom values
    # But we will provide critical new config values.
    initializeGlobalVariables()
    if configurationKey not in globalVariables.keys():
        assignConfigurationValue(configurationKey, configurationValue)
    
    
def getConfigurationValue(configurationKey):
    initializeGlobalVariables()
    if configurationKey in globalVariables.keys():
        return globalVariables[configurationKey]
    else:
        #logEvent ('Configuration Key Not Found:' + configurationKey)
        #raise KeyError('Key Not Found:' + configurationKey)
        return None


def logMessageToFile(currentMessage):
    # The log keeps track of step-by-step things performed by Nanopore Prospector.
    fullLogMessage = str(datetime.now()) + ' : ' + currentMessage + '\n'
    
    if (getConfigurationValue('log_file_location') is None):
        print ('No Logfile exists yet.\nTrying to log this message:\n' + currentMessage)
    else:
        #print ('Logging message: ' + currentMessage) 
        # Append the log with the log entry            
        resultsOutput = open(getConfigurationValue('log_file_location'), 'a')
        resultsOutput.write(fullLogMessage)
        resultsOutput.close()
    
def getBlastSortResourceLocation(resourceName): 
    sortReferencePath = None
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        # TODO: Check the spec file for pyinstaller, make sure it is including all the blast references
        sortReferencePath = join(sys._MEIPASS, resourceName)
        #print('using packaged MEIPASS directory to find barcode.')
    except Exception:
        sortReferencePath = join(abspath('.'),'hla_blast_references/' + resourceName)
        #print('using local path to barcode kit file.')
    return sortReferencePath
    
def getBarcodeFilePath():
    # TODO: I should add some sort of option, in case they are using the 12x barcode kit for some reason
    # For now, I assume they will never do that.
    barcodeFilePath = None
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        barcodeFilePath = join(sys._MEIPASS, 'barcodes_96_bc_kit.txt')
        #print('using packaged MEIPASS directory to find barcode.')
    except Exception:
        barcodeFilePath = join(abspath('.'),'barcodes/barcodes_96_bc_kit.txt')
        #print('using local path to barcode kit file.')

    return barcodeFilePath

def combineBlastReferences(blastReferenceFileList, blastOutputDirectory):
    if(len(blastReferenceFileList) == 0):
        logMessageToFile('Blast Reference File List is empty. I will not sort the reads.')
        return None
        
    blastReferenceRecords = []
    
    for blastRefFile in blastReferenceFileList:
        # Open records.
        print ('Opening a blast reference file:' + str(blastRefFile))
        
        currentReferenceSequences = parseRecords(blastRefFile, 'fasta')
        for record in currentReferenceSequences:
            # add records to combined list
            blastReferenceRecords.append(record)
        
    # write list to a new blast reference file.
    blastReferenceFileName = join(blastOutputDirectory, 'HLA_Blast_Sort_Reference.fasta')
    if not exists(blastOutputDirectory):
        makedirs(blastOutputDirectory)      
        
    logMessageToFile('Writing ' + str(len(blastReferenceRecords)) + ' Blast sort records to file:' + str(blastReferenceFileName))  
    writeRecords(blastReferenceRecords, createOutputFile(blastReferenceFileName), 'fasta')
        
    # return the location of this file.
    return blastReferenceFileName
    
        
def writeConfigurationFile():
    assignIfNotExists('config_file_location',
        join(join(
            expanduser("~"),'nanopore_prospector_temp'),'Nanopore.Prospector.Config.xml'))

    root = ET.Element("config")
    
    for key in globalVariables.keys():
        # Some config values I don't want to store.
        if(key not in [
            'log_file_location'
            ,'input_directory'
            ,'output_directory'
            ]):
            
            #print ('assigning this config value to the file:' + key + ':' + str(globalVariables[key]))
            ET.SubElement(root, key).text = str(globalVariables[key])

    xmlText = ET.tostring(root, encoding='utf8', method='xml')
    prettyXmlText = MD.parseString(xmlText).toprettyxml()
    
    xmlOutput = createOutputFile(globalVariables['config_file_location'])
    xmlOutput.write(prettyXmlText)
    xmlOutput.close()

def loadConfigurationFile():
    assignIfNotExists('config_file_location',
        join(join(
            expanduser("~"),'nanopore_prospector_temp'),'Nanopore.Prospector.Config.xml'))
    
    if not isfile(globalVariables['config_file_location']):
        print ('The config file does not exist yet. I will not load it:\n' + globalVariables['config_file_location'])
        
    else:
        print ('The config file already exists, I will load it:\n' + globalVariables['config_file_location'])
        
        tree = ET.parse(globalVariables['config_file_location'])
        root = tree.getroot()
        
        for child in root:
            assignConfigurationValue(child.tag, child.text)
            
    # TODO: Im assuming we always want an analysis log.
    #assignIfNotExists('logging','0')
    
    assignIfNotExists('demultiplex_reads', 2)    
     
    assignIfNotExists('min_length', 1)
    assignIfNotExists('max_length', 1000000)
    
    assignIfNotExists('min_quality', 0)
    assignIfNotExists('max_quality', 100)
    
    assignIfNotExists('analyze_hla_a', 0)
    assignIfNotExists('analyze_hla_b', 0)
    assignIfNotExists('analyze_hla_c', 0)
    assignIfNotExists('analyze_hla_e', 0)
    assignIfNotExists('analyze_hla_dra', 0)
    assignIfNotExists('analyze_hla_dqa1', 0)
    assignIfNotExists('analyze_hla_drb1', 0)
    
    