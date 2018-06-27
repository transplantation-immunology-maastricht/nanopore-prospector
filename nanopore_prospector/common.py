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

from os.path import split, isdir, isfile, join, expanduser, abspath, exists, splitext
from os import makedirs

from xml.etree import ElementTree as ET
from xml.dom import minidom as MD

from Bio.SeqIO import parse as parseRecords, write as writeRecords

from datetime import datetime

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
    
    