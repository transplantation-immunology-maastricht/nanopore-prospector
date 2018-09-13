#from datetime import datetime
from os import listdir#, makedirs, 
#from os.path import split, isdir, join, splitext
from os.path import isdir, isfile, join, split

#from subprocess import Popen, PIPE, STDOUT, call
#from operator import itemgetter

from Bio.SeqIO import parse as parseReads, write

# copy2 preserves the file characteristics, I don't know if this is good or bad.
from shutil import copy, copy2

#import h5py
#h5py.run_tests()
#from Bio.Blast.NCBIXML import parse as parseNCBIXML

#from io import StringIO

#from punkin_chunker.blast_result import blastResult
#from nanopore_prospector.common import logMessageToFile, getConfigurationValue

#from nit_picker.minion_read_collection import minionReadCollection

#from threading import Thread
#from time import sleep


# The purpose of this script is to look at a fastq file from MinION
# And collect the corresponding fast5 reads. We don't want ALL the fast5 reads, we just want to
# use the fast5 reads from a single fastq file.
# This is used for downstream nanopolish.
def collectReads(fastQReadFile, fast5ReadDir, outputDir):
    #print ('Inside the collect reads method.')
    
    fast5ReadCollection = parseReadDirectory(fast5ReadDir, True)
    shortNameReadCollection = makeShortNames(fast5ReadCollection)
    
    print (str(fast5ReadDir) + ' has ' + str(len(fast5ReadCollection)) + ' files and directories.')

    fastqReadCollection = list(parseReads(fastQReadFile, format='fastq'))    
    fastqReadCount = len(fastqReadCollection)
    print (str(fastQReadFile) + ' has ' + str(fastqReadCount) + ' reads.')
    
    for fastqIndex, fastqRead in enumerate(fastqReadCollection):
        
        # Some progress.
        if (fastqIndex == 0
            or fastqIndex == fastqReadCount - 1
            or (fastqIndex + 1) % 500 == 0):
            print('Collecting fast5 reads for nanopolish: ' 
                  + str(fastqIndex + 1) + '/' + str(fastqReadCount)
                  + ' ' + str(round((100 * (fastqIndex + 1) / fastqReadCount),3)) + '%')

        readDescription = fastqRead.description
        # A fast5 read looks like this:
        # COM0007921_20180105_FAH35044_MN17278_sequencing_run_M180105A_R95_Validatie_AB_126_36732_read_19_ch_303_strand.fast5
        # A description might look like this:
        # f7be2811-7673-4a16-98dc-f910bbdf6622d7e27c0f-fbbf-4dcc-96f7-078f71998f1b runid=81539cea26a947c3a1a44ea2cc07be4ef46299df read=28192 ch=297 start_time=2018-01-05T23:32:54Z
        # I want to pull out the 'read=28192 ch=297' information. This is unique to the read.
        # I search bsed on the text with read and channel number, like this: _19_ch_303_
        
        readNumberBeginIndex = readDescription.rfind(' read=') + 6
        readNumberEndIndex =  readDescription.rfind(' ch=')
        chanNumberBeginIndex = readDescription.rfind(' ch=') + 4
        chanNumberEndIndex =  readDescription.rfind(' start_time=')
        
        #print ('fastq descr:' + readDescription)
        readNumber = readDescription[readNumberBeginIndex:readNumberEndIndex]
        chanNumber = readDescription[chanNumberBeginIndex:chanNumberEndIndex]
        readFindingText = '_' + readNumber + '_ch_' + chanNumber + '_'
        #print ('read looks like this:' + readFindingText)
        
        
        # Look in our fast5 reads to spot the read that corresponds.
        for fast5Index, shortReadName in enumerate(shortNameReadCollection):
            if(shortReadName == readFindingText):
                #print ('I found a read that matches!')
                #print ('short read name:' + shortReadName)
                #print ('long read name:' + fast5ReadCollection[fast5Index])
                # Copy file to output file.
                #print('I will copy the file to the output directory: ' + str(outputDir))
                path, file = split(fast5ReadCollection[fast5Index]) 
                newFileName = join(outputDir, file)
                
                
                #print ('copying file from\n')
                
                copy(fast5ReadCollection[fast5Index], newFileName)
                # We found the read, we can break the fast5 read loop. Should cut time in half roughly.
                break
            
                
        
        

    
def makeShortNames(longReadNames):    
    # The idea is to trim off as much information as we can.
    # The Read # and Channel # are important, these should be unique, I think.
    shortReadNames=[]
    for longReadName in longReadNames:
        beginIndex = longReadName.rfind('_read_') + 5
        endIndex = longReadName.rfind('strand.fast5')
        shortReadName = longReadName[beginIndex:endIndex]
        shortReadNames.append(shortReadName)
        #print ('long name:' + longReadName)
        #print ('short name:' + shortReadName)
    return shortReadNames

def parseReadDirectory(readDirectory, parseSubdirectories):
    # This method is recursive. 
    # We want to make a list of all fast5 files in this directory, and subdirectories
    readFiles = []
    
    filesAndDirectories = sorted(listdir(readDirectory))
    #print (str(readDirectory) + ' has ' + str(len(filesAndDirectories)) + ' files and directories.')
    
    for fileOrDir in filesAndDirectories:
        fullPath = join(readDirectory,fileOrDir)
        if(isdir(fullPath)):
            #print (fullPath + ' is a directory.')
            if(parseSubdirectories):
                subdirFiles = parseReadDirectory(fullPath, parseSubdirectories)
                readFiles = readFiles + subdirFiles
                #print (fullPath + ' has ' + str(len(subdirFiles)) + ' files.')
        elif(isfile(fullPath)):
            #print (fullPath + ' is a file.')
            readFiles.append(fullPath)
        else:
            print('This is not a file or directory:\n' + fullPath)
            raise Exception('File is not a file or directory, Why?:' + fullPath)
    
    
    
    return readFiles