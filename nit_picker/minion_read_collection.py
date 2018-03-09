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


from numpy import mean, random
from Bio.SeqIO import parse, write

from os.path import join

from .read_quality_analyzer import calculateReadQualities
from nanopore_prospector.common import createOutputFile, createScatterPlot

class minionReadCollection:
    # A class defining a collection of minION reads

    def __init__(self, readArray):
        self.readCollection=readArray
        self.readInputFormat=None
        self.readLengths = []
        self.readAvgPhredQualities = []
        self.gene = None
        
        if self.readCollection is None:
            self.readCollection = []
        
        # TODO: I commented this out, i think this will mean sometimes stats are not calculated.
        # TODO: Validate that I am still making read stats for fastq files
        #self.summarizeSimpleReadStats()

    def summarizeSimpleReadStats(self):
        #print('calculating read stats')        
        self.readLengths = []
        self.readAvgPhredQualities = []
        
        for currentRead in self.readCollection:  
            currentSeqLength = len(currentRead)
            self.readLengths.append(currentSeqLength)
            
            if (self.readInputFormat == 'fastq'):
                phredQualities = currentRead.letter_annotations["phred_quality"] 
                currentAvgPhredQuality = mean(phredQualities)            
                self.readAvgPhredQualities.append(currentAvgPhredQuality)   
        
    def concatenate(self, otherCollection): 
        self.readCollection = self.readCollection + otherCollection.readCollection
        self.readLengths = self.readLengths + otherCollection.readLengths
        self.readAvgPhredQualities = self.readAvgPhredQualities + otherCollection.readAvgPhredQualities
        #self.summarizeSimpleReadStats()
        
        
    def outputReadPlotsAndSimpleStats(self, outputDirectory, plotName, sampleID, referenceSequencelocation):
        self.summarizeSimpleReadStats()
        
        # Remove spaces for file names.
        simplePlotName = plotName.replace(' ','_')
        
        # Reads
        readOutputFileName = join(outputDirectory, str(sampleID) + '_' + simplePlotName + '.' + self.readInputFormat)
        allReadOutputFile = createOutputFile(readOutputFileName)
        print ('writing the all read file. the format is ' + self.readInputFormat)
        write(self.readCollection, allReadOutputFile, self.readInputFormat)
        allReadOutputFile.close()
        
        # Simple Statistics  
        simpleStatsOutputFileName = join(outputDirectory, str(sampleID) + '_' + simplePlotName + '.csv')
        self.printSimpleReadStats(simpleStatsOutputFileName)

        # Make up random qualities so we can still make a scatterplot.
        # This is for fasta files, where i want to see read lengths.
        # TODO: I mean, if we don't have read qualities, i could make a bar plot or something instead?
        # This solution is dumb.
        if (self.readInputFormat != 'fastq'):
            randomQualities = random.rand(len(self.readLengths))
            self.readAvgPhredQualities = randomQualities

        # Scatterplot      
        createScatterPlot(plotName
            , self.readLengths
            , self.readAvgPhredQualities 
            , "Read Lengths"
            , "Avg Read Quality(Phred)"
            , join(outputDirectory,  str(sampleID) + '_' + simplePlotName))
        
        # TODO: I really only need to calculate these statistics against the pass reads. 
        # Can i filter that somehow?
        
        # Calculate Statistics against Reference
        if(referenceSequencelocation is not None):
            # TODO: Calculate number of threads somewhere. I should pass it into nit_picker.
            numberThreads = 4
            calculateReadQualities(referenceSequencelocation, readOutputFileName, outputDirectory, simplePlotName, numberThreads)
 
        # Return some simple stats. These are useful downstream.
        # readstats is a 2d array, with lengths and qualities.
        return self.readLengths,self.readAvgPhredQualities 
   
    def printSimpleReadStats(self, outputFileLocation):
        print ('Printing simple statistics to this location:\n' + str(outputFileLocation))
        simpleStatOutputFile = createOutputFile(outputFileLocation)
        
        simpleStatOutputFile.write('Read_ID,Length,Avg_Phred_Quality\n')
        for index, currentRead in enumerate(self.readCollection):    
            
            if(self.readInputFormat == 'fastq'):
                readQual = self.readAvgPhredQualities[index]
            else:
                readQual = 0  
                
            readLength = self.readLengths[index]    
                    
            simpleStatOutputFile.write(str(currentRead.id) + ',' 
                + str(readLength) + ','
                + str(readQual) + '\n')
            
        simpleStatOutputFile.close()
    
        
def createCollectionFromReadFile(readFile):
    global readInputFormat
    #print ('loading reads from:' + readFile)
    
    # Determine Input File Type
    if (".fasta" == readFile[-6:] or ".fa" == readFile[-3:]):
        readInputFormat = "fasta"
        #raise Exception('Fasta files are not supported.  You can try to add this functionality if you want.')
    elif (".fastq"== readFile[-6:] or ".fq" == readFile[-3:]):
        readInputFormat = "fastq"
    else:
        print('I expect a .fasta or .fastq format for read input. Alternatively, specify a directory containing read inputs. Please check your input.')
        raise Exception('Bad Read Input Format')
    
    newCollectionObject = minionReadCollection(list(parse(readFile, readInputFormat)))
    newCollectionObject.readInputFormat = readInputFormat
    
    return newCollectionObject    
