# This file is part of allele_wrangler.
#
# allele_wrangler is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# allele_wrangler is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with allele_wrangler. If not, see <http://www.gnu.org/licenses/>.


from sys import exc_info

from os.path import split, join, isdir, isfile
from os import makedirs, system

from Bio.Seq import Seq
from Bio.SeqIO import write, parse
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.AlignIO import read
from Bio.Align import AlignInfo
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Sequencing.Applications import BwaIndexCommandline

from random import shuffle

from sklearn.cluster import KMeans
from numpy import asarray

#from subprocess import Popen, PIPE, STDOUT

from pysam import AlignmentFile

from .alignment_info import AlignmentInfo, AlignmentColumn

from nit_picker.read_quality_analyzer import calculateTotalCoverage


class AlleleWrangler():   
    
    def __init__(self, readsFile, outputDir, referenceFile, numIters, numThreads, splitHeterozygotes):
         
        print ('Setting up the Allele Wrangler...')
        self.readInput          = readsFile
        self.inputDirectory = False # Is the input a directory? If not, it is a single file.
        self.outputRootDirectory = outputDir
        self.referenceSequenceFileName     = referenceFile
        self.totalIterations      = numIters
        self.numberThreads         = numThreads
        self.splitHeterozygotes = splitHeterozygotes
        
        # Determine if input is a file, or a directory
        if (isfile(self.readInput)):
            #print ('Read input is a file that exists.')
            pass
            
            # Determine Input File Type
            if (".fasta" == self.readInput[-6:] or ".fa" == self.readInput[-3:]):
                self.readInputFormat = "fasta"
            elif (".fastq"== self.readInput[-6:] or ".fq" == self.readInput[-3:]):
                self.readInputFormat = "fastq"
            else:
                print('I expect a .fasta or .fastq format for read input. Alternatively, specify a directory containing read inputs. Please check your input.')
                raise Exception('Bad Read Input Format')
                
        elif (isdir(self.readInput)):
            #print ('Read input is a directory that exists.')
            self.inputDirectory = True
            
        else :
            print('I expect a .fasta or .fastq format for read input. Alternatively, specify a directory containing read inputs. Please check your input.')
            raise Exception('Bad Read Input Format')
 
    # Start here.
    def analyzeReads(self):
        # This method returns a dictionary of coverageResults.
        # Key is the location of the fasta consensus sequence.
        # Value is the coverage. (# reads in this alignment)        
        print ('Beginning analysis of reads.')        

        coverageResults = {}
        
        # Create a LogFile
        self.wranglerLog = createOutputFile(join(self.outputRootDirectory,'analysis_summary.txt'))
        
        self.wranglerLog.write('Read Input:' + self.readInput + '\n')
        self.wranglerLog.write('Result directory:' + str(self.outputRootDirectory) + '\n')
        self.wranglerLog.write('Number Threads:' + str(self.numberThreads) + '\n')

        readCount = len(list(parse(self.readInput, self.readInputFormat)))
        
        if (readCount > 5):
    
            # Was a reference sequence provided?
            if(self.referenceSequenceFileName is None):
                print('No Reference was provided.  No problem, I can generate my own reference.')
                self.createFirstGuessReferenceFromReads()            
    
            else:
                print('I was provided this reference filename:' + self.referenceSequenceFileName)
                pass
                # TODO I just commented this out, i already know i have a reference sequence.
                #consensusSequenceFileName = self.wrangleRecurser(None, 1)            
                #self.summarizeWrangling(consensusSequenceFileName)
    
            self.currentIteration = 1
            consensusSequenceFileName = self.analysisRecurser(self.referenceSequenceFileName)

            coverageResults = self.summarizeAnalysis(consensusSequenceFileName, self.splitHeterozygotes)
        else:
            print ('Skipping this read file. Not enough reads to assemble:' + str(self.readInput))
            self.wranglerLog.write('Skipping this read file. Not enough reads to assemble:' + str(self.readInput))
        
        self.wranglerLog.close()

        return coverageResults

    def summarizeAnalysis(self, finalConsensusSequence, splitHeterozygotes):  
        # TODO: I shouldnt need to pass in the splitHeterozygotes, since it is a class variable.
        print('\n\nSummarizing Analysis Results.') 

        try:
            coverageResults = {}
            if(splitHeterozygotes):
                # TODO: Im passing in a lot of self variables here. I should be able to remove those.
                
                self.heterozygousDirectory = join(self.outputRootDirectory, 'HeterozygousAlignment')                        
                heterozygousAlignedReadCount = self.alignReads(finalConsensusSequence,self.readInput,self.heterozygousDirectory, False)
                
                # Find heterozygous positions
                coverageResults = self.phaseHeterozygousReads()
               
            else:
                # write consensus to file. 
                # copy alignment reference and call it "final" 
                self.finalAlignmentDirectory = join(self.outputRootDirectory, 'FinalAlignment')                        
                alignedReadCount = self.alignReads(finalConsensusSequence,self.readInput,self.finalAlignmentDirectory, False)
                
                finalConsensusFilename = join(self.outputRootDirectory, 'AssembledConsensus.fasta')
                finalConsensusSequence = list(parse(finalConsensusSequence, 'fasta'))[0]
                #finalConsensusSequence.id = 'Assembled_Consensus_Sequence'
                sequenceWriter = createOutputFile(finalConsensusFilename)
                write([finalConsensusSequence], sequenceWriter, 'fasta')
                sequenceWriter.close()
                
                coverageResults[finalConsensusFilename] = alignedReadCount

                pass
            
            # TODO Return a value here. It should be a dictionary.
            # Key is the location of the fasta consensus sequence.
            # Value is the coverage. (# reads in this alignment)
            return coverageResults
                        
        except Exception:
            print ('Exception performing the summarize alignment')                  
            raise 
        

    
    def phaseHeterozygousReads(self):
    # TODO: Should this method accept a cluster count?
    # That will break some things. What things?
    # This method is only called from this file, in the summarizeAnalysis method.
    
        print('Splitting reads by heterozygous positions')
        
        # Get a list of reads for later.
        parsedReads = list(parse(self.readInput, self.readInputFormat))
        
        heterozygousConsensusDirectory = join(self.outputRootDirectory,'HeterozygousAlignment')

        # Open the bam file
        print ('opening final alignment_bamfile')
        bamfile = AlignmentFile(join(heterozygousConsensusDirectory,'alignment.bam'), 'rb')  
        
        # Load up the Alignment Reference file, we'll need it.
        alignmentReferenceFileName = join(heterozygousConsensusDirectory,'AlignmentReference.fasta')
        alignmentRef = list(parse(alignmentReferenceFileName, 'fasta'))[0]
     
        # get list of AlignedReads
        print ('Making a list of Aligned Reads.')
        readIDs = []
        for read in parsedReads:
            if not read.id in readIDs:
                readIDs.append(read.id)
        readIDs.sort()

        # Heterozygous base list
        heterozygousBasesSummaryFile = createOutputFile(join(heterozygousConsensusDirectory,'HeterozygousBases.txt'))         
        heterozygousBasesSummaryFile.write('List of Heterozygous Bases:')

        # get list of Heterozygous Positions
        # TODO: I suppose I don't need to align 100% of reads to determine heterozygosity.
        # Maybe this would speed up if i use a smaller alignment, or stop the loop after X reads
        print('Getting a list of Heterozygous Positions:')
        heterozygousPositions = []
        pileupIterator = bamfile.pileup(alignmentRef.id)
        for pileupColumn in pileupIterator:
            readCount = 0
            matchCount = 0
            mismatchCount = 0
            
            referenceBase = alignmentRef.seq[pileupColumn.pos].upper()
            
            # Iterate the Reads at this position           
            for pileupRead in pileupColumn.pileups:
                readCount += 1
                # indels
                if(pileupRead.is_del == 1 or pileupRead.indel > 0):
                    mismatchCount += 1                
                else:    
                    currentBase = pileupRead.alignment.query_sequence[pileupRead.query_position].upper()                    

                    if(currentBase == referenceBase):
                        matchCount += 1
                    else:
                        mismatchCount += 1
                   
            matchProportion = (1.0 * matchCount / readCount)
            #print ('Position ' + str(pileupColumn.pos) + ', Coverage ' + str(pileupColumn.n) + ', Match/Mismatch : ' + str(matchCount) + '/' + str(mismatchCount))
            #print ('Match Percentage ' + str(matchProportion))
            
            # TODO: Should accepted match proprtion be a commandline parameter?
            # if > 75% of bases match, this is not a heterzygous position
            if(matchProportion > .60):
                pass
                #print ('Position ' + str(pileupColumn.pos) + ', Coverage ' + str(pileupColumn.n) + ', Match/Mismatch : ' + str(matchCount) + '/' + str(mismatchCount))
                #print ('This position does not look heterozygous.')
            # If coverage is very low, we should not use this position
            elif ((1.0 * pileupColumn.n / readCount) < .25):
                pass
            else:
                #print ('HETEROZYGOUS Position ' + str(pileupColumn.pos) + ', Coverage ' + str(pileupColumn.n) + ', Match/Mismatch : ' + str(matchCount) + '/' + str(mismatchCount))
                heterozygousBasesSummaryFile.write (str(pileupColumn.pos) + ', Coverage ' + str(pileupColumn.n) + ', Match/Mismatch : ' + str(matchCount) + '/' + str(mismatchCount) + '\n')
                heterozygousPositions.append(pileupColumn.pos)
                
        heterozygousBasesSummaryFile.close()
            #print ('Pileup Column # ' + str(pileupIterator))

        print('Calculating read distance arrays:')            
        # I'm making this distance array. In this array, a 0 represents a Match.  a 1 represents indels or substitutions.
        # This way I can calculate "distance" in an arbitrary number of dimensions
        # Distance is a euclidian way to represent how far away a read is from the consensus,
        # based on the heterozygous positions.  Each heterozygous position is a "dimension" in this space
        distanceArrays = {}
        for readID in readIDs:
            distanceArrays[readID] = list([0] * len(heterozygousPositions))

        pileupIterator = bamfile.pileup(alignmentRef.id)
        for pileupColumn in pileupIterator:
            currentColumn = pileupColumn.pos
            
            # Only do this if the column number exists in our list of heterozygous positions
            if currentColumn in heterozygousPositions:
                
                heterozygousPositionIndex = heterozygousPositions.index(currentColumn)
                
                referenceBase = alignmentRef.seq[currentColumn].upper()
                for pileupRead in pileupColumn.pileups:
                    readID = pileupRead.alignment.query_name
                    
                    #print('Pos:' + str(currentColumn) + ', Refbase:' + str(referenceBase) + ', Read:' + str(readID))
                    
                    # In this model, the distance is either 0 or 1. This was intentional but
                    # Maybe we can tune the algorithm using these distances.
                    
                    if(pileupRead.is_del == 1):
                        distanceArrays[readID][heterozygousPositionIndex] = 1
                    elif(pileupRead.indel > 0):
                        distanceArrays[readID][heterozygousPositionIndex] = 1    
                    else:   
                        currentBase = pileupRead.alignment.query_sequence[pileupRead.query_position].upper()  
                        if(currentBase == referenceBase):
                            #print('Assinging Match. Column=' + str(currentColumn) + ', CurrentBase:' + str(currentBase) + ', HeterozygousPosIndex=' + str(heterozygousPositionIndex))
                            distanceArrays[readID][heterozygousPositionIndex] = 0
                        else:
                            distanceArrays[readID][heterozygousPositionIndex] = 1

        readIDs1, readIDs2 = self.clusterReads(distanceArrays)
        
        # Write group reads to output files  
        
        strand1OutputFileLocation = join(join(self.outputRootDirectory, 'Strand1ClusteredReads'), 'Strand1Reads.' + self.readInputFormat)
        strand2OutputFileLocation = join(join(self.outputRootDirectory, 'Strand2ClusteredReads'), 'Strand2Reads.' + self.readInputFormat)
        
        strand1OutputFile = createOutputFile(strand1OutputFileLocation)
        strand2OutputFile = createOutputFile(strand2OutputFileLocation)
        
        # Loop parsed reads, sort by read cluster
        for read in parsedReads:
            
            readFound = False
            for readID in readIDs1:
                if(readID in read.id):
                    write([read], strand1OutputFile, self.readInputFormat)
                    readFound = True
                    break
                
            # Don't enter this loop if we already found the read. Save a bit of time.
            if not readFound:
                for readID in readIDs2:
                    if(readID in read.id):
                        write([read], strand2OutputFile, self.readInputFormat)
                        break
        
        strand1OutputFile.close()
        strand2OutputFile.close()
 
        # Dictionary of results to return. Key is location of the consensus sequence.
        # Value is the # of reads represented in this consensus alignment.
        coverageResults = {}
 
        # Assemble those 2 output files
        strand1Wrangler = AlleleWrangler(
            strand1OutputFileLocation
            , join(self.outputRootDirectory, 'Strand1Alignment')
            , join(self.heterozygousDirectory, 'AlignmentReference.fasta')
            , 6
            , self.numberThreads
            , False)        
        strand1CoverageResults = strand1Wrangler.analyzeReads()
        
        strand2Wrangler = AlleleWrangler(
            strand2OutputFileLocation
            , join(self.outputRootDirectory, 'Strand2Alignment')
            , join(self.heterozygousDirectory, 'AlignmentReference.fasta')
            , 6
            , self.numberThreads
            , False)
        strand2CoverageResults = strand2Wrangler.analyzeReads()
        
        # Merge the dictionaries of coverage values ane return them.
        for key in strand1CoverageResults.keys():
            coverageResults[key] = strand1CoverageResults[key]
        for key in strand2CoverageResults.keys():
            coverageResults[key] = strand2CoverageResults[key]
        
        print ('Done Phasing Reads.')
        return coverageResults

        
    def clusterReads(self, dataInput):
        readIDs = list(dataInput.keys())
        readDistanceArrays=asarray([dataInput[k] for k in readIDs if k in dataInput])
        
        print('I will attempt to cluster and separate heterozygous reads.')
        print('There are ' + str(len(readDistanceArrays)) + ' distance arrays.')
        print('They are of length ' + str(len(readDistanceArrays[0])))
        
        # TODO: I can choose number of jobs/threads.
        kmeansObject = KMeans(n_clusters=2 )
        predictedLabels = kmeansObject.fit_predict(readDistanceArrays)
        
        print ('Done Clustering. Here are the predicted labels:\n' + str(predictedLabels))
        
        # Get read IDs from each cluster and populate arrays
        strand1ReadIDs = []
        strand2ReadIDs = []
        
        for labelIndex, label in enumerate(predictedLabels):
            #print ('Read:' + str(readIDs[labelIndex]) + ' , Label: ' + str(label) + ' Data:\n' + str(readDistanceArrays[labelIndex]))
            #Arbitrarily, label 0 = strand 1 and label 1 = strand 2
            if(label == 0):
                strand1ReadIDs.append(readIDs[labelIndex])
            elif(label == 1):
                strand2ReadIDs.append(readIDs[labelIndex])
            else:
                raise Exception('Unknown Cluster Label:' + str(label) + '. Perhaps you ended up with extra clusters.')

        return strand1ReadIDs, strand2ReadIDs
            
    # Input = location of Reference Sequence
    # Output = Location of the Consensus Sequence    
    def analysisRecurser(self, currentReferenceSequence):
        print('\n\nAttempting a read alignment and consensus polish, Iteration ' + str(self.currentIteration))

        try:
            
            # Was a consensus sequence provided?
            if(currentReferenceSequence is None):
                print('No Reference was provided.  No problem, I will use the first read as a reference.')
                currentReferenceSequence = self.createFirstGuessReferenceFromReads()            
            else:
                print('A Reference sequence was provided.')
                #self.prepareReferenceInput() 
                #Copy the Reference to the input folder
                #consensusSequence = self.openConsensusSequence()            
                       
            alignmentSubdir = join(self.outputRootDirectory,'alignments')
            if not isdir(alignmentSubdir):
                makedirs(alignmentSubdir)
            currentIterationSubdirectory = join(alignmentSubdir,'iter_'+ str(self.currentIteration))
            if not isdir(currentIterationSubdirectory):
                makedirs(currentIterationSubdirectory)
            
            # Write a bit about this iteration to the Log.

            self.wranglerLog.write('\nIteration # (' + str(self.currentIteration) + '/' + str(self.totalIterations) + ')\n')
            self.wranglerLog.write('Iteration Alignment Directory:' + currentIterationSubdirectory + '\n')
            self.wranglerLog.write('Reference Filename:' + str(currentReferenceSequence) + '\n')


            self.alignReads(currentReferenceSequence,self.readInput,currentIterationSubdirectory, True)
            self.currentAlignmentInfo = self.analyzeAlignment(currentIterationSubdirectory)
            
            # If we want more iterations, I should Recurse and try again.
            if (int(self.totalIterations) > int(self.currentIteration)):
                # On the next iteration, we want to use the new consensus sequence 
                # as the reference. 
                print('I am on iteration (' + str(self.currentIteration) + '/' + str(self.totalIterations) + ') I will continue...')
                newReferenceSequenceFileName = join(currentIterationSubdirectory, 'Consensus.fasta')
                self.currentIteration += 1
                
                # Return the consensus from one layer deeper.

                return self.analysisRecurser(newReferenceSequenceFileName)
                
                
            else:
                print('That was the last iteration (' + str(self.currentIteration) + '/' + str(self.totalIterations) + '), I am done now.')
                
                # Return the consensus
                return join(currentIterationSubdirectory, 'Consensus.fasta')


        except Exception:
            print ('Exception encountered in analyzeReads()')                  
            print (exc_info()[0])
            print (exc_info()[1])
            print (exc_info()[2])        
            raise 
        
    
    def alignReads(self, referenceLocation, readFileLocation, alignmentOutputDirectory, useReadSubset):
        # Perform BW Alignment.  Align all reads against the Reference.
        # This method returns the # of reads that aligned to this reference.
        print('\nStep 1.) Aligning reads against the reference.')
        
        if not isdir(alignmentOutputDirectory):
            makedirs(alignmentOutputDirectory)
        
        # Part 1 Index the Reference        
        try:
            # Copy the reference sequence to the alignment directory. This is a complicated way to do it.
            newReferenceLocation = join(alignmentOutputDirectory,'AlignmentReference.fasta')
            refSequence = list(parse(referenceLocation, 'fasta'))[0]
            refSequence.id = 'AlignmentReference'
            sequenceWriter = createOutputFile(newReferenceLocation)
            write([refSequence], sequenceWriter, 'fasta')
            sequenceWriter.close()
                        
            # Index The Reference
            referenceIndexName = newReferenceLocation.replace('.fasta','.mmi')        
            cmd = ('minimap2 -d ' + referenceIndexName + ' ' + newReferenceLocation)
            system(cmd)
            #asdf
            #indexCmd = BwaIndexCommandline(infile=newReferenceLocation, algorithm="bwtsw")
            #indexCmd()

            
        except Exception:
            print ('Exception indexing alignment reference. Is bwa installed? folder writing permission issue?')                  
            raise
        
        # TODO: Make this a commandline parameter.  Lower = faster. Higher = more accurate consensus correction
        alignmentReadSubsetCount = 150
        try:
            if useReadSubset:
                # load Reads
                parsedReads = list(parse(readFileLocation, self.readInputFormat))
                
                # If there aren't enough reads for this
                if (len(parsedReads) < alignmentReadSubsetCount):
                    alignmentReadSubsetCount = len(parsedReads)
                
                # choose random subset
                randomIndexes = list(range(0, len(parsedReads)))
                shuffle(randomIndexes)                
                sampledReads = []
                for i in range(0,alignmentReadSubsetCount):
                    sampledReads.append(parsedReads[randomIndexes[i]])
                
                # write random reads to alignment directory
                # Reassign the reads we'll use downstream
                readFileLocation = join(alignmentOutputDirectory, 'ReadSample.fasta')        
                readSampleWriter = createOutputFile(readFileLocation)          
                   
                write(sampledReads, readSampleWriter, 'fasta')
                readSampleWriter.close()

            else:
                # We'll use the whole read file.
                pass
                        
        except Exception:
            print ('Exception selecting a read subset.')                  
            raise
        
        # Part 2 Align
        try:
            # TODO: How can i put this into biopython?  Pipelines are hard.
            # align | sam->bam | sort
            
            #tempAlignmentName = join(alignmentOutputDirectory,'alignment')
            alignmentOutputName = join(alignmentOutputDirectory,'alignment.bam')
            #bwaMemArgs = "-t " + str(self.numberThreads) + " -x ont2d"
            #cmd = ("bwa mem " + 
            #    bwaMemArgs + " " +  
            #    newReferenceLocation + " " +
            #    readFileLocation + 
            #    " | samtools view  -Sb - | samtools sort -o "
            #    + alignmentOutputName)
            #print ('alignment command:\n' + cmd)
            #system(cmd)
            
            # TODO: Honeslty the -ax map-ont settings sometimes allow for some messed up alignments.
            # It allows "secondary" alignments int he bam, which are not quite accurate
            # Full of SNPs and makes the alignment a bit bogus.
            # I am doing some tests, it might make more sense to use the asm5 or asm10 settings instead of the ONT settings.            
            cmd = ("minimap2 -ax map-ont " + 
                referenceIndexName + " " +  
                readFileLocation + 
                " | samtools view -b | samtools sort -o "
                + alignmentOutputName)
        #print ('alignment command:\n' + cmd)
            system(cmd)
            
            
            #alignmentOutputName = tempAlignmentName + '.bam'
            
        except Exception:
            print ('Exception aligning reads against reference. Are bwa and samtools installed?')                  
            raise 
        
        # Part 3 Index Alignment
        try:
            cmd = ("samtools index " + alignmentOutputName)
            #print ('alignment index command:\n' + cmd)
            system(cmd)
            #print ('index command:\n' + cmd)
        except Exception:
            print ('Exception indexing alignment reference. Is bwa installed?')                  
            raise 
        
        alignedReadCount = calculateTotalCoverage(alignmentOutputName)
        return alignedReadCount





    def analyzeAlignment(self, alignmentOutputDirectory):
        print ('\nStep 2.) Parse the alignment and create a new consensus sequence.')
        
        # Load up the Alignment Reference file, we'll need it.
        alignmentReferenceFileName = join(alignmentOutputDirectory,'AlignmentReference.fasta')
        alignmentRef = list(parse(alignmentReferenceFileName, 'fasta'))[0]
        
        # Count the reads in the input file
        totalReadCount = len(list(parse(self.readInput, self.readInputFormat)))
        #self.readInputFormat
        #self.readInput
                
        # We generate a new consensus sequence from the alignment results.
        newConsensusSequence = ""
        
        # Open the bam file
        bamfile = AlignmentFile(join(alignmentOutputDirectory,'alignment.bam'), 'rb')  
        
        # Open alignment analysis text file
        alignmentSummaryFile = createOutputFile(join(alignmentOutputDirectory,'AlignmentSummary.csv')) 
        alignmentSummaryFile.write('Ref_Position,Ref_Base,Reference_Adjustment,Aligned_Count,Unaligned_Count,Match_Count,Mismatch_Count,In_Count,Del_Count,A_Count,G_Count,C_Count,T_Count\n')
        
        # A smaller log. I will provide human-readable descriptions of the
        # bases that were adjusted in the new consensus sequence.
        # TODO: Provide surrounding sequence as well, maybe it's a repeat region....
        # Acutally NAH, I want to just put it in the wrangler log. 
        #adjustedBasesSummaryFile = createOutputFile(join(alignmentOutputDirectory,'AdjustedBases.txt')) 
        
        # Todo: I should keep a more structured array of info for these alignments.
        # Store this info into an object
        #class columnStats():
        alignmentInfo = AlignmentInfo()
        
        # Keep a running total of adjustments made to the reference.
        # If this total is 0, then theoretically the consensus matches the alignment reference, and we're done.
        totalSequenceAdjustments = 0
        
        # Iterate the reference sequence column by column.
        pileupIterator = bamfile.pileup(alignmentRef.id)
        
        for pileupColumn in pileupIterator:
            
            currentAlignmentColumn = AlignmentColumn()
            #columnResults = None
           # columnResults.name='ll'
            #
            """referencePosition = 0
            referenceBase = ''
            referenceAdjustment = '?'
            alignedCount = 0
            unalignedCount = 0
            matchCount = 0
            mismatchCount = 0
            inCount = 0
            delCount = 0
            aCount = 0
            gCount = 0
            cCount = 0
            tCount = 0"""
            
            currentAlignmentColumn.referencePosition = pileupColumn.reference_pos
            currentAlignmentColumn.referenceBase = alignmentRef[pileupColumn.reference_pos].upper()
            currentAlignmentColumn.alignedCount = pileupColumn.nsegments
            currentAlignmentColumn.unalignedCount = totalReadCount - currentAlignmentColumn.alignedCount
            
            # Iterate the Reads at this position           
            for pileupRead in pileupColumn.pileups:
                
                # If this read is a deletion
                if(pileupRead.is_del == 1):
                    currentAlignmentColumn.delCount += 1
                # else if this read is an insertion
                elif(pileupRead.indel > 0):
                    
                    #print ('INSERTION DETECTED, INDEL=' + str(pileupRead.indel))  
                    currentAlignmentColumn.inCount += 1                   
                # Else if it is a refskip (TODO What does this mean? no read aligned? Count these?)
                elif(pileupRead.is_refskip):
                    print('This read is a refskip, i dont know what that means:' + pileupRead.alignment.query_name)
                    raise Exception('This read is a refskip, i dont know what that means:' + pileupRead.alignment.query_name)
                # else this means we have a base aligned at this position for this read.
                else:    
                    currentBase = pileupRead.alignment.query_sequence[pileupRead.query_position].upper()                    
                    #print('Reference,Current:' + referenceBase + ',' + currentBase)
                    #print('Curr')
                    if(currentBase == currentAlignmentColumn.referenceBase):
                        currentAlignmentColumn.matchCount += 1
                    else:
                        currentAlignmentColumn.mismatchCount += 1
                   
                # Count the nucleotide 
                if (currentBase == 'A'):
                    currentAlignmentColumn.aCount += 1
                elif (currentBase == 'G'):
                    currentAlignmentColumn.gCount += 1
                elif (currentBase == 'C'):
                    currentAlignmentColumn.cCount += 1
                elif (currentBase == 'T'):
                    currentAlignmentColumn.tCount += 1
                else:
                    print('Unknown Base found in Alignment at position ' + str(currentAlignmentColumn.referencePosition) + ':' + currentBase)
                    raise Exception('Unknown Base in Alignment')
                
                
                # TODO: What if the query insertion sequence is longer than one base?
                # Maybe I can only adjust one base per iteration, is that okay? Probably for the Best, actually..
                # Don't worry bout it for now.
            
            # Calculate highest frequency base
            # I hope this algorithm makes sense, probably there is a smarter way to do it.
            if(currentAlignmentColumn.aCount >= currentAlignmentColumn.gCount and currentAlignmentColumn.aCount >= currentAlignmentColumn.cCount and currentAlignmentColumn.aCount >= currentAlignmentColumn.tCount):
                mostFrequentBase = 'A'
                mostFrequentBaseCount = currentAlignmentColumn.aCount
            elif(currentAlignmentColumn.gCount >= currentAlignmentColumn.cCount and currentAlignmentColumn.gCount >= currentAlignmentColumn.tCount):
                mostFrequentBase = 'G'
                mostFrequentBaseCount = currentAlignmentColumn.gCount
            elif(currentAlignmentColumn.cCount >= currentAlignmentColumn.tCount):
                mostFrequentBase = 'C'
                mostFrequentBaseCount = currentAlignmentColumn.cCount
            else:
                mostFrequentBase = 'T'
                mostFrequentBaseCount = currentAlignmentColumn.tCount


            
            # Add the next base to the new consensus sequence            
            if (currentAlignmentColumn.matchCount >= currentAlignmentColumn.mismatchCount and currentAlignmentColumn.matchCount >= currentAlignmentColumn.inCount and currentAlignmentColumn.matchCount >= currentAlignmentColumn.delCount):
                # Aligned bases match the reference, add reference base to the consensus.
                referenceAdjustment='-'
                newConsensusSequence += currentAlignmentColumn.referenceBase
                
            elif (currentAlignmentColumn.inCount >= currentAlignmentColumn.mismatchCount and currentAlignmentColumn.inCount >= currentAlignmentColumn.delCount):
                # Aligned bases show an insertion.
                # Add the Reference Base and the Insertion Base to the consensus.  
                totalSequenceAdjustments += 1 
                referenceAdjustment='I'  
                newConsensusSequence += currentAlignmentColumn.referenceBase + mostFrequentBase         
                
                self.wranglerLog.write(str(currentAlignmentColumn.referencePosition) + ':Insertion' +
                    '\n(' + str(currentAlignmentColumn.inCount) + '/' + str(currentAlignmentColumn.alignedCount) + ') = ' + str((100.0 * currentAlignmentColumn.inCount) / currentAlignmentColumn.alignedCount) + '% of aligned reads'
                    '\n(' + currentAlignmentColumn.referenceBase + ' > ' + currentAlignmentColumn.referenceBase + mostFrequentBase + ')' +
                    '\n')
                
                #TODO: I need to insert multiple bases, if that is waht the alignment suggests.

            elif (currentAlignmentColumn.delCount >= currentAlignmentColumn.mismatchCount):
                # Reads show a deletion.
                # Don't add anything to the consensus.
                totalSequenceAdjustments += 1
                referenceAdjustment='D'
                
                self.wranglerLog.write(str(currentAlignmentColumn.referencePosition) + ':Deletion' +
                    '\n(' + str(currentAlignmentColumn.delCount) + '/' + str(currentAlignmentColumn.alignedCount) + ') = ' + str((100.0 * currentAlignmentColumn.delCount) / currentAlignmentColumn.alignedCount) + '% of aligned reads'
                    '\n(' + currentAlignmentColumn.referenceBase + ' > _)' +
                    '\n')
                
            else:
                # Mismatch base.
                # Add the highest read count base to the reference.
                # It might actually be the same base as the reference,
                # Because this just means there are more mismatches than matches.
                # Problematic base, at least we'll notice here.
                # TODO: What to do with highly heterozygous Positions?
                # I should report those that look particularly heterozygous, somewhere.
                newConsensusSequence += mostFrequentBase 
                totalSequenceAdjustments += 1     
                referenceAdjustment='M'   
                
                self.wranglerLog.write(str(currentAlignmentColumn.referencePosition) + ':Mismatch' +
                    '\n(' + str(mostFrequentBaseCount) + '/' + str(currentAlignmentColumn.alignedCount) + ') = ' + str((100.0 * mostFrequentBaseCount) / currentAlignmentColumn.alignedCount) + '% of aligned reads'
                    '\n(' + currentAlignmentColumn.referenceBase + ' > ' + mostFrequentBase + ')' +
                    '\n')
              

            # Write a line to the alignment Summary 
            alignmentSummaryFile.write(str(currentAlignmentColumn.referencePosition) + 
                ',' + str(currentAlignmentColumn.referenceBase) +
                ',' + str(referenceAdjustment) + 
                ',' + str(currentAlignmentColumn.alignedCount) + 
                ',' + str(currentAlignmentColumn.unalignedCount) + 
                ',' + str(currentAlignmentColumn.matchCount) + 
                ',' + str(currentAlignmentColumn.mismatchCount) + 
                ',' + str(currentAlignmentColumn.inCount) + 
                ',' + str(currentAlignmentColumn.delCount) + 
                ',' + str(currentAlignmentColumn.aCount) + 
                ',' + str(currentAlignmentColumn.gCount) + 
                ',' + str(currentAlignmentColumn.cCount) + 
                ',' + str(currentAlignmentColumn.tCount) +
                '\n')
            
            alignmentInfo.alignmentColumns.append(currentAlignmentColumn)
            
        print('\nTotal Sequence Adjustments:' + str(totalSequenceAdjustments) + ' (How many bases the consensus differs from the reference.)\n')    
        
        # Write the newly constructed consensus sequence.
        currentConsensusSequenceFileName = join(alignmentOutputDirectory, 'Consensus.fasta')        
        consensusWriter = createOutputFile(currentConsensusSequenceFileName)          
           
        # TODO: How to i give this a better name? Can I find a gene guess or something?
        sequenceID = "Consensus_Sequence"

        write([SeqRecord(Seq(newConsensusSequence,
            IUPAC.unambiguous_dna),
            id=sequenceID, description="") ], consensusWriter, 'fasta')
        consensusWriter.close()
            
        self.wranglerLog.write('Total Sequence Adjustments:' + str(totalSequenceAdjustments) + '\n')
            
        # Close Summary Files
        alignmentSummaryFile.close()
        #adjustedBasesSummaryFile.close()
        
        return alignmentInfo
        
        
        #return totalSequenceAdjustments
   
    def createFirstGuessReferenceFromReads(self):   
        #TODO: I should make this a commandline parameter. More = MSA takes longer. Less = worse reference
        msaReadCount = 4
        
        print ('I choose ' + str(msaReadCount) + ' random reads.'
            + '\nThese are aligned to form a rough initial consensus sequence. Here:'
            + '\n' + join(self.outputRootDirectory,'Initial_Reference')
            + '\nPerforming ClustalO Multiple Sequence Alignment Now...')
        try:            
            # Load Reads from File

            parsedReads = list(parse(self.readInput, self.readInputFormat))            
            referenceSequence = None

            
            # Reference Directory
            referenceDirectory = join(self.outputRootDirectory,'Initial_Reference')
            if not isdir(referenceDirectory):
                makedirs(referenceDirectory)
                        
            if (len(parsedReads) > msaReadCount):
                

                # Select a subset of reads for Multiple SequneceAlignment. Randomly, i guess.
                randomIndexes = list(range(0, len(parsedReads)))
                shuffle(randomIndexes)                
                rawClustalReads = []
                for i in range(0,msaReadCount):
                    rawClustalReads.append(parsedReads[randomIndexes[i]])
              
                rawClustalReadsFilename = join(referenceDirectory, 'MSARaw.fasta')                
                rawClustalReadsFileWriter = createOutputFile(rawClustalReadsFilename)        
                write(rawClustalReads, rawClustalReadsFileWriter, 'fasta')
                rawClustalReadsFileWriter.close()
            
                #Perform Clustal MSA
                clustalOAlignmentOutputFileName = join(referenceDirectory, 'clustalOAlignment.fasta')
                clustalOCommandLine = ClustalOmegaCommandline(infile=rawClustalReadsFilename, outfile=clustalOAlignmentOutputFileName, verbose=True, auto=True, force=True, threads=int(self.numberThreads))
                clustalOCommandLine()                
        
                # Calculate consensus 
                # A dumb consensus has lots of ambiguous nucleotides.  We'll polish those out later.
                alignmentType = 'fasta'    
                alignmentObject = read(clustalOAlignmentOutputFileName, alignmentType)           
                alignmentSummaryInfo = AlignInfo.SummaryInfo(alignmentObject)                
                dumbConsensus = alignmentSummaryInfo.dumb_consensus(threshold=.5)
                
                referenceSequence = SeqRecord(Seq(str(dumbConsensus) , IUPAC.IUPACUnambiguousDNA),
                    id='Initial_Consensus',
                    description='Initial_Consensus')

                
            # Else
            else:
                # Select the first read, use it as the reference. It's something.
                #referenceSequence = parsedReads[0]
                # You know what? we should just give up. There aren't enough reads to assemble.
                #raise Exception('Not enough reads to continue.')
                referenceSequence = SeqRecord(Seq('' , IUPAC.IUPACUnambiguousDNA),
                    id='Initial_Consensus',
                    description='Initial_Consensus')
                        
             
            #Write reference to file
            self.referenceSequenceFileName = join(referenceDirectory, 'FirstGuessReference.fasta')            
            firstGuessRefFileWriter = createOutputFile(self.referenceSequenceFileName)        
            write([referenceSequence], firstGuessRefFileWriter, 'fasta')

            firstGuessRefFileWriter.close()
            
            return self.referenceSequenceFileName
       
       
            print ('Done making initial consensus sequence.')
      
                                    
                                     
        except Exception:
            print ('Exception encountered in createFirstGuessReferenceFromReads()') 
            print (exc_info()[0])
            print (exc_info()[1])
            print (exc_info()[2]) 
            raise    

# This method is a directory-safe way to open up a write file.
# TODO: Don't ihave this method in common somewhere?
# I should delete it here. It's not hurtin nothin but i should delete this method for consistency. Common?
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not isdir(tempDir):
        print('Making Directory:' + tempDir)
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput
