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
# from Bio.Sequencing.Applications import BwaIndexCommandline

from random import shuffle

from sklearn.cluster import KMeans
from numpy import asarray

#from subprocess import Popen, PIPE, STDOUT

from pysam import AlignmentFile

from .alignment_info import AlignmentInfo, AlignmentColumn

#from nit_picker.read_quality_analyzer import calculateTotalCoverage
from nanopore_prospector.common import alignReads


class AlleleWrangler():   
    
    def __init__(self, readsFile, outputDir, referenceFile, numIters, numThreads, splitHeterozygotes, snps ):
         
        print ('Setting up the Allele Wrangler...')
        self.readInput          = readsFile
        self.inputDirectory = False # Is the input a directory? If not, it is a single file.
        self.outputRootDirectory = outputDir
        self.referenceSequenceFileName     = referenceFile
        self.totalIterations      = numIters
        self.numberThreads         = numThreads
        self.splitHeterozygotes = splitHeterozygotes
        self.snps = snps
        
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

                # i should only do this if we were NOT provided a reference. I don't need to iterate on a known reference. Ok.

                self.currentIteration = 1
                # TODO I think this consensusSequenceFileName is a wasted variable, shouldn't I just re-store it to the self.referenceSequenceFileName ? I think that's better.
                consensusSequenceFileName = self.analysisRecurser(self.referenceSequenceFileName)
    
            else:
                print('I was provided this reference filename:' + self.referenceSequenceFileName)
                consensusSequenceFileName = self.referenceSequenceFileName
                #TODO: there is something funny here and i cant remember why i did this. I'm passing the reference sequence as the consensus sequence.
                # That's not correct and I think it is because I don't want to iterate on a reference sequence
                # Maybe I need to...to generate a consensus sequence from this informaiton.

                pass

            coverageResults = self.summarizeAnalysis(consensusSequenceFileName)
        else:
            print ('Skipping this read file. Not enough reads to assemble:' + str(self.readInput))
            self.wranglerLog.write('Skipping this read file. Not enough reads to assemble:' + str(self.readInput))
        
        self.wranglerLog.close()

        return coverageResults

    def summarizeAnalysis(self, finalConsensusSequence):
        print('\n\nSummarizing Analysis Results.')

        try:
            coverageResults = {}
            if(self.splitHeterozygotes):

                self.heterozygousDirectory = join(self.outputRootDirectory, 'HeterozygousAlignment')
                #alignReads(referenceLocation, readFileLocation, alignmentOutputDirectory, useReadSubset, numberThreads, excludeShortAlignments):
                heterozygousAlignedReadCount = alignReads(finalConsensusSequence,self.readInput,self.heterozygousDirectory, False)
                
                # Find heterozygous positions
                coverageResults = self.phaseHeterozygousReads()
               
            else:
                # write consensus to file. 
                # copy alignment reference and call it "final" 
                self.finalAlignmentDirectory = join(self.outputRootDirectory, 'FinalAlignment')                        
                alignedReadCount = alignReads(finalConsensusSequence,self.readInput,self.finalAlignmentDirectory, False)
                
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
        heterozygousBasesSummaryFile = createOutputFile(join(heterozygousConsensusDirectory, 'HeterozygousBases.txt'))
        heterozygousBasesSummaryFile.write('List of Heterozygous Bases (0-based):\n')


        if (self.snps is not None and len(self.snps) > 0):
            # A string of SNPs was passed in, I don't need to calculate them myself.
            # TODO: I could write alignment stats here, like I do when i self-calculate the hetero positions.
            # This is just a simple list of 0-based positions.
            for snp in self.snps:
                heterozygousBasesSummaryFile.write(str(snp) + '\n')

        else:

            # get list of Heterozygous Positions
            # TODO: I suppose I don't need to align 100% of reads to determine heterozygosity.
            # Maybe this would speed up if i use a smaller alignment, or stop the loop after X reads
            print('Getting a list of Heterozygous Positions:')
            self.snps = []
            pileupIterator = bamfile.pileup(alignmentRef.id)
            for pileupColumn in pileupIterator:
                readCount = 0
                matchCount = 0
                mismatchCount = 0
                insCount = 0
                delCount = 0

                # dictionary of base counts.



                referenceBase = alignmentRef.seq[pileupColumn.pos].upper()

                # Iterate the Reads at this position. Each read at each position is either:
                # ins, Del, match, mismatch.

                #TODO: is it possible to exclude secondary/supplemetnary in the pileups method?  No.
                for pileupRead in pileupColumn.pileups:

                    #TODO: Important. Filter secondary / supplementary reads. This is causing problems, these secondary reads are FULL of snps.
                    # Difficulty: these parameters are on an aligned segment.

                    alignedSegmentObject = pileupRead.alignment

                    if(False):
                        pass
                    elif(alignedSegmentObject.is_secondary):
                        #print ('Secondary read at Position ' + str(pileupColumn.pos))
                        pass
                    elif(alignedSegmentObject.is_supplementary):
                        #print ('Supplementary read at Position ' + str(pileupColumn.pos))
                        pass


                    # Just trying some things, not sure what these mean.
                    #elif (alignedSegmentObject.is_unmapped):
                    #    print('UNMAPPED READ!!!!!!!!!!!!!!!!!!! what does that mean?')
                    #elif (alignedSegmentObject.is_qcfail):
                    #    print('This read was a QC failure. What does that mean?????????????')

                    else:
                        readCount += 1
                        # indels
                        if(pileupRead.is_del == 1):
                            delCount += 1
                        elif(pileupRead.indel > 0):
                            insCount += 1
                        else:
                            currentBase = pileupRead.alignment.query_sequence[pileupRead.query_position].upper()

                            if(currentBase == referenceBase):
                                matchCount += 1
                            else:
                                mismatchCount += 1


                    # This is a cheap way to stop analysis early. I will only analyze the first 250 reads.
                    # Potential problem: are these reads sorted somehow? Maybe my numbers are biased by only looking at the
                    # first reads
                    # Todo: This is another parameter that can be tuned. Add to inputs? Maybe.
                    maxAnalyzedReadCounts = 1000

                    if(readCount > maxAnalyzedReadCounts):
                        break

                matchProportion =      (1.0 * matchCount / readCount)
                insertionProportion =  (1.0 * insCount / readCount)
                deletionProportion =   (1.0 * delCount / readCount)
                mismatchProportion =   (1.0 * mismatchCount / readCount)

                #print ('Position ' + str(pileupColumn.pos) + ', Coverage ' + str(pileupColumn.n) + ', Match/Mismatch : ' + str(matchCount) + '/' + str(mismatchCount))
                #print ('Match Percentage ' + str(matchProportion))

                # TODO: Should accepted match proprtion be a commandline parameter?
                # if > 75% of bases match, this is not a heterzygous position
                baseCutoff = .70

                if(matchProportion > baseCutoff or insertionProportion > baseCutoff or deletionProportion > baseCutoff):
                    pass
                    #print ('Position ' + str(pileupColumn.pos) + ', Coverage ' + str(pileupColumn.n) + ', Deletion/Insertion/Match/Mismatch : ' + str(delCount) + '/' + str(insCount) + '/' + str(matchCount) + '/' + str(mismatchCount))
                    #print ('This position does not look heterozygous.')


                # If coverage is very low, we should not use this position
                # This logic is flawed, i think this is never working.
                elif ((1.0 * pileupColumn.n / readCount) < .25):
                    pass

                elif (mismatchProportion > baseCutoff):
                    pass

                # These are the hardcoded values I used for the DRA analysis. Cheating.

                # # I want to write a condition where we don't use the position if it's not clearly polymorphic.
                # #elif (False):
                # #    pass
                # # If the mismatch proportion is too high, what happens? What if there are 2 different bases that are mismatched, like if both my alleles have a different snp from reference. I'll miss that right now.
                #

                # # TEMP, this is very temporary. This is specific to a reference.
                # # TODO : Fix these hard coded values.
                # TODO: I don't really need this code, this is to ignore regions of my DRA reference.
                # Instead, I can pass in a list of 1-based polymorphic positions to sort based on those. A "whitelist" instead of a "blacklist"
                # # In a perfect world....I could tell what positions are heterozygous, but I can't.



                # # I can tell if this sequence is a homopolymer though, but looking at the bases around it.....But that's not the correct thing to do.
                # # I can keep this logic but make it a parameter. Big deletion regions are hard to analyze so I'm just ignoring them for now.
                # elif(5890 <= pileupColumn.pos <= 5970):
                #     print('WARNING: I am skipping analysis on a region using hardcoded values, check this in allele_wrangler.')
                #     pass
                # elif (6203 <= pileupColumn.pos <= 6212):
                #     print('WARNING: I am skipping analysis on a region using hardcoded values, check this in allele_wrangler.')
                #     pass
                # # Big String of A's
                # elif (774 <= pileupColumn.pos <= 796):
                #     print('WARNING: I am skipping analysis on a region using hardcoded values, check this in allele_wrangler.')
                #     pass
                # #Known homopolymer positions....this is terrible programming.
                # # I could at least pass these in ad ignored positions....
                # elif (pileupColumn.pos in (403,430, 1479, 1510, 1683,
                #         1991, 1996, 1997, 2003, 2009, 2093, 2100, 2133, 2134, 2191,
                #         2262, 2289, 2294, 2342, 2449, 2450, 2524, 2647, 2663, 2732,
                #         2895, 2902, 3113, 3114, 3180, 3197, 3362, 3396, 3453, 3542,
                #         3551, 3665, 3832, 3903, 3953, 4108, 4109, 4400, 4639, 4698,
                #         4703, 4769, 4785, 4786, 4828, 4878, 5084, 5301, 5302, 5449,
                #         5575, 5597, 6155, 6279, 6280, 6314, 6375, 6376, 6712, 6755,
                #         6790, 7084, 7631, 7718, 7769, 7971, 7978, 8132, 8133, 8134,
                #         8314, 8315, 8352, 8476, 8477, 8478, 8642, 8650, 8651, 8652,
                #         8653, 8654, 8655, 8656, 8657, 8698, 8725, 8753, 8759
                #         )):
                #     print('WARNING: I am skipping analysis on a region using hardcoded values, check this in allele_wrangler.')
                #     pass


                else:
                    #heterozygousBasesSummaryFile.write (str(pileupColumn.pos) + ', Coverage ' + str(pileupColumn.n) + ', Deletion/Insertion/Match/Mismatch : ' + str(delCount) + '/' + str(insCount) + '/' + str(matchCount) + '/' + str(mismatchCount) + '\n')
                    heterozygousBasesSummaryFile.write(str(pileupColumn.pos) + ', Coverage ' + str(
                        pileupColumn.n) + ', Deletion/Insertion/Match/Mismatch : ' + str(delCount) + '/' + str(
                        insCount) + '/' + str(matchCount) + '/' + str(mismatchCount)
                        + ' : ' + str(round(deletionProportion,2)) + '/'
                        + str(round(insertionProportion, 2)) + '/'
                        + str(round(matchProportion, 2)) + '/'
                        + str(round(mismatchProportion, 2))
                        + '\n')
                    self.snps.append(pileupColumn.pos)





        heterozygousBasesSummaryFile.close()
            #print ('Pileup Column # ' + str(pileupIterator))

        print('Calculating read distance arrays:')            
        # I'm making this distance array. In this array, a 0 represents a Match.  a 1 represents indels or substitutions.
        # This way I can calculate "distance" in an arbitrary number of dimensions
        # Distance is a euclidian way to represent how far away a read is from the consensus,
        # based on the heterozygous positions.  Each heterozygous position is a "dimension" in this space
        distanceArrays = {}
        for readID in readIDs:
            # TODO: A Bug! Initializing this list as 0s will bias the results.
            # TODO: Pileupcolumn loop is not hitting each read. Only...half sometimes. Some reads are not analyzed.
            # Why? SPOTTED IT! bamfile.pileup has a default to maximum read depth of 8000

            #distanceArrays[readID] = list([999] * len(self.snps))
            distanceArrays[readID] = list([0] * len(self.snps))


        # I spotted the bug!!! pileup defaults to maximum 8000 read depth. That's bad!.
        pileupIterator = bamfile.pileup(alignmentRef.id,max_depth=99999999)
        #pileupIterator = bamfile.pileup(alignmentRef.id)
        for pileupColumn in pileupIterator:
            currentColumn = pileupColumn.pos
            
            # Only do this if the column number exists in our list of heterozygous positions
            if currentColumn in self.snps:
                
                heterozygousPositionIndex = self.snps.index(currentColumn)
                currentAnalyzedReadCount = 0 # A debugging variable, i dont think I actually use this count.
                
                referenceBase = alignmentRef.seq[currentColumn].upper()
                for pileupRead in pileupColumn.pileups:
                    currentAnalyzedReadCount += 1
                    readID = pileupRead.alignment.query_name
                    
                    #print('Pos:' + str(currentColumn) + ', Refbase:' + str(referenceBase) + ', Read:' + str(readID))
                    
                    # In this model, the distance is either 0 or 1. This was intentional but
                    # Maybe we can tune the algorithm using these distances.
                    # This could actually be tuned to do the heterozygous split using ONLY snps.
                    # TODO: if we're having problems splitting based on homopolymers check this spot.
                    # Maybe, I want to count indels as 0, no distance.
                    # TODO: Something to try: indels are -1. SNPS are 1. Match = 0
                    # Maybe that would help the sorting?
                    # TODO: Newest idea. Default to 0. 1 is match, -1 is indels. -1 is mismatches. I think that's it.
                    
                    if(pileupRead.is_del == 1):
                        distanceArrays[readID][heterozygousPositionIndex] = -1
                    elif(pileupRead.indel > 0):
                        distanceArrays[readID][heterozygousPositionIndex] = -1
                    else:   
                        currentBase = pileupRead.alignment.query_sequence[pileupRead.query_position].upper()  
                        if(currentBase == referenceBase):
                            #print('Assinging Match. Column=' + str(currentColumn) + ', CurrentBase:' + str(currentBase) + ', HeterozygousPosIndex=' + str(heterozygousPositionIndex))
                            distanceArrays[readID][heterozygousPositionIndex] = 1
                        else:
                            distanceArrays[readID][heterozygousPositionIndex] = -1

                print('At position ' + str(heterozygousPositionIndex + 1) + ' I analyzed ' + str(currentAnalyzedReadCount) + ' reads.')

        self.printDistanceArrays(distanceArrays, join(self.heterozygousDirectory, 'DistanceArrays.csv'))

        # TODO: Im making 3 clusters. that worked. I need to make a parameter for cluster count.
        clusteredReadIDs = self.clusterReads(distanceArrays, 2)

        # Dictionary of results to return. Key is location of the consensus sequence.
        # Value is the # of reads represented in this consensus alignment.
        coverageResults = {}

        for zeroBasedClusterIndex, readCluster in enumerate(clusteredReadIDs):
            # I want to call the Strand (1 and 2), not Strand (0 and 1).
            clusterIndex = zeroBasedClusterIndex + 1

            clusteredReadIDs = readCluster.keys()

            clusterOutputDir = join(self.outputRootDirectory, 'Strand' + str(clusterIndex) + 'ClusteredReads')

            distanceArrayFileName = join(clusterOutputDir, 'Strand' + str(clusterIndex) + 'DistanceArrays.csv')
            self.printDistanceArrays(readCluster, distanceArrayFileName)

            readOutputFileName = join(clusterOutputDir, 'Strand' + str(clusterIndex) + 'Reads.' + self.readInputFormat)
            readOutputFile = createOutputFile(readOutputFileName)

            # Loop parsed reads, grab reads belonging to this cluster.
            # FYI it looks like each input is clustered in the output, i haven't found a missing read yet. I should still check.
            for readObject in parsedReads:

                #print ('ReadClusterIndex=' + str(zeroBasedClusterIndex))
                #print ('AllReadID=' + str(readObject.id))

                for clusteredReadID in clusteredReadIDs:
                    #print ('clusteredReadID=' + str(clusteredReadID))

                    if (readObject.id == clusteredReadID):
                        write([readObject], readOutputFile, self.readInputFormat)
                        break

            readOutputFile.close()

            currentWranglerObject = AlleleWrangler(
                readOutputFileName
                , join(self.outputRootDirectory, 'Strand' + str(clusterIndex) + 'Alignment')
                , join(self.heterozygousDirectory, 'AlignmentReference.fasta')
                , 6
                , self.numberThreads
                , False
                , self.snps)
            currentCoverageResults = currentWranglerObject.analyzeReads()

            # Merge the dictionaries of coverage values ane return them.
            for key in currentCoverageResults.keys():
                coverageResults[key] = currentCoverageResults[key]

        print ('Done Phasing Reads.')
        return coverageResults


    def printDistanceArrays(self, distanceArrays, fileLocation):
        distanceArrayOutputFile = createOutputFile(fileLocation)
        readIDs = distanceArrays.keys()

        for readID in readIDs:
            distanceArrayOutputFile.write(readID + ',')
            for distanceValue in distanceArrays[readID]:
                distanceArrayOutputFile.write(str(distanceValue) + ',')
            distanceArrayOutputFile.write('\n')

    def clusterReads(self, dataInput, numClusters):
        readIDs = list(dataInput.keys())
        readDistanceArrays=asarray([dataInput[k] for k in readIDs if k in dataInput])
        
        print('I will attempt to cluster and separate heterozygous reads.')
        print('There are ' + str(len(readDistanceArrays)) + ' distance arrays.')
        print('They are of length ' + str(len(readDistanceArrays[0])))
        
        # TODO: I can choose number of jobs/threads.
        kmeansObject = KMeans(n_clusters= numClusters)
        predictedLabels = kmeansObject.fit_predict(readDistanceArrays)
        
        print ('Done Clustering. Here are the predicted labels:\n' + str(predictedLabels))
        
        # Get read IDs from each cluster and populate arrays
        # Store them as their own dictionary instead of arrays. I hope that didn't break much.

        #strand1ReadIDs = {}
        #strand2ReadIDs = {}

        clusteredReadIDs = []
        for i in range(0,numClusters):
            clusteredReadIDs.append({})


        # Each read has been assigned a label, for which cluster it should appear in. Loop through them.
        for labelIndex, label in enumerate(predictedLabels):
            readID = readIDs[labelIndex]

            #print ('LabelIndex=' + str(labelIndex))
            #print ('ReadID=' + str(readID))
            #print ('Length of clusteredReadIDs:' + str(len(clusteredReadIDs)))

            # Loop through my clusters, assign each read to one of the clusters.
            labelFound = False
            for clusterIndex in range (0,numClusters):
                #print ('ClusterIndex=' + str(clusterIndex))
                if(label == clusterIndex):
                    labelFound = True
                #    print ('Label Found.')
                    #print ()

                    clusteredReadIDs[clusterIndex][readID] = dataInput[readID]

            if(not labelFound):
                raise Exception('This read ( ' + str(readID) + ' ) was not clustered successfully, why did clustering break? It has a label of: ' + str(label))


        return clusteredReadIDs

        #return strand1ReadIDs, strand2ReadIDs
            
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


            alignReads(currentReferenceSequence,self.readInput,currentIterationSubdirectory, True)
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
