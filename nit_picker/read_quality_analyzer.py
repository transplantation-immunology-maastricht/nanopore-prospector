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


from os import makedirs, system
#from .minion_read_collection import minionReadCollection
#from Bio.Seq import Seq
from Bio.SeqIO import parse, write
#from os.path import join

from os.path import join, isdir

#from random import shuffle

#from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Sequencing.Applications import BwaIndexCommandline

#from .quality_statistics import QualityStatistics

from nanopore_prospector.common import createOutputFile, getReadFileType

from pysam import AlignmentFile


class QualityStatistics:
# This class is confusing because I'm using it for both:
# Read quality analysis, relative to a reference, and 
# allele polymorphism. How are different alleles polymorphic relative to the reference?
    
    
    def __init__(self, referenceSequence):
        self.referenceSequence = referenceSequence
        self.referenceLength = len(referenceSequence)
        
        # Storing the new consensus as a list of characters. So i can modify it.
        #self.newConsensus = list(referenceSequence)
        self.newConsensus = referenceSequence
        
        self.alignedReadCountsByPosition   = [0]*self.referenceLength
        # For this application, I assume that a read will always be a match, mismatch,insertion, or deletion. Never both. An indel is not a match.
        self.matchReadCountsByPosition     = [0]*self.referenceLength
        self.insertionReadCountsByPosition = [0]*self.referenceLength
        self.deletionReadCountsByPosition  = [0]*self.referenceLength
        self.mismatchReadCountsByPosition  = [0]*self.referenceLength
        self.homopolymerSizesByPosition    = [0]*self.referenceLength
        self.segmentCountByPosition        = [0]*self.referenceLength
        self.supplementaryReadCountsByPosition = [0]*self.referenceLength
        self.snpCallsByPosition = [None]*self.referenceLength
        
        self.ACountsByPosition = [0]*self.referenceLength
        self.GCountsByPosition = [0]*self.referenceLength
        self.CCountsByPosition = [0]*self.referenceLength
        self.TCountsByPosition = [0]*self.referenceLength
        
        # I want to keep track of what alleles / reads have what polymorphisms. 
        # This is a dictionary. 
        # Key = allele id.
        # Value = a dictionary.
            # key = position, from the alignment reference
            # value = alternate base(s) from reference. 
            # only have a key/value pair if there is a deviation from reference.        
        self.alleleSpecificPolymorphisms = {}

        self.totalReadCount = 0

    def storeReadHomopolymerLength(self, referencePosition, referenceBase, expectedLength, observedLength):
        # Store the statistic of a single read, how much we expect, how much we observed.
        pass

    def getReadHomopolymerLengths(self, referenceBase, expectedLength):
        # For an expected homopolymer length, what was found in the reads?
        # return an array of all homopolymer read lengths for an expected homopolymer.
        pass

    def writePolymorphicAlleles(self, alignmentSummaryOutputFilename):

        polymorphicAllelesOutputFile = createOutputFile(alignmentSummaryOutputFilename)

        # This output file has positionZeroBased# in the column headers.
        # The allele names are in the row headers.
        distinctAlleleNames = sorted(self.alleleSpecificPolymorphisms.keys())
        distinctPolymorphicPositions = []
        for alleleName in distinctAlleleNames:
            for positionZeroBased in self.alleleSpecificPolymorphisms[alleleName].keys():
                if positionZeroBased not in distinctPolymorphicPositions:
                    distinctPolymorphicPositions.append(positionZeroBased)
        distinctPolymorphicPositions = sorted(distinctPolymorphicPositions)
        
        # print header, these are the 1-based positions to match IGV
        # The first column is empty, that column is for allele/read ids.
        polymorphicAllelesOutputFile.write('Reference Position 1-based')
        for positionZeroBased in distinctPolymorphicPositions:
            polymorphicAllelesOutputFile.write(',' + str(positionZeroBased + 1))
        polymorphicAllelesOutputFile.write('\n')
        
        # Print the reference bases underneath the position.
        polymorphicAllelesOutputFile.write('Reference Base')
        for positionZeroBased in distinctPolymorphicPositions:
            polymorphicAllelesOutputFile.write(',' + str(self.referenceSequence[positionZeroBased]))
        polymorphicAllelesOutputFile.write('\n')
            
        # print a line for each "allele"
        for alleleName in distinctAlleleNames:
            polymorphicAllelesOutputFile.write(alleleName)
            for positionZeroBased in distinctPolymorphicPositions:
                # If this read has a polymorphism at this position
                if(positionZeroBased in self.alleleSpecificPolymorphisms[alleleName].keys()):
                    polymorphicAllelesOutputFile.write(',' + str(self.alleleSpecificPolymorphisms[alleleName][positionZeroBased]))
                else:
                    polymorphicAllelesOutputFile.write(',')
                 
            polymorphicAllelesOutputFile.write('\n')
            
        polymorphicAllelesOutputFile.close()

    def writeLDMatrix(self, ldMatrixFilename):
        
        # The level of linkage disequilibrium between A and B can be quantified by the coefficient of linkage disequilibrium {\displaystyle D_{AB}} D_{AB}, which is defined as
        # D_{AB}=p_{AB}-p_{A}p_{B},
        # https://en.wikipedia.org/wiki/Linkage_disequilibrium
        # In the case {D_{{AB}}=0, said to be in linkage equilibrium
        print ('Writing a Linkage Disequilibrium Matrix')

        # Quick way to identify the polymorphic positions
        # TODO I already did this in writePolymorphicAlleles, maybe I can store it and reuse.
        distinctAlleleNames = sorted(self.alleleSpecificPolymorphisms.keys())
        distinctPolymorphicPositions = []
        for alleleName in distinctAlleleNames:
            for positionZeroBased in self.alleleSpecificPolymorphisms[alleleName].keys():
                if positionZeroBased not in distinctPolymorphicPositions:
                    distinctPolymorphicPositions.append(positionZeroBased)
        distinctPolymorphicPositions = sorted(distinctPolymorphicPositions)

        #ldMatrixUnscaled=[len(distinctPolymorphicPositions) , len(distinctPolymorphicPositions)]
        ldMatrixUnscaled = [[0 for x in range(len(distinctPolymorphicPositions))] for y in range(len(distinctPolymorphicPositions))] 
        
        #Maximum value for scaling?
        maxLDValue = 0      
        
        # Loop through snps in 2d array, calculate LD positions, write position in first column.
        for snpYIndex, snpYPosition in enumerate(distinctPolymorphicPositions):
                    
            # Loop through snps, write LD calculation in corresponding column.
            for snpXIndex, snpXPosition in enumerate(distinctPolymorphicPositions):
                if(snpYPosition == snpXPosition):
                    ldMatrixUnscaled[snpXIndex][snpYIndex] = 0
                    #ldMatrixOutputFile.write(',' + str(0))
                    
                # Checking this allows me to calculate the LD in both directions at the same time.
                # To reduce redundancy.
                elif(snpYPosition > snpXPosition):
                    # Calculate LD.
                    # Coefficient of Linkage Disequilibrium D{AB}
                    # D{AB} = p{AB} - p{A} * p{B}
                    # Define the allele as sequences that matches the reference.
                    # Any insert/delete/mismatch is treated as the "other" allele.
                    
                    # p{A}
                    pA = 1.0 * self.matchReadCountsByPosition[snpYPosition] / self.alignedReadCountsByPosition[snpYPosition]
                    
                    # p{B}
                    pB = 1.0 * self.matchReadCountsByPosition[snpXPosition] / self.alignedReadCountsByPosition[snpXPosition]
                    
                    # p{AB}
                    # Must look at specific alleles to calculate this. 
                    # Loop thru alleles and count if BOTH of them match.
                    matchBothCount = 0
                    for alleleName in self.alleleSpecificPolymorphisms.keys():
                        # alleleSpecificPolymorphisms only tracks snps. 
                        # If there is no dictionary entry for either,
                        # that means the allele matches the reference at this position.
                        if (snpYPosition not in self.alleleSpecificPolymorphisms[alleleName] and snpXPosition not in self.alleleSpecificPolymorphisms[alleleName] ):
                            matchBothCount += 1
                    # TODO: Am i biasing this by dividing by aligned read counts at X position?
                    # Maybe. Only matters for snps with low coverage.....
                    pAB = 1.0 * matchBothCount / self.alignedReadCountsByPosition[snpXPosition]
                    
                    # Negative LD values mean that they're linked, but in the opposite direction.
                    # So it is fair to use a negative linkage disequilbrium.
                    # The allele in LD calculations is defined as "matching the reference" or not.
                    # any insert, delete, or mismatch is the "other allele"
                    ldValue = abs(pAB - (pA * pB))                   
                    ldMatrixUnscaled[snpXIndex][snpYIndex] = ldValue
                    ldMatrixUnscaled[snpYIndex][snpXIndex] = ldValue
                    
                    #Maximum value for scaling?
                    if (ldValue > maxLDValue):
                        maxLDValue = ldValue
                  
                else:
                    # Don't need to calculate the LD a second time.
                    pass
                    
            #ldMatrixOutputFile.write('\n')
         
        # Scale the LD values.
        ldMatrixScaled = [[0 for x in range(len(distinctPolymorphicPositions))] for y in range(len(distinctPolymorphicPositions))]
        #ldMatrixScaled=[len(distinctPolymorphicPositions) , len(distinctPolymorphicPositions)]

        for snpYIndex, snpYPosition in enumerate(distinctPolymorphicPositions):
            for snpXIndex, snpXPosition in enumerate(distinctPolymorphicPositions):
                ldMatrixScaled[snpXIndex][snpYIndex] = ldMatrixUnscaled[snpXIndex][snpYIndex] / maxLDValue

        # Write calculated LD values to matrix file.   
        ldMatrixOutputFile = createOutputFile(ldMatrixFilename)   
        
        # Loop through snps, write their positions on the first row.  Skip one.
        # Use 1-based positions.        
        for positionZeroBased in distinctPolymorphicPositions:
            ldMatrixOutputFile.write(',' + str(positionZeroBased + 1))
        ldMatrixOutputFile.write('\n')

        for snpYIndex, snpYPosition in enumerate(distinctPolymorphicPositions):
            ldMatrixOutputFile.write(str(snpYPosition + 1))        
        
            for snpXIndex, snpXPosition in enumerate(distinctPolymorphicPositions):
                ldMatrixOutputFile.write(',' + str(round(ldMatrixScaled[snpXIndex][snpYIndex], 2)))
                
            ldMatrixOutputFile.write('\n')        
            
        ldMatrixOutputFile.close()


    def writeAlignmentSummaryHeader(self, alignmentSummaryOutputFile):
        return alignmentSummaryOutputFile.write(
            'Reference_Position_1_based' 
            + ',Reference_Base' 
            + ',SNP_base'
            + ',Aligned_Reads' 
            + ',Segment_Count'
            + ',Supplementary_Read_Almts'
            + ',Is_Polymorphic'
            + ',A_Count'
            + ',G_Count'
            + ',C_Count'
            + ',T_Count' 
            + ',Match_Count' 
            + ',Match_Percent' 
            + ',Insertion_Count' 
            + ',Insertion_Percent' 
            + ',Deletion_Count' 
            + ',Deletion_Percent' 
            + ',Mismatch_Count' 
            + ',Mismatch_Percent' 
            + ',Homopolymer_Size\n')

    def callSNP(self, referencePositionZeroBased, referenceBase):
        # Assign an alternate base if there is a snp at this position.
        # Do nothing if there is no SNP here.
        
        # Declaring a constant here so it's easier to change. .90 is arbitrary
        # TODO: Should this be a commandline param? Yes.
        # .9999999 is useful for identifying any position that might be polymorphic sometimes.
        # Like if we are searching for a single snp in consensus sequences.
        # .90 is probably more useful for minion reads
        #polymorphicBaseCutoff = .99999999
        
        # This method is not perfect, we're calling some basecalling errors as SNPS and homopolymers.
        
        polymorphicBaseCutoff = .80



        alignedReads = self.alignedReadCountsByPosition[referencePositionZeroBased]
            
        if(alignedReads == 0):
            # this is a hack to avoid divide by zero errors. Dividing zero by 1 should give a zero.
            matchPercent = 0
            #alignedReads = 1            
        else:
            matchPercent = self.matchReadCountsByPosition[referencePositionZeroBased] / alignedReads

        # Is there an insertion? Insertion count > match count * cutoff
        if(self.insertionReadCountsByPosition[referencePositionZeroBased]  > self.matchReadCountsByPosition[referencePositionZeroBased] * polymorphicBaseCutoff):
            # TODO: Calculate the inserted text based on insertion length. For now i will call this SNP an I.
            self.snpCallsByPosition[referencePositionZeroBased] = 'I'
            
            # TODO Change the consensus sequence for this insertion?
            
        # Is there a deletion? Insertion count > match count * cutoff
        if(self.deletionReadCountsByPosition[referencePositionZeroBased]  > self.matchReadCountsByPosition[referencePositionZeroBased] * polymorphicBaseCutoff):
            self.snpCallsByPosition[referencePositionZeroBased] = '-'

            # TODO Delete this base from the consensus sequence
     
        # If this position is polymorphic, also print this to the polymorphic positions output file.
        if( self.alignedReadCountsByPosition[referencePositionZeroBased] > 0 and matchPercent < polymorphicBaseCutoff):
            # Not an indel, this is probably a mismatch.

            maxBaseCount = max(
                self.ACountsByPosition[referencePositionZeroBased]
                ,self.GCountsByPosition[referencePositionZeroBased]
                ,self.CCountsByPosition[referencePositionZeroBased]
                ,self.TCountsByPosition[referencePositionZeroBased]                
                )
            
            if (self.ACountsByPosition[referencePositionZeroBased] == maxBaseCount):
                self.snpCallsByPosition[referencePositionZeroBased] = 'A'
            elif (self.GCountsByPosition[referencePositionZeroBased] == maxBaseCount):
                self.snpCallsByPosition[referencePositionZeroBased] = 'G'
            elif (self.CCountsByPosition[referencePositionZeroBased] == maxBaseCount):
                self.snpCallsByPosition[referencePositionZeroBased] = 'C'
            elif (self.TCountsByPosition[referencePositionZeroBased] == maxBaseCount):
                self.snpCallsByPosition[referencePositionZeroBased] = 'T'
            else:
                raise Exception('I made a mistake calculating the maximum base count for SNP detection.')
            
            
            # TODO THis is a very rough way to get a new consensus sequence with SNPs. Fix this logic.
            
            # If we're not in the deletion regions...
            # TODO: Feed this as parameters somehow. I just don't want to analyze the regions with the complicated repeats.
            if 4550 <= referencePositionZeroBased <= 4695:
                # We are in the complicated delete region, lets not assign this snp.
                pass
            
            elif 4909 <= referencePositionZeroBased <= 4918:
                # This is the shorter deletion.
                pass
            else:
            
                # manipulating strings because they are immutable
                #tempid = self.newConsensus.id
                tempConsensus = list(self.newConsensus)
                tempConsensus[referencePositionZeroBased] = self.snpCallsByPosition[referencePositionZeroBased]
                self.newConsensus = "".join(tempConsensus)
                #self.newConsensus.id = 'GeneratedConsensus'
                # Put the SNP in our new consensus sequence
                #self.newConsensus[referencePositionZeroBased] = self.snpCallsByPosition[referencePositionZeroBased]
                
            #self.snpCallsByPosition[referencePositionZeroBased] = '!'
    
        
        pass

    def writeSNPs(self, snpSummaryOutputFilename):
        snpSummaryOutputFile = createOutputFile(snpSummaryOutputFilename)
        
        row1Text = 'Pos,'
        row2Text = 'Ref,'
        row3Text = 'SNP,'

        for referencePositionZeroBased, referenceBase in enumerate(self.referenceSequence):
            if(self.snpCallsByPosition[referencePositionZeroBased] is not None):
                row1Text += str(referencePositionZeroBased + 1) + ','
                row2Text += referenceBase + ','
                row3Text += self.snpCallsByPosition[referencePositionZeroBased] + ','
            else:
                pass

        snpSummaryOutputFile.write(row1Text + '\n')
        snpSummaryOutputFile.write(row2Text + '\n')
        snpSummaryOutputFile.write(row3Text + '\n')
        snpSummaryOutputFile.close()

    def writeAlignmentSummary(self, alignmentSummaryOutputFilename):
        alignmentSummaryOutputFile = createOutputFile(alignmentSummaryOutputFilename)
        polymorphicPositionsOutputFile = createOutputFile(alignmentSummaryOutputFilename.replace('.csv','.PolymorphicPositions.csv'))
        
        newConsensusOutputFile = createOutputFile(alignmentSummaryOutputFilename.replace('.csv','.NewConsensus.fasta'))        
        
      
        #createOutputFile(outputfileName)
        self.writeAlignmentSummaryHeader(alignmentSummaryOutputFile)
        self.writeAlignmentSummaryHeader(polymorphicPositionsOutputFile)


        
        for referencePositionZeroBased, referenceBase in enumerate(self.referenceSequence):
            
            self.callSNP(referencePositionZeroBased, referenceBase)

            #print('analyzing base # ' + str(referencePositionZeroBased))
            #print ('the coverage should be:' + str(self.alignedReadCountsByPosition[referencePositionZeroBased]))
            #print ('the reference base should be:' + str(referenceBase))
            #alignedReads = self.alignedReadCountsByPosition[referencePositionZeroBased]
            
            #if(alignedReads == 0):
                # this is a hack to avoid divide by zero errors. Dividing zero by 1 should give a zero.
            #    matchPercent = 0
                #alignedReads = 1
                
            #else:
            #    matchPercent = self.matchReadCountsByPosition[referencePositionZeroBased] / alignedReads
            
            self.writeAlignmentSummaryLine(alignmentSummaryOutputFile, referencePositionZeroBased, referenceBase, (self.snpCallsByPosition[referencePositionZeroBased] is not None))
            
            # If this position is polymorphic, also print this to the polymorphic positions output file.
            #if( self.alignedReadCountsByPosition[referencePositionZeroBased] > 0 and matchPercent < polymorphicBaseCutoff):
            #    self.writeAlignmentSummaryLine(polymorphicPositionsOutputFile, referencePositionZeroBased, referenceBase, (matchPercent < polymorphicBaseCutoff))
                
            if(self.snpCallsByPosition[referencePositionZeroBased] is not None):
                self.writeAlignmentSummaryLine(polymorphicPositionsOutputFile, referencePositionZeroBased, referenceBase, True)
    
        alignmentSummaryOutputFile.close()
        polymorphicPositionsOutputFile.close()

        # Write new consensus
        #write([self.newConsensus], newConsensusOutputFile, 'fasta')  
        newConsensusOutputFile.write('>GeneratedConsensusSequence\n')
        newConsensusOutputFile.write(self.newConsensus + '\n')
        newConsensusOutputFile.close()
                
        self.writeSNPs(alignmentSummaryOutputFilename.replace('.csv','.SNPs.csv'))

    

    def writeAlignmentSummaryLine(self, outputFile, referencePositionZeroBased, referenceBase, isPolymorphic):
        alignedReads = self.alignedReadCountsByPosition[referencePositionZeroBased]
        if(self.alignedReadCountsByPosition[referencePositionZeroBased] == 0):
            # this is a hack to avoid divide by zero errors. Dividing zero by 1 should give a zero.
            matchPercent = 0
            alignedReads = 1
                
        else:
            matchPercent = self.matchReadCountsByPosition[referencePositionZeroBased] / alignedReads
                    
        outputFile.write(
            str(referencePositionZeroBased + 1) +
            ',' + str(referenceBase) +
            ',' + str(self.snpCallsByPosition[referencePositionZeroBased]) + 
            # What is the difference? pysam reports an aligned read count. I'm also counting them myself.
            ',' + str(self.alignedReadCountsByPosition[referencePositionZeroBased]) + 
            ',' + str(self.segmentCountByPosition[referencePositionZeroBased]) + 
            ',' + str(self.supplementaryReadCountsByPosition[referencePositionZeroBased]) + 
            ',' + ('YES' if isPolymorphic else '') +
            ',' + str(self.ACountsByPosition[referencePositionZeroBased]) + 
            ',' + str(self.GCountsByPosition[referencePositionZeroBased]) + 
            ',' + str(self.CCountsByPosition[referencePositionZeroBased]) + 
            ',' + str(self.TCountsByPosition[referencePositionZeroBased]) + 
            ',' + str(self.matchReadCountsByPosition[referencePositionZeroBased]) +
            ',' + str(matchPercent) +
            ',' + str(self.insertionReadCountsByPosition[referencePositionZeroBased]) +
            ',' + str(self.insertionReadCountsByPosition[referencePositionZeroBased] / alignedReads) +
            ',' + str(self.deletionReadCountsByPosition[referencePositionZeroBased]) +
            ',' + str(self.deletionReadCountsByPosition[referencePositionZeroBased] / alignedReads) +
            ',' + str(self.mismatchReadCountsByPosition[referencePositionZeroBased]) +
            ',' + str(self.mismatchReadCountsByPosition[referencePositionZeroBased] / alignedReads) +
            ',' + str(self.homopolymerSizesByPosition[referencePositionZeroBased]) +
            '\n' 
            )
 
    def storeAlleleSpecificPolymorphism(self, alleleID, baseIndexZeroBased, substituteBases ):
        #print ('allele ' + alleleID + ' has a substitute base at position ' + str(baseIndexZeroBased + 1) + ':' + substituteBases)
        
        # If we haven't seen this allele before
        if alleleID not in self.alleleSpecificPolymorphisms.keys():
            # Set up a dictionary of it's polymorphisms.
            self.alleleSpecificPolymorphisms[alleleID] = {}
            
        # Store this polymorphic base in this allele's dictionary.
        self.alleleSpecificPolymorphisms[alleleID][baseIndexZeroBased] = substituteBases
        

def calculateTotalCoverage(alignmentOutputName):
    # The number of reads that are involved in this alignment
    bamfile = AlignmentFile(alignmentOutputName, 'rb')
    alignedReads = list(bamfile.fetch())
    return len(alignedReads)
 
 
def calculateReadQualities(referenceSequencelocation, readFileLocation, outputDirectory, sampleID, numberThreads):
    print ('Calculating Read Quality Statistics.')
    print ('I have this reference File:\n' + str(referenceSequencelocation))
    
    alignmentOutputDir = join(outputDirectory, sampleID + '_ReadQualityAnalysis')    
    # Copy reference to outputdir
    referenceSequenceTargetLocation = join(alignmentOutputDir, 'AlignmentReference.fasta')
    # Store all the stats in a handy object.
    alignmentRef = list(parse(referenceSequencelocation, 'fasta'))[0]
    qualStats = QualityStatistics(alignmentRef.seq)
    #qualStats.referenceSequence = parse(referenceSequencelocation, 'fasta')
    write(alignmentRef, createOutputFile(referenceSequenceTargetLocation), 'fasta')  
    
    # Align Reads against Reference
    alignReads(referenceSequenceTargetLocation, readFileLocation, alignmentOutputDir, numberThreads)
    
    # Open Alignment
    alignmentRef = list(parse(referenceSequenceTargetLocation, 'fasta'))[0]
    bamfile = AlignmentFile(join(alignmentOutputDir,'alignment.bam'), 'rb')

    # Declare variables to store the information.
    qualStats.totalReadCount = len(list(parse(readFileLocation, getReadFileType(readFileLocation))))
    
    # Initialize the read dictionary, so we don't lose track of any of these read / alleles
    for readObject in list(parse(readFileLocation, getReadFileType(readFileLocation))):
        if readObject.id not in qualStats.alleleSpecificPolymorphisms.keys():
            qualStats.alleleSpecificPolymorphisms[readObject.id] = {}

    # Loop alignment:     
    # Iterate the reference sequence column by column.
    pileupIterator = bamfile.pileup(alignmentRef.id)
    
    for pileupColumn in pileupIterator:

        referencePositionZeroBased = pileupColumn.reference_pos  
        referenceBase = str(alignmentRef[referencePositionZeroBased])
        # think this is wrong. It seems like nsegments 
        # aligned reads + deletions.
        qualStats.segmentCountByPosition[referencePositionZeroBased] = pileupColumn.nsegments      
        
        # Progress bar...
        if (referencePositionZeroBased%250 == 0):
            print ('Checking aligned reads at genomic position: (' + str(referencePositionZeroBased) + '/' + str(qualStats.referenceLength) + ')')
            print (str(pileupColumn.nsegments) + ' segments found at that position.')
        
        #print ('at 0-based pos:' + str(referencePositionZeroBased) + ' i see a ' + str(alignmentRef[referencePositionZeroBased]) + ' with  this many reads:' + str(alignedCount))
        
        #qualStats.alignedReadCountsByPosition[referencePositionZeroBased] = alignedCount

        # Search for homopolymers in the reference sequence.
        # If a homopolymer begins at this position, store the LENGTH of the homopolymer.
        # TODO: This is a slight bug. I'm only measuring homopolymers where there are reads aligned. 
        # Homopolymers outside aligned regions are ignored. That's probably fine.
                
        # The most 3' base in the reference. A homopolymer cannot begin at the rightmost base.
        if(referencePositionZeroBased == qualStats.referenceLength - 1): 
            pass 
        else:
            
            # Check the most 5' base in the reference. Don't look to the left, we will go out of bounds.
            if(referencePositionZeroBased == 0): 
                pass  
            
            # If the current base matches the base to the left, a homopolymer cannot begin here.
            elif(referenceBase == str(alignmentRef[referencePositionZeroBased - 1])):
                pass # Nothing to do.
            
            else:
            
                # Check the bases to the right, to determine homopolymer length. (don't exceed length).                
                homopolymerLengthIterator = 1
                while (referencePositionZeroBased + homopolymerLengthIterator < qualStats.referenceLength):
                    if (referenceBase == str(alignmentRef[referencePositionZeroBased + homopolymerLengthIterator])):
                        
                        # Record homopolymer length.
                        qualStats.homopolymerSizesByPosition[referencePositionZeroBased] = homopolymerLengthIterator + 1
                        homopolymerLengthIterator += 1
                        
                    else:
                        break
        
        
        qualStats.alignedReadCountsByPosition[referencePositionZeroBased] = 0
        # Iterate the Reads at this position      
        # We use the pileup method to iterate through column-by-column.     
        for pileupRead in pileupColumn.pileups:
            
            
            
            # Supplementary alignments are part of the minimap algorithm. I dont want to use them for quality analysis.
            if(pileupRead.alignment.is_supplementary):
                qualStats.supplementaryReadCountsByPosition[referencePositionZeroBased] += 1
            
                pass

            # If this read is a deletion
            elif(pileupRead.is_del == 1):
                qualStats.deletionReadCountsByPosition[referencePositionZeroBased] += 1
                qualStats.storeAlleleSpecificPolymorphism(pileupRead.alignment.query_name, referencePositionZeroBased, '-' )
            # else if this read is an insertion
            elif(pileupRead.indel > 0):
                # This apparently indicates an insertion AFTER this base. The query text at this position is the reference
                # , and there are inserted bases after this position
                # Indexing the query sequence is tricky.
                # I am also recording the correct base, because I hope that is more clear, and doesn't just look like a base substitution.
                qualStats.insertionReadCountsByPosition[referencePositionZeroBased] += 1
                
                insertedBasesCount = pileupRead.indel
                insertedBases = pileupRead.alignment.query_sequence[pileupRead.query_position:pileupRead.query_position+(insertedBasesCount+1)].upper() 
                qualStats.storeAlleleSpecificPolymorphism(pileupRead.alignment.query_name, referencePositionZeroBased, insertedBases )
                
                queryPosition = int(pileupRead.query_position)
                # Because....Insertions still have aligned bases....which can match the reference
                # TODO: SHould I actually be checking for insertions in the bases BEFORE this position?
                # This could be a major misunderstanding on my part....
                currentBase = pileupRead.alignment.query_sequence[queryPosition].upper() 
                
                if(currentBase.upper() == referenceBase.upper()):
                    qualStats.matchReadCountsByPosition[referencePositionZeroBased] += 1
                else:
                    qualStats.mismatchReadCountsByPosition[referencePositionZeroBased] += 1
                    qualStats.storeAlleleSpecificPolymorphism(pileupRead.alignment.query_name, referencePositionZeroBased, currentBase.upper() )
                
                if(currentBase =='A'):
                    qualStats.ACountsByPosition[referencePositionZeroBased] += 1
                elif(currentBase =='G'):
                    qualStats.GCountsByPosition[referencePositionZeroBased] += 1
                elif(currentBase =='C'):
                    qualStats.CCountsByPosition[referencePositionZeroBased] += 1
                elif(currentBase =='T'):
                    qualStats.TCountsByPosition[referencePositionZeroBased] += 1
                
                
                
            # Else if it is a refskip (TODO: What does this mean? no read aligned? Count these?) Doesn't happen often.
            elif(pileupRead.is_refskip):
                print('This read is a refskip, i dont know what that means:' + pileupRead.alignment.query_name)
                raise Exception('This read is a refskip, i dont know what that means:' + pileupRead.alignment.query_name)
            
            # else this means we have a base aligned at this position for this read.
            # It is either a match or a mismatch.
            else:    

                qualStats.alignedReadCountsByPosition[referencePositionZeroBased] += 1
                # TODO: Sometimes, there is no query position. What if this is 'None'?
                # This is because in the case of indels, there is no query position.
                
                queryPosition = int(pileupRead.query_position)
                
                currentBase = pileupRead.alignment.query_sequence[queryPosition].upper()                    

                if(currentBase.upper() == referenceBase.upper()):
                    qualStats.matchReadCountsByPosition[referencePositionZeroBased] += 1
                else:
                    qualStats.mismatchReadCountsByPosition[referencePositionZeroBased] += 1
                    qualStats.storeAlleleSpecificPolymorphism(pileupRead.alignment.query_name, referencePositionZeroBased, currentBase.upper() )
                    
                if(currentBase =='A'):
                    qualStats.ACountsByPosition[referencePositionZeroBased] += 1
                elif(currentBase =='G'):
                    qualStats.GCountsByPosition[referencePositionZeroBased] += 1
                elif(currentBase =='C'):
                    qualStats.CCountsByPosition[referencePositionZeroBased] += 1
                elif(currentBase =='T'):
                    qualStats.TCountsByPosition[referencePositionZeroBased] += 1

            # TODO: I just put a false in this if statement to "comment out" this section.
            # Figure out the homopolymers here!!!! polish.
            if (False and qualStats.homopolymerSizesByPosition[referencePositionZeroBased] >= 3):
            #if pileupcolumn.pos == 14693:
                print ('Homopolymer. (1)-based position:' + str(referencePositionZeroBased + 1) + ', Base:' + str(referenceBase) + ', Length:' + str(qualStats.homopolymerSizesByPosition[referencePositionZeroBased]))
                # TODO: Im working on homopolymers here.
                readID = pileupRead.alignment.query_name
                
                
                readHomopolymerSequence = pileupRead.alignment.query_sequence[queryPosition - 5: queryPosition + 5].upper() 
                #isDeletion = (pileupRead.is_del == 1)
                
                # TODO: I left off here. This print statement does not work.
                print('Deletion / Insertion: ' + str(pileupRead.is_del) + ' / ' ) 
                #qualStats.deletionReadCountsByPosition[referencePositionZeroBased] += 1
                # else if this read is an insertion
                # elif(pileupRead.indel > 0):
                
                
                
                
                print (readID + ' : ' + str(readHomopolymerSequence))
      
    qualStats.writeAlignmentSummary(join(alignmentOutputDir, 'AlignmentSummary.csv'))  
    qualStats.writePolymorphicAlleles(join(alignmentOutputDir, 'AlleleSpecificPolymorphisms.csv'))
    # TODO: I commented this part out because it takes X amount of time.
    #qualStats.writeLDMatrix(join(alignmentOutputDir, 'LDMatrix.csv'))
      
      
"""
        # What do we do if the reference is a homopolymer. Minimum length is 3.
        if (qualStats.homopolymerSizesByPosition[referencePositionZeroBased] >= 3):
            #print ('Homopolymer. (1)-based position:' + str(referencePositionZeroBased + 1) + ', Base:' + str(referenceBase) + ', Length:' + str(qualStats.homopolymerSizesByPosition[referencePositionZeroBased]))

            # Fetch the reads from this specific homopolymer region.
            #iter = bamfile.fetch(alignmentRef.id, referencePositionZeroBased , referencePositionZeroBased + qualStats.homopolymerSizesByPosition[referencePositionZeroBased])
            
            iter = bamfile.fetch(alignmentRef.id, referencePositionZeroBased , referencePositionZeroBased + 10)
            alignedReadList = list(iter)
            
            for alignedRead in alignedReadList:
                print ('Homopolymer. (1)-based position:' + str(referencePositionZeroBased + 1) + ', Base:' + str(referenceBase) + ', Length:' + str(qualStats.homopolymerSizesByPosition[referencePositionZeroBased]))
                print ('Aligned Read:' + str(alignedRead.query_name))
                
                #alignedPairs = alignedRead.get_aligned_pairs()
                #print ('AlignedPairs:' + str(alignedPairs))
                
                print ('QueryAlignment:' + str(alignedRead.query_alignment_sequence))
                print ('Query Begin/End' + str(alignedRead.query_alignment_start) + ',' + str(alignedRead.query_alignment_end))
                print ('Reference Length vs Query Length:' + str(alignedRead.reference_length) + ' vs.' + str(alignedRead.query_length))
            #for x in iter:
                
                
                #print (str(x))
            #    pass
"""

            # How long is the homopolymer reported by the read?
            
            # Check : Every position of the read matches every position of the reference homopolymer.
            
            # Check: What if there are mismatches? We can't really use this bit of data for homopolymer analysis i think.
            
            
            
              
            
    # Calculate max coverage. Probably nothing is 100% coverage, so i should scale my graphs wisely.
    # Make graphs
    # Coverage over length of reference.
    # %match over length of reference. Percent and Q Score.
    # Expected homopolymer count vs observed homopolymer count.
    
    # Can I output the observed quality of each read?
    
    
    



# Perform minimap2 Alignment.  Align all reads against the Reference.
# TODO: Wasn't this supposed to be in a common methods somewhere?  Maybe not. Check if this method is duplicated.
# 
def alignReads(referenceLocation, readFileLocation, alignmentOutputDirectory, numberThreads):
    print('\nAnalyzing read quality.')
    print('\nStep 1.) Aligning reads against the reference.')
    
    if not isdir(alignmentOutputDirectory):
        makedirs(alignmentOutputDirectory)
    
    # Part 1 Index the Reference  
    # Actually i don't need to do this,    
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
        
    except Exception:
        print ('Exception indexing alignment reference. Is bwa installed? folder writing permission issue?')                  
        raise
    
    
    # Part 2 Align
    try:
        alignmentOutputName = join(alignmentOutputDirectory,'alignment.bam')

        # Parameters depend on "read type"
        # if i have a fastq, i assume these are minion reads.
        # if i have a fasta, i assume these are consensus sequences.
        if(getReadFileType(readFileLocation) == 'fasta'):
            minimapParams = '-ax asm5'
        elif(getReadFileType(readFileLocation) == 'fastq'):
            print ('attempting the consensus parameters instead of ONT settings.')
            minimapParams = '-ax map-ont'
        else:
            raise Exception('Unknown read file type....')

        cmd = ("minimap2 " + minimapParams + " " + 
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


