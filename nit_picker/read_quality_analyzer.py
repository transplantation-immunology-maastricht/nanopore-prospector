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
        
        self.alignedReadCountsByPosition   = [0]*self.referenceLength
        # For this application, I assume that a read will always be a match, mismatch,insertion, or deletion. Never both. An indel is not a match.
        self.matchReadCountsByPosition     = [0]*self.referenceLength
        self.insertionReadCountsByPosition = [0]*self.referenceLength
        self.deletionReadCountsByPosition  = [0]*self.referenceLength
        self.mismatchReadCountsByPosition  = [0]*self.referenceLength
        self.homopolymerSizesByPosition    = [0]*self.referenceLength
        
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


    def writeAlignmentSummaryHeader(self, alignmentSummaryOutputFile):
        return alignmentSummaryOutputFile.write(
            'Reference_Position_1_based' 
            + ',Reference_Base' 
            + ',Aligned_Reads' 
            + ',Is_Polymorphic' 
            + ',Match_Count' 
            + ',Match_Percent' 
            + ',Insertion_Count' 
            + ',Insertion_Percent' 
            + ',Deletion_Count' 
            + ',Deletion_Percent' 
            + ',Mismatch_Count' 
            + ',Mismatch_Percent' 
            + ',Homopolymer_Size\n')

    def writeAlignmentSummary(self, alignmentSummaryOutputFilename):
        alignmentSummaryOutputFile = createOutputFile(alignmentSummaryOutputFilename)
        polymorphicPositionsOutputFile = createOutputFile(alignmentSummaryOutputFilename.replace('.csv','.PolymorphicPositions.csv'))
      
        #createOutputFile(outputfileName)
        self.writeAlignmentSummaryHeader(alignmentSummaryOutputFile)
        self.writeAlignmentSummaryHeader(polymorphicPositionsOutputFile)

        # TODO: I want to print an alignment summary for ONLY the positions that look polymorphic. 
        # I might already have code for this, I am already lookin for the heterozygous positions.
        # Maybe it can be reused.


        # Declaring a constant here so it's easier to change. .90 is arbitrary
        # TODO: Should this be a commandline param? Yes.
        # .9999999 is useful for identifying any position that might be polymorphic sometimes.
        # .90 is probably more useful for minion reads
        polymorphicBaseCutoff = .99999999
        #polymorphicBaseCutoff = .90
        
        for referencePositionZeroBased, referenceBase in enumerate(self.referenceSequence):
            
            #print('analyzing base # ' + str(referencePositionZeroBased))
            #print ('the coverage should be:' + str(self.alignedReadCountsByPosition[referencePositionZeroBased]))
            #print ('the reference base should be:' + str(referenceBase))
            alignedReads = self.alignedReadCountsByPosition[referencePositionZeroBased]
            
            if(alignedReads == 0):
                # this is a hack to avoid divide by zero errors. Dividing zero by 1 should give a zero.
                matchPercent = 0
                #alignedReads = 1
                
            else:
                matchPercent = self.matchReadCountsByPosition[referencePositionZeroBased] / alignedReads
            
            self.writeAlignmentSummaryLine(alignmentSummaryOutputFile, referencePositionZeroBased, referenceBase, (matchPercent < polymorphicBaseCutoff))
            
            # If this position is polymorphic, also print this to the polymorphic positions output file.
            if( self.alignedReadCountsByPosition[referencePositionZeroBased] > 0 and matchPercent < polymorphicBaseCutoff):
                self.writeAlignmentSummaryLine(polymorphicPositionsOutputFile, referencePositionZeroBased, referenceBase, (matchPercent < polymorphicBaseCutoff))
        
        alignmentSummaryOutputFile.close()
        polymorphicPositionsOutputFile.close()

    

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
            # What is the difference? pysam reports an aligned read count. I'm also counting them myself.
            ',' + str(self.alignedReadCountsByPosition[referencePositionZeroBased]) + 
            ',' + ('YES' if isPolymorphic else '') +
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
    # TODO: Check where this pileup iterator starts. Are there reads mapped before or after the reference begins/ends?
    # I think the iterator will only look at where reads are mapped on the reference. No mapped reads = not in the pileup columns.
    for pileupColumn in pileupIterator:

        referencePositionZeroBased = pileupColumn.reference_pos  
        referenceBase = str(alignmentRef[referencePositionZeroBased])
        alignedCount = pileupColumn.nsegments      
        
        # Progress bar...
        if (referencePositionZeroBased%250 == 0):
            print ('Checking aligned reads at genomic position: (' + str(referencePositionZeroBased) + '/' + str(qualStats.referenceLength) + ')')
            print (str(alignedCount) + ' reads found at that position.')
        
        #print ('at 0-based pos:' + str(referencePositionZeroBased) + ' i see a ' + str(alignmentRef[referencePositionZeroBased]) + ' with  this many reads:' + str(alignedCount))
        
        qualStats.alignedReadCountsByPosition[referencePositionZeroBased] = alignedCount

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
        
        # Iterate the Reads at this position      
        # We use the pileup method to iterate through column-by-column.     
        for pileupRead in pileupColumn.pileups:
            
            qualStats.alignedReadCountsByPosition[referencePositionZeroBased] = alignedCount

            # If this read is a deletion
            if(pileupRead.is_del == 1):
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
            # Else if it is a refskip (TODO: What does this mean? no read aligned? Count these?) Doesn't happen often.
            elif(pileupRead.is_refskip):
                print('This read is a refskip, i dont know what that means:' + pileupRead.alignment.query_name)
                raise Exception('This read is a refskip, i dont know what that means:' + pileupRead.alignment.query_name)
            
            # else this means we have a base aligned at this position for this read.
            # It is either a match or a mismatch.
            else:    
                
                # TODO: Sometimes, there is no query position. What if this is 'None'?
                # This is because in the case of indels, there is no query position.
                queryPosition = int(pileupRead.query_position)
                
                currentBase = pileupRead.alignment.query_sequence[queryPosition].upper()                    

                if(currentBase.upper() == referenceBase.upper()):
                    qualStats.matchReadCountsByPosition[referencePositionZeroBased] += 1
                else:
                    qualStats.mismatchReadCountsByPosition[referencePositionZeroBased] += 1
                    qualStats.storeAlleleSpecificPolymorphism(pileupRead.alignment.query_name, referencePositionZeroBased, currentBase.upper() )
                    
            

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
    
    
    



# Perform BW Alignment.  Align all reads against the Reference.
# TODO: Wasn't this supposed to be in a common methods somewhere?  Maybe not. Check if this method is duplicated.
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


