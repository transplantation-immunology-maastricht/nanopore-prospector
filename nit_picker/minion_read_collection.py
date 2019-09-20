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

# TODO: I wonder if i can merge the ShannonEntropy code into this file.

from nanopore_prospector.common import createOutputFile, createScatterPlot, alignReads, getReadFileType
from Bio.SeqIO import parse, write
from numpy import mean, random, amin, amax, median
from os.path import join

from pysam import AlignmentFile

from math import log10, pow

class MinionReadCollection:
    # A class defining a collection of minION reads

    def __init__(self, readArray):
        self.readCollection = readArray
        self.readInputFormat = None

        # These arrays represent all the reads in a batch.
        self.readLengths = []
        self.readAvgReportedPhredQualities = []

        self.readFileLocation = None


        self.gene = None

        if self.readCollection is None:
            self.readCollection = []

        # TODO: I commented this out, i think this will mean sometimes stats are not calculated.
        # TODO: Validate that I am still making read stats for fastq files
        # self.summarizeSimpleReadStats()

    # This method evolved from the old constructor for the QualityStatistics class. It used to accept a reference sequence (String, i think)
    # now I made a method to handle that behavior.
    # If we have a reference sequence than we can calculate all sorts of quality statistics.
    def setReferenceSequence(self, referenceSequence):
        self.referenceSequence = referenceSequence

        if(self.referenceSequence is not None):
            self.referenceLength = len(referenceSequence)

            # Storing the new consensus as a list of characters. So i can modify it.
            # self.newConsensus = list(referenceSequence)
            self.newConsensus = referenceSequence

            self.alignedReadCountsByPosition = [0] * self.referenceLength
            # For this application, I assume that a read will always be a match, mismatch,insertion, or deletion. Never both. An indel is not a match.
            self.matchReadCountsByPosition = [0] * self.referenceLength
            self.insertionReadCountsByPosition = [0] * self.referenceLength
            self.deletionReadCountsByPosition = [0] * self.referenceLength
            self.mismatchReadCountsByPosition = [0] * self.referenceLength
            self.homopolymerSizesByPosition = [0] * self.referenceLength
            self.segmentCountByPosition = [0] * self.referenceLength
            self.supplementaryReadCountsByPosition = [0] * self.referenceLength
            self.snpCallsByPosition = [None] * self.referenceLength

            self.ACountsByPosition = [0] * self.referenceLength
            self.GCountsByPosition = [0] * self.referenceLength
            self.CCountsByPosition = [0] * self.referenceLength
            self.TCountsByPosition = [0] * self.referenceLength

            # I want to keep track of what alleles / reads have what polymorphisms.
            # This is a dictionary.
            # Key = allele id.
            # Value = a dictionary.
                # key = position, from the alignment reference
                # value = alternate base(s) from reference.
                # only have a key/value pair if there is a deviation from reference.
            self.alleleSpecificPolymorphisms = {}

            # Keep track of match percent for each allele/read. This is storage for OBSERVED read qualities, not reported in the reads.
            # Should only have entries for the reads that mapped to the reference.
            # This is a dictionary.
            # Key = allele/read id.
            # matchPercent.
            self.readMatchPercentages = {}

            self.totalAlignedReadCount = 0
            self.alignedReadLengths = []
            self.alignedReadPhredScores = []
            self.alignedReadQualityScores = []

            # Store a big pile of homopolymer data.
            # (Reference Pos, Reference Base, Reference Homopolymer Length, Read Homopolymer Length)
            self.homopolymerTuples = []
        else:
            print('No reference sequence provided.')


    def summarizeSimpleReadStats(self):
        # print('calculating read stats')
        self.readLengths = []
        self.readAvgReportedPhredQualities = []

        for currentRead in self.readCollection:
            currentSeqLength = len(currentRead)
            self.readLengths.append(currentSeqLength)

            if (self.readInputFormat == 'fastq'):
                phredQualities = currentRead.letter_annotations["phred_quality"]
                currentAvgPhredQuality = mean(phredQualities)
                self.readAvgReportedPhredQualities.append(currentAvgPhredQuality)

    def concatenate(self, otherCollection):
        self.readCollection = self.readCollection + otherCollection.readCollection
        self.readLengths = self.readLengths + otherCollection.readLengths
        self.readAvgReportedPhredQualities = self.readAvgReportedPhredQualities + otherCollection.readAvgReportedPhredQualities
        # self.summarizeSimpleReadStats()

    # TODO: Can i remove referenceSequenceLocation and use the "self" value? Probably.
    def outputReadPlotsAndSimpleStats(self, outputDirectory, plotName, sampleID, referenceSequenceLocation, calculateReadStatistics, minAlignmentLength, minimumSnpPctCutoff):
        self.summarizeSimpleReadStats()

        # Remove spaces for file names.
        simplePlotName = plotName.replace(' ', '_')

        self.readFileLocation = join(outputDirectory, str(sampleID) + '_' + simplePlotName + '.' + self.readInputFormat)
        allReadOutputFile = createOutputFile(self.readFileLocation)
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
            self.readAvgReportedPhredQualities = randomQualities

        # Scatterplot
        createScatterPlot(plotName
                          , self.readLengths
                          , self.readAvgReportedPhredQualities
                          , "Read Lengths"
                          , "Avg Read Quality(Phred)"
                          , join(outputDirectory, str(sampleID) + '_' + simplePlotName + '.png'))

        # TODO: I really only need to calculate these statistics against the pass reads.
        # Can i filter that somehow?

        # Calculate Statistics against Reference. This may crash if there is no reference sequence.
        if (calculateReadStatistics):


            # TODO: Calculate number of threads somewhere. I should pass it into nit_picker.
            numberThreads = 4
            self.calculateReadQualities(referenceSequenceLocation, outputDirectory, simplePlotName,
                                   numberThreads, minAlignmentLength, minimumSnpPctCutoff)

        # Return some simple stats. These are useful downstream.
        # readstats is a 2d array, with lengths and qualities.
        return self.readLengths, self.readAvgReportedPhredQualities

    def printSimpleReadStats(self, outputFileLocation):
        print ('Printing simple statistics to this location:\n' + str(outputFileLocation))
        simpleStatOutputFile = createOutputFile(outputFileLocation)

        simpleStatOutputFile.write('Read_ID,Length,Avg_Phred_Quality\n')
        for index, currentRead in enumerate(self.readCollection):

            if (self.readInputFormat == 'fastq'):
                readQual = self.readAvgReportedPhredQualities[index]
            else:
                readQual = 0

            readLength = self.readLengths[index]

            simpleStatOutputFile.write(str(currentRead.id) + ','
                                       + str(readLength) + ','
                                       + str(readQual) + '\n')

        simpleStatOutputFile.close()

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
                if (positionZeroBased in self.alleleSpecificPolymorphisms[alleleName].keys()):
                    polymorphicAllelesOutputFile.write(
                        ',' + str(self.alleleSpecificPolymorphisms[alleleName][positionZeroBased]))
                else:
                    polymorphicAllelesOutputFile.write(',')

            polymorphicAllelesOutputFile.write('\n')

        polymorphicAllelesOutputFile.close()

    def writeLDMatrix(self, outputDirectory):
        print ('Writing a Linkage Disequilibrium Matrix')

        # Quick way to identify the polymorphic positions
        # TODO I already did this in writePolymorphicAlleles, maybe I can store it and reuse.
        distinctAlignedSeqNames = sorted(self.alleleSpecificPolymorphisms.keys())
        distinctPolymorphicPositions = []
        for alleleName in distinctAlignedSeqNames:
            for positionZeroBased in self.alleleSpecificPolymorphisms[alleleName].keys():
                if positionZeroBased not in distinctPolymorphicPositions:
                    distinctPolymorphicPositions.append(positionZeroBased)
        distinctPolymorphicPositions = sorted(distinctPolymorphicPositions)

        # Create LD Matrix, and write it out.
        asymmetricLdMatrix, symmetricLdMatrix = self.calculateLinkDiseqMatrices(distinctPolymorphicPositions, distinctAlignedSeqNames)
        asymmetricLdMatrixOutputFile = createOutputFile(join(outputDirectory, 'AsymmetricLDMatrix.csv'))
        symmetricLdMatrixOutputFile = createOutputFile(join(outputDirectory, 'SymmetricLDMatrix.csv'))
        # Use 1-based positions.
        for positionZeroBased in distinctPolymorphicPositions:
            asymmetricLdMatrixOutputFile.write(',' + str(positionZeroBased + 1))
            symmetricLdMatrixOutputFile.write(',' + str(positionZeroBased + 1))
        asymmetricLdMatrixOutputFile.write('\n')
        symmetricLdMatrixOutputFile.write('\n')


        for snpYIndex, snpYPosition in enumerate(distinctPolymorphicPositions):
            asymmetricLdMatrixOutputFile.write(str(snpYPosition + 1))
            symmetricLdMatrixOutputFile.write(str(snpYPosition + 1))

            for snpXIndex, snpXPosition in enumerate(distinctPolymorphicPositions):
                asymmetricLdMatrixOutputFile.write(',' + str(round(asymmetricLdMatrix[snpXIndex][snpYIndex], 2)))
                symmetricLdMatrixOutputFile.write(',' + str(round(symmetricLdMatrix[snpXIndex][snpYIndex], 2)))

            asymmetricLdMatrixOutputFile.write('\n')
            symmetricLdMatrixOutputFile.write('\n')

        asymmetricLdMatrixOutputFile.close()
        symmetricLdMatrixOutputFile.close()


    def calculateLinkDiseqMatrices(self, distinctPolymorphicPositions, distinctAlignedSeqNames):
        # I'm referring to this paper:
        # https://www.genetics.org/content/198/1/321
        # I'm switching to a asymmetric LD constant.

        asymmetricLDMatrix = [[0 for x in range(len(distinctPolymorphicPositions))] for y in
                            range(len(distinctPolymorphicPositions))]

        symmetricWnMatrix = [[0 for x in range(len(distinctPolymorphicPositions))] for y in
                              range(len(distinctPolymorphicPositions))]

        # Step 1 is to calculate a list of possible snps/alleles for each position.
        # This structure is so I can loop over possible values of alleles at a position, and so I can compare phased SNPs.
        uniqueAllelesPerPosition = {}
        allelesPerPositionPerAlignedSeq = {}
        for snpIndex, snpPosition in enumerate(distinctPolymorphicPositions):
            if (snpIndex % 100 == 0):
                print('Step 1. calculating Unique SNPs at position # (' + str(snpIndex) + '/' + str(len(distinctPolymorphicPositions)) + ')')

            uniqueAllelesPerPosition[snpPosition] = []
            allelesPerPositionPerAlignedSeq[snpPosition] = {}

            for alignedSeqName in distinctAlignedSeqNames:

                if (snpPosition in self.alleleSpecificPolymorphisms[alignedSeqName].keys()):
                    currentSNP = self.alleleSpecificPolymorphisms[alignedSeqName][snpPosition]
                else:
                    # If the SNP doesn't exist, I assume it's the reference value.
                    # Todo: This might not always be true. It might be that there isn't a sequence aligned at this position and we don't see the SNP.
                    currentSNP = self.referenceSequence[snpPosition]

                allelesPerPositionPerAlignedSeq[snpPosition][alignedSeqName] = currentSNP
                if (currentSNP not in uniqueAllelesPerPosition[snpPosition]):
                    uniqueAllelesPerPosition[snpPosition].append(currentSNP)

            #print('At 1-based position ' + str(snpPosition + 1) + ' I can see ' + str(len(uniqueSNPsPerPosition[snpPosition])) + ' SNPS:' + str(uniqueSNPsPerPosition[snpPosition]))
            #print(str)


        # Step 2... Calculate Allele and Haplotype Frequencies.
        frequenciesPerPositionPerAllele = {}
        haplotypeFrequenciesPerPositionsPerAlleles = {}
        for snpAIndex, snpAPosition in enumerate(distinctPolymorphicPositions):

            # Progress bar...
            if(snpAIndex % 100 == 0):
                print('Step 2. calculating SNP & Haplotype Frequencies for SNP# (' + str(snpAIndex) + '/' + str(len(distinctPolymorphicPositions)) + ')')

            # Allele Frequencies...
            frequenciesPerPositionPerAllele[snpAPosition] = {}
            for currentSNPValueI in uniqueAllelesPerPosition[snpAPosition]:
                alleleFrequencyAI = 0
                for alignedSeqName in distinctAlignedSeqNames:
                    if allelesPerPositionPerAlignedSeq[snpAPosition][alignedSeqName] == currentSNPValueI:
                        alleleFrequencyAI += 1
                frequenciesPerPositionPerAllele[snpAPosition][currentSNPValueI] = alleleFrequencyAI / len(distinctAlignedSeqNames)

            # Haplotype Frequencies...
            haplotypeFrequenciesPerPositionsPerAlleles[snpAPosition] = {}
            for snpBIndex, snpBPosition in enumerate(distinctPolymorphicPositions):
                haplotypeFrequenciesPerPositionsPerAlleles[snpAPosition][snpBPosition] = {}
                for currentSNPValueI in uniqueAllelesPerPosition[snpAPosition]:
                    haplotypeFrequenciesPerPositionsPerAlleles[snpAPosition][snpBPosition][currentSNPValueI] = {}
                    for currentSNPValueJ in uniqueAllelesPerPosition[snpBPosition]:
                        haplotypeFrequency = 0
                        for alignedSeqName in distinctAlignedSeqNames:
                            if (allelesPerPositionPerAlignedSeq[snpAPosition][alignedSeqName] == currentSNPValueI
                                and allelesPerPositionPerAlignedSeq[snpBPosition][alignedSeqName] == currentSNPValueJ):
                                haplotypeFrequency += 1
                        haplotypeFrequenciesPerPositionsPerAlleles[snpAPosition][snpBPosition][currentSNPValueI][currentSNPValueJ] = haplotypeFrequency / len(distinctAlignedSeqNames)



        # Step 3... Loop through snps in 2d array, calculate LD positions, write position in first column.
        for snpAIndex, snpAPosition in enumerate(distinctPolymorphicPositions):

            # Progress bar...
            if(snpAIndex % 100 == 0):
                print('Step 3. calculating asymmetric linkage constants for SNP# (' + str(snpAIndex) + '/' + str(len(distinctPolymorphicPositions)) + ')')

            # Loop through snps, write LD calculation in corresponding column.
            for snpBIndex, snpBPosition in enumerate(distinctPolymorphicPositions):
                if (snpBPosition == snpAPosition):
                    asymmetricLDMatrix[snpAIndex][snpBIndex] = 0

                else:
                    # Asymmetric Linkage Disequlibrium
                    # W^2 (a/b) = [ Sum(i) Sum(j) ( D^2(ij) / p(Bj) ) ] / (1-FA)

                    # Calculate FA, the homozygosity constant?
                    # 1 - [Sum(i) (PAi^2)]

                    linkageSum = 0
                    homozygositySum = 0
                    #loop each possible SNP value
                    for currentSNPValueI in uniqueAllelesPerPosition[snpAPosition]:

                        homozygositySum += (frequenciesPerPositionPerAllele[snpAPosition][currentSNPValueI]) ** 2

                        for currentSNPValueJ in uniqueAllelesPerPosition[snpBPosition]:
                            # calculate the D^2 value, add it to the sum
                            # ( D^2(ij) / p(Bj) )

                            # Dij = fij - pAi * pBj
                            DStatisticSquared = (haplotypeFrequenciesPerPositionsPerAlleles[snpAPosition][snpBPosition][currentSNPValueI][currentSNPValueJ]
                                - (frequenciesPerPositionPerAllele[snpAPosition][currentSNPValueI]
                                * frequenciesPerPositionPerAllele[snpBPosition][currentSNPValueJ])
                                )**2
                            # The squared D stat is divided by PBj and added to the sum
                            linkageSum += (DStatisticSquared / frequenciesPerPositionPerAllele[snpBPosition][currentSNPValueJ])


                    heterozygosityConstant = 1 - homozygositySum

                    # Divide Sum by (1/FA)
                    ldValue = linkageSum / (heterozygosityConstant)
                    asymmetricLDMatrix[snpAIndex][snpBIndex] = ldValue

        # Step 4... Loop through snps in 2d array, calculate symmetric Wn values.
        for snpAIndex, snpAPosition in enumerate(distinctPolymorphicPositions):

            # Progress bar...
            if (snpAIndex % 100 == 0):
                print('Step 4. calculating symmetric Wn linkage constants for SNP# (' + str(snpAIndex) + '/' + str(
                    len(distinctPolymorphicPositions)) + ')')

            # Loop through snps, write LD calculation in corresponding column.
            for snpBIndex, snpBPosition in enumerate(distinctPolymorphicPositions):
                if (snpBPosition == snpAPosition):
                    symmetricWnMatrix[snpAIndex][snpBIndex] = 0

                # Symmsetric, to reduce redundancy we don't need to calculate every pair in both directions.
                elif( snpBPosition > snpAPosition):
                    #Wn2 is the multi-alleleic extension of the r2 measure.
                    #Wn2 = [ Sum(i) Sum(j) ( D^2(ij) / (pAi * pBj))  ] / min((kA - 1, kB - 1))

                    linkageSum = 0
                    # loop each possible SNP value
                    for currentSNPValueI in uniqueAllelesPerPosition[snpAPosition]:


                        for currentSNPValueJ in uniqueAllelesPerPosition[snpBPosition]:
                            # calculate the D^2 value, add it to the sum
                            #(D ^ 2(ij) / (pAi * pBj))

                            # Dij = fij - pAi * pBj
                            DStatisticSquared = (haplotypeFrequenciesPerPositionsPerAlleles[snpAPosition][snpBPosition][
                                     currentSNPValueI][currentSNPValueJ]
                                 - (frequenciesPerPositionPerAllele[snpAPosition][currentSNPValueI]
                                    * frequenciesPerPositionPerAllele[snpBPosition][currentSNPValueJ])
                                 ) ** 2
                            # The squared D stat is divided by PBj and added to the sum
                            linkageSum += (DStatisticSquared /
                                (frequenciesPerPositionPerAllele[snpBPosition][currentSNPValueJ]
                                *frequenciesPerPositionPerAllele[snpAPosition][currentSNPValueI])
                                )

                    # The entire sum is divided by min((kA - 1, kB - 1))
                    # Those are degrees of freedom. The kA and kB are the numbers of alleles at each locus.
                    degreesOfFreedom = min(
                        len(uniqueAllelesPerPosition[snpAPosition]) - 1
                        , len(uniqueAllelesPerPosition[snpBPosition]) - 1
                    )

                    ldValue = linkageSum / degreesOfFreedom

                    symmetricWnMatrix[snpAIndex][snpBIndex] = ldValue
                    symmetricWnMatrix[snpBIndex][snpAIndex] = ldValue

                else:
                    pass
                    #print('snpBPosition(' + str(snpBPosition) + ') <= snpAPosition(' + str(snpAPosition) + '). This is fine, im just not calculating it again.')


        return asymmetricLDMatrix, symmetricWnMatrix


    def calculateDPrimeLDMatrix(self, distinctPolymorphicPositions):
        # TODO: I know this method is faulty. It treats each SNP as bi-allelic. Either it's a perfect match, or not.
        # TODO: I'm deprecating this method and switching to the asymmetric LD measure matrix.
        # I'm switching to a asymmetric LD constant, look in the other method.
        # ldMatrixUnscaled=[len(distinctPolymorphicPositions) , len(distinctPolymorphicPositions)]
        # The level of linkage disequilibrium between A and B can be quantified by the coefficient of linkage disequilibrium {\displaystyle D_{AB}} D_{AB}, which is defined as
        # D_{AB}=p_{AB}-p_{A}p_{B},
        # https://en.wikipedia.org/wiki/Linkage_disequilibrium
        # In the case {D_{{AB}}=0, said to be in linkage equilibrium
        ldMatrixUnscaled = [[0 for x in range(len(distinctPolymorphicPositions))] for y in
                            range(len(distinctPolymorphicPositions))]
        # Maximum value for scaling?
        maxLDValue = 0
        # Loop through snps in 2d array, calculate LD positions, write position in first column.
        for snpYIndex, snpYPosition in enumerate(distinctPolymorphicPositions):

            # Loop through snps, write LD calculation in corresponding column.
            for snpXIndex, snpXPosition in enumerate(distinctPolymorphicPositions):
                if (snpYPosition == snpXPosition):
                    ldMatrixUnscaled[snpXIndex][snpYIndex] = 0
                    # ldMatrixOutputFile.write(',' + str(0))

                # Checking this allows me to calculate the LD in both directions at the same time.
                # To reduce redundancy.
                elif (snpYPosition > snpXPosition):
                    # Calculate LD.
                    # Coefficient of Linkage Disequilibrium D{AB}
                    # D{AB} = p{AB} - p{A} * p{B}
                    # Define the allele as sequences that matches the reference.
                    # Any insert/delete/mismatch is treated as the "other" allele.

                    # TODO: This is a multi-allelic loci. Should I instead be calculating Wn, instead of D?
                    # This way, it's either match, or everything else.
                    # A match is an allele, a deletion is an allele.

                    # p{A}
                    pA = 1.0 * self.matchReadCountsByPosition[snpYPosition] / self.alignedReadCountsByPosition[
                        snpYPosition]

                    # p{B}
                    pB = 1.0 * self.matchReadCountsByPosition[snpXPosition] / self.alignedReadCountsByPosition[
                        snpXPosition]

                    # p{AB}
                    # Must look at specific alleles to calculate this.
                    # Loop thru alleles and count if BOTH of them match.
                    matchBothCount = 0
                    for alleleName in self.alleleSpecificPolymorphisms.keys():
                        # alleleSpecificPolymorphisms only tracks snps.
                        # If there is no dictionary entry for either,
                        # that means the allele matches the reference at this position.
                        if (snpYPosition not in self.alleleSpecificPolymorphisms[alleleName] and snpXPosition not in
                                self.alleleSpecificPolymorphisms[alleleName]):
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

                    # Maximum value for scaling?
                    if (ldValue > maxLDValue):
                        maxLDValue = ldValue

                else:
                    # Don't need to calculate the LD a second time.
                    pass

            # ldMatrixOutputFile.write('\n')
        # Scale the LD values.
        ldMatrixScaled = [[0 for x in range(len(distinctPolymorphicPositions))] for y in
                          range(len(distinctPolymorphicPositions))]
        # ldMatrixScaled=[len(distinctPolymorphicPositions) , len(distinctPolymorphicPositions)]
        for snpYIndex, snpYPosition in enumerate(distinctPolymorphicPositions):
            for snpXIndex, snpXPosition in enumerate(distinctPolymorphicPositions):
                ldMatrixScaled[snpXIndex][snpYIndex] = ldMatrixUnscaled[snpXIndex][snpYIndex] / maxLDValue
        return ldMatrixScaled

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

    def createAlignmentQualityGraph(self, alignmentQualityFileName):
        # This is a graph of %of reads that match at any position across the reference seqeunce.
        # Problems: The graph Y axis scales, so the distribution always looks the same
        # Problem: Low percentages at the homopolymers, Often times only 20% of reads are binding there.
        # Because the reads suggest shorter homopolymers than in the reference.

        xValues = [0] * self.referenceLength
        yValues = [0] * self.referenceLength

        for referencePositionZeroBased in range(0, self.referenceLength):
            # print('checking position' + str(referencePositionZeroBased))
            # Y = Percent Match
            if(self.segmentCountByPosition[referencePositionZeroBased] != 0):
                yValues[referencePositionZeroBased] = 100 * self.matchReadCountsByPosition[referencePositionZeroBased] / \
                    self.segmentCountByPosition[referencePositionZeroBased]
            else:
                yValues[referencePositionZeroBased] = 0
            # X = 0-based reference position.
            xValues[referencePositionZeroBased] = referencePositionZeroBased

        createScatterPlot('Alignment Quality', xValues, yValues, 'Reference Position', 'Percent Read Match',
                          alignmentQualityFileName)

    def createSlidingScaleAlignmentQualityGraph(self, alignmentQualityFileName):
        # Assuming a window size of 5.
        # I think odd numbers are easier to math.
        windowSize = 5
        halfWindowSize = int(windowSize / 2)  # The int() function should always round DOWN, which is what i want here.
        print ('this half window size should be 2:' + str(halfWindowSize))

        xValues = []
        yValues = []

        for referencePositionZeroBased in range(0, self.referenceLength):
            # print('checking position' + str(referencePositionZeroBased))
            if ((1 + referencePositionZeroBased + halfWindowSize) % windowSize == 0):
                # print('calculating at position' + str(referencePositionZeroBased))

                # Initialize with the match and alignedread counts at the currentPosition.
                currentMatchCount = self.matchReadCountsByPosition[referencePositionZeroBased]
                currentAlignedSectionCount = self.segmentCountByPosition[referencePositionZeroBased]
                # currentValueCount = 1

                # Collect all the values within the window size. We can iterate using the halfWindowSize.
                for windowIterator in range(1, halfWindowSize + 1):
                    # Value to the left...
                    currentMatchCount += self.matchReadCountsByPosition[referencePositionZeroBased - windowIterator]
                    currentAlignedSectionCount += self.segmentCountByPosition[
                        referencePositionZeroBased - windowIterator]
                    # currentValueCount += 1

                    # Value to the right. Must ensure we don't walk off the end of the arrays.
                    if (referencePositionZeroBased + windowIterator < len(self.matchReadCountsByPosition)):
                        currentMatchCount += self.matchReadCountsByPosition[referencePositionZeroBased + windowIterator]
                        currentAlignedSectionCount += self.segmentCountByPosition[
                            referencePositionZeroBased + windowIterator]
                        # currentValueCount += 1

                # Y = Percent Match
                if(currentAlignedSectionCount == 0):
                    yValues.append(0)
                else:
                    yValues.append(100 * currentMatchCount / currentAlignedSectionCount)
                # X = 0-based reference position.
                xValues.append(referencePositionZeroBased)

        createScatterPlot('Alignment Quality', xValues, yValues, 'Reference Position', 'Percent Read Match',
                          alignmentQualityFileName)

    def callSNP(self, referencePositionZeroBased, referenceBase, minimumSnpPctCutoff):
        # Assign an alternate base if there is a snp at this position.
        # Do nothing if there is no SNP here.
        # TODO: I made minimumSnpPctCutoff a parameter. What will break?
        # expectedMatchingBaseAccuracy is the inverse of the minimumSnpPctCutoff. It's a measure of "accuracy", not a
        # minimum polymorphism amount.
        if minimumSnpPctCutoff is None:
            minimumSnpPctCutoff = .00001
        expectedMatchingBaseAccuracy = (1-minimumSnpPctCutoff)

        alignedReads = self.alignedReadCountsByPosition[referencePositionZeroBased]

        if (alignedReads == 0):
            # this is a hack to avoid divide by zero errors. Dividing zero by 1 should give a zero.
            matchPercent = 0
            # alignedReads = 1
        else:
            matchPercent = self.matchReadCountsByPosition[referencePositionZeroBased] / alignedReads

        # Bad: Is there an insertion? Insertion count > match count * accuracy
        # Better Logic: Is there a deletion? Insertion count > aligned count * snpCutoff
        #if (self.insertionReadCountsByPosition[referencePositionZeroBased] > self.matchReadCountsByPosition[referencePositionZeroBased] * expectedMatchingBaseAccuracy):
        if (self.insertionReadCountsByPosition[referencePositionZeroBased] > self.alignedReadCountsByPosition[referencePositionZeroBased] * minimumSnpPctCutoff):
            # TODO: Calculate the inserted text based on insertion length. For now i will call this SNP an I.
            # TODO: I'm disabling this for now. It's not so useful to just place an I in the sequence.
            # Since I don't know what the insertion should be, I skip for now.
            # TODO: If I'm trying to change the consensus, I should figure out what text to change in the sequence. What if there are different inserted sequences?
            # TODO Change the consensus sequence for this insertion?
            self.snpCallsByPosition[referencePositionZeroBased] = 'I'
            #pass


        # TODO: I found a bug in how I detect deletions. It happens when I use .99999 as a expected accuracy.
        # TODO: When there are more matching sequences than deletions detected, I do not detect the SNP.
        # TODO: Solution: change the deletion detection logic.
        # TODO: Probably the same for insertions.
        # Flawed Logic: Is there a deletion? Deletion count > match count * accuracy
        # Better Logic: Is there a deletion? Deletion count > aligned count * snpCutoff
        #if (self.deletionReadCountsByPosition[referencePositionZeroBased] > self.matchReadCountsByPosition[referencePositionZeroBased] * expectedMatchingBaseAccuracy):
        if (self.deletionReadCountsByPosition[referencePositionZeroBased] > self.alignedReadCountsByPosition[referencePositionZeroBased] * minimumSnpPctCutoff):
            self.snpCallsByPosition[referencePositionZeroBased] = '-'

            # TODO Delete this base from the consensus sequence

        # If this position is polymorphic, also print this to the polymorphic positions output file.
        if (self.alignedReadCountsByPosition[referencePositionZeroBased] > 0 and matchPercent < expectedMatchingBaseAccuracy):
            # Not an indel, this is probably a mismatch.

            maxBaseCount = max(
                self.ACountsByPosition[referencePositionZeroBased]
                , self.GCountsByPosition[referencePositionZeroBased]
                , self.CCountsByPosition[referencePositionZeroBased]
                , self.TCountsByPosition[referencePositionZeroBased]
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

            # TODO: Find a way to selectively ignore regions. My snp caller doesn't work in deletions.
            # if 4550 <= referencePositionZeroBased <= 4695:
            #    # We are in the complicated delete region, lets not assign this snp.
            #    pass

            # elif 4909 <= referencePositionZeroBased <= 4918:
            #    # This is the shorter deletion.
            #    pass
            # else:

            # manipulating strings because they are immutable
            # tempid = self.newConsensus.id
            # TODO This is a very rough way to get a new consensus sequence with SNPs. Fix this logic. The indexes are probably being messed up.
            tempConsensus = list(self.newConsensus)
            tempConsensus[referencePositionZeroBased] = self.snpCallsByPosition[referencePositionZeroBased]
            self.newConsensus = "".join(tempConsensus)
            # self.newConsensus.id = 'GeneratedConsensus'
            # Put the SNP in our new consensus sequence
            # self.newConsensus[referencePositionZeroBased] = self.snpCallsByPosition[referencePositionZeroBased]

            # self.snpCallsByPosition[referencePositionZeroBased] = '!'

        pass

    def writeSNPs(self, snpSummaryOutputFilename):
        snpSummaryOutputFile = createOutputFile(snpSummaryOutputFilename)

        row1Text = 'Pos,'
        row2Text = 'Ref,'
        row3Text = 'SNP,'

        for referencePositionZeroBased, referenceBase in enumerate(self.referenceSequence):
            if (self.snpCallsByPosition[referencePositionZeroBased] is not None):
                row1Text += str(referencePositionZeroBased + 1) + ','
                row2Text += referenceBase + ','
                row3Text += self.snpCallsByPosition[referencePositionZeroBased] + ','
            else:
                pass

        snpSummaryOutputFile.write(row1Text + '\n')
        snpSummaryOutputFile.write(row2Text + '\n')
        snpSummaryOutputFile.write(row3Text + '\n')
        snpSummaryOutputFile.close()

    def writeAlignmentSummary(self, alignmentSummaryOutputFilename, minimumSnpPctCutoff):
        alignmentSummaryOutputFile = createOutputFile(alignmentSummaryOutputFilename)
        polymorphicPositionsOutputFile = createOutputFile(
            alignmentSummaryOutputFilename.replace('.csv', '.PolymorphicPositions.csv'))

        newConsensusOutputFile = createOutputFile(alignmentSummaryOutputFilename.replace('.csv', '.NewConsensus.fasta'))

        # createOutputFile(outputfileName)
        self.writeAlignmentSummaryHeader(alignmentSummaryOutputFile)
        self.writeAlignmentSummaryHeader(polymorphicPositionsOutputFile)

        for referencePositionZeroBased, referenceBase in enumerate(self.referenceSequence):

            self.callSNP(referencePositionZeroBased, referenceBase, minimumSnpPctCutoff)

            # print('analyzing base # ' + str(referencePositionZeroBased))
            # print ('the coverage should be:' + str(self.alignedReadCountsByPosition[referencePositionZeroBased]))
            # print ('the reference base should be:' + str(referenceBase))
            # alignedReads = self.alignedReadCountsByPosition[referencePositionZeroBased]

            # if(alignedReads == 0):
            # this is a hack to avoid divide by zero errors. Dividing zero by 1 should give a zero.
            #    matchPercent = 0
            # alignedReads = 1

            # else:
            #    matchPercent = self.matchReadCountsByPosition[referencePositionZeroBased] / alignedReads

            self.writeAlignmentSummaryLine(alignmentSummaryOutputFile, referencePositionZeroBased, referenceBase,
                                           (self.snpCallsByPosition[referencePositionZeroBased] is not None))

            # If this position is polymorphic, also print this to the polymorphic positions output file.
            # if( self.alignedReadCountsByPosition[referencePositionZeroBased] > 0 and matchPercent < polymorphicBaseCutoff):
            #    self.writeAlignmentSummaryLine(polymorphicPositionsOutputFile, referencePositionZeroBased, referenceBase, (matchPercent < polymorphicBaseCutoff))

            if (self.snpCallsByPosition[referencePositionZeroBased] is not None):
                self.writeAlignmentSummaryLine(polymorphicPositionsOutputFile, referencePositionZeroBased,
                                               referenceBase, True)

        alignmentSummaryOutputFile.close()
        polymorphicPositionsOutputFile.close()

        # Write new consensus
        # TODO: It seems like sometimes self.newConsensus is a String, sometimes its a Seq object. Not sure why.
        # write([self.newConsensus], newConsensusOutputFile, 'fasta')
        newConsensusOutputFile.write('>GeneratedConsensusSequence\n')
        newConsensusOutputFile.write(str(self.newConsensus) + '\n')
        newConsensusOutputFile.close()

        self.writeSNPs(alignmentSummaryOutputFilename.replace('.csv', '.SNPs.csv'))

    def writeAlignmentSummaryLine(self, outputFile, referencePositionZeroBased, referenceBase, isPolymorphic):
        # I suspect that "segmentCount" is more accurate than the number of aligned reads.
        # If i calculate match percentage based on alignedReads, sometimes it's >100%.
        # Because insertions don't count as alignedReads?
        alignedReads = self.alignedReadCountsByPosition[referencePositionZeroBased]
        segmentCount = self.segmentCountByPosition[referencePositionZeroBased]
        if (self.alignedReadCountsByPosition[referencePositionZeroBased] == 0):
            # this is a hack to avoid divide by zero errors. Dividing zero by 1 should give a zero.
            matchPercent = 0
            alignedReads = 1

        else:
            # matchPercent = self.matchReadCountsByPosition[referencePositionZeroBased] / alignedReads
            matchPercent = self.matchReadCountsByPosition[referencePositionZeroBased] / segmentCount

        #print('Reference0based: ' + str(referencePositionZeroBased))
        #print('Debug Segment Count: ' + str(segmentCount))

        if(segmentCount == 0):
            insertionPercent = 0
            deletionPercent = 0
            mismatchPercent = 0
        else:
            insertionPercent = self.insertionReadCountsByPosition[referencePositionZeroBased] / segmentCount
            deletionPercent = self.deletionReadCountsByPosition[referencePositionZeroBased] / segmentCount
            mismatchPercent = self.mismatchReadCountsByPosition[referencePositionZeroBased] / segmentCount

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
            ',' + str(insertionPercent) +
            ',' + str(self.deletionReadCountsByPosition[referencePositionZeroBased]) +
            ',' + str(deletionPercent) +
            ',' + str(self.mismatchReadCountsByPosition[referencePositionZeroBased]) +
            ',' + str(mismatchPercent) +
            ',' + str(self.homopolymerSizesByPosition[referencePositionZeroBased]) +
            '\n'
        )

    def storeAlleleSpecificPolymorphism(self, alleleID, baseIndexZeroBased, substituteBases):
        # print ('allele ' + alleleID + ' has a substitute base at position ' + str(baseIndexZeroBased + 1) + ':' + substituteBases)

        # If we haven't seen this allele before
        if alleleID not in self.alleleSpecificPolymorphisms.keys():
            # Set up a dictionary of it's polymorphisms.
            self.alleleSpecificPolymorphisms[alleleID] = {}

        # Store this polymorphic base in this allele's dictionary.
        self.alleleSpecificPolymorphisms[alleleID][baseIndexZeroBased] = substituteBases


    # TODO: this method is a big mess, try to break it up a bit. Can I move it within QualityStatistics? No, Into MinionReadCollection
    # TODO: The minAlignmentLength is never used! Am i accidentally using a hardcoded value in there?
    #def calculateReadQualities(self, referenceSequencelocation, outputDirectory, sampleID, numberThreads):
    def calculateReadQualities(self, referenceSequencelocation, outputDirectory, sampleID, numberThreads, minAlignmentLength, minimumSnpPctCutoff):

        print ('Calculating Read Quality Statistics.')
        print ('I have this reference File:\n' + str(referenceSequencelocation))

        alignmentOutputDir = join(outputDirectory, sampleID + '_ReadQualityAnalysis')
        # Copy reference to outputdir
        referenceSequenceTargetLocation = join(alignmentOutputDir, 'AlignmentReference.fasta')
        # Try to do this with all reference sequences, not just the first.
        alignmentRef = list(parse(referenceSequencelocation, 'fasta'))[0]
        #alignmentRef = list(parse(referenceSequencelocation, 'fasta'))

        # Here I'm getting rid of the QualityStatistcs objects and using a MinIONReadCOllection insted.
        # qualStats = QualityStatistics(alignmentRef.seq)
        # I think Im inside my MinionReadCollection object, i think it's not necessary to make a new one.
        # I still need to set the reference file, because that initializes all the lists I need.
        #currentReadCollection = createCollectionFromReadFile(self.readFileLocation)
        self.setReferenceSequence(alignmentRef.seq)

        write(alignmentRef, createOutputFile(referenceSequenceTargetLocation), 'fasta')

        # Align Reads against Reference
        print('Attempting the alignment, the read file location should be something:' + str(self.readFileLocation))
        # TODO: I changed this to use a read subset to make it faster, that might cause problems later.
        # TODO: I'm passing a "True" as minimum alignmentLength? Fix this, duh. Actually it's not being used. It's an extra argument right now.
        #alignReads(referenceSequenceTargetLocation, self.readFileLocation, alignmentOutputDir, True, numberThreads, True)
        alignReads(referenceSequenceTargetLocation, self.readFileLocation, alignmentOutputDir, False, numberThreads, minAlignmentLength)


        # Open Alignment
        alignmentRef = list(parse(referenceSequenceTargetLocation, 'fasta'))[0]
        bamfile = AlignmentFile(join(alignmentOutputDir, 'alignment.bam'), 'rb')

        # Declare variables to store the information.
        # Commenting this out because....we already have a totalReadCount, I can calculate from list of reads.
        # What will this break? Maybe it breaks the HLA Allele Analysis...
        # self.totalReadCount = len(list(parse(readFileLocation, getReadFileType(readFileLocation))))

        # Initialize the read dictionary, so we don't lose track of any of these read / alleles
        # TODO: not sure if i should also initialize readMatchPercentages. Because I only want those if the read is mapped.
        for readObject in list(parse(self.readFileLocation, getReadFileType(self.readFileLocation))):
            if readObject.id not in self.alleleSpecificPolymorphisms.keys():
                self.alleleSpecificPolymorphisms[readObject.id] = {}

        # Loop alignment:
        # Iterate the reference sequence base by base (column by column)
        # Calculate percentages for each base.
        pileupIterator = bamfile.pileup(alignmentRef.id)
        for pileupColumn in pileupIterator:

            referencePositionZeroBased = pileupColumn.reference_pos
            referenceBase = str(alignmentRef[referencePositionZeroBased])
            # think this is wrong. It seems like nsegments
            # aligned reads + deletions.
            self.segmentCountByPosition[referencePositionZeroBased] = pileupColumn.nsegments

            # Progress bar...
            if (referencePositionZeroBased % 250 == 0):
                print ('Checking aligned reads at genomic position: (' + str(referencePositionZeroBased) + '/' + str(
                    self.referenceLength) + ')')
                print (str(pileupColumn.nsegments) + ' segments found at that position.')

            # Search for homopolymers in the reference sequence.
            # If a homopolymer begins at this position, store the LENGTH of the homopolymer.
            # TODO: This is a slight bug. I'm only measuring homopolymers where there are reads aligned.
            # Big problem with low read counts. Or is it? I think small gaps are fine. Ok no problem.
            # Homopolymers outside aligned regions are ignored. That's probably fine.
            # The most 3' base in the reference. A homopolymer cannot begin at the rightmost base.
            if (referencePositionZeroBased == self.referenceLength - 1):
                pass
            else:

                # Check the most 5' base in the reference. Don't look to the left, we will go out of bounds.
                if (referencePositionZeroBased == 0):
                    pass  # Is this logic right? seems like we skip the leftmost base
                # If the current base matches the base to the left, a homopolymer cannot begin here.
                elif (referenceBase == str(alignmentRef[referencePositionZeroBased - 1])):
                    pass
                else:
                    # Check the bases to the right, to determine homopolymer length. (don't exceed length).
                    homopolymerLengthIterator = 1
                    while (referencePositionZeroBased + homopolymerLengthIterator < self.referenceLength):
                        if (referenceBase == str(alignmentRef[referencePositionZeroBased + homopolymerLengthIterator])):

                            # Record homopolymer length.
                            self.homopolymerSizesByPosition[
                                referencePositionZeroBased] = homopolymerLengthIterator + 1
                            homopolymerLengthIterator += 1

                        else:
                            break

            self.alignedReadCountsByPosition[referencePositionZeroBased] = 0
            # Iterate the Reads at this position
            # We use the pileup method to iterate through column-by-column.
            for pileupRead in pileupColumn.pileups:

                # Supplementary alignments are part of the minimap algorithm. I dont want to use them for quality analysis.
                if (pileupRead.alignment.is_supplementary):
                    self.supplementaryReadCountsByPosition[referencePositionZeroBased] += 1
                    pass

                # If this read is a deletion
                elif (pileupRead.is_del == 1):
                    # TODO: Found a bug involving several base deletions.
                    # A sequence with several deleted bases should only be counted as a single deletion.
                    # The deletion is located at the very first base. The read does NOT count as a deletion if the
                    # previous base is also a deletion.

                    # Solution: Get the read id, check previous base.
                    # If previous base is also a deletion, do nothing.
                    #print('Deletion, RefPos(0based)=' + str(referencePositionZeroBased) + ', Read=' + str(pileupRead.alignment.query_name))
                    alignedPairs = pileupRead.alignment.get_aligned_pairs()

                    # If we are in the "middle" of a deletion, then there is also a deletion at the base to the left.
                    alignedReadPositionsJustBeforeDeletion = [pair[0] for pair in alignedPairs if pair[1] == referencePositionZeroBased - 1][0]
                    if(alignedReadPositionsJustBeforeDeletion is None):
                        # Previous base is a deletion. We should NOT count this one.
                        #print ('NOT a real deletion\n')
                        pass
                    else:
                        #print( 'REAL deletion.\n')
                        self.deletionReadCountsByPosition[referencePositionZeroBased] += 1
                        self.storeAlleleSpecificPolymorphism(pileupRead.alignment.query_name, referencePositionZeroBased, '-')


                # else if this read is an insertion
                elif (pileupRead.indel > 0):
                    # This apparently indicates an insertion AFTER this base. The query text at this position is the reference
                    # , and there are inserted bases after this position
                    # Indexing the query sequence is tricky.
                    # I am also recording the correct base, because I hope that is more clear, and doesn't just look like a base substitution.
                    self.insertionReadCountsByPosition[referencePositionZeroBased] += 1

                    insertedBasesCount = pileupRead.indel
                    insertedBases = pileupRead.alignment.query_sequence[
                                    pileupRead.query_position:pileupRead.query_position + (insertedBasesCount + 1)].upper()
                    self.storeAlleleSpecificPolymorphism(pileupRead.alignment.query_name,
                                                                          referencePositionZeroBased,
                                                                          insertedBases)  # Get it?

                    queryPosition = int(pileupRead.query_position)
                    # Because....Insertions still have aligned bases....which can match the reference
                    # TODO: SHould I actually be checking for insertions in the bases BEFORE this position?
                    # This could be a major misunderstanding on my part....
                    currentBase = pileupRead.alignment.query_sequence[queryPosition].upper()

                    if (currentBase.upper() == referenceBase.upper()):
                        self.matchReadCountsByPosition[referencePositionZeroBased] += 1
                    else:
                        self.mismatchReadCountsByPosition[referencePositionZeroBased] += 1
                        self.storeAlleleSpecificPolymorphism(pileupRead.alignment.query_name,
                                                                              referencePositionZeroBased,
                                                                              currentBase.upper())

                    if (currentBase == 'A'):
                        self.ACountsByPosition[referencePositionZeroBased] += 1
                    elif (currentBase == 'G'):
                        self.GCountsByPosition[referencePositionZeroBased] += 1
                    elif (currentBase == 'C'):
                        self.CCountsByPosition[referencePositionZeroBased] += 1
                    elif (currentBase == 'T'):
                        self.TCountsByPosition[referencePositionZeroBased] += 1



                # Else if it is a refskip TODO: What does this mean? no read aligned? Count these? Doesn't happen often.
                elif (pileupRead.is_refskip):
                    print('This read is a refskip, i dont know what that means:' + pileupRead.alignment.query_name)
                    raise Exception(
                        'This read is a refskip, i dont know what that means:' + pileupRead.alignment.query_name)

                # else this means we have a base aligned at this position for this read.
                # It is either a match or a mismatch.
                else:

                    self.alignedReadCountsByPosition[referencePositionZeroBased] += 1
                    # TODO: Sometimes, there is no query position. What if this is 'None'?
                    # This is because in the case of indels, there is no query position.

                    queryPosition = int(pileupRead.query_position)

                    currentBase = pileupRead.alignment.query_sequence[queryPosition].upper()

                    if (currentBase.upper() == referenceBase.upper()):
                        self.matchReadCountsByPosition[referencePositionZeroBased] += 1
                    else:
                        self.mismatchReadCountsByPosition[referencePositionZeroBased] += 1
                        self.storeAlleleSpecificPolymorphism(pileupRead.alignment.query_name,
                                                                              referencePositionZeroBased,
                                                                              currentBase.upper())

                    if (currentBase == 'A'):
                        self.ACountsByPosition[referencePositionZeroBased] += 1
                    elif (currentBase == 'G'):
                        self.GCountsByPosition[referencePositionZeroBased] += 1
                    elif (currentBase == 'C'):
                        self.CCountsByPosition[referencePositionZeroBased] += 1
                    elif (currentBase == 'T'):
                        self.TCountsByPosition[referencePositionZeroBased] += 1

        # TODO: Wait on this until the end? Maybe it's more useful after i calculate read qualities. Maybe not...
        self.writeAlignmentSummary(join(alignmentOutputDir, 'AlignmentSummary.csv'), minimumSnpPctCutoff)
        self.writePolymorphicAlleles(join(alignmentOutputDir, 'AlleleSpecificPolymorphisms.csv'))
        self.createAlignmentQualityGraph(join(alignmentOutputDir, 'AlignmentQuality.png'))
        self.createSlidingScaleAlignmentQualityGraph(
            join(alignmentOutputDir, 'SlidingScaleAlignmentQuality.png'))

        # TODO: I commented this part out because it takes X^2 amount of time.
        # TODO: Put this as an optional parameter. Right now I'm just not doing it.
        #self.writeLDMatrix(alignmentOutputDir)

        # Loop in the opposite direction, read-by-read then position by position.
        # That loop is more clear for analyzing individual read quality
        # I can use fetch() to get the reads.
        alignedReads = bamfile.fetch()


        self.totalAlignedReadCount = bamfile.count()

        for index, read in enumerate(alignedReads):
            queryName = read.query_name
            #print('Looking at aligned read#' + str(index) + '/' + str(self.totalAlignedReadCount) + ':' + str(queryName))
            alignedPairs = read.get_aligned_pairs()
            alignedSequence = read.query_sequence

            # Couple of Flags, because I only want to calculate qualities within the aligned region.
            BeforeReference = True
            AfterReference = False

            correctAlignedBaseScore = 0
            # Loop through each aligned base to calculate a alignment score.
            for alignedPairIndex, alignedPair in enumerate(alignedPairs):

                readPos = alignedPair[0]
                referencePos = alignedPair[1]

                if (BeforeReference):
                    # Record nothing until we are inside the reference.
                    if (referencePos is None):
                        pass
                    elif (referencePos >= read.reference_start):
                        BeforeReference = False

                if (not BeforeReference and not AfterReference):
                    # Inside the aligned region now.
                    # Calculate read and reference bases.
                    if (readPos == None):
                        readSeq = '-'
                    else:
                        readSeq = alignedSequence[int(readPos)]
                    if (referencePos == None):
                        refSeq = '-'
                    else:
                        refSeq = alignmentRef[int(referencePos)]

                    # Increment/Decrement the score.
                    # TODO: Pass in Alignment penalties, the metrics should maybe change. Maybe. This is my opinion.
                    matchScore = 1
                    mismatchScore = 0
                    deletionScore = 0
                    insertionScore = -1

                    if (readSeq == refSeq):
                        # print ('Pair(read,ref) ' + str(alignedPair) + ' looks like a match: ' + refSeq)
                        correctAlignedBaseScore += matchScore
                    elif (readSeq == '-'):
                        # print ('Pair(read,ref) ' + str(alignedPair) + ' looks like a deletion: ' + alignmentRef[int(referencePos)])
                        correctAlignedBaseScore += deletionScore
                    elif (refSeq == '-'):
                        # print ('Pair(read,ref) ' + str(alignedPair) + ' looks like an insertion: ' + alignedSequence[int(readPos)])
                        correctAlignedBaseScore += insertionScore
                    else:
                        # print ('Pair(read,ref) ' + str(alignedPair) + ' looks like a mismatch: ' + refSeq + '->' + readSeq)
                        correctAlignedBaseScore += mismatchScore

                    # If we're outside the reference then we should set a flag. No need to calculate outside the region.
                    if (referencePos is None):
                        pass
                    elif (referencePos >= read.reference_end):
                        AfterReference = True
                    if (readPos is None):
                        pass
                    elif (readPos > read.query_alignment_end):
                        AfterReference = True

            # Read Score = MatchScore / Alignment Length.
            # We're only looking at aligned portion, so we use reference_length.
            readScore = correctAlignedBaseScore / read.reference_length

            # I have seen with ion torrent the reads are perfect. a readScore of 1 breaks the log10 calculation.
            # I give an arbitrary high read score.
            if (readScore == 1):
                #print('READ SCORE GREATER THAN 1: ' + str(read))
                readScore = .999999


            # Pretend this readScore is a correctly mapped percentage. we can calculate a Q(phred) score for each read.
            # Q = -10 log(10) P
            # Where P is the basecalling error probability. (P = 1 - readScore)
            phredScore = (-10) * (log10(1 - readScore))

            #print ('This read has a score of:' + str(readScore))
            #print ('Phred Score:' + str(phredScore))

            # But for the Scatter plot, I want quality vs total read length (query_length). A bit dishonest but probably informative.
            self.alignedReadLengths.append(read.query_length)
            self.alignedReadPhredScores.append(phredScore)
            self.alignedReadQualityScores.append(readScore)


            #Loop thru the reference to record homopolymers in the reads.
            for referencePos in range(0,len(self.referenceSequence)):

                homopolymerLength = self.homopolymerSizesByPosition[referencePos]
                # TODO 3 is arbitrary for min length, I think 4 or 3 are important. Change this parameter for interesting results.
                if (homopolymerLength >= 3):

                    homopolymerBase = self.referenceSequence[referencePos]

                    #print('Homopolymer(' + homopolymerBase + ') of length ' + str(self.homopolymerSizesByPosition[referencePos]) + ' detected at position ' + str(referencePos))

                    # What if I just get the aligned pairs for this homopolymer?
                    # Select from Aligned Pairs where referencePosition (alignedPair[1]) is mapped to the homo region, meaning:
                    # between (referencePos : referencePos + homopolymerLength)
                    alignedHomopolymerReadPositions = [t[0] for t in alignedPairs
                        if t[1] is not None and (referencePos <= t[1] < referencePos + homopolymerLength)]

                    # Loop through the homopolymer aligned pairs and count the bases that match.
                    # TODO: I think there's a problem with this, None of my homopolymers are longer than the reference.
                    # TODO: Some should be. Am I detecting the insertions properly?

                    #print('I found ' + str(len(alignedHomopolymerReadPositions)) + ' aligned read positions:' + str(alignedHomopolymerReadPositions))
                    baseCountMatch = 0




                    # Check the aligned bases, count the bases that match.
                    for alignedReadPosition in alignedHomopolymerReadPositions:
                        # if the readPosition is None, that means there is a "deletion"
                        # so we don't count it.
                        if(alignedReadPosition is not None):
                            readSeq = alignedSequence[int(alignedReadPosition)]
                            #print('ReadSeq=' + readSeq)
                            if ((readSeq) == homopolymerBase):
                                baseCountMatch += 1


                    # Check the base BEFORE the aligned region, these bases are not aligned with the homopolymer region.
                    if(len(alignedHomopolymerReadPositions) > 0 and alignedHomopolymerReadPositions[0] is not None):
                        insertedBaseCount = 1
                        previousBase = str(alignedSequence[int(alignedHomopolymerReadPositions[0])-insertedBaseCount])

                        while(previousBase == homopolymerBase):
                            insertedBaseCount += 1
                            previousBase = str(alignedSequence[int(alignedHomopolymerReadPositions[0]) - insertedBaseCount])
                            baseCountMatch += 1



                    # TODO: I...hope the inserted bases are always to the left of the homopolymer sequences.
                    # May not be true if the aligner is different.
                    # I'm not checking for matching homopolymer inserted bases to the right of the region, but this
                    # is a potential bug.

                    #print (read.query_name + ' has this many aligned bases: ' + str(baseCountMatch) + '\n')

                    # Store this homopolymer in our list of tuples.
                    # (Reference Pos, Reference Base, Reference Homopolymer Length, Read Homopolymer Length)
                    if(baseCountMatch != 0):
                        self.homopolymerTuples.append((referencePos, homopolymerBase, homopolymerLength, baseCountMatch))

        self.printHomopolymerTuples(alignmentOutputDir)



        # TODO: Integrate with the code in nit-picker. I feel like there is duplicate logic.
        # Instead of just printing this, I could call the "output read stats and diagrams" or whatever method that is in nitpicker.
        createScatterPlot('Calculated Mapped-Read Qualities'
                          , self.alignedReadLengths
                          , self.alignedReadPhredScores
                          , "Read Lengths"
                          , "Avg Read Quality(Phred)"
                          , join(alignmentOutputDir, 'CalculatedMappedReadQualities.png'))


    def printHomopolymerTuples(self, outputDirectory):
        tupleOutputFile = createOutputFile(join(outputDirectory, 'ReadHomopolymer.csv'))
        tupleOutputFile.write('Reference_Pos,Reference_Base,Homopolymer_Ref_Length,Homopolymer_Read_Length\n')
        for homopolymerTuple in self.homopolymerTuples:
            tupleOutputFile.write(str(homopolymerTuple[0]) + ',' + str(homopolymerTuple[1]) + ',' + str(homopolymerTuple[2]) + ',' + str(homopolymerTuple[3]) + '\n')

        if (len(self.homopolymerTuples) == 0):
            maxHomopolymerReferenceLength = 0
            maxHomopolymerReadLength = 0
        else:

            maxHomopolymerReferenceLength = max([t[2] for t in self.homopolymerTuples])
            maxHomopolymerReadLength = max([t[3] for t in self.homopolymerTuples])

        #print ('Max Length =' + str(maxHomopolymerLength))

        # TODO: Make a file for each base/length combo, print them to the same directory.
        # TODO: Summarize the information. CSV LIKE
        # Make a CSV like: Base, Length, Percent3, Percent4, Percent5, Percent6
        tupleSummaryOutputFile = createOutputFile(join(outputDirectory, 'HomopolymerSummary.csv'))
        tupleSummaryOutputFile.write('Base,Reference_Length,Count,Avg_Observed_Length,Median_Observed_Length')
        # Loop header to maximum observed homopolymer size.
        for i in range(1,maxHomopolymerReadLength + 1):
            tupleSummaryOutputFile.write(',%' + str(i))
        tupleSummaryOutputFile.write('\n')

        # Sort first by length, then by Base. 3 is the minimum homopolymer length we care about, maybe that changes.
        for referenceHomoLength in range(3, maxHomopolymerReferenceLength + 1):
            for nucleotide in ['A','G','C','T']:
                # Get list of homopolymer lengths that mapped to the current nucleotide and length
                currentHomopolymerData = [t[3] for t in self.homopolymerTuples if (t[1]==nucleotide and t[2]==referenceHomoLength)]
                currentHomopolymerData = list(map(int,currentHomopolymerData))

                if(len(currentHomopolymerData) > 0):
                    #averageObservedLength = round(sum(float(d) for d in currentHomopolymerData) / len(currentHomopolymerData) , 3)
                    #medianObseredLength = median(currengHomopolymerData)
                    averageObservedLength = round(mean(currentHomopolymerData) , 3)
                    medianObservedLength = median(currentHomopolymerData)
                else:
                    averageObservedLength = 0
                    medianObservedLength = 0

                #print('Summarizing homopolymers for length ' + str(referenceHomoLength) + ':' + nucleotide)
                #print('Data:' + str(currentHomopolymerData))
                #print('Average:' + str(averageObservedLength))

                tupleSummaryOutputFile.write(str(nucleotide) + ',' + str(referenceHomoLength) + ',' + str(len(currentHomopolymerData)) + ',' + str(averageObservedLength) + ',' + str(medianObservedLength))
                # Now we loop through observed read homopolymer lengths.
                for readHomoLength in range(1, maxHomopolymerReadLength + 1):
                    #print ('Looking for reads of length: ' + str(readHomoLength))

                    if(len(currentHomopolymerData) > 0):
                        percentObserved = round(currentHomopolymerData.count(readHomoLength) / len (currentHomopolymerData) , 3)
                    else:
                        percentObserved = 0
                    tupleSummaryOutputFile.write(',' + str(percentObserved))

                tupleSummaryOutputFile.write('\n')













def createCollectionFromReadFile(readFile):
    global readInputFormat
    # print ('loading reads from:' + readFile)

    # Determine Input File Type
    if (".fasta" == readFile[-6:] or ".fa" == readFile[-3:]):
        readInputFormat = "fasta"
        # raise Exception('Fasta files are not supported.  You can try to add this functionality if you want.')
    elif (".fastq" == readFile[-6:] or ".fq" == readFile[-3:]):
        readInputFormat = "fastq"
    else:
        print(
            'I expect a .fasta or .fastq format for read input. Alternatively, specify a directory containing read inputs. Please check your input.')
        raise Exception('Bad Read Input Format')

    newCollectionObject = MinionReadCollection(list(parse(readFile, readInputFormat)))
    newCollectionObject.readInputFormat = readInputFormat
    newCollectionObject.readFileLocation = readFile

    return newCollectionObject


# I made a bug here.
# TODO: readStats is a dictionary with stats for pass, reject lenght, and reject quality reads.
# alignedReadCount, meanAlignedReadLength, and meanCalculatedQuality are not split up in this way...
# To fix: allow MinION ReadCollection to store rejected(len and qual) reads. That way NitPicker doesn't have to split them up.
# That way, I can put these methods inside the class where they belong.
def writeReadStats(readStats, outputDirectory, alignedReadCount, meanAlignedReadLength, meanCalculatedQuality ):
    outputFileName = join(outputDirectory, 'Read_Summary.txt')
    outputFile = createOutputFile(outputFileName)

    for key in readStats.keys():
        outputFile.write(writeReadStatsSub(key, readStats[key], alignedReadCount, meanAlignedReadLength, meanCalculatedQuality ))
    outputFile.close()


def writeReadStatsSub(readType, readStats, alignedReadCount, meanAlignedReadLength, meanCalculatedQuality ):
    statsSubText = ''

    if (len(readStats) > 0):
        readLengths = readStats[0]
        readPhredQualities = readStats[1]
        avgReportedPhredQuality = mean(readPhredQualities)
        # Phred quality -> Accuracy
        # P = 10 ^ (-Q / 10)
        # Where P is probability of basecalling error.
        avgReportedReadAccuracy = 1 - pow(10 , (-1) * avgReportedPhredQuality / 10  )

        statsSubText += (readType + '\n'
            + 'Total Read Count        :' + str(len(readLengths)) + '\n'
            + 'Minimum Length          :' + str(int(amin(readLengths))) + '\n'
            + 'Maximum Length          :' + str(int(amax(readLengths))) + '\n'
            + 'Mean Length             :' + str(mean(readLengths)) + '\n'
            + 'Mean Reported Phred     :' + str(mean(readPhredQualities)) + '\n'
            + 'Mean Reported Accuracy  :' + str(avgReportedReadAccuracy) + '\n'
            + 'Aligned Read Count      :' + str(alignedReadCount) + '\n'
            + 'Mean Aligned Length     :' + str(meanAlignedReadLength) + '\n'
            + 'Mean Calculated Quality :' + str(meanCalculatedQuality) + '\n'
            + '\n'
            )

    else:
        statsSubText += (readType + '\n'
                         + 'Total Read Count :0\n\n'
                         )

    return statsSubText
