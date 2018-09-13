# This file is part of HLA-allele-analysis.
#
# HLA-allele-analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HLA-allele-analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with HLA-allele-analysis. If not, see <http://www.gnu.org/licenses/>.

import xml.etree.ElementTree
import sys
import os
from os import listdir, makedirs
from os.path import isfile, join, split
from Parse_IMGT_HLA.HLA_Allele import HLA_Allele
from Parse_IMGT_HLA.HLA_Allele_Group import HLA_Allele_Group
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo
import traceback
import subprocess
from _ast import Or

#Class I Genes
classIGenes = ['A','B','C']
justBForTesting = ['B']
justEForTesting = ['E']

genesForAnalysis = classIGenes
#genesForAnalysis = justEForTesting
#genesForAnalysis = justBForTesting

# This method is a directory-safe way to open up a write file.
# TODO Put this in a common methods section somehow.
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not os.path.isdir(tempDir):
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput

# Print a list of alleles in Fasta format.
def printFasta(alleleList, fileName, printAPDSequences, printIntronSequences, printFullLenMinusUTRs):
    outputFile = createOutputFile(fileName)     
    for currentAllele in alleleList:    
        
        currentHeader = str(currentAllele.getFastaHeader(printAPDSequences, printIntronSequences, printFullLenMinusUTRs))
        currentSequence = str(currentAllele.getFastaSequence(printAPDSequences, printIntronSequences, printFullLenMinusUTRs))
        
        if(len(currentSequence) > 0):
            outputFile.write('>' + currentHeader + '\n')
            outputFile.write(currentSequence + '\n')
        
    outputFile.close()

# Load the XML file, and parse to extract the HLA alleles.
def createAlleleList(inputFileName, outputDirectory):
    #document root is an "alleles" node.  It has a series of "allele" nodes underneath it.  Nothing else, according to the schema.
    print('Loading input XML file...')
    documentRoot = xml.etree.ElementTree.parse(inputFileName).getroot()
    print(str(len(documentRoot)) + ' alleles found.')

    alleleFullList = []
    # Loop through allele nodes and create allele objects
    for index, alleleNode in enumerate(documentRoot):

        #if (index % 1000 == 0):
        #    This goes quickly, don't need to print a progress bar.
        #    print ('Allele#:' + str(index)) 
        #    pass
        currentAllele = HLA_Allele()
        currentAllele.outputDirectory = outputDirectory
        currentAllele.parseAlleleNomenclature(alleleNode)
        #print ('Assigning output directory:' + outputDirectory)
        alleleFullList.append(currentAllele)
    
    return alleleFullList

# This is a short method used for sorting alleles into allele groups.
def getGroupIndex(alleleGroups, Gene, AlleleGroup):
    for i in range(0,len(alleleGroups)):
        currentAlleleGroup = alleleGroups[i]
        if (currentAlleleGroup.AlleleGroup == AlleleGroup and currentAlleleGroup.Gene == Gene):
            return i
    return -1

def getGeneIndex(alleleGenes, Gene):
    for i in range(0,len(alleleGenes)):
        currentAlleleGene = alleleGenes[i]
        if (currentAlleleGene.Gene == Gene):
            return i
    return -1

# Split the big list of alleles into allele groups.
# In this instance, an Allele Group is a combination of the Gene and the first field of HLA nomenclature.
def getAlleleGroups(alleleFullList):
    print ('Splitting Alleles into Groups.')
    alleleGroups = []   
    for allele in alleleFullList:
        groupIndex = getGroupIndex(alleleGroups, allele.geneName, allele.alleleGroup)

        # No, the group doesn't already exist
        if (groupIndex == -1):

            currentAlleleGroup = HLA_Allele_Group()
            currentAlleleGroup.Alleles = []
            #print('Group should be empty:' + str(len(currentAlleleGroup.Alleles)))
            currentAlleleGroup.Gene = allele.geneName
            currentAlleleGroup.AlleleGroup = allele.alleleGroup
            currentAlleleGroup.Alleles.append(allele)
            currentAlleleGroup.FileName = 'HLA_' + allele.geneName + '_' + allele.alleleGroup + '.fasta'
            alleleGroups.append(currentAlleleGroup)

        # Yes, this group exists already.  Add this to the list.
        else:
            alleleGroups[groupIndex].Alleles.append(allele)

    return alleleGroups

def getAlleleGenes(alleleFullList):
    print ('Splitting Alleles into Genes.')
    alleleGenes = []   
    for allele in alleleFullList:
        geneIndex = getGeneIndex(alleleGenes, allele.geneName)

        # No, the group doesn't already exist
        if (geneIndex == -1):

            currentAlleleGroup = HLA_Allele_Group()
            currentAlleleGroup.Alleles = []
            #print('Group should be empty:' + str(len(currentAlleleGroup.Alleles)))
            currentAlleleGroup.Gene = allele.geneName
            currentAlleleGroup.Alleles.append(allele)
            currentAlleleGroup.FileName = 'HLA_' + allele.geneName + '.fasta'
            alleleGenes.append(currentAlleleGroup)

        # Yes, this group exists already.  Add this to the list.
        else:
            alleleGenes[geneIndex].Alleles.append(allele)

    return alleleGenes

def reverseComplementList(alleleFullList):
    reversedAlleles = []   
    for allele in alleleFullList:
        reversedAlleles.append(allele.reverseComplement())
    return reversedAlleles

def forwardComplementList(alleleFullList):
    forwardAlleles = []   
    for allele in alleleFullList:
        forwardAlleles.append(allele.forwardComplement())
    return forwardAlleles

def reverseList(alleleFullList):
    reversedAlleles = []   
    for allele in alleleFullList:
        reversedAlleles.append(allele.reverse())
    return reversedAlleles


def performClustalWAlignmentsForGroupwiseReference(outputDirectory, alleleFullList, runAlignments):
    print ('Performing clustalW Alignments and finding consensus sequences')

    # Create the output directories for clustalw.
    clustalwOutputDirectory = join(outputDirectory, 'ClustalwAlignmentsAPD')
    if not os.path.isdir(clustalwOutputDirectory):
        os.mkdir(clustalwOutputDirectory)

    clustalwConsensusOutputDirectory = clustalwOutputDirectory.replace('Alignments','Consensus')
    if not os.path.isdir(clustalwConsensusOutputDirectory):
        os.mkdir(clustalwConsensusOutputDirectory)
    #clustalwAlignmentScriptFile = createOutputFile(clustalwAlignmentScriptFileName)


    alleleGroups = getAlleleGroups(alleleFullList)
    for index, alleleGroup in enumerate(alleleGroups):
        
        print('(' + str(index + 1) + '/' + str(len(alleleGroups)) + '): HLA-' + alleleGroup.Gene + '*' + alleleGroup.AlleleGroup)

        if (True):
        #if (alleleGroup.Gene in genesForAnalysis):
            outputGroupFileName = join(outputDirectory, 
                join('AlleleGroupsAPD',alleleGroup.FileName))

            clustalwAlignmentOutputFileName = outputGroupFileName.replace(
                '.fasta','.aln').replace('/AlleleGroupsAPD/','/ClustalwAlignmentsAPD/')

            clustalwConsensusOutputFileName = outputGroupFileName.replace('/AlleleGroupsAPD/','/ClustalwConsensusAPD/')
     
            # If the alignment does not already exist
            if not (os.path.isfile(clustalwAlignmentOutputFileName)):

                # if there is more than one allele in the group
                if (len(alleleGroup.Alleles) > 1):
                    print (str(len(alleleGroup.Alleles)) + ' Alleles Found.')                

                    clustalwCommandline = ClustalwCommandline("clustalw", infile=outputGroupFileName, outfile=clustalwAlignmentOutputFileName)
                    print ('ClustalW Alignment Commandline:\n' + str(clustalwCommandline))

                    if (runAlignments):

                       # print ('Performing Clustalw Alignment...')
                        #clustalwAlignmentScriptFile.write(str(clustalwCommandline) + '\n') 
                        #Perform the alignment
                        clustalwCommandline()

                        if (os.path.isfile(clustalwAlignmentOutputFileName)):  
                            # If consensus does not exist yet
                            if not (os.path.isfile(clustalwConsensusOutputFileName)):  
                                #Perform the consensus
                                alignmentType = 'clustal'

                                align = AlignIO.read(clustalwAlignmentOutputFileName, alignmentType)
                            
                                #print ('Consensus FileName = ' + clustalwConsensusOutputFileName)
                            
                                summary_align = AlignInfo.SummaryInfo(align)
                                
                                dumb_consensus = summary_align.dumb_consensus()
                                #print('LengthDumbConsensus:' + str(len(dumb_consensus)))
                                gap_consensus = summary_align.gap_consensus()
                                #print('LengthGapConsensus:' + str(len(gap_consensus)))
                                #print ('Consensus=' + str(gap_consensus))

                                # Print Consensus to fasta.

                                # I can cheat and just create an HLA_Allele object, and print that.
                                currentAllele = HLA_Allele()
                                # I think I'll use the dumb consensus.  The only difference is that a gap consensus allows gaps.
                                currentAllele.APDSequence = str(dumb_consensus)
                                currentAllele.alleleName = os.path.basename(clustalwConsensusOutputFileName).replace('.fasta','')
                                currentAllele.outputDirectory = outputDirectory
                                #print ('Consensus2=' + currentAllele.APDSequence)
                                printFasta([currentAllele], clustalwConsensusOutputFileName, True, False, False)

                                pass
                            else:
                                print ('Consensus file ' + clustalwConsensusOutputFileName + ' already exists.  Moving on...')
                        else:
                            print ('Cannot find alignment file after completing alignment:' + clustalwAlignmentOutputFileName) 
                            #raise Exception('Cannot find alignment file after completing alignment:' + clustalwAlignmentOutputFileName) 
                            pass

                    else:
                        print ('Not running Alignments because you told me not to.')

                # There is only one allele in this group.
                else:
                    print ('Only one allele found')
                    currentGene = alleleGroup.Alleles[0].geneName

                    # Only class 1.
                    #if (currentGene in genesForAnalysis):
                    currentAllele = HLA_Allele()
                    currentAllele.sequence = alleleGroup.Alleles[0].sequence
                    currentAllele.alleleName = os.path.basename(clustalwConsensusOutputFileName).replace('.fasta','')
                    printFasta([currentAllele], clustalwConsensusOutputFileName, False, False, False)
          

            else:
                print ('Alignment file ' + clustalwAlignmentOutputFileName + ' already exists.  Moving on...')

        else:
            print ('Skipping alignment, because this gene isnt included in genesForAnalysis')


def createBlastReferences(alleleList, outputDirectory):
    print ('Creating a list of BLAST reference files.')
    #outputAPDFileName = join(outputDirectory, 'HLA_Alleles_APD.fasta')
    #printFasta(alleleFullList, outputAPDFileName, True)
    
    blastRefOutputDir = join(outputDirectory, 'blast_references')
    
    # Iterate each gene.
    for gene in getAlleleGenes(alleleList):
        
        blastReferenceFileName = join(blastRefOutputDir, 'HLA_' + gene.Gene + '_BlastReference.fasta')
        print ('Writing to this file:' + blastReferenceFileName)
        
        # Alleles to print in this gene's reference
        geneBlastAlleles = []
               
        for group in getAlleleGroups(alleleList):
            #print (gene)
            print (group.Gene)
            print (group.AlleleGroup)
            
            
            
            if (group.Gene == gene.Gene):
                # foreach alllele in the group
                for allele in group.Alleles:
                    # If allele has a full length sequence
                    if allele.APDSequence() is not None:
                        geneBlastAlleles.append(allele)
                        break
                        
        printFasta(geneBlastAlleles, blastReferenceFileName, True, False, False)
                
        
        
        
        
        
# I already made a method like this.
#def getDistinctGenes(alleleList):    
#    distinctGenes = []
#    for allele in alleleList:
#        if allele.geneName not in distinctGenes:
#            #print ('new gene found:' + allele.geneName)
#            distinctGenes.append(allele.geneName)
#    return distinctGenes
    

        
def createAlleleSequenceOutputFiles(alleleFullList, outputDirectory):
    # Main allele output file
    outputAllAllelesFileName = join(outputDirectory, 'HLA_Alleles.fasta')
    printFasta(alleleFullList, outputAllAllelesFileName, False, False, False)
    
    # Antigen Presenting Domain allele output file
    outputAPDFileName = join(outputDirectory, 'HLA_Alleles_APD.fasta')
    printFasta(alleleFullList, outputAPDFileName, True, False, False)



    # Full Length Allele output filea
    outputFullLengthAllelesFileName = join(outputDirectory, 'HLA_Alleles_Full_Length.fasta')
    fullLengthAlleles=[]
    for currentAllele in alleleFullList:
        header = currentAllele.getFastaHeader(False, False, False)
        #if(
        #    '(3UTR, 5UTR, EX_1, EX_2, EX_3, EX_4, EX_5, EX_6, EX_7, EX_8, IN_1, IN_2, IN_3, IN_4, IN_5, IN_6, IN_7)' in header or
        #    '(3UTR, 5UTR, EX_1, EX_2, EX_3, EX_4, EX_5, EX_6, EX_7, IN_1, IN_2, IN_3, IN_4, IN_5, IN_6)' in header
        #   ):
        if(
            # TODO: If it has 2 UTR sequences, then this is a valid sequence. 
            # I'm sure this logic is wrong sometimes.
            '3UTR' in header and '5UTR' in header
            ):
            fullLengthAlleles.append(currentAllele)
    printFasta(fullLengthAlleles, outputFullLengthAllelesFileName, False, False, False)

    # Intron sequences aligned. We were using this for evolutionary phylo trees.
    intronOutputFileName = join(outputDirectory, 'HLA_Allele_Intron_Sequences.fasta')
    printFasta(fullLengthAlleles, intronOutputFileName, False, True, False)

    # Output file full length, without the UTRs.
    fullLenMinusUTRFileName = join(outputDirectory, 'HLA_Allele_FullLen_Minus_UTRs.fasta')
    printFasta(fullLengthAlleles, fullLenMinusUTRFileName, False, False, True)

    
# Output Allele information in Fasta format.
def printAlleleGroupsAndInfo(inputFileName, alleleFullList, outputDirectory):
    print ('Creating a fasta reference for all HLA Alleles:' + join(outputDirectory, 'HLA_Alleles.fasta') )

    # Allele Information Output File
    alleleInfoOutputFilename = join(outputDirectory, 'AlleleInfo.txt')
    alleleInfoOutputFile = createOutputFile(alleleInfoOutputFilename)      

 
    # Reverses and Complements
    # TODO: This might be useful in the future. For now lets' skip it.
    #outputAPDRevComFileName = join(outputDirectory, 'HLA_Alleles_APD_RevCom.fasta')
    #printFasta(reverseComplementList(alleleFullList), outputAPDRevComFileName, True)
    #outputAPDForComFileName = join(outputDirectory, 'HLA_Alleles_APD_ForCom.fasta')
    #printFasta(forwardComplementList(alleleFullList), outputAPDForComFileName, True)
    #outputAPDRevFileName = join(outputDirectory, 'HLA_Alleles_APD_Rev.fasta')
    #printFasta(reverseList(alleleFullList), outputAPDRevFileName, True)

    alleleInfoOutputFile.write('Input File:' + inputFileName + '\n')
    alleleInfoOutputFile.write('Output Directory:' + outputDirectory + '\n\n')

    # Print outputfiles and info for each allele group.
    print ('Generating output files for each HLA Allele Group')
    alleleGroups = getAlleleGroups(alleleFullList)
    for index, alleleGroup in enumerate(alleleGroups):
        
        print('(' + str(index + 1) + '/' + str(len(alleleGroups)) + '): HLA-' + alleleGroup.Gene + '*' + alleleGroup.AlleleGroup)
        alleleInfoOutputFile.write('HLA-' + alleleGroup.Gene + '*' 
            + alleleGroup.AlleleGroup + ' contains ' + str(len(alleleGroup.Alleles))
            + ' alleles.\n')

        outputGroupFileName = join(outputDirectory, 
            join('AlleleGroupsAPD',alleleGroup.FileName))

        alleleInfoOutputFile.write('Sorted Group fasta: ' + outputGroupFileName + '\n')
        #alleleInfoOutputFile.write('Alignment: ' + clustalwAlignmentOutputFileName + '\n')
        #alleleInfoOutputFile.write('Consensus: ' + clustalwConsensusOutputFileName + '\n')
        
        # if there is more than one allele in the group
        if (len(alleleGroup.Alleles) > 1):
            print (str(len(alleleGroup.Alleles)) + ' Alleles Found.')

            # Print allele group to a fasta file
            printFasta(alleleGroup.Alleles, outputGroupFileName, True, False, False)

   

        # There is only one allele in this group.
        else:
            print ('Only one allele found')

            printFasta([alleleGroup.Alleles[0]], outputGroupFileName, True, False, False)

    alleleInfoOutputFile.close()

# Combine the consensus sequences into an HLA Groupwise Reference
def combineGroupConsensusIntoReference(outputDirectory):
    print("Combining Group Consensus\' into a Groupwise Reference")
    #consensusList = ['ClustalwConsensus', 'MuscleConsensus']
    consensusList = ['ClustalwConsensusAPD']

    for consensusName in consensusList:
        consensusOutputSubDir = join(outputDirectory, consensusName)

        referenceOutputFileName = join(
            outputDirectory
            , consensusName + '.HLA.Groupwise.Reference.fasta')

        try:
            sequenceList = []
            fileNames = [f for f in listdir(consensusOutputSubDir) if isfile(join(consensusOutputSubDir, f))]
            fileNames = sorted(fileNames)
            for fileName in fileNames:
                if('.fasta' in fileName): 
                    records = SeqIO.parse( join(consensusOutputSubDir, fileName) , "fasta")
                    for index, record in enumerate(records):
                        currentAllele = HLA_Allele()
                        currentAllele.sequence = str(record.seq)
                        currentAllele.alleleName = str(record.id)
                        currentAllele.outputDirectory = outputDirectory
                        sequenceList.append(currentAllele)


            #I should check if this printFasta boolean parameter is right.  These should be just the APD Sequences right?  Should
            #No maybe not, because I'm feeding them the consensus sequence directly into the .sequence allele.
            printFasta(sequenceList, referenceOutputFileName, False, False, False)

            
        except Exception:
            print ('Unexpected problem when creating HLA Reference for ' + consensusName)
            print (sys.exc_info()[0])
            print (sys.exc_info()[1])
            print (sys.exc_info()[2])


def generateIntron2Consensus(alleleFullList, outputDirectory):
    # TODO: This method does not seem to work anymore. I am not assigning the in2Sequence anywhere.
    # Do I need this code anymore? Why would I want to simulate an Intron 2 consensus sequence?
    for featureName in ['Intron 2']:
        shortFeatureName = featureName.replace(' ', '')
        
        #Im deciding to quit here.  Late enough.  I want to fix this method tomorrow.
        
        print ('Creating a ' + featureName + ' Reference:' + join( join(outputDirectory,shortFeatureName + 'References'), 'HLA_Intron2.fasta') )
        
        intron2Alleles = []
    
        for allele in alleleFullList:
            #TODO fix featuresInFullSequence.  Might work this way.
            if('Intron 2' in allele.featuresInFullSequence):
         
                currentIntron2Allele = allele.copy()
                #TODO I don't know if i'm still gonna use in2Sequence.
                currentIntron2Allele.sequence = allele.in2Sequence
                intron2Alleles.append(currentIntron2Allele)
    
        # Intron 2 output file, for analyizing *just* the intron 2
        outputIn2FileName = join(join(outputDirectory,'Intron2References'), 'HLA_Intron2.fasta')
        printFasta(intron2Alleles, outputIn2FileName, False, False, False)
    
    
        # Print outputfiles and info for each allele group.
        print ('Generating output files for each HLA Allele Group')
        alleleGroups = getAlleleGroups(intron2Alleles)
        alleleGenes = getAlleleGenes(intron2Alleles)
        combinedAlleleGroups = alleleGroups + alleleGenes
        for index, alleleGroup in enumerate(combinedAlleleGroups):
            
            print('(' + str(index + 1) + '/' + str(len(combinedAlleleGroups)) + '): ' + alleleGroup.FileName)
    
            outputGroupFileName = join(outputDirectory, 
                join('Intron2References',alleleGroup.FileName))
    
            clustalwAlignmentOutputFileName = outputGroupFileName.replace('.fasta','.aln')
            clustalwConsensusOutputFileName = outputGroupFileName.replace('.fasta','.consensus.fasta')
            # if there is more than one allele in the group
            if (len(alleleGroup.Alleles) > 1):
                print (str(len(alleleGroup.Alleles)) + ' Alleles Found.  Writing to file: ' + outputGroupFileName)
    
                # Print allele group to a fasta file
                # So this should actually be a false, I don't want to use the APD sequence here.
                printFasta(alleleGroup.Alleles, outputGroupFileName, False, False, False)
    
                if (not os.path.isfile(clustalwAlignmentOutputFileName)): 
                    clustalwCommandline = ClustalwCommandline("clustalw", infile=outputGroupFileName, outfile=clustalwAlignmentOutputFileName)
                    print ('Performing  ClustalW Alignment : \n' + str(clustalwCommandline))
    
                    #Perform the alignment
                    clustalwCommandline()
        
                    # sanity check to make sure it exists.
                    if (os.path.isfile(clustalwAlignmentOutputFileName)):  
                        # If consensus does not exist yet
                        if not (os.path.isfile(clustalwConsensusOutputFileName)):  
                            #Find consensus
                            alignmentType = 'clustal'    
                            align = AlignIO.read(clustalwAlignmentOutputFileName, alignmentType)
                        
                            print ('Consensus FileName = ' + clustalwConsensusOutputFileName)
                        
                            #print('Alignment:' + str(align))
                            summary_align = AlignInfo.SummaryInfo(align)
    
                            dumb_consensus = summary_align.dumb_consensus()
                            #print('LengthDumbConsensus:' + str(len(dumb_consensus)))
                            gap_consensus = summary_align.gap_consensus()
                            #print('LengthGapConsensus:' + str(len(gap_consensus)))
        
                            # Print Consensus to fasta.    
                            # I can cheat and just create an HLA_Allele object, and print that.
                            currentAllele = HLA_Allele()
                            currentAllele.APDSequence = str(dumb_consensus)
                            currentAllele.alleleName = os.path.basename(clustalwConsensusOutputFileName).replace('.fasta','')
                            currentAllele.outputDirectory = outputDirectory
                            #print ('Consensus2=' + currentAllele.APDSequence)
                            printFasta([currentAllele], clustalwConsensusOutputFileName, True, False, False)
        
                            pass
                        else:
                            print ('Consensus file ' + clustalwConsensusOutputFileName + ' already exists.  Moving on...')
                    else:
                        print ('Cannot find alignment file after completing alignment:' + clustalwAlignmentOutputFileName) 
                        #raise Exception('Cannot find alignment file after completing alignment:' + clustalwAlignmentOutputFileName) 
                        pass
    
                else:
                    print('This alignment file ' + clustalwAlignmentOutputFileName + ' already exists.  Moving on...')
                #else:
                #    print ('Not running Alignments because you told me not to.')   
    
            # There is only one allele in this group.
            else:
                print ('Only one allele found. Writing to file: ' + outputGroupFileName)
    
                #writing it out twice, that's kind of silly but whatever.
                printFasta([alleleGroup.Alleles[0]], outputGroupFileName, True, False, False)
                printFasta([alleleGroup.Alleles[0]], clustalwConsensusOutputFileName, True, False, False)
    
        #alleleInfoOutputFile.close()
        
        
def getShortAlleleNameForHeader(featureName):
        shortFeatureName = str(featureName.replace('Simulated ','SIM-').replace('\' ','').replace('Intron','IN').replace('Exon','EX').replace(' ','_'))
        return shortFeatureName


def createFeatureReferences(alleleFullList, outputDirectory):
    featuresList = ['5\' UTR', '3\' UTR',
                    'Exon 1', 'Exon 2', 'Exon 3', 'Exon 4', 'Exon 5', 'Exon 6', 'Exon 7', 'Exon 8',
                    'Intron 1', 'Intron 2', 'Intron 3', 'Intron 4', 'Intron 5', 'Intron 6', 'Intron 7'
                    ]
    for featureName in featuresList:
        shortFeatureName = str(
            featureName.replace('Simulated ', 'SIM-').replace('\' ', '').replace('Intron ', 'IN').replace('Exon ',
                                                                                                          'EX').replace(
                ' ', '_'))

        featureAlleles = []
        featureReferenceOutputDirectory = join(outputDirectory, shortFeatureName + '_sequences')

        print ('Creating a ' + featureName + ' Reference:' + join(featureReferenceOutputDirectory,
                                                                  'HLA_' + shortFeatureName + '.fasta'))

        for allele in alleleFullList:
            if (featureName in allele.featuresInFullSequence):
                currentFeatureAllele = allele.copy()
                currentFeatureAllele.sequence = allele.featuresInFullSequence[featureName]
                # Update the features that are in this sequence
                currentFeatureAllele.featuresInFullSequence = {}
                currentFeatureAllele.featuresInFullSequence[featureName] = allele.featuresInFullSequence[featureName]
                featureAlleles.append(currentFeatureAllele)

        outputFeatureFileName = join(featureReferenceOutputDirectory, 'HLA_' + shortFeatureName + '.fasta')
        printFasta(featureAlleles, outputFeatureFileName, False, False, False)

        # Print outputfiles and info for each allele group.
        print ('Generating output files for each HLA Allele Group')
        alleleGroups = getAlleleGroups(featureAlleles)
        alleleGenes = getAlleleGenes(featureAlleles)
        combinedAlleleGroups = alleleGroups + alleleGenes
        for index, alleleGroup in enumerate(combinedAlleleGroups):
            print('(' + str(index + 1) + '/' + str(len(combinedAlleleGroups)) + '): ' + alleleGroup.FileName)

            outputGroupFileName = join(featureReferenceOutputDirectory, alleleGroup.FileName)

            printFasta(alleleGroup.Alleles, outputGroupFileName, False, False, False)

# This was formarly a main method, no longer because it is a part of prospector now.
#if __name__=='__main__':
def parseImgtHla(inputFileName, outputDirectory):
    try:
        #First arg is the input file.  Second arg is the output directory.
        #inputFileName = sys.argv[1]
        #outputDirectory = sys.argv[2]
        print('*** Generating a Fasta reference file from a IMGT HLA XML. ***')
        print('Input:' + inputFileName + '\nOutput:' + outputDirectory)
        print('Just a second...')

        # TODO use makedirs here, it's better. Probably not necessary though.
        if not os.path.isdir(outputDirectory):
            os.mkdir(outputDirectory)

        # alleleList is all alleles read from the HLA XML file.
        alleleList = createAlleleList(inputFileName, outputDirectory)
        
        # This is a generic method to just make output files for every gene feature.
        # This shoudl overlap with the generateIntron2Consensus method, so I think I'll need to trim that method 
        # so as not to redo efforts. 
        createFeatureReferences(alleleList, outputDirectory)
        
        createBlastReferences(alleleList, outputDirectory)
        
        # Main sequence output files.
        createAlleleSequenceOutputFiles(alleleList, outputDirectory)
        
        # Output many allele references
        # inputFileName is just for some info to print in the allele info file.
        printAlleleGroupsAndInfo(inputFileName, alleleList, outputDirectory)

        # TODO: I decided not to use the false generated intron 2 sequences. This is now commented out.
        # Does this have repercussions?        
        #generateIntron2Consensus(alleleList, outputDirectory)

        # Align and find consensus for allele groups.
        # This takes a long time.        
        #performClustalWAlignmentsForGroupwiseReference(outputDirectory, alleleList, True)
        
        # Make a groupwise reference consensus.        
        #combineGroupConsensusIntoReference(outputDirectory)

        print('Done.  Ben did a great job.')

    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print ('Unexpected problem during execution:')
        print (sys.exc_info()[1])
        raise
