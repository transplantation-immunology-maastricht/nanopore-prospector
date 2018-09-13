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

import os
import time
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO

def seqForwardComplement(s):
    origSeq = Seq(s, IUPAC.ambiguous_dna)
    return str(origSeq.complement())

def seqReverseComplement(s):    
    origSeq = Seq(s, IUPAC.ambiguous_dna)
    return str(origSeq.reverse_complement())

def seqReverse(s):
    #Extended Slice Syntax.  Not sure what that means but step size is -1
    return s[::-1]

# Check if a string is a number.  
# Another dumb method, I can import something for this.
def isNumber(s):
    # TODO: there's gotta be a numpy method or something, this is pointless.
    try:
        float(s)
        return True
    except ValueError:
        return False
    

def readFirstSeqFromFasta(fileName):
    try:

        fastaHandler = open(fileName, "rU")
        parsedFile = SeqIO.parse(fastaHandler, "fasta")
        firstRecord = parsedFile.next()
        
        return str(firstRecord.seq)

    except Exception as e:

        print('EXCEPTION READING SEQUENCE FROM FASTA:' + str(e))




class HLA_Allele:
    def __init__(self):
    
        self.alleleName = ''
        self.allelePrefix = ''
        self.geneName = ''
        self.nomenclatureGroupCount = 0
        self.alleleGroup = ''
        self.specificProtein = ''
        self.synonymousSubstitution = ''
        self.noncodingChanges = ''
        self.expressionSuffix = ''
        self.sequence = ''
        self.featuresInFullSequence = {}
        self.featuresInAPDSequence = {}
        self.geneFilter = ''
        self.outputDirectory = ''


    # TODO: Pass these variables in a smarter way, two booleans doesn't make sense.
    def getFastaHeader(self, printAPDSequences, printIntronSequences, printFullLenMinusUTRs):
        currentFastaHeader = self.alleleName + ' ('


        if(printAPDSequences):
            if('Exon 2' in self.featuresInFullSequence):
                currentFastaHeader += ('EX2, ')
            if('Intron 2' in self.featuresInFullSequence):
                currentFastaHeader += ('IN2, ')
            if('Exon 3' in self.featuresInFullSequence):
                currentFastaHeader += ('EX3, ')

        elif(printIntronSequences):
            #if('5\' UTR' in self.featuresInFullSequence):
            #    currentFastaHeader += ('5UTR, ')
            if('Intron 1' in self.featuresInFullSequence):
                currentFastaHeader += ('IN1, ')
            if('Intron 2' in self.featuresInFullSequence):
                currentFastaHeader += ('IN2, ')
            if('Intron 3' in self.featuresInFullSequence):
                currentFastaHeader += ('IN3, ')
            if('Intron 4' in self.featuresInFullSequence):
                currentFastaHeader += ('IN4, ')
            if('Intron 5' in self.featuresInFullSequence):
                currentFastaHeader += ('IN5, ')
            if('Intron 6' in self.featuresInFullSequence):
                currentFastaHeader += ('IN6, ')
            if('Intron 7' in self.featuresInFullSequence):
                currentFastaHeader += ('IN7, ')
            #if('3\' UTR' in self.featuresInFullSequence):
            #    currentFastaHeader += ('3UTR, ')

        elif(printFullLenMinusUTRs):
            if('Exon 1' in self.featuresInFullSequence):
                currentFastaHeader += ('EX1, ')
            if('Intron 1' in self.featuresInFullSequence):
                currentFastaHeader += ('IN1, ')
            if('Exon 2' in self.featuresInFullSequence):
                currentFastaHeader += ('EX2, ')
            if('Intron 2' in self.featuresInFullSequence):
                currentFastaHeader += ('IN2, ')
            if('Exon 3' in self.featuresInFullSequence):
                currentFastaHeader += ('EX3, ')
            if('Intron 3' in self.featuresInFullSequence):
                currentFastaHeader += ('IN3, ')
            if('Exon 4' in self.featuresInFullSequence):
                currentFastaHeader += ('EX4, ')
            if('Intron 4' in self.featuresInFullSequence):
                currentFastaHeader += ('IN4, ')
            if('Exon 5' in self.featuresInFullSequence):
                currentFastaHeader += ('EX5, ')
            if('Intron 5' in self.featuresInFullSequence):
                currentFastaHeader += ('IN5, ')
            if('Exon 6' in self.featuresInFullSequence):
                currentFastaHeader += ('EX6, ')
            if('Intron 6' in self.featuresInFullSequence):
                currentFastaHeader += ('IN6, ')
            if('Exon 7' in self.featuresInFullSequence):
                currentFastaHeader += ('EX7, ')
            if('Intron 7' in self.featuresInFullSequence):
                currentFastaHeader += ('IN7, ')
            if('Exon 8' in self.featuresInFullSequence):
                currentFastaHeader += ('EX8, ')

        else:
            if('5\' UTR' in self.featuresInFullSequence):
                currentFastaHeader += ('5UTR, ')
            if('Exon 1' in self.featuresInFullSequence):
                currentFastaHeader += ('EX1, ')
            if('Intron 1' in self.featuresInFullSequence):
                currentFastaHeader += ('IN1, ')
            if('Exon 2' in self.featuresInFullSequence):
                currentFastaHeader += ('EX2, ')
            if('Intron 2' in self.featuresInFullSequence):
                currentFastaHeader += ('IN2, ')
            if('Exon 3' in self.featuresInFullSequence):
                currentFastaHeader += ('EX3, ')
            if('Intron 3' in self.featuresInFullSequence):
                currentFastaHeader += ('IN3, ')
            if('Exon 4' in self.featuresInFullSequence):
                currentFastaHeader += ('EX4, ')
            if('Intron 4' in self.featuresInFullSequence):
                currentFastaHeader += ('IN4, ')
            if('Exon 5' in self.featuresInFullSequence):
                currentFastaHeader += ('EX5, ')
            if('Intron 5' in self.featuresInFullSequence):
                currentFastaHeader += ('IN5, ')
            if('Exon 6' in self.featuresInFullSequence):
                currentFastaHeader += ('EX6, ')
            if('Intron 6' in self.featuresInFullSequence):
                currentFastaHeader += ('IN6, ')
            if('Exon 7' in self.featuresInFullSequence):
                currentFastaHeader += ('EX7, ')
            if('Intron 7' in self.featuresInFullSequence):
                currentFastaHeader += ('IN7, ')
            if('Exon 8' in self.featuresInFullSequence):
                currentFastaHeader += ('EX8, ')
            if('3\' UTR' in self.featuresInFullSequence):
                currentFastaHeader += ('3UTR, ')

        # Don't mess with sorting. Just check for each feature, it's easy enough.
        # if(not printAPDSequences):
        #     if('5\' UTR' in self.featuresInFullSequence):
        #         currentFastaHeader += ('5UTR, ')
        #     if('Exon 1' in self.featuresInFullSequence):
        #         currentFastaHeader += ('EX1, ')
        #     if('Intron 1' in self.featuresInFullSequence):
        #         currentFastaHeader += ('IN1, ')
        #
        # if('Exon 2' in self.featuresInFullSequence):
        #     currentFastaHeader += ('EX2, ')
        # if('Intron 2' in self.featuresInFullSequence):
        #     currentFastaHeader += ('IN2, ')
        # if('Exon 3' in self.featuresInFullSequence):
        #     currentFastaHeader += ('EX3, ')
        #
        # if(not printAPDSequences):
        #     if('Intron 3' in self.featuresInFullSequence):
        #         currentFastaHeader += ('IN3, ')
        #     if('Exon 4' in self.featuresInFullSequence):
        #         currentFastaHeader += ('EX4, ')
        #     if('Intron 4' in self.featuresInFullSequence):
        #         currentFastaHeader += ('IN4, ')
        #     if('Exon 5' in self.featuresInFullSequence):
        #         currentFastaHeader += ('EX5, ')
        #     if('Intron 5' in self.featuresInFullSequence):
        #         currentFastaHeader += ('IN5, ')
        #     if('Exon 6' in self.featuresInFullSequence):
        #         currentFastaHeader += ('EX6, ')
        #     if('Intron 6' in self.featuresInFullSequence):
        #         currentFastaHeader += ('IN6, ')
        #     if('Exon 7' in self.featuresInFullSequence):
        #         currentFastaHeader += ('EX7, ')
        #     if('Intron 7' in self.featuresInFullSequence):
        #         currentFastaHeader += ('IN7, ')
        #     if('Exon 8' in self.featuresInFullSequence):
        #         currentFastaHeader += ('EX8, ')
        #     if('3\' UTR' in self.featuresInFullSequence):
        #         currentFastaHeader += ('3UTR, ')
            
        currentFastaHeader = currentFastaHeader[0:len(currentFastaHeader)-2] + ')'
        
        #if(printAPDSequences):
        #    currentFeatures = self.featuresInAPDSequence.copy()
            
        #else:
        #    currentFeatures = self.featuresInFullSequence.copy()
           
        # TODO: Should I sort the feature keys? Pros and Cons. 
        #for featureKey in sorted(currentFeatures.keys()):
        #for featureKey in currentFeatures.keys():
        #    shortFeatureName = str(featureKey.replace('Simulated ','SIM-').replace('\' ','').replace('Intron','IN').replace('Exon','EX').replace(' ','_')) 
        #    currentFastaHeader += (shortFeatureName + ', ')
            
        #currentFastaHeader = currentFastaHeader[0:len(currentFastaHeader)-2] + ')'

        return currentFastaHeader
    
             
    
    def getFastaSequence(self, printAPDSequences, printIntronSequences, printFullLenMinusUTRs):
        if(printAPDSequences):
            return self.APDSequence()
        elif(printIntronSequences):
            return self.intronSequence()
        elif(printFullLenMinusUTRs):
            return self.fullLengthMinusUTRSequence()
        else:
            return self.sequence  


    def fullLengthMinusUTRSequence(self):
        fullLenMinusUTRSeq = ''

        if ('Exon 1' in self.featuresInFullSequence):
            fullLenMinusUTRSeq += str(self.featuresInFullSequence['Exon 1'])
        if ('Intron 1' in self.featuresInFullSequence):
            fullLenMinusUTRSeq += str(self.featuresInFullSequence['Intron 1'])
        if ('Exon 2' in self.featuresInFullSequence):
            fullLenMinusUTRSeq += str(self.featuresInFullSequence['Exon 2'])
        if ('Intron 2' in self.featuresInFullSequence):
            fullLenMinusUTRSeq += str(self.featuresInFullSequence['Intron 2'])
        if ('Exon 3' in self.featuresInFullSequence):
            fullLenMinusUTRSeq += str(self.featuresInFullSequence['Exon 3'])
        if ('Intron 3' in self.featuresInFullSequence):
            fullLenMinusUTRSeq += str(self.featuresInFullSequence['Intron 3'])
        if ('Exon 4' in self.featuresInFullSequence):
            fullLenMinusUTRSeq += str(self.featuresInFullSequence['Exon 4'])
        if ('Intron 4' in self.featuresInFullSequence):
            fullLenMinusUTRSeq += str(self.featuresInFullSequence['Intron 4'])
        if ('Exon 5' in self.featuresInFullSequence):
            fullLenMinusUTRSeq += str(self.featuresInFullSequence['Exon 5'])
        if ('Intron 5' in self.featuresInFullSequence):
            fullLenMinusUTRSeq += str(self.featuresInFullSequence['Intron 5'])
        if ('Exon 6' in self.featuresInFullSequence):
            fullLenMinusUTRSeq += str(self.featuresInFullSequence['Exon 6'])
        if ('Intron 6' in self.featuresInFullSequence):
            fullLenMinusUTRSeq += str(self.featuresInFullSequence['Intron 6'])
        if ('Exon 7' in self.featuresInFullSequence):
            fullLenMinusUTRSeq += str(self.featuresInFullSequence['Exon 7'])
        if ('Intron 7' in self.featuresInFullSequence):
            fullLenMinusUTRSeq += str(self.featuresInFullSequence['Intron 7'])
        if ('Exon 8' in self.featuresInFullSequence):
            fullLenMinusUTRSeq += str(self.featuresInFullSequence['Exon 8'])
        if ('Intron 8' in self.featuresInFullSequence):
            fullLenMinusUTRSeq += str(self.featuresInFullSequence['Intron 8'])
        if ('Exon 9' in self.featuresInFullSequence):
            fullLenMinusUTRSeq += str(self.featuresInFullSequence['Exon 9'])
        if ('Intron 9' in self.featuresInFullSequence):
            fullLenMinusUTRSeq += str(self.featuresInFullSequence['Intron 9'])

        return fullLenMinusUTRSeq


    def APDSequence(self):
        
        # TODO: I wonder if i should just do exon 2, for class II genes.
        # Lets see how well sorting is working.
             
        currentAPDSequence = ''
        
        if('Exon 2' in self.featuresInAPDSequence):
            currentAPDSequence += str(self.featuresInAPDSequence['Exon 2'])            
        if('Intron 2' in self.featuresInAPDSequence):
            currentAPDSequence += str(self.featuresInAPDSequence['Intron 2'])            
        if('Exon 3' in self.featuresInAPDSequence):
            currentAPDSequence += str(self.featuresInAPDSequence['Exon 3'])
            
            
                
        #elif('Exon 2' in self.featuresInAPDSequence and
        #   'Simulated Intron 2' in self.featuresInAPDSequence and
        #   'Exon 3' in self.featuresInAPDSequence):
        #    currentAPDSequence = (
        #        str(self.featuresInAPDSequence['Exon 2']) +
        #        str(self.featuresInAPDSequence['Simulated Intron 2']) +
        #        str(self.featuresInAPDSequence['Exon 3']))
        
        #else:
        #    print('I cannot construct a reliable APD sequence for this allele.' + self.alleleName)
            #for featureKey in self.featuresInAPDSequence.keys():
            #    currentAPDSequence += str(self.featuresInAPDSequence[featureKey])
        #    return ''
        
        return currentAPDSequence


    def intronSequence(self):
        currentIntronSequence = ''

        # How many introns do we expect lol?
        # TODO Put this in a loop, this is silly.

        #if ('5\' UTR' in self.featuresInFullSequence):
        #    currentIntronSequence += str(self.featuresInFullSequence['5\' UTR'])
        if ('Intron 1' in self.featuresInFullSequence):
            currentIntronSequence += str(self.featuresInFullSequence['Intron 1'])
        if ('Intron 2' in self.featuresInFullSequence):
            currentIntronSequence += str(self.featuresInFullSequence['Intron 2'])
        if ('Intron 3' in self.featuresInFullSequence):
            currentIntronSequence += str(self.featuresInFullSequence['Intron 3'])
        if ('Intron 4' in self.featuresInFullSequence):
            currentIntronSequence += str(self.featuresInFullSequence['Intron 4'])
        if ('Intron 5' in self.featuresInFullSequence):
            currentIntronSequence += str(self.featuresInFullSequence['Intron 5'])
        if ('Intron 6' in self.featuresInFullSequence):
            currentIntronSequence += str(self.featuresInFullSequence['Intron 6'])
        if ('Intron 7' in self.featuresInFullSequence):
            currentIntronSequence += str(self.featuresInFullSequence['Intron 7'])
        if ('Intron 8' in self.featuresInFullSequence):
            currentIntronSequence += str(self.featuresInFullSequence['Intron 8'])
        if ('Intron 9' in self.featuresInFullSequence):
            currentIntronSequence += str(self.featuresInFullSequence['Intron 9'])
        #if ('3\' UTR' in self.featuresInFullSequence):
        #    currentIntronSequence += str(self.featuresInFullSequence['3\' UTR'])

        # elif('Exon 2' in self.featuresInAPDSequence and
        #   'Simulated Intron 2' in self.featuresInAPDSequence and
        #   'Exon 3' in self.featuresInAPDSequence):
        #    currentAPDSequence = (
        #        str(self.featuresInAPDSequence['Exon 2']) +
        #        str(self.featuresInAPDSequence['Simulated Intron 2']) +
        #        str(self.featuresInAPDSequence['Exon 3']))

        # else:
        #    print('I cannot construct a reliable APD sequence for this allele.' + self.alleleName)
        # for featureKey in self.featuresInAPDSequence.keys():
        #    currentAPDSequence += str(self.featuresInAPDSequence[featureKey])
        #    return ''

        return currentIntronSequence



    def copy(self):
        # What a dumb method. I'm better than this.
        copyAllele = HLA_Allele()
        copyAllele.allelePrefix = self.allelePrefix
        copyAllele.geneName = self.geneName
        copyAllele.nomenclatureGroupCount = self.nomenclatureGroupCount
        copyAllele.alleleGroup = self.alleleGroup
        copyAllele.specificProtein = self.specificProtein
        copyAllele.synonymousSubstitution = self.synonymousSubstitution
        copyAllele.noncodingChanges = self.noncodingChanges
        copyAllele.expressionSuffix = self.expressionSuffix
        copyAllele.featuresInFullSequence = self.featuresInFullSequence.copy()
        copyAllele.featuresInAPDSequence = self.featuresInAPDSequence.copy()
        copyAllele.geneFilter = self.geneFilter
        copyAllele.alleleName = self.alleleName
        #copyAllele.APDSequence = self.APDSequence
        copyAllele.sequence = self.sequence
        #copyAllele.in2Sequence = self.in2Sequence
        copyAllele.outputDirectory = self.outputDirectory
        return copyAllele

    def reverseComplement(self):     
        newAllele = self.copy()
        newAllele.sequence = seqReverseComplement(self.sequence)
        
        for featureKey in self.featuresInFullSequence.keys():
            newAllele.featuresInFullSequence[featureKey] = seqReverseComplement(self.featuresInFullSequence[featureKey])
        
        for featureKey in self.featuresInAPDSequence.keys():
            newAllele.featuresInAPDSequence[featureKey] = seqReverseComplement(self.featuresInAPDSequence[featureKey])
  
        return newAllele

    def forwardComplement(self):     
        newAllele = self.copy()
        newAllele.sequence = seqForwardComplement(self.sequence)
        
        for featureKey in self.featuresInFullSequence.keys():
            newAllele.featuresInFullSequence[featureKey] = seqForwardComplement(self.featuresInFullSequence[featureKey])
        
        for featureKey in self.featuresInAPDSequence.keys():
            newAllele.featuresInAPDSequence[featureKey] = seqForwardComplement(self.featuresInAPDSequence[featureKey])
  
        return newAllele

    def reverse(self):     
        newAllele = self.copy()
        newAllele.sequence = seqReverse(self.sequence)
        
        for featureKey in self.featuresInFullSequence.keys():
            newAllele.featuresInFullSequence[featureKey] = seqReverse(self.featuresInFullSequence[featureKey])
        
        for featureKey in self.featuresInAPDSequence.keys():
            newAllele.featuresInAPDSequence[featureKey] = seqReverse(self.featuresInAPDSequence[featureKey])
  
        return newAllele


    def parseAlleleNomenclatureFiltered(self, alleleNode, geneFilter):
        self.geneFilter = geneFilter
        self.parseAlleleNomenclature(alleleNode)

    def parseAlleleNomenclature(self, alleleNode):

        self.alleleName = alleleNode.get('name')
     
        # Split the HLA allele apart.  See this as a reference:
        # http://hla.alleles.org/nomenclature/naming.html
        alleleSplit = str.split(self.alleleName, '*')
        
        #HLA
        self.allelePrefix  = str.split(alleleSplit[0],'-')[0]
        #A
        self.geneName      = str.split(alleleSplit[0],'-')[1]

        #01:01:01:01
        nomenclatureFields = str.split(alleleSplit[1],':')

        fieldCount = len(nomenclatureFields)
        self.nomenclatureGroupCount = fieldCount

        if (fieldCount == 0):
            print ('For allele ' + self.alleleName + ' no nomenclature fields were found.  This might be a problem.')
        if (fieldCount > 0):
            self.alleleGroup = nomenclatureFields[0]
        if (fieldCount > 1):
            self.specificProtein = nomenclatureFields[1]
        if (fieldCount > 2):
            self.synonymousSubstitution = nomenclatureFields[2]
        if (fieldCount > 3):
            #Testing for suffix at the end of an allele:
            #Is the last character a digit?
            lastToken = nomenclatureFields[3] 
            if (isNumber(lastToken[ len(lastToken)-1 : len(lastToken) ])):
                self.noncodingChanges = lastToken
            else:
                self.noncodingChanges = lastToken[ 0 : len(lastToken)-1 ]
                self.expressionSuffix = lastToken[ len(lastToken)-1 : len(lastToken) ]               
        if (fieldCount > 4):
            print ('For allele ' + self.alleleName + 
                ' I found more than 4 allele groups.  This seems like a nomenclature problem.')
        
        # each 'allele' node should have a 'sequence' underneath it.  I think just one.
        sequenceList = alleleNode.findall('{http://hla.alleles.org/xml}sequence')
        
        # Is it really just one?
        if (len(sequenceList) != 1):
            # In the case of deleted alleles, the allele still has an ID, but no sequence node.
            if (len(sequenceList) == 0):
                print( 'No sequences found for this allele name. This is OK if the allele has been deleted:' 
                    + str(self.alleleName))
            else:
                print( 'Error: More than one sequence node found, thats strange:' + str(len(sequenceList)))
                raise Exception('More than one sequence node found for an allele.' + str(self.alleleName))

        else:            
            for sequenceNode in sequenceList:
                nucSequenceList = sequenceNode.findall('{http://hla.alleles.org/xml}nucsequence')

                # Is it really just one nucleotide sequence?
                if (len(nucSequenceList) != 1):
                    print ('Error: More than one nucsequence node found: '  + str(len(nucSequenceList)))
                    raise Exception('Error: More than one nucsequence node on ' + str(self.alleleName) + ' found: '  + str(len(nucSequenceList)))
                else:
                    self.sequence = nucSequenceList[0].text
                    
                    # This is simple.  Don't filter on gene name, no point to it really.
                    if (len(self.geneFilter)==0 or self.geneName in self.geneFilter):

                        # Calculate APD sequence, which is the sequence including exon 2, intron 2, exon 3.
                        # The Antigen Presenting Domain 
                        featureList = sequenceNode.findall('{http://hla.alleles.org/xml}feature')

                        print('Allele:' + self.alleleName)

                        print('Features found: ' + str(len(featureList)))


                        for featureNode in featureList:                        
                            featureName = featureNode.get('name')
                            featureType = featureNode.get('featuretype')

                            # Record what features are present
                            # This excludes the translated protein feature, we just want UTRs,Introns,Exons
                            if (featureType in ['Intron','Exon','UTR']):
                                #featureShortName = (featureName
                                #    .replace('Exon ', 'X')
                                #    .replace('Intron ', 'I')
                                #   .replace('5\' UTR', '5UTR')
                                #    .replace('3\' UTR', '3UTR'))
                                #self.featuresInFullSequence += featureShortName + ' '

                                coordinate = featureNode.findall('{http://hla.alleles.org/xml}SequenceCoordinates')[0]
                                
                                # Subtract one from the start index because 
                                # IMGT XML uses 1-based indexing
                                # Python string indexing uses 0 based.        
                                featureStart = int(coordinate.get('start')) - 1
                                featureEnd = int(coordinate.get('end'))
                                
                                print('Storing this feature:' + str(featureName) + ' which is located here:(' + str(featureStart) + ':' + str(featureEnd) + ')')
                              
                                
                                featureSequence = self.sequence[featureStart:featureEnd]
                                
                                print('sequence=' + featureSequence)
                                
                                self.featuresInFullSequence[featureName] = featureSequence
                     
                  
                        if('Exon 2' in self.featuresInFullSequence and
                           'Intron 2' in self.featuresInFullSequence and
                           'Exon 3' in self.featuresInFullSequence):
                     
                            print('Ex2-In2-Ex3 Found!')
        
                            self.featuresInAPDSequence['Exon 2']   = self.featuresInFullSequence['Exon 2']
                            self.featuresInAPDSequence['Intron 2'] = self.featuresInFullSequence['Intron 2']
                            self.featuresInAPDSequence['Exon 3']   = self.featuresInFullSequence['Exon 3']
                            
                            print('Lengths Ex2:In2:Ex3:Total ' + str(len(self.featuresInAPDSequence['Exon 2'])) 
                                  + ':' + str(len(self.featuresInAPDSequence['Intron 2'])) 
                                  + ':' + str(len(self.featuresInAPDSequence['Exon 3']))
                                  + ':' + str(len(self.APDSequence()))) 

                        #elif (ex2Start != 0 and ex2End != 0 and 
                        #    ex3Start != 0 and ex3End != 0): 
                        elif('Exon 2' in self.featuresInFullSequence and
                           'Exon 3' in self.featuresInFullSequence):
                            
                            print('Ex2 and Ex3 found for allele ' + self.alleleName + ', but no intron between them. This sequence does not have a valid Antigen Presentation sequence.')
                            
                            # I am removing the logic to use a generated intron 2.
                            #print('Ex2 and Ex3 found, but no intron between them.  I will replace In2 with a simulated intron.')
                            #intronSequence = self.loadConsensusIntron2Sequence()
                            #self.featuresInAPDSequence['Exon 2']   = self.featuresInFullSequence['Exon 2']
                            #self.featuresInAPDSequence['Simulated Intron 2'] = intronSequence
                            #self.featuresInAPDSequence['Exon 3']   = self.featuresInFullSequence['Exon 3']
                            #self.featuresInAPDSequence = 'X2 SIM-I2 X3'
                            #print('Lengths Ex2:Sim-In2:Ex3:Total ' + str(len(self.featuresInAPDSequence['Exon 2'])) 
                                #+ ':' + str(len(self.featuresInAPDSequence['Simulated Intron 2'])) 
                                #+ ':' + str(len(self.featuresInAPDSequence['Exon 3']))
                                #+ ':' + str(len(self.APDSequence()))) 


                        else:
                            # TODO  Add some more logic in here.  What can be done?
                            #print('Problem with these coordinates:' + str(ex2Start) + ',' + str(ex2End) + ',' + 
                            #    str(in2Start) + ',' + str(in2End) + ',' + str(ex3Start) + ',' + str(ex3End))
                            #raise Exception('Cannot find Exons 2 and 3 for allele ' + self.alleleName) 
                            print('Cannot find Ex2-In2-Ex3 for allele ' + self.alleleName + 'I give up.') 
                            #self.APDSequence = self.sequence
                            #self.featuresInAPDSequence = self.featuresInFullSequence
                            #print('Lengths Total :' + str(len(self.APDSequence()))) 


                    # Not filtering on gene name anymore.
                    else:
                        print('Not parsing the feature information, because this gene isn\'t in the gene filter.')

    def loadConsensusIntron2Sequence(self):
        #print('Loading an Intron 2 consensus.')
        
        intronSequence = ''
        
        groupwiseIntron2ReferenceFilename = os.path.join('Intron2References', 'HLA_' + self.geneName + '_' + self.alleleGroup + '.consensus.fasta')
        genewiseIntron2ReferenceFilename = os.path.join('Intron2References', 'HLA_' + self.geneName + '.consensus.fasta')

        groupwiseIntron2ReferenceFilenameFull = os.path.join(self.outputDirectory, groupwiseIntron2ReferenceFilename)
        genewiseIntron2ReferenceFilenameFull = os.path.join(self.outputDirectory, genewiseIntron2ReferenceFilename)

        #print('Groupwise Intron 2 Reference Filename:' + groupwiseIntron2ReferenceFilenameFull)
        #print('Genewise Intron 2 Reference Filename:' + genewiseIntron2ReferenceFilenameFull)
        
        #If the groupwise file exists, use that.
        if (os.path.isfile(groupwiseIntron2ReferenceFilenameFull)): 
            print('Groupwise intron 2 reference found.  Using: ' + groupwiseIntron2ReferenceFilename)
            intronSequence = readFirstSeqFromFasta(groupwiseIntron2ReferenceFilenameFull)
        #Otherwise, use the genewise reference.  It's a little less specific than the groupwise consensus
        elif(os.path.isfile(genewiseIntron2ReferenceFilenameFull)):
            print('Groupwise intron 2 reference is missing, using the genewise reference: ' + genewiseIntron2ReferenceFilename)
            intronSequence = readFirstSeqFromFasta(genewiseIntron2ReferenceFilenameFull)
        #I don't know what to do here, I'll use a bunch of ambiguous nucleotides.  Who has a better idea?
        else:
            print('Suitable reference was not found.  I will fill with ambiguous nucleotides.')
            intronSequence = ('NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN' + 
                            'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN' + 
                            'NNNNNNNNNNNNNNNNNNNNNNNN')

        return intronSequence 
