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


class blastResult:
    def __init__(self):
        self.NomenclatureFields = []
        self.BlastScore = 0
        self.Gene = None
        self.ForwardMatch = True
        self.readRecord = None
    
    def processBlastRecord(self, blastRecord):
        #self.readRecord = blastRecord
        blastRecordFrame = blastRecord.alignments[0].hsps[0].frame
        self.ForwardMatch = (blastRecordFrame[0] == blastRecordFrame[1])
        self.parseNomenclatureLine(blastRecord.descriptions[0].title)
        self.BlastScore = blastRecord.descriptions[0].score
        #self.ForwardMatch = matchIsForwardDirection
    
    def parseNomenclatureLine(self, alleleNameLine):
        # "gnl|BL_ORD_ID|10 HLA-A*31:01:02:01 (X2 I2 X3)"
        #print('Parsing Nomenclature Line:' + alleleNameLine)
        nomenclatureFields = []
        
        alleleNameToken = alleleNameLine.split()[1]
        #print('Allele Name:' + alleleNameToken)
        
        hyphenLocation = alleleNameToken.find('-')
        asteriskLocation = alleleNameToken.find('*')
        
        #print('The Hyphen is at position:' + str(hyphenLocation))        
        #print('The Asterisk is at position:' + str(asteriskLocation))
        
        self.Gene = alleleNameToken[hyphenLocation+1:asteriskLocation]
        #print('Gene Found:' + self.Gene)
        
        nomenclatureFieldString = alleleNameToken[asteriskLocation+1:]        
        #print('Nomenclature String:' + nomenclatureFieldString)
        
        nomenclatureFieldTokens = nomenclatureFieldString.split(':')     
        for nomenclatureFieldToken in nomenclatureFieldTokens:
            nomenclatureFields.append(nomenclatureFieldToken)
        
        self.NomenclatureFields = nomenclatureFields
 
