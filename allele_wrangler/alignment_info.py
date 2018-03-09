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

class AlignmentInfo():  
    
    def __init__(self):
        self.sequence = ''
        self.alignmentColumns = []

class AlignmentColumn():   
    
    def __init__(self):

        self.referencePosition = 0
        self.referenceBase = ''
        self.referenceAdjustment = '?'
        self.alignedCount = 0
        self.unalignedCount = 0
        self.matchCount = 0
        self.mismatchCount = 0
        self.inCount = 0
        self.delCount = 0
        self.aCount = 0
        self.gCount = 0
        self.cCount = 0
        self.tCount = 0
            
       