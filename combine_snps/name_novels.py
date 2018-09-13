from os.path import splitext, isdir, split, join
from os import makedirs, system, listdir

from Bio.Seq import Seq
from Bio.SeqIO import write, parse


def nameNovels(imgtFileName, draFileName, outputDirectory):  

    print ('I am naming my novel alleles.')

    # Define Allele Call Association File
    alleleCallOutputFile = createOutputFile(join(outputDirectory, 'SequenceAlleleCalls.csv'))
    alleleCallOutputFile.write('DRA_Sequence,Allele_Call\n')
    
    # Define Novel Sequence Output File
    novelAlleleOutputFile = createOutputFile(join(outputDirectory, 'FoundAlleles.fasta'))

    # Open IMGT files
    referenceSequences = list(parse(imgtFileName, 'fasta'))
    
    # Open DRA Sequences
    draSequences = list(parse(draFileName, 'fasta'))
    
    # For naming our novel sequences.
    novelSequenceIndex = 1
    
    # For sequence in DRA
    for draRecord in draSequences:
        
        currentDRASeq = draRecord.seq
        currentDRAID = draRecord.id
        
        print('Checking DRA sequence...:' + currentDRAID)
        
        #print ('I have this many ref sequences:' + str(len(referenceSequences)))
 
        sequenceWasAlleleCalled = False
        # For sequence in IMGT
        for referenceRecord in referenceSequences:
            
            currentRefSeq = referenceRecord.seq
            currentRefID = referenceRecord.id
        
            # Does sequence Match exactly?  Allele Call.
            if(currentRefSeq == currentDRASeq):
                sequenceWasAlleleCalled = True
                
                print('It matches: ' + currentRefID)
                
                draRecord.id = currentDRAID + '_' + currentRefID
                alleleCallOutputFile.write(str(currentDRAID) + ',' + str(currentRefID) + '\n')
                
                
                
                break
            
        # If not? Give it a name (new_1...)
        if (not sequenceWasAlleleCalled):
            novelAlleleName = 'new_' + str(novelSequenceIndex)
            novelSequenceIndex += 1
            
            print('Novel allele: ' + novelAlleleName)
            
            # "Allele Call" the novel allele.
            alleleCallOutputFile.write(str(currentDRAID) + ',' + str(novelAlleleName) + '\n')
            
            # Add this novel allele to reference sequences
            draRecord.id = novelAlleleName
            draRecord.description = ''
            referenceSequences.append(draRecord)
            
            
            
        # TODO: What snps are different than DRA 01:01:01:01

    write(referenceSequences, novelAlleleOutputFile, 'fasta')

    alleleCallOutputFile.close()
    novelAlleleOutputFile.close()
    

    
    
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not isdir(tempDir):
        print('Making Directory:' + tempDir)
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput     

