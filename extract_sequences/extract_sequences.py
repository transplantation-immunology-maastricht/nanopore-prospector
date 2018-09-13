from os.path import splitext, isdir, split, join
from os import makedirs, system

from Bio.SeqIO import parse, write

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from pysam import AlignmentFile

from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA


# TODO: Put this logic in nanopore prospector. Why not just work out of there from now on
def extractGeneSequencesFromLongerSequences(fullLengthSequenceFile, geneReferenceFastaFile):  
    # This method extracts a specific sequence from a larger amplicon.
    # For example, we are extracting the HLA-DRA coding sequence from our DRA amplicon.
    # I added functionality for fastq files, we can use this to extract the "interesting" part of reads.

    referenceSequence = parse(geneReferenceFastaFile, 'fasta')
    referenceSequenceLength = len(next(referenceSequence).seq)

    #print ('trying to split this path:' + fullLengthSequenceFile)
    # Get output directory.
    # TODO Pass this in instead of calculating it.
    outputDirPath, shortSequenceFileName = split(fullLengthSequenceFile)
    
    # Sequence file read type, is it fasta or fastq?
    sequenceType = getReadFileType(shortSequenceFileName)
    
    # Align sequences against gene reference.
   # alignReads(geneReferenceFastaFile, outputFileName, join(outputDirPath, outputFileName + '_alignment'), 4)
    almtOutputDir = join(outputDirPath, shortSequenceFileName + '_alignment')
    alignReads(geneReferenceFastaFile, fullLengthSequenceFile, almtOutputDir, 4)
   
    # Open read file and keep a list of the sequences. Dictionary?
    fullLengthSequenceReads = list(parse(fullLengthSequenceFile, sequenceType))
    shortSequenceReads = []

    rejectSequenceReads = []
        
    # Open the alignment.
    bamfile = AlignmentFile(join(almtOutputDir,'alignment.bam'), 'rb')
    
    # Open the Reference
    alignmentRef = list(parse(geneReferenceFastaFile, 'fasta'))[0]
                    
    # Loop sequences in alignment. 
    reads = bamfile.fetch()
    for alignedRead in reads:
        
        readID = alignedRead.query_name

        # Nevermind i can just pull the aligned sequence. This is slick.
        readSequence = alignedRead.query_alignment_sequence

        # Phred Scores
        if(alignedRead.query_alignment_qualities is not None):
            readQualities = alignedRead.query_alignment_qualities
        else:
            print('This aligned read has no quality scores, that means it is fasta. I think no problem. Delete this warning if it works ok')
            readQualities = None
            pass

        # Read Begin/End.
        # Reference Begin/End
        # IsSecondary
        # IsSupplementary

        isSupplementary = alignedRead.is_supplementary
        isSecondary = alignedRead.is_secondary
        isReverse = alignedRead.is_reverse
        readBeginIndex = alignedRead.query_alignment_start
        readEndIndex = alignedRead.query_alignment_end
        referenceBeginIndex = alignedRead.reference_start
        referenceEndIndex = alignedRead.reference_end


        # Adjust the sequences we already read.
        # Why do I need this loop? I think I don't.  What was ben thinking, I'm already in a loop.
        #for fullLengthSequenceRead in fullLengthSequenceReads:
        #    if fullLengthSequenceRead.id == readID:
                #print ('FOund a read match.')
                #shortRead=fullLengthSequenceRead.copy()

                # Problem - I don't have quality scores here.

        newReadID = (readID + '_REF(' + str(referenceBeginIndex) + ',' + str(referenceEndIndex)
            + ')_READ(' + str(readBeginIndex) + ',' + str(readEndIndex)
            + ')_REVERSE=' + str(isReverse)
            + ')_SECONDARY=' + str(isSecondary)
            + '_SUPPLEMENTARY=' + str(isSupplementary))

        shortRead = SeqRecord(Seq(readSequence,
            IUPACUnambiguousDNA),
            id=newReadID,
            description='',
            letter_annotations={'phred_quality':readQualities}
            )

        # Had a problem with short reads, only have partial alignments. If the read is shorter than, say
        # 90% of the reference lenght, I should reject it.
        # Note: I think these shorter reads are supplementary alignments from minimap2. It is aligning sequences in strange places
        # I can deal with that, I'm recording some extra information from the indices and properties.
        if(len(readSequence) > (.90) * referenceSequenceLength):
            shortSequenceReads.append(shortRead)
        else:
            rejectSequenceReads.append(shortRead)

    outputFileName = fullLengthSequenceFile.replace('.' + sequenceType, '.short.' + sequenceType)
    sequenceWriter = createOutputFile(outputFileName)
    write(shortSequenceReads, sequenceWriter, sequenceType)
    sequenceWriter.close()
    alignReads(geneReferenceFastaFile, outputFileName, join(outputDirPath, outputFileName+'_alignment'), 4)

    outputFileName = fullLengthSequenceFile.replace('.' + sequenceType, '.reject.' + sequenceType)
    rejectSequenceWriter = createOutputFile(outputFileName)
    write(rejectSequenceReads, rejectSequenceWriter, sequenceType)
    rejectSequenceWriter.close()
    alignReads(geneReferenceFastaFile, outputFileName, join(outputDirPath, outputFileName + '_alignment'), 4)

    
def getReadFileType(fileNameFull):
    # This method is for deciding if a file is fasta or fastq.
    # It can be extended for other file types.
    fileName, fileExtension = splitext(fileNameFull)
    
    if(fileExtension == '.fasta' or fileExtension == '.fa'):
        return 'fasta'
    elif(fileExtension == '.fastq' or fileExtension == '.fq'):
        return 'fastq'
    else:
        raise Exception('I do not understand this file type: ' + str(fileNameFull))
    

def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not isdir(tempDir):
        print('Making Directory:' + tempDir)
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput     

# Perform minimap2 Alignment.  Align all reads against the Reference.
# TODO: Wasn't this supposed to be in a common methods somewhere?  Maybe not. Check if this method is duplicated.
# 
def alignReads(referenceLocation, readFileLocation, alignmentOutputDirectory, numberThreads):
    print('\nAligning Reads.')
    #print('\nStep 1.) Aligning reads against the reference.')
    
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
    
    
    # TODO: I should mess with the alignment parameters, because sometimes sequence
    # at the 3' end is being trimmed from alignments.
    # homopolymers make an insertion/deletion, and the aligner prefers to end the alignment,
    # rather than inserting a deletion. 
       
    # Part 2 Align
    try:
        alignmentOutputName = join(alignmentOutputDirectory,'alignment.bam')

        # Parameters depend on "read type"
        # if i have a fastq, i assume these are minion reads.
        # if i have a fasta, i assume these are consensus sequences.
        if(getReadFileType(readFileLocation) == 'fasta'):
            minimapParams = '--secondary=no -ax asm5'
        elif(getReadFileType(readFileLocation) == 'fastq'):
            minimapParams = '--secondary=no -ax map-ont'
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

