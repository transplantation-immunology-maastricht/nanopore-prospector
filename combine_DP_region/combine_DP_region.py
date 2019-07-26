#from datetime import datetime
#from os import listdir#, makedirs,
#from os.path import split, isdir, join, splitext
#from os.path import isdir, isfile, join, split

#from subprocess import Popen, PIPE, STDOUT, call
#from operator import itemgetter

from Bio.SeqIO import parse as parseSequences, write as writeSequences
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from os.path import split, join, isdir
from os import makedirs

from nanopore_prospector.common import createOutputFile

# copy2 preserves the file characteristics, I don't know if this is good or bad.
#from shutil import copy, copy2



def combineDPSequences(inputFasta, outputDir):
    print('Let us begin combining the DP sequences.')
    # Create output directory.
    if not isdir(outputDir):
        makedirs(outputDir)

    # Load unaligned.fasta
    unalignedSequences = list(parseSequences(inputFasta,format='fasta'))
    print('I loaded ' + str(len(unalignedSequences)) + ' unaligned, unreversed sequences...')
    for sequenceIndex, sequence in enumerate(unalignedSequences):
        print(str(sequence.id) + ' has length ' + str(len(sequence.seq)))

    # Reverse Complement the DPA sequence.
    print('Reverse Complementing the DPA1 sequence...')
    oldId = unalignedSequences[1].id
    unalignedSequences[1] = unalignedSequences[1].reverse_complement()
    unalignedSequences[1].id = oldId + '_revcom'
    unalignedSequences[1].description = ''
    # Write Reverse-Complement fasta
    # TODO: This will break with multiple .fasta in the file name or if the extension is different than .fasta
    inputDirPath, shortInputFileName = split(inputFasta)
    shortRevComFileName = shortInputFileName.replace('.fasta','_dparevcom.fasta')
    revcomFileName = join(outputDir, shortRevComFileName)
    writeSequences(unalignedSequences, revcomFileName, format='fasta')

    # Clustalo MSA.
    print('Starting Clustalo multiple sequence alingnment...')
    alignedFileName = revcomFileName.replace('_dparevcom.fasta','_dparevcom_aligned.fasta')
    clustalomega_cline = ClustalOmegaCommandline(infile=revcomFileName, outfile=alignedFileName, verbose=True, auto=True, threads=8, force=True, wrap=9999999)
    #print(clustalomega_cline)
    clustalomega_cline()# I don't need to execute this during testing.

    # Load aligned fasta. Store Aligned sequences for reference, A, B, and Promoter.
    alignedSequences = list(parseSequences(alignedFileName, format='fasta'))
    print('I loaded ' + str(len(alignedSequences)) + ' aligned sequences...')
    for sequenceIndex, sequence in enumerate(alignedSequences):
        print(str(sequence.id) + ' has length ' + str(len(sequence.seq)))

    alignedReferenceSequence = alignedSequences[0]
    alignedDpaRevComSequence = alignedSequences[1]
    alignedDpbSequence       = alignedSequences[2]
    alignedPromoterSequence  = alignedSequences[3]


    # Declaring a state variable to keep track of where in the alignment we are.
    # Valid options = BeforeAlignment, DPA, DPAPromoter, Promoter, DPBPromoter, DPB, AfterAlignment
    alignmentState = 'BeforeAlignment'

    consensusSequence = ''
    # Iterate base by base to build the consensus sequence.
    for baseIndex, alignedReferenceBase in enumerate(alignedReferenceSequence):
        # I don't use "elif" for this structure on purpose. Or else I skip processing the base when i transition state.
        if (alignmentState == 'BeforeAlignment'):
            # Outside/Before the alignedRegion
            # Just use dashes.
            if(alignedDpaRevComSequence[baseIndex] != '-'):
                alignmentState = 'DPA'
                print('DPA region starts at 0-based index ' + str(baseIndex))
            else:
                consensusSequence += '-'

        if (alignmentState == 'DPA'):
            # While in the DPA region:
            # Only use the DPA sequence. Never use the reference. Keep Insertions.
            if(alignedPromoterSequence[baseIndex] != '-'):
                alignmentState = 'DPAPromoter'
                print('DPAPromoter region starts at 0-based index ' + str(baseIndex))
            else:
                consensusSequence += alignedDpaRevComSequence[baseIndex]

        if (alignmentState == 'DPAPromoter'):
            # While in A/Prom overlap region
            # Check each base A vs Prom. Panic if there is a mismatch!
            if(alignedDpaRevComSequence[baseIndex] != alignedPromoterSequence[baseIndex]):
                # "mismatch" is fine if the DPA has a '-', that means we switch states.
                if(alignedDpaRevComSequence[baseIndex] == '-'):
                    alignmentState = 'Promoter'
                    print('Promoter region starts at 0-based index ' + str(baseIndex))
                else:
                    raise Exception('Dpa(' + str(alignedDpaRevComSequence[baseIndex])
                        + ') and promoter(' + str(alignedPromoterSequence[baseIndex])
                        + ') have a different base at 0-based index ' + str(baseIndex))
            else:
                # If they match, use that base
                consensusSequence += alignedDpaRevComSequence[baseIndex]

        if (alignmentState == 'Promoter'):
            # While in the Promoter Region:
            # Only use the promoter sequence, never the reference.Keep Insertions.
            if(alignedDpbSequence[baseIndex] != '-'):
                alignmentState = 'DPBPromoter'
                print('DPBPromoter region starts at 0-based index ' + str(baseIndex))
            else:
                consensusSequence += alignedPromoterSequence[baseIndex]

        if (alignmentState == 'DPBPromoter'):
            # While in Prom /DPB overlap
            # Check each base B vs Promoter. Panic if there is a mismatch.
            # If they match, use that base.
            if(alignedDpbSequence[baseIndex] != alignedPromoterSequence[baseIndex]):
                # "mismatch" is fine if the Promoter has a '-', that means we switch states.
                if(alignedPromoterSequence[baseIndex] == '-'):
                    alignmentState = 'DPB'
                    print('DPB region starts at 0-based index ' + str(baseIndex))
                else:
                    raise Exception('Dpb(' + str(alignedDpbSequence[baseIndex])
                        + ') and promoter(' + str(alignedPromoterSequence[baseIndex])
                        + ') have a different base at 0-based index ' + str(baseIndex))
            else:
                # If they match, use that base
                consensusSequence += alignedDpbSequence[baseIndex]

        if (alignmentState == 'DPB'):
            # While in the DPB region
            # Only use the DPB sequence, never the reference.Keep Insertions.
            # This also includes the "AfterAlignment" region, but i dont need to detect it. Just use the DPB sequence, which might be a "-"
            consensusSequence += alignedDpbSequence[baseIndex]


    print('The consensus sequence has a length ' + str(len(consensusSequence)))
    # Output: Consensus aligned with the Reference, A, B, Prom.
    alignedSequencesWithConsensus = alignedSequences.copy()
    alignedSequencesWithConsensus.append(SeqRecord(Seq(consensusSequence))) #cheap way to make a new sequence object
    alignedSequencesWithConsensus[4].id = 'ConsensusSequence'
    alignedSequencesWithConsensus[4].description = ''
    #alignedSequencesWithConsensus[4].seq = Seq(consensusSequence)
    alignedOutputWithConsensusFileName = shortInputFileName.replace('.fasta','_aligned_withconsensus.fasta')
    alignedOutputWithConsensusFileName = join(outputDir, alignedOutputWithConsensusFileName)
    writeSequences(alignedSequencesWithConsensus, alignedOutputWithConsensusFileName, format='fasta')

    # Output: Just Consensus sequence, without the inserted sequences etc.\
    consensusOutputFileName = shortInputFileName.replace('.fasta', '_Consensus.fasta')
    consensusOutputFileName = join(outputDir, consensusOutputFileName)
    consensusOutputFile = createOutputFile(consensusOutputFileName)
    cleanConsensus = consensusSequence.replace('-','')
    # Manual. Just to put the whole sequence on one line. I hate writing sequences in blocks.
    consensusOutputFile.write('>' + shortInputFileName.replace('.fasta', '_Consensus\n'))
    consensusOutputFile.write(cleanConsensus)
    consensusOutputFile.close()


    # Output: Just Consensus sequence, without the inserted sequences etc.\
    # But this time give it a decent name.
    consensusOutputFileName = shortInputFileName.replace('.fasta', '_ConsensusWithAlleles.fasta')
    consensusOutputFileName = join(outputDir, consensusOutputFileName)
    consensusOutputFile = createOutputFile(consensusOutputFileName)
    cleanConsensus = consensusSequence.replace('-','')
    # Manual. Just to put the whole sequence on one line. I hate writing sequences in blocks.

    dpaAlleleCall = alignedSequencesWithConsensus[1].id.replace('_revcom','')
    dpbAlleleCall = alignedSequencesWithConsensus[2].id
    promoterAlleleCall = alignedSequencesWithConsensus[3].id
    consensusOutputFile.write('>' + promoterAlleleCall + '_' + dpaAlleleCall + '_' + dpbAlleleCall + '\n')
    consensusOutputFile.write(cleanConsensus)
    consensusOutputFile.close()







