#from datetime import datetime
#from os import listdir#, makedirs,
#from os.path import split, isdir, join, splitext
#from os.path import isdir, isfile, join, split

#from subprocess import Popen, PIPE, STDOUT, call
#from operator import itemgetter

from Bio.SeqIO import parse as parseSequences, write as writeSequences
from Bio.Align.Applications import ClustalOmegaCommandline

from os.path import split, join, isdir
from os import makedirs

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
    clustalomega_cline = ClustalOmegaCommandline(infile=revcomFileName, outfile=alignedFileName, verbose=True, auto=True, threads=8, force=True)
    #print(clustalomega_cline)
    clustalomega_cline()# I don't need to execute this during testing.

    # Load aligned fasta. Store Aligned sequences for reference, A, B, and Promoter.
    alignedSequences = list(parseSequences(inputFasta, format='fasta'))
    print('I loaded ' + str(len(alignedSequences)) + ' aligned sequences...')
    for sequenceIndex, sequence in enumerate(alignedSequences):
        print(str(sequence.id) + ' has length ' + str(len(sequence.seq)))

    # Iterate base by base.
    # Outside/Before the alignedRegion
    # Just use dashes.
    # While in the DPA region:
    # Only use the DPA sequence. Never use the reference. Keep Insertions.
    # While in A/Prom overlap region
    # Check each base A vs Prom. Panic if there is a mismatch!
    # If they match, use that base.
    # While in the Promoter Region:
    # Only use the promoter sequence, never the reference.Keep Insertions.
    # While in Prom /DPB overlap
    # Check each base B vs Promoter. Panic if there is a mismatch.
    # If they match, use that base.
    # While in the DPB region
    # Only use the DPB sequence, never the reference.Keep Insertions.
    # Outside the aligned region:
    # Just use dashes.
    # Output: Consensus aligned with the Reference, A, B, Prom.
    # Output: Just Consensus sequence, without the inserted sequences etc.




