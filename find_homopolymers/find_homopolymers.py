from os.path import join, basename, split, isdir
from os import makedirs

from shutil import copyfile
#from os.path import splitext, isdir, split, join
#from os import makedirs, system

from Bio.SeqIO import parse, write

#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord

#from pysam import AlignmentFile

#from Bio.Alphabet.IUPAC import IUPACUnambiguousDNA




def findHomopolymers(inputFile, outputDirectory):
    print('Looking for homopolymers in ' + str(inputFile))

    # Open fasta file
    sequenceRecords = list(parse(inputFile, 'fasta'))
    print('I found ' + str(len(sequenceRecords)) + ' records in ' + inputFile)

    # Create output files.
    inputFileShortName = basename(inputFile)

    # Name, Position, Base, Length, Previous Base, Next Base
    homopolymerDetailOutputFile = createOutputFile(
        join(outputDirectory,
        inputFileShortName.replace('.fasta', 'HomopolymerDetails.csv')))
    homopolymerDetailOutputFile.write('AlleleName,HP_Position,HP_Base,HP_Length,Prev_Base,Next_Base\n')

    # Allele Name, HP_Count, FeatureLength, Homopolymer Density, 3_count,4_count,5_count,6_count,7_count,8plus_count
    homopolymerSummaryOutputFile = createOutputFile(
        join(outputDirectory,
        inputFileShortName.replace('.fasta', 'HomopolymerSummary.csv')))
    homopolymerSummaryOutputFile.write('AlleleName,HP_Count,Feature_Length,Density,3_count,4_count,5_count,6_count,7_count,8plus_count\n')

    # Copy the file to the output directory.
    sequenceFastaOutputFileName = join(outputDirectory, inputFileShortName)
    #write(sequenceRecords, createOutputFile(sequenceFastaOutputFileName),'fasta')
    copyfile(inputFile, sequenceFastaOutputFileName)

    # loop entries in fasta file.
    for currentSequenceRecord in sequenceRecords:
        alleleName = currentSequenceRecord.name
        sequence = currentSequenceRecord.seq
        homopolymerRecords = []

        print('Searching for homopolymers in ' + str(alleleName))

        # loop through nucleotide sequence. Probably a manual while loop so i can mess with index
        nucleotideIndex = 0
        while nucleotideIndex < len(sequence):
            currentBase = str(sequence[nucleotideIndex])
            #print('Checking base:' + str(nucleotideIndex) + ':' + str(sequence[nucleotideIndex]))

            # Detect homopolymers to the right of this base.  Store them.
            homopolymerIndex = 1
            # Walk forward, checking for homopolymers.
            # Loop condition so I don't walk off the end of the sequence
            while (nucleotideIndex + homopolymerIndex) < len(sequence):
                if(currentBase == str(sequence[nucleotideIndex + homopolymerIndex])):
                    homopolymerIndex = homopolymerIndex + 1
                else:
                    break

            # If my math is right, homopolymerIndex is currently equal to the homopolymer length.
            # 3 could be tuned here, TODO pass minimum homopolymer size as a commandline parameter.
            if(homopolymerIndex > 2):
                previousBase = str(sequence[nucleotideIndex - 1]) if (nucleotideIndex > 0) else '?'
                nextBase = str(sequence[nucleotideIndex + homopolymerIndex]) if ((nucleotideIndex + homopolymerIndex) < len(sequence)) else '?'

                # Position, Base, Length, Previous Base, Next Base
                # Should I store 1 or 0 based indices? I think 1 based in this case.
                homopolymerRecord = [(nucleotideIndex + 1), currentBase, homopolymerIndex, previousBase, nextBase]
                homopolymerRecords.append(homopolymerRecord)

            # Skip bases, i don't need to check the other bases within the same homopolymer
            # homopolymerIndex will always be at least 1, so I can just use this as my loop iterator condition. Slick.
            nucleotideIndex = nucleotideIndex + homopolymerIndex

        # Specific data points for summaries and drawing charts.
        count3 = 0
        count4 = 0
        count5 = 0
        count6 = 0
        count7 = 0
        count8plus = 0

        # loop through the stored homopolymer informations.
        for homopolymerRecord in homopolymerRecords:
            # Print entries to detail file.
            # Name, Position, Base, Length, Previous Base, Next Base
            homopolymerDetailOutputFile.write(str(alleleName)
                + ',' + str(homopolymerRecord[0])
                + ',' + str(homopolymerRecord[1])
                + ',' + str(homopolymerRecord[2])
                + ',' + str(homopolymerRecord[3])
                + ',' + str(homopolymerRecord[4])
                + '\n')


            if(int(homopolymerRecord[2]) == 3):
                count3 = count3 + 1
            elif(int(homopolymerRecord[2]) == 4):
                count4 = count4 + 1
            elif (int(homopolymerRecord[2]) == 5):
                count5 = count5 + 1
            elif (int(homopolymerRecord[2]) == 6):
                count6 = count6 + 1
            elif (int(homopolymerRecord[2]) == 7):
                count7 = count7 + 1
            elif (int(homopolymerRecord[2]) > 8):
                count8plus = count8plus + 1



        # Print an entry in the summary file.
        # Allele Name, HP_Count, FeatureLength, Homopolymer Density, 3_count, 4_count,5_count,6_count,7_count,8plus_count
        hpDensity = len(homopolymerRecords) / len(sequence)
        homopolymerSummaryOutputFile.write(alleleName
            + ',' + str(len(homopolymerRecords))
            + ',' + str(len(sequence))
            + ',' + str(hpDensity)
            + ',' + str(count3)
            + ',' + str(count4)
            + ',' + str(count5)
            + ',' + str(count6)
            + ',' + str(count7)
            + ',' + str(count8plus)
            + '\n')
    # close files
    homopolymerDetailOutputFile.close()
    homopolymerSummaryOutputFile.close()


def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not isdir(tempDir):
        print('Making Directory:' + tempDir)
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput