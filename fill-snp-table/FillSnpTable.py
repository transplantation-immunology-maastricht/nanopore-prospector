import csv
from os.path import split, isdir

# This method is a directory-safe way to open up a write file.
def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not isdir(tempDir):
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput

# I missed a SNP, so with these inputs I can re-run the fill_snp program.
rawSNPDataFile='/home/minion/MinIONData/2018.DRA.Basecalling/Combined_SNP_Arrays/CombinedSnps.1463.csv'
SNPOutputFile='/home/minion/MinIONData/2018.DRA.Basecalling/Combined_SNP_Arrays/CombinedSnps.1463.clean.csv'

# Input values, these should be commandline parameters.
# RawSNPs comes from nanopore prospector.
#rawSNPDataFile='/home/minion/MinIONData/2018.DRA.Basecalling/AllelesAgainstIMGTReference/RawDRASnps.csv'
#SNPOutputFile='/home/minion/MinIONData/2018.DRA.Basecalling/AllelesAgainstIMGTReference/RawDRASnpsClean.csv'

# This is how i run this program to find the "normal" imgt snps for my table.
#rawSNPDataFile='/home/minion/MinIONData/2018.DRA.Basecalling/AllelesAgainstIMGTReference/SNP_Analysis/All_ReadQualityAnalysis/AlleleSpecificPolymorphisms.csv'
#SNPOutputFile='/home/minion/MinIONData/2018.DRA.Basecalling/AllelesAgainstIMGTReference/DRAReferenceSNPS.csv'


# Load file, it's a csv.
with open(rawSNPDataFile) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0

    # A list for each Allele, containing an array of SNPs
    correctedSNPs = []

    for row in csv_reader:

        #Remember column 0 is the allele name.

        # First line is the indexes of the SNPs
        print('line' + str(line_count) + str(row))
        if line_count == 0:
            SNPIndexes = row

        # 2nd line is the reference Base.
        elif line_count == 1:
            referenceBases = row
        else:

            # I don't actually need to store the allele name.
            alleleName = row[0]

            # Loop each column, if the entry is blank, set it to the reference base.
            for columnIndex in range(1,len(row)):
                if(len(row[columnIndex]) == 0):
                    row[columnIndex] = referenceBases[columnIndex]

            # Store the new allele with the filled SNPs.
            correctedSNPs.append(row)

            print('new row:'+ str(row))

        line_count += 1
    print(f'Processed {line_count} lines.')






# Write my new csv to the output file.
    csvOutputFile = createOutputFile(SNPOutputFile)

    csvOutputFile.write(','.join(SNPIndexes) +  '\n')

    csvOutputFile.write(','.join(referenceBases) +  '\n')

    # Loop each allele, print table of SNPS
    for alleleIndex in range(0, len(correctedSNPs)):
        csvOutputFile.write(','.join(correctedSNPs[alleleIndex]) +  '\n')

    csvOutputFile.close()

