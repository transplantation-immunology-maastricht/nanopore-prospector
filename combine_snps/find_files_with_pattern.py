from os.path import splitext, isdir, split, join, relpath
from os import makedirs, system, listdir, walk

from Bio.Seq import Seq
from Bio.SeqIO import write, parse

from shutil import copyfile



# TODO: this method should actually be called "combine fasta files" or something. I want to have a method to just put files with a simliar pattern in the same directory.
def combineFastaFiles(inputDirectory, outputDirectory, pattern):
    # This method extracts a specific sequence from a larger amplicon.
    # For example, we are extracting the HLA-DRA coding sequence from our DRA amplicon.

    seqObjects = []

    print('I am inside extract sequences method.')

    for root, subdirs, files in walk(inputDirectory):
        # print('walking in this directory:' + root)
        for fileName in files:
            # print('checking this file:' + fileName)
            fullPath = join(root, fileName)
            relativePath = relpath(fullPath, inputDirectory)
            # print('full path:' + fullPath)
            # print('The relative path:' + relativePath)

            # If filename has the pattern in it

            if (pattern in relativePath):
                print('This is the file we are looking for: ' + relativePath)

                for record in parse(fullPath, "fasta"):
                    record.id = relativePath.replace('/', '_').replace('\\', '_')
                    record.description = ''

                    seqObjects.append(record)

    # Print sequences to file.
    outputFileName = join(outputDirectory, 'CombinedSequences.fasta')
    outputFile = createOutputFile(outputFileName)
    write(seqObjects, outputFile, 'fasta')
    outputFile.close()


def collectFilesWithPattern(inputDirectory, outputDirectory, pattern):
    # This method extracts a specific sequence from a larger amplicon.
    # For example, we are extracting the HLA-DRA coding sequence from our DRA amplicon.

    seqObjects = []

    makedirs(outputDirectory)

    print('I am inside extract sequences method.')

    for root, subdirs, files in walk(inputDirectory):
        # print('walking in this directory:' + root)
        for fileName in files:
            # print('checking this file:' + fileName)
            fullPath = join(root, fileName)
            relativePath = relpath(fullPath, inputDirectory)
            # print('full path:' + fullPath)
            # print('The relative path:' + relativePath)

            # If filename has the pattern in it

            if (pattern in relativePath):
                print('This is the file we are looking for: ' + relativePath)
                newpathinfo = relativePath.replace('/', '_').replace('\\', '_')
                newFileName = join(outputDirectory, fileName.replace(pattern, (newpathinfo + pattern)))
                print('new file name:' + newFileName)


                copyfile(fullPath, newFileName)
                #for record in parse(fullPath, "fasta"):
                #    record.id = relativePath.replace('/', '_').replace('\\', '_')
                #    record.description = ''

                #    seqObjects.append(record)

    # Print sequences to file.
    #outputFileName = join(outputDirectory, 'CombinedSequences.fasta')
    #outputFile = createOutputFile(outputFileName)
    #write(seqObjects, outputFile, 'fasta')
    #outputFile.close()


def createOutputFile(outputfileName):
    tempDir, tempFilename = split(outputfileName)
    if not isdir(tempDir):
        print('Making Directory:' + tempDir)
        makedirs(tempDir)
    resultsOutput = open(outputfileName, 'w')
    return resultsOutput     

