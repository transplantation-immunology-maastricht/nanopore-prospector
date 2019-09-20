# # This file is part of nanopore_prospector.
# #
# # nanopore_prospector is free software: you can redistribute it and/or modify
# # it under the terms of the GNU Lesser General Public License as published by
# # the Free Software Foundation, either version 3 of the License, or
# # (at your option) any later version.
# #
# # nanopore_prospector is distributed in the hope that it will be useful,
# # but WITHOUT ANY WARRANTY; without even the implied warranty of
# # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# # GNU Lesser General Public License for more details.
# #
# # You should have received a copy of the GNU Lesser General Public License
# # along with nanopore_prospector. If not, see <http://www.gnu.org/licenses/>.
#
#
# import sys
# from getopt import getopt, GetoptError
# from os.path import exists, isfile, isdir
# from os import makedirs
#
# from punkin_chunker.punkin_chunker import sortDirectory
#
#
# SoftwareVersion = "punkin-chunker Version 1.0"
#
# # TODO Fix usage
# def usage():
#     print("usage:\n" +
#     "\tThis script is written for python 2.7.11\n" +
#     "\tDo this instruction because it is wrong right now.\n" +
#     "\t\tOR\n" +
#     "\t\tSubdirectories containing .fast5 reads (Barcoded Reads, ex. /BCO1)\n\n" +
#     "\tThe output directory will be filled with .fasta and .fastq reads, sorted by subfolders.\n\n" +
#
#     "\tOptions:\n" +
#     "\t-i\t--idir   \tInput Directory (required)\n" +
#     "\t-o\t--odir   \tOutput Directory (required)\n" +
#     "\t-m\t--minlen \tMinimum Read Length Filter\n" +
#     "\t-M\t--maxlen \tMaximum Read Length Filter\n" +
#     "\t-h\t--help   \tPrint this message\n" +
#     "\t-r\t--rundate\tSequencing Date or other information which will be included in the extract filename.\n" +
#
#     "\n\tSee README.MD for instructions on how to set up an anaconda environment for this script\n"
#     )
#     # Usage: you must provide a read input and output directory.
#     # reads are fastq or a directory with fastq.
#     # The rest are optional.
#
# # Read Commandline Arguments.  Return true if everything looks okay for read extraction.
# def readArgs():
#     # Default to None.  So I can easily check if they were not passed in.
#
#     global readInput
#     global outputResultDirectory
#     global sampleID
#     global sortReference
#     global threadCount
#
#     readInput                = None
#     outputResultDirectory    = None
#     sampleID                 = None
#     sortReference            = None
#     threadCount              = None
#
#     # https://www.tutorialspoint.com/python/python_command_line_arguments.htm
#     try:
#         opts, args = getopt(sys.argv[1:]
#             ,"hvo:r:R:s:t:"
#             ,["help", "version", "outputdir=", "reads=", "reference=", "sampleid=", "threads="])
#
#         for opt, arg in opts:
#
#             if opt in ('-h', '--help'):
#                 print (SoftwareVersion)
#                 usage()
#                 return False
#
#             elif opt in ('-v', '--version'):
#                 print (SoftwareVersion)
#                 return False
#
#             elif opt in ("-o", "--outputdir"):
#                 outputResultDirectory = arg
#             elif opt in ("-r", "--reads"):
#                 readInput = arg
#
#             elif opt in ("-R", "--reference"):
#                 sortReference = arg
#
#             elif opt in ("-s", "--sampleid"):
#                 sampleID = arg
#
#             elif opt in ("-t", "--threads"):
#                 threadCount = arg
#
#             else:
#                 print('Unknown Commandline Option:' + str(opt) + ':' + str(arg))
#                 raise Exception('Unknown Commandline Option:' + str(opt) + ':' + str(arg))
#
#         if(len(sys.argv) < 3):
#             print ('I don\'t think you have enough arguments.\n')
#             usage()
#             return False
#
#     except GetoptError, errorMessage:
#         print ('Something seems wrong with your commandline parameters.')
#         print (errorMessage)
#         usage()
#         return False
#
#
#     # Sanity Checks
#
#     # default 1 thread
#     if threadCount is None:
#         threadCount = 1
#
#
#     # This output directory should exist
#     if not exists(outputResultDirectory):
#         makedirs(outputResultDirectory)
#
#     # TODO: trim out spaces and special characters from the sample id.
#     # Check on all the values.
#
#
#     return True
#
#
# if __name__=='__main__':
#     try:
#         if(readArgs()):
#             print ('args look good, lets sort the reads.')
#
#             if (isfile(readInput)):
#                 print ('Read input is a file that exists.')
#             elif (isdir(readInput)):
#                 print ('Read input is a directory that exists.')
#                 sortDirectory(readInput, outputResultDirectory, sortReference, threadCount)
#             else :
#                 print ('I don\'t understand the read input specified, it is not a file or directory:' + readInput)
#                 raise Exception('Bad Read Input Format')
#
#         #BlastMinionReadsAgainstGroupwiseReference()
#         #BlastMinionReadsAgainstAPDRef()
#         print('Done.  Yay.')
#
#     except Exception:
#         # Top Level exception handling like a pro.
#         # This is not really doing anything.
#         print 'Unexpected problem during execution:'
#         print sys.exc_info()[1]
#         raise
