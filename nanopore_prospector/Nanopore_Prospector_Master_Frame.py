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

import sys
from datetime import datetime

from os import makedirs, listdir, pardir
from os.path import exists, split, isdir, isfile, join, abspath, basename, normpath

from numpy import mean, amin, amax

from Bio.SeqIO import parse as parseReads, write

from tkinter import Tk, filedialog, messagebox, Frame, StringVar, Label, Button, Entry, Scrollbar, Text, Toplevel
from tkinter.constants import HORIZONTAL, BOTH, W, BOTTOM, X, Y, RIGHT, NONE, DISABLED, END, NORMAL

# TODO: Put these in try/catch statements, because I will need this to run within pyinstaller. For example:
#Try something like this to fix the imports:
#try:
#    from AlleleSubMainGui import AlleleSubMainGui
#    from AlleleSubCommon import loadConfigurationFile
#except:
#    from saddlebags.AlleleSubMainGui import AlleleSubMainGui
#    from saddlebags.AlleleSubCommon import loadConfigurationFile

from nanopore_prospector.analysis_options_input_form import AnalysisOptionsInputForm
from punkin_chunker.punkin_chunker import sortDirectory
from nit_picker.nit_picker import prepareReads
from allele_wrangler.allele_wrangler import AlleleWrangler
from nanopore_prospector.common import logMessageToFile, writeConfigurationFile, loadConfigurationFile, getConfigurationValue, assignConfigurationValue, getBlastSortResourceLocation, getBarcodeFilePath, combineBlastReferences, createOutputFile


# TODO: Split the GUI stuff from the Analysis stuff.

class NanoporeProspectorMasterFrame(Frame):
    def __init__(self, root):

        Frame.__init__(self, root)
        root.title("Nanopore Prospector")

        self.parent = root        
        self.initialize()

    # Initialize GUI elements
    def initialize(self):
        # Load the configuration file.
        loadConfigurationFile()
        
        # This is the directory the python executable is running from.
        FileAndPath = abspath(__file__)
        self.idir, self.ifile = split(FileAndPath)

        # GUI options
        #self.label_opt = {'fill': Tkconstants.BOTH, 'padx': 10, 'pady': 10}
        self.button_opt = {'fill': BOTH, 'padx': 50, 'pady': 15}
        self.dir_opt = {'initialdir': self.idir, 'mustexist': True, 
            'parent': self, 'title': 'Choose a directory'}
        
        # A frame for choosing the instructions
        self.instructionsFrame = self.makeDemultiplexInstructionsFrame()
        self.instructionsFrame.pack()
        
        # Frame forinput/output directories
        self.chooseDirectoryFrame = self.makeChooseIODirectoriesFrame()
        self.chooseDirectoryFrame.pack()

        # Frame for Analysis Options
        self.analysisButtonsFrame = self.makeAnalysisButtonsFrame()
        self.analysisButtonsFrame.pack()
        
        # Frame for Analysis Log
        self.logFrame = self.makeAnalysisLogFrame()
        self.logFrame.pack()
        
        # Frame for the More Information Panel
        self.moreInfoFrame = self.makeMoreInfoFrame()
        self.moreInfoFrame.pack()
        
        # Assemble all the frames and get ready.
        self.pack()
        
    def recordAnalysisStep(self, currentMessage):
        # The log keeps track of step-by-step things performed by Nanopore Prospector.
        fullLogMessage = str(datetime.now()) + ' : ' + currentMessage + '\n'
        
        # Append the GUI with the Log Entry 
        # TODO: Move to the bottom of the log. The GUI is not updating in real time.
        # Alternatively, reverse the log and put new stuff on top
        self.submOutputGuiObject.insert(END, fullLogMessage) 
        
        logMessageToFile(currentMessage)
            
        # This repaints the window, we can see what's going on.
        self.update()  
        
    def reportResults(self, currentMessage):
        # This Report is intended to Summarize the prospector / allele calling results.
        fullReportMessage = currentMessage + '\n'
             
        if (getConfigurationValue('report_file_location') is None):
            print ('The report file name does not exist yet.\nTrying to report this message:\n' + currentMessage)
        else:
            #print ('Logging message: ' + currentMessage) 
            # Append the log with the log entry            
            resultsOutput = open(getConfigurationValue('report_file_location'), 'a')
            resultsOutput.write(fullReportMessage)
            resultsOutput.close()
    
    def makeDemultiplexInstructionsFrame(self):
        instructionsFrame = Frame(self)

        self.instructionText = StringVar()
        self.instructionText.set('Use this interface to look at a group of \nMinION Reads containing HLA amplicon sequences:')
        Label(instructionsFrame, width=80, height=5, textvariable=self.instructionText).pack()
        return instructionsFrame
      
    def makeChooseIODirectoriesFrame(self):
        chooseDirectoryFrame = Frame(self)        

        self.chooseInputButton = Button(chooseDirectoryFrame, text='Choose Input Directory ', command=self.chooseReadInputDirectory)
        self.chooseInputButton.grid(row=0, column=0, sticky=W)
        self.chooseOutputButton = Button(chooseDirectoryFrame, text='Choose Output Directory', command=self.chooseReadOutputDirectory)
        self.chooseOutputButton.grid(row=1, column=0, sticky=W)
        
        self.inputDirectoryText = StringVar()
        self.inputDirectoryText.set('Where is your MinION Read Directory?')
        Entry(chooseDirectoryFrame, width=60, textvariable=self.inputDirectoryText).grid(row=0, column=1)
                
        self.outputDirectoryText = StringVar()
        self.outputDirectoryText.set('Output Directory')
        Entry(chooseDirectoryFrame, width=60, textvariable=self.outputDirectoryText).grid(row=1, column=1)
        
        return chooseDirectoryFrame
        
    def makeAnalysisButtonsFrame(self):
        analysisButtonsFrame = Frame(self) 
        self.analysisOptionsButton = Button(analysisButtonsFrame, text='Analysis Options', command=self.specifyReadStatOptions)
        self.analysisOptionsButton.grid(row=0, column=0)
        self.prepareReadsButton = Button(analysisButtonsFrame, text='Demultiplex + Prepare Reads', command=self.constructInitialReadStats)
        self.prepareReadsButton.grid(row=1, column=0)
        self.runAnalysisButton = Button(analysisButtonsFrame, text='Run Full Analysis', command=self.runFullAnalysis)
        self.runAnalysisButton.grid(row=2, column=0)
        return analysisButtonsFrame
        
    def makeAnalysisLogFrame(self):
        logFrame = Frame(self)
        
        # TODO: Did i set this log label up correctly?
        self.logLocationText = StringVar()
        self.logLocationText.set('Choose a read directory to start logging...')
        Label(logFrame, width=80, height=5, textvariable=self.logLocationText).pack()
       

        self.submOutputXScrollbar = Scrollbar(logFrame, orient=HORIZONTAL)
        self.submOutputXScrollbar.pack(side=BOTTOM, fill=X)

        self.submOutputYScrollbar = Scrollbar(logFrame)
        self.submOutputYScrollbar.pack(side=RIGHT, fill=Y)

        self.submOutputGuiObject = Text(
            logFrame, width=100, height=10, wrap=NONE
            , xscrollcommand=self.submOutputXScrollbar.set
            , yscrollcommand=self.submOutputYScrollbar.set
        )

        self.submOutputXScrollbar.config(command=self.submOutputGuiObject.xview)
        self.submOutputYScrollbar.config(command=self.submOutputGuiObject.yview) 

        self.submOutputGuiObject.pack() 
        
        return logFrame
    
    def makeMoreInfoFrame(self):    
        moreInfoFrame = Frame(self)  
        self.howToUseButton = Button(moreInfoFrame, text='How to use this tool', command=self.howToUse)
        self.howToUseButton.grid(row=0, column=0)
        self.citationButton = Button(moreInfoFrame, text='Contacting or Citing MUMC', command=self.contactInformation)
        self.citationButton.grid(row=0, column=1)
        self.saddlebagsButton = Button(moreInfoFrame, text='SaddleBags - A (Novel) Allele Submission Tool', command=self.launchSaddleBags)
        self.saddlebagsButton.grid(row=0, column=2)
        return moreInfoFrame
     
    def launchSaddleBags(self):
        #messagebox.showinfo('Need to open allelesub tool','Popup a window containing the allele submission tool please.')
        # TODO: This works but the interface is messed up.
        # I think saddlebags is using "self" to assign variables when it shouldnt. They should assign to the frame instead?
        # Problem is, i should be making a toplevel object.  I can fix this.
        # See how i did it in specifyReadStatOptions
        # TODO This actually doesn't work at all, I renamed the saddlebags classes.
        # This should be a toplevel object, not a Tk object.  I think.
        #saddleBagsRoot = Tk()
        #AlleleGui(saddleBagsRoot).pack()
        #saddleBagsRoot.mainloop()
        raise Exception('Implement Saddlebags Interface, this code does not work right now.')
 
        
    def howToUse(self):
        # TODO: This howToUse method is pretty weak.    
        messagebox.showinfo('Select a directory containing reads.\n'
            + 'Do some analysis on the reads.'
            + 'Ben:Fill in better how-to-use instructions.'
            )
        
    def contactInformation(self):
        # This method should list contact information for MUMC, and a link to the github page.  
        messagebox.showinfo('Contact Information',
            'This software was created at\n'
            + 'Maastricht University Medical Center\n'
            + 'Transplantation Immunology\n'
            + 'Tissue Typing Laboratory.\n'
            + 'by Ben Matern:\n'
            + 'ben.matern@mumc.nl\n\n'
            
            + 'Please send Ben your bioinformatics\n'
            + 'and data related questions.\n\n'
            
            + 'all other inquiries can be directed\n'
            + 'to Marcel Tilanus:\n'
            + 'm.tilanus@mumc.nl\n\n'

            )

  
    def specifyReadStatOptions(self):
        print('Specifying ReadStat Options.....')
    
        self.disableGUI()
        
        readStatOptionsRoot = Toplevel()
        readStatOptionsRoot.bind("<Destroy>", self.enableGUI)
        AnalysisOptionsInputForm(readStatOptionsRoot, self).pack()
        
        # Set the X and the Y Position of the options window, so it is nearby.  
        readStatOptionsRoot.update()        
        windowXpos = str(self.parent.winfo_geometry().split('+')[1])
        windowYpos = str(self.parent.winfo_geometry().split('+')[2])
        newGeometry = (str(readStatOptionsRoot.winfo_width()) + 'x' 
            + str(readStatOptionsRoot.winfo_height()) + '+' 
            + str(windowXpos) + '+' 
            + str(windowYpos))
        readStatOptionsRoot.geometry(newGeometry)
        
        readStatOptionsRoot.mainloop()

    def enableGUI(self, event=None):
        self.toggleGUI(True)  
        
    def disableGUI(self):
        self.toggleGUI(False)   
        
    def toggleGUI(self, isEnabled): 
        #print ('Toggling GUI Widgets:' + str(isEnabled))
         
        newState = (NORMAL if (isEnabled) else DISABLED)
        
        # Choosing the widgets individually, this makes the most sense I think.
        self.chooseInputButton.config(state=newState) 
        self.chooseOutputButton.config(state=newState) 
        
        self.prepareReadsButton.config(state=newState)
        self.analysisOptionsButton.config(state=newState) 
        self.runAnalysisButton.config(state=newState) 
        
        self.howToUseButton.config(state=newState)
        self.citationButton.config(state=newState)
        self.saddlebagsButton.config(state=newState)   
        
        # This repaints the window, we can see what's going on.
        self.update()      
        
    def runFullAnalysis(self):
        #print('Starting a full analysis.....')
        self.recordAnalysisStep('Starting a full analysis...')
        
        # Set the input and output directories
        self.assignIODirectories()
        
        # Start the Analysis report   
        self.reportResults('HLA Analysis') 
        self.reportResults('Start: ' + str(datetime.now()))         
        
        # 1) Nitpicker to prepare reads
        self.constructInitialReadStats()  
          
        # 2) Sort reads by HLA Locus
        self.sortByLocus()        
        
        # 3) Assemble each locus file
        self.assembleSortedReads()     
           
        # 4) Allele Call results if possible
        self.alleleCall()        
        
        # 5) Summarize Results
        self.summarizeResults()
        
    def assignIODirectories(self):
        assignConfigurationValue('input_directory', self.inputDirectoryText.get())
        assignConfigurationValue('output_directory', self.outputDirectoryText.get())
        
        assignConfigurationValue('log_file_location', join(self.outputDirectoryText.get(), 'nanopore-prospector-log.txt'))
        assignConfigurationValue('report_file_location', join(self.outputDirectoryText.get(), 'AnalysisReport.txt'))
        
        #self.consensusSequenceFileName = join(self.outputDirectoryText.get(), 'Consensus.fasta')
        
        self.recordAnalysisStep('Initialized Nanopore Prospector Log File')
        self.recordAnalysisStep('Read Input Directory:' + getConfigurationValue('input_directory'))
        self.recordAnalysisStep('Analysis Output Directory:' + getConfigurationValue('output_directory'))
        
        # TODO Create the results output folder here.   
        

    def reportReadStats(self, preparedReadResults):
        self.reportResults('\nRead Statistics:\n')
        #all_reads, length_rejected, quality_rejected, pass_reads,
        if ('all_reads' in preparedReadResults.keys()):
            readLengths = preparedReadResults['all_reads'][0]
            readQualities = preparedReadResults['all_reads'][1]
            self.reportResults('All Reads')
            self.reportResults(str(len(readLengths)) + ' reads analyzed.')
            self.reportResults('Length  (Min/Avg/Max): ' + str(int(amin(readLengths))) + '/' + str(int(mean(readLengths))) + '/' + str(int(amax(readLengths))))
            self.reportResults('Q-Score (Min/Avg/Max): ' + str(int(amin(readQualities))) + '/' + str(int(mean(readQualities))) + '/' + str(int(amax(readQualities))) + '\n')
            
            if ('pass_reads' in preparedReadResults.keys()):
                readLengths = preparedReadResults['pass_reads'][0]
                readQualities = preparedReadResults['pass_reads'][1]
                self.reportResults('"Pass" Reads')
                self.reportResults(str(len(readLengths)) + ' reads have the correct length and quality.')
                self.reportResults('Length  (Min/Avg/Max): ' + str(int(amin(readLengths))) + '/' + str(int(mean(readLengths))) + '/' + str(int(amax(readLengths))))
                self.reportResults('Q-Score (Min/Avg/Max): ' + str(int(amin(readQualities))) + '/' + str(int(mean(readQualities))) + '/' + str(int(amax(readQualities))) + '\n')
                
            if ('length_rejected' in preparedReadResults.keys()):
                readLengths = preparedReadResults['length_rejected'][0]
                readQualities = preparedReadResults['length_rejected'][1]
                self.reportResults('Length-Reject Reads')
                self.reportResults(str(len(readLengths)) + ' reads were rejected for improper length.')
                self.reportResults('Length  (Min/Avg/Max): ' + str(int(amin(readLengths))) + '/' + str(int(mean(readLengths))) + '/' + str(int(amax(readLengths))))
                self.reportResults('Q-Score (Min/Avg/Max): ' + str(int(amin(readQualities))) + '/' + str(int(mean(readQualities))) + '/' + str(int(amax(readQualities))) + '\n')
                
            if ('quality_rejected' in preparedReadResults.keys()):
                readLengths = preparedReadResults['quality_rejected'][0]
                readQualities = preparedReadResults['quality_rejected'][1]
                self.reportResults('Quality-Reject Reads')
                self.reportResults(str(len(readLengths)) + ' reads were rejected for improper quality.')
                self.reportResults('Length  (Min/Avg/Max): ' + str(int(amin(readLengths))) + '/' + str(int(mean(readLengths))) + '/' + str(int(amax(readLengths))))
                self.reportResults('Q-Score (Min/Avg/Max): ' + str(int(amin(readQualities))) + '/' + str(int(mean(readQualities))) + '/' + str(int(amax(readQualities))) + '\n')
                
                
        else:
            raise Exception ('There is no all_reads category in the read results. This is unexpected, why?')

    def constructInitialReadStats(self):
        self.recordAnalysisStep('Step 1.) Calculating initial read stats')
        self.disableGUI()
        
        # Set the input and output directories
        self.assignIODirectories()
        writeConfigurationFile()
   
        # Run nit-picker for output directory
        # TODO: fix these parameters, especially sample ID. 
        # TODO: RUn this in a thread, so the GUI can update.
        
        preparedReadsOutputDirectory = join(getConfigurationValue('output_directory'), '1_prepared_reads')
        
        #if(self.demultiplexReads):            
        #    sampleID = 'READS'
        #else:
        #    sampleID = 'READS'
        
        # TODO: I am always using the leaf directory name. Does this work for multiplexed and demultiplexed samples?
        #self.inputDirectoryText.set(currentInputDirectory)
        #parentDir = abspath(join(currentInputDirectory, os.pardir))
        leafDirName = basename(normpath(self.inputDirectoryText.get()))
        sampleID = leafDirName
        # suggestedOutputDirectory = join(parentDir,leafDirName + '_analysis')

        if not exists(preparedReadsOutputDirectory):
            makedirs(preparedReadsOutputDirectory)

        # preparedReadResults is a dictionary. 
        # The result will have 2 arrays
        # readstats is a 2d array, with lengths and qualities.
        preparedReadResults = None
        
        # If demultiplex option is set to 1, demultiplex is "on".
        if (getConfigurationValue('demultiplex_reads') == '1'):
            
            barcodeFilePath = getBarcodeFilePath()
      
            preparedReadResults = prepareReads(getConfigurationValue('input_directory')
                , preparedReadsOutputDirectory
                , sampleID
                , barcodeFilePath
                , None # Reference File, we don't have one.
                , int(getConfigurationValue('min_length'))
                , int(getConfigurationValue('max_length'))
                , int(getConfigurationValue('min_quality'))
                , int(getConfigurationValue('max_quality')))
        else:
            preparedReadResults = prepareReads(getConfigurationValue('input_directory')
                , preparedReadsOutputDirectory
                , sampleID
                , None # No Barcoding file.
                , None # No Allele Reference
                , int(getConfigurationValue('min_length'))
                , int(getConfigurationValue('max_length'))
                , int(getConfigurationValue('min_quality'))
                , int(getConfigurationValue('max_quality')))
            
        self.reportReadStats(preparedReadResults)
        self.recordAnalysisStep('Done calculating initial read stats')
        self.enableGUI()  
        
    def sortByLocus(self):
        self.recordAnalysisStep('Step 2.) Sort reads by HLA Locus')
        self.disableGUI()
        
        # TODO: What if it's demultiplexed?
        # I should sort each barcode.
        # Make a loop, find barcodes, sort each one in a subrdirectory
        # What if i use some other sample ID, other than READS
        # TODO: Look for fastq files dynamically in output directory.        
        #sampleID = 'READS'        
        # TODO: I need to rethink how to handle sampleID when I am demultiplexing.
        # I can handle one sample one folder.
        
        # TODO: threadCount should be a parameter somewhere.
        # Specified in the options, probably.
        threadCount = 4
        
        #preparedReadsInputFile = join(join(self.outputDirectoryText.get(), '1_prepared_reads'),sampleID + '_Pass.fastq')
        preparedReadsInput = join(self.outputDirectoryText.get(), '1_prepared_reads')
        #sortedReadsOutputDirectory = join(join(self.outputDirectoryText.get(), '2_sorted_reads'), sampleID)
        sortedReadsOutputDirectory = join(self.outputDirectoryText.get(), '2_sorted_reads')
        if not exists(sortedReadsOutputDirectory):
            makedirs(sortedReadsOutputDirectory)
            
        # TODO: Make more than one gene references. Split by locus.
        # Fine. That's what I'm doing.
        
        # Make a list of the reference files I need.
        referenceFileList = []
        
        print('the configuration value for analyze_hla_a is:' + str(getConfigurationValue('analyze_hla_a')))
        
        if(str(getConfigurationValue('analyze_hla_a')) == '1'):
            referenceFileList.append(getBlastSortResourceLocation('HLA_A_BlastReference.fasta'))
        if(str(getConfigurationValue('analyze_hla_b')) == '1'):
            referenceFileList.append(getBlastSortResourceLocation('HLA_B_BlastReference.fasta'))
        if(str(getConfigurationValue('analyze_hla_c')) == '1'):
            referenceFileList.append(getBlastSortResourceLocation('HLA_C_BlastReference.fasta'))
        if(str(getConfigurationValue('analyze_hla_e')) == '1'):
            referenceFileList.append(getBlastSortResourceLocation('HLA_E_BlastReference.fasta'))
           
        if(str(getConfigurationValue('analyze_hla_dra')) == '1'):
            referenceFileList.append(getBlastSortResourceLocation('HLA_DRA_BlastReference.fasta'))    
        if(str(getConfigurationValue('analyze_hla_dqa1')) == '1'):
            referenceFileList.append(getBlastSortResourceLocation('HLA_DQA1_BlastReference.fasta'))    
        if(str(getConfigurationValue('analyze_hla_dqb1')) == '1'):
            referenceFileList.append(getBlastSortResourceLocation('HLA_DQB1_BlastReference.fasta'))    
        if(str(getConfigurationValue('analyze_hla_drb1')) == '1'):
            referenceFileList.append(getBlastSortResourceLocation('HLA_DRB1_BlastReference.fasta'))
        
        print('I found this many blast references:' + str(len(referenceFileList)))
        
        sortReferencePath = combineBlastReferences(referenceFileList, join(self.outputDirectoryText.get(), 'blast_sort_reference'))      
        
        # the key is the name of the file analyzed.
        # The value is a list of minion_read_collections, pertaining to what gene they sorted to.          
        sortResults = sortDirectory(preparedReadsInput, sortedReadsOutputDirectory, sortReferencePath, threadCount)
        
        # Report the BLAST sorting results in the results summary.
        self.reportResults('Read Sorting:\n')
        # Loop through the analyzed read files (probably just one here)
        for analyzedReadResult in sortResults.keys():
            
            #print ('Looking at this read collection key:' + analyzedReadResult)
            blastGeneList = sortResults[analyzedReadResult]
            #print ('Looking at this list of sorted read groups\n:' + str(blastGeneList))
            
            self.reportResults(analyzedReadResult + ':')
            
            # loop through each entry in this 
            for genewiseSortedReadGroup in blastGeneList:
            
                if (genewiseSortedReadGroup.gene is None):
                    self.reportResults('Unsorted: ' + str(len(genewiseSortedReadGroup.readCollection)) + ' reads.')
                else:
                    self.reportResults('HLA-' + str(genewiseSortedReadGroup.gene) + ': ' + str(len(genewiseSortedReadGroup.readCollection)) + ' reads.')
                

        self.reportResults('')
        self.recordAnalysisStep('Done sorting reads by HLA Locus')
        self.enableGUI()
        
    def assembleSortedReads(self):
        self.recordAnalysisStep('Step 3.) Assemble sorted Reads')
        self.disableGUI()
        
        self.reportResults('Read Assembly:\n')
        
        # TODO: Deal with sample IDs better. 
        # What about demultiplexed samples etc?
        #sampleID = 'Pass_Reads'  
        
        #sortedReadsOutputDirectory = join(join(self.outputDirectoryText.get(), '2_sorted_reads'), sampleID)
        #readAssemblyDirectory = join(join(self.outputDirectoryText.get(), '3_assembled_reads'), sampleID)
        sortedReadsOutputDirectory = join(self.outputDirectoryText.get(), '2_sorted_reads')
        readAssemblyDirectory = join(self.outputDirectoryText.get(), '3_assembled_reads')
        
        # Dictionary.  Key is the location of the generated fasta consensus sequences.
        # Value is the # of reads represented in this consensus.
        readAssemblyResults = {}

        if not exists(readAssemblyDirectory):
            makedirs(readAssemblyDirectory)
        
        # for each barcode folder in this directory
        for sortedReadDirectoryObject in listdir(sortedReadsOutputDirectory):
            
            fullBarcodeDirPath = join(sortedReadsOutputDirectory,sortedReadDirectoryObject)
            # Should only be directories under here.
            # Either barcode directories, or a "pass" folder
            if (isfile(fullBarcodeDirPath)):
                raise Exception('There is a file in the sorted read output directory. I expected only barcode directories:' + str(sortedReadDirectoryObject))
            elif (isdir(fullBarcodeDirPath)):
                # We found a barcode folder. I should assemble the fastq files in here.
                #self.recordAnalysisStep('Assembling reads in this directory:' + fullBarcodeDirPath)
                
                # Do not assemble the "All Reads" directory, waste of time.
                # TODO: How should this work? Make sure this makes sense when I am debarcoding this.
                # TODO: I dont think this is working at all. All_Reads is not used?
                if('All_Reads' not in fullBarcodeDirPath):
                    
                    logMessageToFile('Assembling this directory:' + str(fullBarcodeDirPath))
                    for currentReadFileName in listdir(fullBarcodeDirPath):
                        # Only fastq files
                        # Also, skip assembling the unsorted reads.
                        if(".fastq"== currentReadFileName[-6:] or ".fq" == currentReadFileName[-3:]):
                            
                            currentReadFilePath = join(fullBarcodeDirPath, currentReadFileName)
                            currentAssemblyOutputDirectory = join(join(readAssemblyDirectory, sortedReadDirectoryObject),currentReadFileName)
                            
                            # TODO: This will cause problems. Sometimes I WANT to assemble unsorted reads.
                            if('unsorted_' in currentReadFileName):
                                self.recordAnalysisStep('Skip assembly of these unsorted reads:' + str(currentReadFilePath))
                                
                            else:
                                self.recordAnalysisStep('Assembling these reads:' + str(currentReadFilePath))
                                
                                # TODO Parameters
                                numberIterations = 6
                                numberThreads = 4
                                splitHeterozygotes = True

                                myAlleleWrangler = AlleleWrangler(currentReadFilePath, currentAssemblyOutputDirectory, None, numberIterations, numberThreads, splitHeterozygotes)
                                currentAssemblyResults = myAlleleWrangler.analyzeReads()
                                
                                #Merge into completed results
                                for key in currentAssemblyResults.keys():
                                    readAssemblyResults[key] = currentAssemblyResults[key]
                     
                    else:
                        logMessageToFile('Skipping assembly on this directory:' + str(fullBarcodeDirPath))
                pass
            else:
                raise Exception('This object is neither a file or directory. How mysterious. I expected only barcode directories:' + str(fullBarcodeDirPath))

        # Gather each fasta consensus sequence and write it to a summary.fasta file.
        consensusSequences = []
        
        #self.reportResults('AssemblyResults:\n' + str(readAssemblyResults))
        
        for index, key in enumerate(readAssemblyResults.keys()):
            readCount = readAssemblyResults[key]
            
            self.reportResults('Consensus Sequence: ' + str(key))
            self.reportResults(str(readCount) + ' reads represented in this assembly.\n')
            # Read sequence from assembly fasta            
            currentConsensusSeq = list(parseReads(key, 'fasta'))[0]
            currentConsensusSeq.id = 'Consensus_Sequence_' + str(index + 1)
            
            
            currentConsensusSeq.description = 'Read_Count=' + str(readCount)
                       
            # Add sequence to list
            consensusSequences.append(currentConsensusSeq)

        # Write out one fasta with all consensus sequences.
        combinedConsensusOutputFilename = join(self.outputDirectoryText.get(), 'Consensus_Sequences.fasta')
        sequenceWriter = createOutputFile(combinedConsensusOutputFilename)
        write(consensusSequences, sequenceWriter, 'fasta')
        sequenceWriter.close()
        
        self.recordAnalysisStep('Done assembling sorted Reads')
        self.enableGUI()
        
    def alleleCall(self):
        self.recordAnalysisStep('Step 4.) HLA Allele Calling')
        self.disableGUI()
        
        self.recordAnalysisStep('Just kidding. This feature must be implemented')
        
        self.recordAnalysisStep('Done with HLA Allele Calling')
        self.enableGUI()
        
    def summarizeResults(self):
        self.recordAnalysisStep('Step 5.) Summarize Results')
        self.disableGUI()
        
        self.recordAnalysisStep('Just kidding. This feature must be implemented')
        
        self.recordAnalysisStep('Done Summarizing Results')
        self.enableGUI()

    # chooseReadInputDirectory method is called when the user presses the input directory button
    def chooseReadInputDirectory(self):
        print ('Choosing an input directory.')
        # TODO: What happens when I cancel this? I need to re-enable the GUI.
        self.disableGUI()
        self.setInputDir(filedialog.askdirectory(**self.dir_opt), True)
        self.enableGUI()
        
    # chooseReadInputDirectory method is called when the user presses the input directory button
    
    def chooseReadOutputDirectory(self):        
        print ('Choosing an output directory.')
        self.disableGUI()
        self.setOutputDir(filedialog.askdirectory(**self.dir_opt))
        self.enableGUI()
    
    def setInputDir(self, currentInputDirectory, setOutputDir):
        # setInputDir is a subprocess to be called when an input directory is selected.
        self.inputDirectoryText.set(currentInputDirectory)
        
        #print ('Just set the self.inputDirectory text to this:' + self.inputDirectoryText.get())
        
        # TODO: Is there at least one fastq file in this directory?  It has reads?
        # If not, this might be a directory of pre-demultiplexed reads.
        
        # TODO: Do I need to parse subfolders for any reason? I think not.
        # TODO: In the case of pre-demultiplexed reads, maybe.
        
        # Popup, Should I use the suggested output directory?
        # Include a boolean, in case i don't want to mess with output directory yet.
        if setOutputDir:
            parentDir = abspath(join(currentInputDirectory, pardir))
            leafDirName = basename(normpath(currentInputDirectory))
            suggestedOutputDirectory = join(parentDir,leafDirName + '_analysis')
            queryResult = messagebox.askquestion('Use Suggested Output Directory?'
                , 'Use this output directory?\n' 
                + suggestedOutputDirectory 
                + '\nSelect No to choose your own output directory', icon='warning')
            if queryResult == 'yes':
                self.setOutputDir(suggestedOutputDirectory)
            else:
                self.chooseReadOutputDirectory()
        # If Yes, Do it
        # If no, call chooseReadOutputDirectory 

    def setOutputDir(self, currentOutputDirectory):
        self.outputDirectoryText.set(currentOutputDirectory)
        
        # This output directory should exist
        if not exists(currentOutputDirectory):
            makedirs(currentOutputDirectory)
  
        # Code for input == output?  I guess I really don't care if it's the same dir.
        #if (self.outputDirectoryText.get() == self.inputDirectoryText.get()):
        #    messagebox.showwarning(
        #    "Input = Output",
        #    "The input directory is the same as the output directory." % filename
        #)
        # TODO: If output == input, there are problems and i should just quit.
        #self.initializeLogAndReport()
    
  