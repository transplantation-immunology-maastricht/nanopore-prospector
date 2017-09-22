
import os
import datetime
from os import makedirs#, listdir
from os.path import exists, split, isdir, isfile, join, abspath, basename, normpath#, isfile,
import Tkinter, Tkconstants, tkFileDialog, tkMessageBox
from Tkconstants import *
from Tkinter import *
from AnalysisOptionsInputForm import *
#from serial.win32 import ReadFile

# Import my other submodules
# Code exists in 
# ../{submodule}/src
# Add that path to the sys.path for each submodule
sys.path.insert(0, join(os.pardir,join('nit-picker','src')))
sys.path.insert(0, join(os.pardir,join('saddle-bags','src')))
sys.path.insert(0, join(os.pardir,join('punkin-chunker','src')))
sys.path.insert(0, join(os.pardir,join('allele-wrangler','src')))
# TODO I think these path imports aren't very good.

#from AlleleGui import AlleleGui
from nit_picker import prepareReads
#import nit_picker

#from punkin_chunker import punkin_chunker
from punkin_chunker import sortDirectory
#import punkin_chunker

from AlleleWrangler import AlleleWrangler

class NanoporeProspectorMasterFrame(Tkinter.Frame):
    def __init__(self, root):
        
        #self.useBarcoding = False

        Tkinter.Frame.__init__(self, root)
        root.title("nanopore prospector")

        self.parent = root
        
        self.initialize()

    # Initialize GUI elements
    def initialize(self):
        # TODO: Load from a configuration file.  I should be able to save / load a configuration file
        # Otherwise provide defaults
        # I have good config file code in Saddlebags

        self.logFileName = None
        
        self.demultiplexReads = False
           
        self.minimumReadLength = 0
        self.maximumReadLength = 999999
        
        self.minimumQuality = 0
        self.maximumQuality = 99 

        # This is the directory the python executable is running from.
        FileAndPath = abspath(__file__)
        self.idir, self.ifile = split(FileAndPath)

        # GUI options
        #self.label_opt = {'fill': Tkconstants.BOTH, 'padx': 10, 'pady': 10}
        self.button_opt = {'fill': Tkconstants.BOTH, 'padx': 50, 'pady': 15}
        self.dir_opt = {'initialdir': self.idir, 'mustexist': True, 
            'parent': self, 'title': 'Choose a directory'}
        
        # A frame for choosing the instructions
        self.instructionsFrame = self.makeInstructionsFrame()
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
        
    def logMessage(self, currentMessage):
        fullLogMessage = str(datetime.datetime.now()) + ' : ' + currentMessage + '\n'
        
        # Append the GUI with the Log Entry 
        self.submOutputGuiObject.insert(END, fullLogMessage) 
        
        if (self.logFileName is None):
            print ('No Logfile exists yet.\nTrying to log this message:\n' + currentMessage)
        else:
            #print ('Logging message: ' + currentMessage) 
            # Append the log with the log entry            
            resultsOutput = open(self.logFileName, 'a')
            resultsOutput.write(fullLogMessage)
            resultsOutput.close()
            
        # This repaints the window, we can see what's going on.
        self.update()   
    
    def makeInstructionsFrame(self):
        instructionsFrame = Tkinter.Frame(self)

        self.instructionText = Tkinter.StringVar()
        self.instructionText.set('Use this interface to look at a group of \nMinION Reads containing HLA amplicon sequences:')
        Tkinter.Label(instructionsFrame, width=80, height=5, textvariable=self.instructionText).pack()
        return instructionsFrame
      
    def makeChooseIODirectoriesFrame(self):
        chooseDirectoryFrame = Tkinter.Frame(self)        

        self.chooseInputButton = Tkinter.Button(chooseDirectoryFrame, text='Choose Input Directory ', command=self.chooseReadInputDirectory)
        self.chooseInputButton.grid(row=0, column=0, sticky=Tkconstants.W)
        self.chooseOutputButton = Tkinter.Button(chooseDirectoryFrame, text='Choose Output Directory', command=self.chooseReadOutputDirectory)
        self.chooseOutputButton.grid(row=1, column=0, sticky=Tkconstants.W)
        
        self.inputDirectoryText = Tkinter.StringVar()
        self.inputDirectoryText.set('Where is your MinION Read Directory?')
        Tkinter.Entry(chooseDirectoryFrame, width=60, textvariable=self.inputDirectoryText).grid(row=0, column=1)
                
        self.outputDirectoryText = Tkinter.StringVar()
        self.outputDirectoryText.set('Output Directory')
        Tkinter.Entry(chooseDirectoryFrame, width=60, textvariable=self.outputDirectoryText).grid(row=1, column=1)
        
        return chooseDirectoryFrame
        
    def makeAnalysisButtonsFrame(self):
        analysisButtonsFrame = Tkinter.Frame(self) 
        self.analysisOptionsButton = Tkinter.Button(analysisButtonsFrame, text='Analysis Options', command=self.specifyReadStatOptions)
        self.analysisOptionsButton.grid(row=0, column=0)
        self.prepareReadsButton = Tkinter.Button(analysisButtonsFrame, text='Demultiplex + Prepare Reads', command=self.constructInitialReadStats)
        self.prepareReadsButton.grid(row=1, column=0)
        self.runAnalysisButton = Tkinter.Button(analysisButtonsFrame, text='Run Full Analysis', command=self.runFullAnalysis)
        self.runAnalysisButton.grid(row=2, column=0)
        return analysisButtonsFrame
        
    def makeAnalysisLogFrame(self):
        logFrame = Tkinter.Frame(self)
        
        # TODO: Did i set this log label up correctly?
        self.logLocationText = Tkinter.StringVar()
        self.logLocationText.set('Choose a read directory to start logging...')
        Tkinter.Label(logFrame, width=80, height=5, textvariable=self.logLocationText).pack()
       

        self.submOutputXScrollbar = Scrollbar(logFrame, orient=HORIZONTAL)
        self.submOutputXScrollbar.pack(side=BOTTOM, fill=X)

        self.submOutputYScrollbar = Scrollbar(logFrame)
        self.submOutputYScrollbar.pack(side=RIGHT, fill=Y)

        self.submOutputGuiObject = Tkinter.Text(
            logFrame, width=100, height=10, wrap=NONE
            , xscrollcommand=self.submOutputXScrollbar.set
            , yscrollcommand=self.submOutputYScrollbar.set
        )

        self.submOutputXScrollbar.config(command=self.submOutputGuiObject.xview)
        self.submOutputYScrollbar.config(command=self.submOutputGuiObject.yview) 

        self.submOutputGuiObject.pack() 
        
        return logFrame
    
    def makeMoreInfoFrame(self):    
        moreInfoFrame = Tkinter.Frame(self)  
        self.howToUseButton = Tkinter.Button(moreInfoFrame, text='How to use this tool', command=self.howToUse)
        self.howToUseButton.grid(row=0, column=0)
        self.citationButton = Tkinter.Button(moreInfoFrame, text='Contacting or Citing MUMC', command=self.contactInformation)
        self.citationButton.grid(row=0, column=1)
        self.saddlebagsButton = Tkinter.Button(moreInfoFrame, text='SaddleBags - A (Novel) Allele Submission Tool', command=self.launchSaddleBags)
        self.saddlebagsButton.grid(row=0, column=2)
        return moreInfoFrame
     
    def launchSaddleBags(self):
        #tkMessageBox.showinfo('Need to open allelesub tool','Popup a window containing the allele submission tool please.')
        # TODO: This works but the interface is messed up.
        # I think saddlebags is using "self" to assign variables when it shouldnt. They should assign to the frame instead?
        # Problem is, i should be making a toplevel object.  I can fix this.
        # See how i did it in specifyReadStatOptions
        # TODO This actually doesn't work at all, I renamed the saddlebags classes.
        saddleBagsRoot = Tkinter.Tk()
        AlleleGui(saddleBagsRoot).pack()
        saddleBagsRoot.mainloop()
        
    # This method should popup some instruction text in a wee window.
    # This should be explicit on how to use the tool.    
    # TODO: This howToUse method is pretty weak.    
    def howToUse(self):
        tkMessageBox.showinfo('Select a directory containing reads.\n'
            + 'Do some analysis on the reads.'
            + 'Ben:Fill in better how-to-use instructions.'
            )
        
    def contactInformation(self):
        # This method should list contact information for MUMC, and a link to the github page.  
        tkMessageBox.showinfo('Contact Information',
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


    # TODO: all of these analysis buttons.  Do something....  
    def specifyReadStatOptions(self):
        print('Specifying ReadStat Options.....')
    
        self.disableGUI()
        
        readStatOptionsRoot = Tkinter.Toplevel()
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
        
        #print ('the children:' + str(readStatOptionsRoot.children))
        #print ('it is a type of:' + str(type(readStatOptionsRoot.children)))
        childGUIObject = readStatOptionsRoot.children[list(readStatOptionsRoot.children)[0]]
        #print ('found child object:' + str(type(childGUIObject)))
        
        # Set initial options
        if(self.demultiplexReads):
            childGUIObject.chooseDemultiplexIntVar.set(1)
        else: 
            childGUIObject.chooseDemultiplexIntVar.set(2)
            
        childGUIObject.inputMinLength.set(self.minimumReadLength)
        childGUIObject.inputMaxLength.set(self.maximumReadLength)
        childGUIObject.inputMinQuality.set(self.minimumQuality)
        childGUIObject.inputMaxQuality.set(self.maximumQuality)

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
        self.logMessage('Starting a full analysis.....')
        
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
        
        
        
        
    def constructInitialReadStats(self):
        self.logMessage('Step 1.) Calculating initial read stats')
        self.disableGUI()

        # Run nit-picker for output directory
        # TODO: fix these parameters, especially sample ID. 
        # TODO: RUn this in a thread, so the GUI can update.

        preparedReadsOutputDirectory = join(self.outputDirectoryText.get(), '1_prepared_reads')
        
        
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

        if (self.demultiplexReads):
            try:
                # PyInstaller creates a temp folder and stores path in _MEIPASS
                barcodeFilePath = join(sys._MEIPASS, 'barcodes_96_bc_kit.txt')
                #print('using packaged MEIPASS directory to find barcode.')
            except Exception:
                barcodeFilePath = join(os.path.abspath('.'),'../nit-picker/barcodes/barcodes_96_bc_kit.txt')
                #print('using local path to barcode kit file.')

            prepareReads(self.inputDirectoryText.get()
                , preparedReadsOutputDirectory
                , sampleID
                , barcodeFilePath
                , self.minimumReadLength
                , self.maximumReadLength
                , self.minimumQuality
                , self.maximumQuality )
        else:
            prepareReads(self.inputDirectoryText.get()
                , preparedReadsOutputDirectory
                , sampleID
                , None
                , self.minimumReadLength
                , self.maximumReadLength
                , self.minimumQuality
                , self.maximumQuality )

        # TODO: Popup information about what we just did.        
        # Are we in windows? 
        # If yes, open explorer with the read stats directory
        # If no, popup a dialog with the directory, user can do it themselves.

        self.logMessage('Done calculating initial read stats')
        self.enableGUI()  
        
    def sortByLocus(self):
        self.logMessage('Step 2.) Sort reads by HLA Locus')
        self.disableGUI()
        
        # TODO: What if it's demultiplexed?
        # I should sort each barcode.
        # Make a loop, find barcodes, sort each one in a subrdirectory
        # What if i use some other sample ID, other than READS
        # TODO: Look for fastq files dynamically in output directory.        
        sampleID = 'READS'        
        
        # TODO: threadCount should be a parameter somewhere.
        # Specified in the options, probably.
        threadCount = 2
        
        #preparedReadsInputFile = join(join(self.outputDirectoryText.get(), '1_prepared_reads'),sampleID + '_Pass.fastq')
        preparedReadsInput = join(self.outputDirectoryText.get(), '1_prepared_reads')
        sortedReadsOutputDirectory = join(join(self.outputDirectoryText.get(), '2_sorted_reads'), sampleID)
        if not exists(sortedReadsOutputDirectory):
            makedirs(sortedReadsOutputDirectory)
            
            
        try:
            # PyInstaller creates a temp folder and stores path in _MEIPASS
            sortReferencePath = join(sys._MEIPASS, 'HLA_ClassI_GeneRef.fasta')
            #print('using packaged MEIPASS directory to find barcode.')
        except Exception:
            sortReferencePath = join(os.path.abspath('.'),'../punkin-chunker/inputData/HLA_ClassI_GeneRef.fasta')
            #print('using local path to barcode kit file.')
            
        sortDirectory(preparedReadsInput, sortedReadsOutputDirectory, sortReferencePath, threadCount)
            
        
        
        
        self.logMessage('Done sorting reads by HLA Locus')
        self.enableGUI()
        
    def assembleSortedReads(self):
        self.logMessage('Step 3.) Assemble sorted Reads')
        self.disableGUI()
        
        # TODO: Deal with sample IDs better. 
        sampleID = 'READS'  
        
        sortedReadsOutputDirectory = join(join(self.outputDirectoryText.get(), '2_sorted_reads'), sampleID)
        readAssemblyDirectory = join(join(self.outputDirectoryText.get(), '3_assembled_reads'), sampleID)
        
        if not exists(readAssemblyDirectory):
            makedirs(readAssemblyDirectory)
        
        # for each barcode folder in this directory
        for sortedReadDirectoryObject in os.listdir(sortedReadsOutputDirectory):
            
            fullBarcodeDirPath = join(sortedReadsOutputDirectory,sortedReadDirectoryObject)
            # Should only be directories under here.
            # Either barcode directories, or a "pass" folder
            if (isfile(fullBarcodeDirPath)):
                raise Exception('There is a file in the sorted read output directory. I expected only barcode directories:' + str(sortedReadDirectoryObject))
            elif (isdir(fullBarcodeDirPath)):
                # We found a barcode folder. I should assemble the fastq files in here.
                #self.logMessage('Assembling reads in this directory:' + fullBarcodeDirPath)
                for currentReadFileName in os.listdir(fullBarcodeDirPath):
                    # Only fastq files
                    if(".fastq"== currentReadFileName[-6:] or ".fq" == currentReadFileName[-3:]):
                        
                        currentReadFilePath = join(fullBarcodeDirPath, currentReadFileName)
                        currentAssemblyOutputDirectory = join(join(readAssemblyDirectory, sortedReadDirectoryObject),currentReadFileName)
                        
                        self.logMessage('Assembling these reads:' + str(currentReadFilePath))
                        
                        # TODO Parameters
                        numberIterations = 6
                        numberThreads = 4
                        splitHeterozygotes = True
                        
                        # TODO: fix these parameters.
                        # I should pass in iteration numbers etc.
                        myAlleleWrangler = AlleleWrangler(currentReadFilePath, currentAssemblyOutputDirectory, None, numberIterations, numberThreads, splitHeterozygotes)
                        myAlleleWrangler.analyzeReads()

                pass
            else:
                raise Exception('This object is neither a file or directory. How mysterious. I expected only barcode directories:' + str(fullBarcodeDirPath))

        
        
        self.logMessage('Done assembling sorted Reads')
        self.enableGUI()
        
    def alleleCall(self):
        self.logMessage('Step 4.) HLA Allele Calling')
        self.disableGUI()
        
        self.logMessage('Just kidding. This feature must be implemented')
        
        self.logMessage('Done with HLA Allele Calling')
        self.enableGUI()
        
    def summarizeResults(self):
        self.logMessage('Step 5.) Summarize Results')
        self.disableGUI()
        
        self.logMessage('Just kidding. This feature must be implemented')
        
        self.logMessage('Done Summarizing Results')
        self.enableGUI()

    # chooseReadInputDirectory method is called when the user presses the input directory button
    def chooseReadInputDirectory(self):
        print ('Choosing an input directory.')
        #self.enableInterface(False)
        self.disableGUI()
        self.setInputDir(tkFileDialog.askdirectory(**self.dir_opt))
        #print('my file dialog returned this value:')
        #self.enableInterface(True)        
        self.enableGUI()
        
    # chooseReadInputDirectory method is called when the user presses the input directory button
    def chooseReadOutputDirectory(self):        
        print ('Choosing an output directory.')
        self.disableGUI()
        #self.enableInterface(False)        
        self.setOutputDir(tkFileDialog.askdirectory(**self.dir_opt))
        #self.enableInterface(True)
        self.enableGUI()
        

    # setInputDir is a subprocess to be called when an input directory is selected.
    def setInputDir(self, currentInputDirectory):
        self.inputDirectoryText.set(currentInputDirectory)
        
        print ('Just set the self.inputDirectory text to this:' + self.inputDirectoryText.get())
        
        # TODO: Is there at least one fastq file in this directory?  
        
        # TODO: Do I need to parse subfolders for any reason? I think not.
        
        # Popup, Should I use the suggested output directory?
        parentDir = abspath(join(currentInputDirectory, os.pardir))
        leafDirName = basename(normpath(currentInputDirectory))
        suggestedOutputDirectory = join(parentDir,leafDirName + '_analysis')
        queryResult = tkMessageBox.askquestion('Use Suggested Output Directory?'
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
        #    tkMessageBox.showwarning(
        #    "Input = Output",
        #    "The input directory is the same as the output directory." % filename
        #)
        # TODO: If output == input, there are problems and i should just quit.
        self.initializeLog()
    
    
    def initializeLog(self):
        self.logFileName = join(self.outputDirectoryText.get(), 'nanopore-prospector-log.txt')
        #self.logFile = createOutputFile(logFileName)
        self.logMessage('Initialized Nanopore Prospector Log File')
        self.logMessage('Read Input Directory:' + self.inputDirectoryText.get())
        self.logMessage('Analysis Output Directory:' + self.outputDirectoryText.get())
        # TODO: Put an analysis log file in the output directory.
