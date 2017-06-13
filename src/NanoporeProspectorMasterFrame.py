
import os
import datetime
from os import makedirs#, listdir
from os.path import exists, split, isdir, join, abspath, basename, normpath#, isfile,
import Tkinter, Tkconstants, tkFileDialog, tkMessageBox
from Tkconstants import *
from Tkinter import *

# Import my other submodules
# Code exists in 
# ../{submodule}/src
# Add that path to the sys.path for each submodule
sys.path.insert(0, join(os.pardir,join('nit-picker','src')))
sys.path.insert(0, join(os.pardir,join('saddle-bags','src')))

from AlleleGui import AlleleGui
#from nit_picker import *
import nit_picker

# This method is a directory-safe way to open up a write file.
# TODO: Do i need this method?  
#def createOutputFile(outputfileName):
#    tempDir, tempFilename = split(outputfileName)
#    if not isdir(tempDir):
#        print('Making Directory:' + tempDir)
#        makedirs(tempDir)
#    resultsOutput = open(outputfileName, 'w')
#    return resultsOutput

class NanoporeProspectorMasterFrame(Tkinter.Frame):
    def __init__(self, root):

        Tkinter.Frame.__init__(self, root)
        root.title("nanopore prospector")

        self.parent = root
        
        self.initialize()

    # Initialize GUI elements
    def initialize(self):
        # TODO: Load from a configuration file.  I should be able to save / load a configuration file
        # Otherwise provide defaults

        self.logFileName = None

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
    
    def makeInstructionsFrame(self):
        instructionsFrame = Tkinter.Frame(self)

        self.instructionText = Tkinter.StringVar()
        self.instructionText.set('Use this interface to look at a group of \nMinION Reads containing HLA amplicon sequences:')
        Tkinter.Label(instructionsFrame, width=80, height=5, textvariable=self.instructionText).pack()
        return instructionsFrame
        
        
    def makeChooseIODirectoriesFrame(self):
        chooseDirectoryFrame = Tkinter.Frame(self)        

        Tkinter.Button(chooseDirectoryFrame, text='Choose Input Directory ', command=self.chooseReadInputDirectory).grid(row=0, column=0, sticky=Tkconstants.W)
        Tkinter.Button(chooseDirectoryFrame, text='Choose Output Directory', command=self.chooseReadOutputDirectory).grid(row=1, column=0, sticky=Tkconstants.W)
        
        self.inputDirectoryText = Tkinter.StringVar()
        self.inputDirectoryText.set('Where is your MinION Read Directory?')
        Tkinter.Entry(chooseDirectoryFrame, width=60, textvariable=self.inputDirectoryText).grid(row=0, column=1)
                
        self.outputDirectoryText = Tkinter.StringVar()
        self.outputDirectoryText.set('Output Directory')
        Tkinter.Entry(chooseDirectoryFrame, width=60, textvariable=self.outputDirectoryText).grid(row=1, column=1)
        
        return chooseDirectoryFrame
        
        
    def makeAnalysisButtonsFrame(self):
        analysisButtonsFrame = Tkinter.Frame(self) 
        Tkinter.Button(analysisButtonsFrame, text='Initial Readstats', command=self.constructInitialReadStats).grid(row=0, column=0)
        Tkinter.Button(analysisButtonsFrame, text='Readstat Options', command=self.specifyReadStatOptions).grid(row=0, column=1)
        Tkinter.Button(analysisButtonsFrame, text='Run Full Analysis', command=self.runFullAnalysis).grid(row=1, column=0)
        Tkinter.Button(analysisButtonsFrame, text='Analysis Options', command=self.specifyAnalysisOptions).grid(row=1, column=1)
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
        Tkinter.Button(moreInfoFrame, text='How to use this tool', command=self.howToUse).grid(row=0, column=0)
        Tkinter.Button(moreInfoFrame, text='Contacting or Citing MUMC', command=self.contactInformation).grid(row=0, column=1)
        Tkinter.Button(moreInfoFrame, text='SaddleBags - A (Novel) Allele Submission Tool', command=self.launchSaddleBags).grid(row=0, column=2)
        return moreInfoFrame
        


        # Create a frame so people can choose the read directory
        #self.featureInputFrame = Tkinter.Frame(self)
        
        
    def launchSaddleBags(self):
        #tkMessageBox.showinfo('Need to open allelesub tool','Popup a window containing the allele submission tool please.')
        # TODO: This works but the interface is messed up.
        # I think saddlebags is using "self" to assign variables when it shouldnt. They should assign to the frame instead?
        saddleBagsRoot = Tkinter.Tk()
        AlleleGui(saddleBagsRoot).pack()
        saddleBagsRoot.mainloop()
        
    # This method should popup some instruction text in a wee window.
    # This should be explicit on how to use the tool.    
    def howToUse(self):
        tkMessageBox.showinfo('How to use this tool',
            'This software is to be used to create an\n'
            + 'EMBL-formatted submission document,\n'
            + 'which specifies a (novel) HLA allele.\n\n'       
                       
            + 'This tool requires you to submit a\n'
            + 'full length HLA allele, including\n'
            + '5\' and 3\' UTRs.\n\n'
            
            + 'Use capital letters for exons,\n'
            + 'lowercase for introns & UTRs.\n\n'
            
            + 'Push the "Example Sequence" button to see a small example of'
            + ' a formatted sequence.\n'
            + 'Sequences should follow this pattern:\n'
            + '5\'utr EX1 int1 EX2 ... EX{X} 3\'utr\n\n'
            
            + 'To use this tool:\n'
            + '1.) Fill in a Sample ID, Gene Name, and Allele.'
            + ' This text will be included in the submission.\n'
            + '2.) Paste your formatted sequence in the\n'
            + 'Annotated Sequence text area.\n'
            + '3.) Push \"Generate an EMBL submission\" button'
            + ' to generate a submission.\n'
            + '4.) Push the "Save the submission" button'
            + ' to store the submission on your computer.\nYou can submit this file to EMBL.\n\n'
            
            + 'All spaces, tabs, and newlines are'
            + ' removed before the nucleotide sequence is translated.'
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
            
            + 'This code will be hosted at:\n'
            + 'https://github.com/transplantation-\nimmunology/EMBL-HLA-Submission\n'
            + 'You will find more information on\n'
            + 'EMBL\'s data format on that page.'

            )

 
 
 
      
    # TODO: all of these analysis buttons.  Do something....  
    
    def constructInitialReadStats(self):
        self.logMessage('Calculating initial read stats')
        
        self.enableInterface(False)
        
        # Do we have all the options available?
        # What options are there?
        # Threads? No.
        # Barcoding?  How to get this info from the user?
        # Input and Output directories look good?
        # Minimum Read Length / Maximum Read Length
        # MinMax Quality

        # Run nit-picker for output directory
        # TODO: fix these parameters, especially sample ID. 
        # Get the quality filters from somewhere.
        # get the barcode location from somewhwere, or use a different path.
        
        #Quality :  23 <= AvgQ <= 32
        #Length: 2400 <= Length <= 4000
        # TODO: RUn this in a thread, so the GUI can update.

        nitPickerOutputDirectory = join(self.outputDirectoryText.get(), '1_prepared_reads')
        nit_picker.prepareReads(self.inputDirectoryText.get()
            , nitPickerOutputDirectory
            , 'TEMPSAMPLEID'
            , '/home/ben/Github/nanopore-prospector/nit-picker/barcodes/barcodes_96_bc_kit.txt'
            , 2400
            , 4000
            , 23
            , 32 )
            
        
        
        # Popup information about what we just did.
        
        # Are we in windows? 
        #If yes, open explorer with the read stats directory
        # If no, popup a dialog with the directory, user can do it themselves.
        
        self.enableInterface(True)
        
        self.logMessage('Done calculating initial read stats')
        

    def specifyReadStatOptions(self):
        print('Specifying ReadStat Options.....')
        
    def runFullAnalysis(self):
        print('Starting a full analysis.....')
                
    def specifyAnalysisOptions(self):
        print('Specifying Analysis Options.....')

    def enableInterface(self, enable):
        #self.submitButton.config(state='disabled')
         #self.submitButton.config(state='normal')
        if(enable):            
            # TODO: Something in this method, it does nothing right now.; 
            # disable the buttons and stuff, so people dont break stuff
            print('Enabling the user interface')
        else:
            print('Disabling the user interface')

    # chooseReadInputDirectory method is called when the user presses the input directory button
    def chooseReadInputDirectory(self):
        print ('Choosing an input directory.')
        self.enableInterface(False)
        self.setInputDir(tkFileDialog.askdirectory(**self.dir_opt))
        self.enableInterface(True)        
        
    # chooseReadInputDirectory method is called when the user presses the input directory button
    def chooseReadOutputDirectory(self):        
        print ('Choosing an output directory.')
        self.enableInterface(False)        
        self.setOutputDir(tkFileDialog.askdirectory(**self.dir_opt))
        self.enableInterface(True)
        

    # setInputDir is a subprocess to be called when an input directory is selected.
    def setInputDir(self, currentInputDirectory):
        self.inputDirectoryText.set(currentInputDirectory)
        
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
