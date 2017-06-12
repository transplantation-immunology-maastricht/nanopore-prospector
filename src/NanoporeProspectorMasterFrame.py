
import os
from os import listdir
from os.path import isfile, join, split
import Tkinter, Tkconstants, tkFileDialog, tkMessageBox
from Tkconstants import CURRENT


class NanoporeProspectorMasterFrame(Tkinter.Frame):
    def __init__(self, root):


        Tkinter.Frame.__init__(self, root)
        root.title("nanopore prospector")

        self.parent = root

        
        self.initialize()

    # Initialize GUI elements
    def initialize(self):

        
        
        # This is the directory the python executable is running from.
        FileAndPath = os.path.abspath(__file__)
        self.idir, self.ifile = os.path.split(FileAndPath)

        # GUI options
        self.label_opt = {'fill': Tkconstants.BOTH, 'padx': 10, 'pady': 10}
        self.button_opt = {'fill': Tkconstants.BOTH, 'padx': 50, 'pady': 15}
        self.dir_opt = {'initialdir': self.idir, 'mustexist': True, 
            'parent': self, 'title': 'Choose a directory'}
        
        
        
        # A frame for the instructions
        self.instructionsFrame = Tkinter.Frame(self)

        self.instructionText = Tkinter.StringVar()
        self.instructionText.set('Use this interface to look at a group of \nMinION Reads containing HLA amplicon sequences:')
        Tkinter.Label(self.instructionsFrame, width=80, height=5, textvariable=self.instructionText).pack()
        self.instructionsFrame.pack()
        
        # A frame for the input directory
        self.inputDirectoryFrame = Tkinter.Frame(self)
        
        Tkinter.Button(self, 'Choose Input Directory', command=self.chooseReadInputDirectory()).pack(**self.button_opt)
        self.inputDirectoryText = Tkinter.StringVar()
        self.inputDirectoryText.set('C:/folder//directory???')
        Tkinter.Label(self.inputDirectoryFrame, width=80, height=1, textvariable=self.inputDirectoryText).pack()
        self.inputDirectoryFrame.pack()
        
        
        # Assemble all the frames and get ready.
        self.pack()


        # Create a frame so people can choose the read directory
        #self.featureInputFrame = Tkinter.Frame(self)
        
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
        
        self.enableInterface(False)
        #self.submitButton.config(state='disabled')
        self.setInputDir(tkFileDialog.askdirectory(**self.dir_opt))
        self.enableInterface(False)
        #s#elf.submitButton.config(state='normal')

    # setInputDir is a subprocess to be called when an input directory is selected.
    def setInputDir(self, currentFast5Dir):
        self.inputDirectoryText= currentFast5Dir
        #currentFast5Dir = os.path.normpath(currentFast5Dir)

        # The default output directory is called 'extracted_reads'.
        # This folder should be at the same directory tree level as the input folder,
        # To allow the barcoded reads to be in the same directory.
        #outputDirectory = os.path.normpath(os.path.join(
        #    os.path.normpath(currentFast5Dir + os.sep + os.pardir)
         #   ,'extracted_reads'))

        # Open an extractor for the root directory. 
        # This extractor is run whether or not there are subfolders. 
        #self.instructionText.set('Extracting sequences from this root directory:\n'
        #    + currentFast5Dir)         
        #tempDir, tempFileName = os.path.split(currentFast5Dir)  
        #sampleID = self.getSampleID(tempFileName) 
        #if not (sampleID == 'SKIPFOLDER'):      
        #    self.openSubFrame(currentFast5Dir, outputDirectory, sampleID)
        #self.defaultMessage()     
