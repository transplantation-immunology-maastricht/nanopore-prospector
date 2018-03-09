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


from tkinter import Frame, StringVar, Label, Entry, Button, IntVar, Radiobutton, Checkbutton
from tkinter.constants import BOTH, W

from nanopore_prospector.common import assignConfigurationValue, getConfigurationValue, writeConfigurationFile, loadConfigurationFile

class AnalysisOptionsInputForm(Frame):
        
    # Initialize the GUI
    def __init__(self, root, prospectorFormObject):
        Frame.__init__(self, root)
        root.title("Choose Read Stat and Demultiplexing options")
        self.parent = root
        
        # keep a handle on the parent object. I can specify configuration using this handle
        self.prospectorFormObject = prospectorFormObject
        
        # To define the exit behavior.  Save and exit.
        self.parent.protocol('WM_DELETE_WINDOW', self.saveOptions)
        
        # Define the return behavior.  Same as "close window" etc
        root.bind('<Return>', self.returnFunction)
  
        # This window should not be resizeable. I guess.
        self.parent.resizable(width=False, height=False)
     
        #Standard Inputs widths for the form elements
        self.formInputWidth = 10
        self.labelInputWidth = 60
        
        self.makeAnalysisGeneralOptionsFrame()
        self.makeDemultiplexInstructionsFrame()
        self.makeDemultiplexOptionsFrame()
        self.makeLengthAndQualityFilterFrame()
        self.makeSaveOptionsFrame()
        self.makeBlastInstructionsFrame()
        self.makeBlastOptionsFrame()
        
        self.loadOptions()

    def makeAnalysisGeneralOptionsFrame(self):
        self.generalOptionsFrame = Frame(self)
        # TODO: Add parameters for # of threads.
        # Other generic parameters?
        
        # TODO Should I worry  about the # of ram the user is using? I have never ran into the limit before.
        self.generalOptionsFrame.pack()

    def makeDemultiplexInstructionsFrame(self):
        self.instructionsFrame = Frame(self)  
        self.instructionText = StringVar()       
        self.instructionText.set('\nOptions: Read statistics and Demultiplexing.\n')        
        Label(self.instructionsFrame, width=85, height=3, textvariable=self.instructionText).pack()
        self.instructionsFrame.pack()

    def makeDemultiplexOptionsFrame(self):
        # Demultiplexing
        self.demultiplexOptionFrame = Frame(self)        
        self.demultiplexInstrStringVar = StringVar()
        self.demultiplexInstrStringVar.set('Demultiplex the reads, using 96x barcode kit:')
        self.demultiplexInstrLabel = Label(self.demultiplexOptionFrame, width=self.labelInputWidth, height=1, textvariable=self.demultiplexInstrStringVar).grid(row=0, column=0)           
        self.chooseDemultiplexIntVar = IntVar()
        #self.chooseDemultiplexIntVar.set(2)
        Radiobutton(self.demultiplexOptionFrame, text="Yes", variable=self.chooseDemultiplexIntVar, value=1).grid(row=0, column=1)
        Radiobutton(self.demultiplexOptionFrame, text="No", variable=self.chooseDemultiplexIntVar, value=2).grid(row=0, column=2)  
        self.demultiplexOptionFrame.pack()

    def makeLengthAndQualityFilterFrame(self):
        # Length and Quality
        self.lengthParamsFrame = Frame(self)          
        self.lengthParamsStringVar = StringVar()
        self.lengthParamsStringVar.set('Minimum & Maximum Read Length:')
        self.demultiplexInstrLabel = Label(self.lengthParamsFrame, width=self.labelInputWidth, height=1, textvariable=self.lengthParamsStringVar).grid(row=0, column=0) 
        self.inputMinLength = StringVar()
        self.inputMinLengthEntry = Entry(self.lengthParamsFrame, width=self.formInputWidth, textvariable=self.inputMinLength).grid(row=0, column=1)
        self.inputMaxLength = StringVar()
        self.inputMaxLengthEntry = Entry(self.lengthParamsFrame, width=self.formInputWidth, textvariable=self.inputMaxLength).grid(row=0, column=2)
        self.qualityParamsStringVar = StringVar()
        self.qualityParamsStringVar.set('Minimum & Maximum Quality (Q-Score):')
        self.demultiplexInstrLabel = Label(self.lengthParamsFrame, width=self.labelInputWidth, height=1, textvariable=self.qualityParamsStringVar).grid(row=1, column=0)   
        self.inputMinQuality = StringVar()
        self.inputMinLengthEntry = Entry(self.lengthParamsFrame, width=self.formInputWidth, textvariable=self.inputMinQuality).grid(row=1, column=1)
        self.inputMaxQuality = StringVar()
        self.inputMaxLengthEntry = Entry(self.lengthParamsFrame, width=self.formInputWidth, textvariable=self.inputMaxQuality).grid(row=1, column=2)
        self.lengthParamsFrame.pack()
        
    def makeBlastInstructionsFrame(self):
        self.blastInstructionsFrame = Frame(self)  
        self.blastInstructionText = StringVar()       
        self.blastInstructionText.set('\nOptions: Blast Sorting.\n'
            + 'What genes do you want to sort against?\n'
            + '(More genes = more analysis time)')        
        Label(self.blastInstructionsFrame, width=85, height=5, textvariable=self.blastInstructionText).pack()
        self.blastInstructionsFrame.pack()    
        
    def makeBlastOptionsFrame(self):
        self.blastOptionsFrame = Frame(self)  
        
        # TODO: I have been adding genes as we use them. Should I add more?
        # HLA Allele Analysis makes these reference files now. 
        # I can add DMA DMB DOA DOB DPA1 DPB1 DQA1 DQB1 DRA DRB1 DRB3 DRB4 DRB5 
        # And E and F
        
        self.geneAIntVar = IntVar()
        Checkbutton(self.blastOptionsFrame, text="HLA-A", variable=self.geneAIntVar).grid(row=0, sticky=W)
        
        self.geneBIntVar = IntVar()
        Checkbutton(self.blastOptionsFrame, text="HLA-B", variable=self.geneBIntVar).grid(row=1, sticky=W)
        
        self.geneCIntVar = IntVar()
        Checkbutton(self.blastOptionsFrame, text="HLA-C", variable=self.geneCIntVar).grid(row=2, sticky=W)
        
        self.geneEIntVar = IntVar()
        Checkbutton(self.blastOptionsFrame, text="HLA-E", variable=self.geneEIntVar).grid(row=3, sticky=W)
        
        self.geneDRAIntVar = IntVar()
        Checkbutton(self.blastOptionsFrame, text="HLA-DRA", variable=self.geneDRAIntVar).grid(row=4, sticky=W)

        self.geneDQA1IntVar = IntVar()
        Checkbutton(self.blastOptionsFrame, text="HLA-DQA1", variable=self.geneDQA1IntVar).grid(row=5, sticky=W)
        
        self.geneDQB1IntVar = IntVar()
        Checkbutton(self.blastOptionsFrame, text="HLA-DQB1", variable=self.geneDQB1IntVar).grid(row=6, sticky=W)
        
        self.geneDRB1IntVar = IntVar()
        Checkbutton(self.blastOptionsFrame, text="HLA-DRB1", variable=self.geneDRB1IntVar).grid(row=7, sticky=W)

        self.blastOptionsFrame.pack()
                
                
    def makeSaveOptionsFrame(self):  
        # Make a frame for the save options button.
        self.saveOptionsFrame = Frame(self)
        Button(self.saveOptionsFrame, text='Save Options', command=self.saveOptions).grid(row=0, column=0)
        self.saveOptionsFrame.pack()
        
    # I needed a function for the return keypress to latch onto.
    # It is just a wrapper for the saveOptions method.
    def returnFunction(self, event):
        self.saveOptions()
 
 
    def loadOptions(self):
        # Read the configuration file.
        loadConfigurationFile()

        #TODO: Number of threads. Configure this.
        
        if getConfigurationValue('demultiplex_reads') is not None:
            self.chooseDemultiplexIntVar.set(int(getConfigurationValue('demultiplex_reads')))
        
        if getConfigurationValue('min_length') is not None:
            self.inputMinLength.set(getConfigurationValue('min_length'))            
        if getConfigurationValue('max_length') is not None:
            self.inputMaxLength.set(getConfigurationValue('max_length'))
        
        if getConfigurationValue('min_quality') is not None:
            self.inputMinQuality.set(getConfigurationValue('min_quality'))            
        if getConfigurationValue('max_quality') is not None:
            self.inputMaxQuality.set(getConfigurationValue('max_quality'))
            
        if getConfigurationValue('analyze_hla_a') is not None:
            self.geneAIntVar.set(getConfigurationValue('analyze_hla_a'))
        if getConfigurationValue('analyze_hla_b') is not None:
            self.geneBIntVar.set(getConfigurationValue('analyze_hla_b'))
        if getConfigurationValue('analyze_hla_c') is not None:
            self.geneCIntVar.set(getConfigurationValue('analyze_hla_c'))
        if getConfigurationValue('analyze_hla_e') is not None:
            self.geneEIntVar.set(getConfigurationValue('analyze_hla_e'))
            
        if getConfigurationValue('analyze_hla_dra') is not None:
            self.geneDRAIntVar.set(getConfigurationValue('analyze_hla_dra'))
        if getConfigurationValue('analyze_hla_dqa1') is not None:
            self.geneDQA1IntVar.set(getConfigurationValue('analyze_hla_dqa1'))
        if getConfigurationValue('analyze_hla_dqb1') is not None:
            self.geneDQB1IntVar.set(getConfigurationValue('analyze_hla_dqb1'))
        if getConfigurationValue('analyze_hla_drb1') is not None:
            self.geneDRB1IntVar.set(getConfigurationValue('analyze_hla_drb1'))
            
    def saveOptions(self):
        
        print ('Saving Options')
        
        #TODO: Number of threads. Configure this.

        assignConfigurationValue('demultiplex_reads', str(self.chooseDemultiplexIntVar.get()))
           
        assignConfigurationValue('min_length', self.inputMinLength.get())
        assignConfigurationValue('max_length', self.inputMaxLength.get())
        
        assignConfigurationValue('min_quality', self.inputMinQuality.get())
        assignConfigurationValue('max_quality', self.inputMaxQuality.get())
        
        assignConfigurationValue('analyze_hla_a', (self.geneAIntVar.get()))
        assignConfigurationValue('analyze_hla_b', (self.geneBIntVar.get()))
        assignConfigurationValue('analyze_hla_c', (self.geneCIntVar.get()))
        assignConfigurationValue('analyze_hla_e', (self.geneEIntVar.get()))
        assignConfigurationValue('analyze_hla_dra', (self.geneDRAIntVar.get()))
        assignConfigurationValue('analyze_hla_dqa1', (self.geneDQA1IntVar.get()))
        assignConfigurationValue('analyze_hla_dqb1', (self.geneDQB1IntVar.get()))
        assignConfigurationValue('analyze_hla_drb1', (self.geneDRB1IntVar.get()))

        writeConfigurationFile()
        self.parent.destroy() 
        
        
    