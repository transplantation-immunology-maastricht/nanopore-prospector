# This file is part of nanopore-prospector.
#
# nanopore-prospector is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# nanopore-prospector is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with nanopore-prospector. If not, see <http://www.gnu.org/licenses/>.

import os

import Tkinter, Tkconstants
#, tkFileDialog, tkMessageBox
from Tkinter import *


class AnalysisOptionsInputForm(Tkinter.Frame):
        
    # Initialize the GUI
    def __init__(self, root, prospectorFormObject):
        Tkinter.Frame.__init__(self, root)
        root.title("Choose Read Stat and Demultiplexing options")
        self.parent = root
        
        # keep a handle on the parent object. I can specify configuration using this handle
        self.prospectorFormObject = prospectorFormObject

        button_opt = {'fill': Tkconstants.BOTH, 'padx': 35, 'pady': 5}
        
        # To define the exit behavior.  Save and exit.
        self.parent.protocol('WM_DELETE_WINDOW', self.saveOptions)
        
        # Define the return behavior.  Same as "close window" etc
        root.bind('<Return>', self.returnFunction)
  
        # This window should not be resizeable. I guess.
        self.parent.resizable(width=False, height=False)
        
        
        
        #Standard Inputs widths for the form elements
        formInputWidth = 10
        labelInputWidth = 60
        
        # Instructions
        self.instructionsFrame = Tkinter.Frame(self)  
        self.instructionText = Tkinter.StringVar()       
        self.instructionText.set('\nOptions: Read statistics and Demultiplexing.\n')        
        Tkinter.Label(self.instructionsFrame, width=85, height=3, textvariable=self.instructionText).pack()
        self.instructionsFrame.pack()
        
        
        # Demultiplexing
        self.demultiplexOptionFrame = Tkinter.Frame(self)        
        self.demultiplexInstrStringVar = Tkinter.StringVar()
        self.demultiplexInstrStringVar.set('Demultiplex the reads, using 96x barcode kit:')
        self.demultiplexInstrLabel = Tkinter.Label(self.demultiplexOptionFrame, width=labelInputWidth, height=1, textvariable=self.demultiplexInstrStringVar).grid(row=0, column=0)           
        self.chooseDemultiplexIntVar = IntVar()
        #self.chooseDemultiplexIntVar.set(2)
        Radiobutton(self.demultiplexOptionFrame, text="Yes", variable=self.chooseDemultiplexIntVar, value=1).grid(row=0, column=1)
        Radiobutton(self.demultiplexOptionFrame, text="No", variable=self.chooseDemultiplexIntVar, value=2).grid(row=0, column=2)  
        self.demultiplexOptionFrame.pack()
        
        # Length and Quality
        self.lengthParamsFrame = Tkinter.Frame(self)          
        self.lengthParamsStringVar = Tkinter.StringVar()
        self.lengthParamsStringVar.set('Minimum & Maximum Read Length:')
        self.demultiplexInstrLabel = Tkinter.Label(self.lengthParamsFrame, width=labelInputWidth, height=1, textvariable=self.lengthParamsStringVar).grid(row=0, column=0) 
        self.inputMinLength = Tkinter.StringVar()
        #self.inputMinLength.set('0')
        self.inputMinLengthEntry = Tkinter.Entry(self.lengthParamsFrame, width=formInputWidth, textvariable=self.inputMinLength).grid(row=0, column=1)
        self.inputMaxLength = Tkinter.StringVar()
        #self.inputMaxLength.set('999999')
        self.inputMaxLengthEntry = Tkinter.Entry(self.lengthParamsFrame, width=formInputWidth, textvariable=self.inputMaxLength).grid(row=0, column=2)
        self.qualityParamsStringVar = Tkinter.StringVar()
        self.qualityParamsStringVar.set('Minimum & Maximum Quality (Q-Score):')
        self.demultiplexInstrLabel = Tkinter.Label(self.lengthParamsFrame, width=labelInputWidth, height=1, textvariable=self.qualityParamsStringVar).grid(row=1, column=0)   
        self.inputMinQuality = Tkinter.StringVar()
        #self.inputMinQuality.set('0')
        self.inputMinLengthEntry = Tkinter.Entry(self.lengthParamsFrame, width=formInputWidth, textvariable=self.inputMinQuality).grid(row=1, column=1)
        self.inputMaxQuality = Tkinter.StringVar()
        #self.inputMaxQuality.set('99')
        self.inputMaxLengthEntry = Tkinter.Entry(self.lengthParamsFrame, width=formInputWidth, textvariable=self.inputMaxQuality).grid(row=1, column=2)
        self.lengthParamsFrame.pack()
        
        # Make a frame for the save options button.
        self.saveOptionsFrame = Tkinter.Frame(self)
        Tkinter.Button(self.saveOptionsFrame, text='Save Options', command=self.saveOptions).grid(row=0, column=0)
        self.saveOptionsFrame.pack()


    # I needed a function for the return keypress to latch onto.
    # It is just a wrapper for the saveOptions method.
    def returnFunction(self, event):
        self.saveOptions()
 

    def saveOptions(self):
        
        print ('Saving Options')

        if (self.chooseDemultiplexIntVar.get() == 1):
            self.prospectorFormObject.demultiplexReads = True
        else:
            self.prospectorFormObject.demultiplexReads = False
           
        self.prospectorFormObject.minimumReadLength = int(self.inputMinLength.get())
        self.prospectorFormObject.maximumReadLength = int(self.inputMaxLength.get())
        
        self.prospectorFormObject.minimumQuality = int(self.inputMinQuality.get())
        self.prospectorFormObject.maximumQuality = int(self.inputMaxQuality.get())     
        
        self.parent.destroy() 
 
    
    def closeWindow(self):
        self.parent.destroy()        
    