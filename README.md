# nanopore-prospector
GUI Interface for analysis of MinION nanopore reads of HLA amplicon sequences

Note: This project uses submodules to load code from related repos. It is necessary to download submodule code, in addition to the git clone instructions.

Use git clone --recursive https://github.com/transplantation-immunology/nanopore-prospector.git

More info 

https://stackoverflow.com/questions/8090761/pull-using-git-including-submodule

http://www.vogella.com/tutorials/GitSubmodules/article.html

To update submodule code:
git submodule update --recursive --remote

to add a submodule to project
git submodule add -b master [URL to Git repo] 
git submodule init 
git submodule update --remote



## To configure Anaconda
Anaconda uses separate environments to run your programs in.  
Install Anaconda for python 2.7.  
https://www.continuum.io/downloads  
To set up the environment in anaconda:  

Linux/Mac:  
```
conda create --name minionvironment biopython six pycurl pysam
source activate minionvironment  
pip install pyinstaller packaging matplotlib
source deactivate  
```  
Windows:  
```  
conda create --name minionvironment biopython six pycurl pywin32 pysam  
call activate minionvironment && pip install pyinstaller packaging matplotlib && call deactivate  
```



Also, you must install other things:
pywin32 (for installing in windows.  Available in Conda.)
pylab.  This is matplotlib. Pip can install this inside environment.

BLAST
```  
Windows:  
```  
https://www.ncbi.nlm.nih.gov/books/NBK52637/

aligners.  Don't I need to install thes? or no?
clustalo
bwalign


pysam must be upgraded using pip.  
pip install pysam --upgrade