# nanopore-prospector
GUI Interface for analysis of MinION nanopore reads of HLA amplicon sequences

To setup environment:
Install python 3.6
https://askubuntu.com/questions/865554/how-do-i-install-python-3-6-using-apt-get

Create virtual environment
virtualenv -p python3.6 minionvenv

activate virtual environment:
source /home/ben/minionvenv/bin/activate

blast and clustalo.
Blast:
https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download

sudo apt-get install clustalo

numpy:
pip install numpy

biopython
pip install biopython

pylab / matplotlib
pip install matplotlib

skikit learn for kmeans clustering:
pip install -U scikit-learn

pip install scipy

pip install pysam

# Example Uses

Run the main user interface:

bash Run_Nanopore_Prospector.sh

Analyze Read Quality:

bash Run_Analyze_Read_Quality.sh

Run Prospector from command line:

bash Run_Prospector_CL.sh
