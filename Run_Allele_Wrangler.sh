# This file is part of Allele-Wrangler.
#
# Allele-Wrangler is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Allele-Wrangler is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Allele-Wrangler. If not, see <http://www.gnu.org/licenses/>.

# Version 1.0 

# See the file README.MD for how to set up your anaconda environment.

# You should change the following input variables before you run this program:
# The $Rundate variable is included in the name of the output files, 
# You can fill that variable with whatever text you wish to appear there
# Don't use special characters ("*", "/", "\", "$", "?", etc.)


ReadInputFile="/home/minion/MinIONData/2017.EfiAbstractProject/43785/Analysis/blast_results/SortedReads/HLA-C_15.fastq"
ResultsOutputDirectory="/home/minion/MinIONData/43785_C15_wrangler_test"
#NumberIterations="5"

#ReadInputFile="/home/eggs/Workspace/Github/allele-wrangler/data/EvenSmallerReads.fasta"
#ReadInputFile="/home/eggs/Workspace/Github/allele-wrangler/data/ToyReads.fasta"
#ResultsOutputDirectory="/home/eggs/Workspace/wrangled_results_toy"

#ReadInputFile="/home/eggs/Workspace/Data/SampleReads_extracts/SampleReads_BC01_TwoDir_reads.fasta"


#ReadInputFile="/home/eggs/Workspace/Data/HLA-B_07_FewReads.fastq"
#ResultsOutputDirectory="/home/eggs/Workspace/Data/HLA-B_07_FewReads_Assembled"

#ReadInputFile="/home/eggs/Workspace/Data/43785/HLA-B_07_08.fastq"
#ResultsOutputDirectory="/home/eggs/Workspace/Data/HLA-B_07_08_Combined"


NumberIterations="8"
ThreadCount="4"

source activate minionvironment

cd src

python AlleleWrangler_Main.py \
 --reads=$ReadInputFile \
 --outputdir=$ResultsOutputDirectory \
 --iterations=$NumberIterations \
 --threads=$ThreadCount

source deactivate
