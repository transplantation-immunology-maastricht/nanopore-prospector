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

# Version 1.0

# See the file README.MD for how to set up your anaconda environment.

source /home/ben/minionvenv/bin/activate

cd /home/ben/Github/nanopore_prospector


ReadFile="/MinIONData/Diagnostics.DataQualityAnalysis/M190905A/BC13.LocusSorted/BC13/HLA-A_Reads.fastq"
ResultsOutputDirectory="/MinIONData/Diagnostics.DataQualityAnalysis/M190905A/A_HeteroSplit"
ReferenceFile="/MinIONData/Diagnostics.DataQualityAnalysis/M190905A/HLAReference/A_Reference.fasta"

python Nanopore_Prospector_Main.py \
 --reads=$ReadFile \
 --outputdir=$ResultsOutputDirectory \
 --reference=$ReferenceFile \
 --action="heterosplit"