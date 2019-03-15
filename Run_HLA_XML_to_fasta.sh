# This file is part of HLA-allele-analysis.
#
# HLA-allele-analysis is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HLA-allele-analysis is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with HLA-allele-analysis. If not, see <http://www.gnu.org/licenses/>.

source /home/ben/minionvenv/bin/activate

InputFile="/home/ben/Github/nanopore_prospector/Sample_Data/input/hla_xml.3.35.0/hla.xml"
OutputDirectory="/home/ben/Github/nanopore_prospector/Sample_Data/output/hla_alleles"

python Nanopore_Prospector_Main.py \
 --inputfile=$InputFile \
 --outputdirectory=$OutputDirectory \
 --action="hlaalleleanalysis"
