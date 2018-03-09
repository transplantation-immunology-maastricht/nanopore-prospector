# This file is part of nit_picker.
#
# nit_picker is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# nit_picker is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with nit_picker. If not, see <http://www.gnu.org/licenses/>.


# TODO change these to relative paths. Probably easy.

source /home/ben/minionvenv/bin/activate

ReadDirectory="/home/ben/Github/nanopore_prospector/33332_analyze_read_quality/A24_Reads"
ResultsOutputDirectory="/home/ben/Github/nanopore_prospector/33332_analyze_read_quality/33332_A24_nitpicker_analysis"
ReferenceFile="/home/ben/Github/nanopore_prospector/33332_analyze_read_quality/HLA_Reference_33332_A24.fasta"
BarcodeFile="/home/ben/Github/nanopore_prospector/nit_picker/barcodes/barcodes_96_bc_kit.txt"

python nit_picker_main.py \
 --reads=$ReadDirectory \
 --outputdir=$ResultsOutputDirectory \
 --reference=$ReferenceFile \
 --minqual=34 \
 --minlen=3100 \
 --maxlen=3300
