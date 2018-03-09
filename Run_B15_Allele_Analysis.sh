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

# TODO: Giving a reference file instructs nitpicker to calculate read stats.
# make CalculateReadStats a commandline option, which requires a reference sequence.

python nit_picker_main.py \
 --reads="/home/ben/Github/analyze-polymorphic-regions/B15_Alleles/nitpicker_analysis/AllAlleles" \
 --outputdir="/home/ben/Github/analyze-polymorphic-regions/B15_Alleles/nitpicker_analysis/AllAlleles_Analysis" \
 --reference="/home/ben/Github/analyze-polymorphic-regions/B15_Alleles/nitpicker_analysis/B15_01_Reference.fasta"
 
python nit_picker_main.py \
 --reads="/home/ben/Github/analyze-polymorphic-regions/B15_Alleles/nitpicker_analysis/B62" \
 --outputdir="/home/ben/Github/analyze-polymorphic-regions/B15_Alleles/nitpicker_analysis/B62_Analysis" \
 --reference="/home/ben/Github/analyze-polymorphic-regions/B15_Alleles/nitpicker_analysis/B15_01_Reference.fasta"
 
python nit_picker_main.py \
 --reads="/home/ben/Github/analyze-polymorphic-regions/B15_Alleles/nitpicker_analysis/B63" \
 --outputdir="/home/ben/Github/analyze-polymorphic-regions/B15_Alleles/nitpicker_analysis/B63_Analysis" \
 --reference="/home/ben/Github/analyze-polymorphic-regions/B15_Alleles/nitpicker_analysis/B15_01_Reference.fasta"
 
python nit_picker_main.py \
 --reads="/home/ben/Github/analyze-polymorphic-regions/B15_Alleles/nitpicker_analysis/B70" \
 --outputdir="/home/ben/Github/analyze-polymorphic-regions/B15_Alleles/nitpicker_analysis/B70_Analysis" \
 --reference="/home/ben/Github/analyze-polymorphic-regions/B15_Alleles/nitpicker_analysis/B15_01_Reference.fasta"
 
python nit_picker_main.py \
 --reads="/home/ben/Github/analyze-polymorphic-regions/B15_Alleles/nitpicker_analysis/B75" \
 --outputdir="/home/ben/Github/analyze-polymorphic-regions/B15_Alleles/nitpicker_analysis/B75_Analysis" \
 --reference="/home/ben/Github/analyze-polymorphic-regions/B15_Alleles/nitpicker_analysis/B15_01_Reference.fasta"
 
python nit_picker_main.py \
 --reads="/home/ben/Github/analyze-polymorphic-regions/B15_Alleles/nitpicker_analysis/B76" \
 --outputdir="/home/ben/Github/analyze-polymorphic-regions/B15_Alleles/nitpicker_analysis/B76_Analysis" \
 --reference="/home/ben/Github/analyze-polymorphic-regions/B15_Alleles/nitpicker_analysis/B15_01_Reference.fasta"
 