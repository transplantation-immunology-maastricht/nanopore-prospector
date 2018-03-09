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
 --reads="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_All_Features_A.fasta" \
 --outputdir="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/A_cDNA_Polymorphism_Analysis" \
 --reference="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_A_Reference.fasta"
 
python nit_picker_main.py \
 --reads="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_All_Features_B.fasta" \
 --outputdir="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/B_cDNA_Polymorphism_Analysis" \
 --reference="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_B_Reference.fasta"
 
python nit_picker_main.py \
 --reads="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_All_Features_C.fasta" \
 --outputdir="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/C_cDNA_Polymorphism_Analysis" \
 --reference="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_C_Reference.fasta"

python nit_picker_main.py \
 --reads="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_Alleles_Full_Length_A.fasta" \
 --outputdir="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/A_Polymorphism_Analysis" \
 --reference="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_A_Reference.fasta"
 
python nit_picker_main.py \
 --reads="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_Alleles_Full_Length_B.fasta" \
 --outputdir="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/B_Polymorphism_Analysis" \
 --reference="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_B_Reference.fasta"
 
python nit_picker_main.py \
 --reads="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_Alleles_Full_Length_C.fasta" \
 --outputdir="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/C_Polymorphism_Analysis" \
 --reference="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_C_Reference.fasta"
 
