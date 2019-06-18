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

cd /home/ben/Github/nanopore_prospector


python Nanopore_Prospector_Main.py \
 --reads="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_Alleles_Combined_A.fasta" \
 --outputdirectory="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/A_Polymorphism_Analysis_Combined" \
 --reference="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_A_Reference.fasta" \
 --action="preparereads" \
 --snpcutoff=".00000001"

python Nanopore_Prospector_Main.py \
 --reads="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_Alleles_Combined_B.fasta" \
 --outputdirectory="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/B_Polymorphism_Analysis_Combined" \
 --reference="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_B_Reference.fasta" \
 --action="preparereads" \
 --snpcutoff=".00000001"

python Nanopore_Prospector_Main.py \
 --reads="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_Alleles_Combined_C.fasta" \
 --outputdirectory="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/C_Polymorphism_Analysis_Combined" \
 --reference="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_C_Reference.fasta" \
 --action="preparereads" \
 --snpcutoff=".00000001"

python Nanopore_Prospector_Main.py \
 --reads="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_Alleles_Full_Length_A.fasta" \
 --outputdirectory="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/A_Polymorphism_Analysis_Full_Length" \
 --reference="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_A_Reference.fasta" \
 --action="preparereads" \
 --snpcutoff=".00000001"

python Nanopore_Prospector_Main.py \
 --reads="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_Alleles_Full_Length_B.fasta" \
 --outputdirectory="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/B_Polymorphism_Analysis_Full_Length" \
 --reference="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_B_Reference.fasta" \
 --action="preparereads" \
 --snpcutoff=".00000001"

python Nanopore_Prospector_Main.py \
 --reads="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_Alleles_Full_Length_C.fasta" \
 --outputdirectory="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/C_Polymorphism_Analysis_Full_Length" \
 --reference="/home/ben/Github/nanopore_prospector/HLA_Polymorphism_Analysis/HLA_C_Reference.fasta" \
 --action="preparereads" \
 --snpcutoff=".00000001"


