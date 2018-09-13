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

source activate minionvironment

python CalculatePolymorphism.py \
 --input="/home/ben/MUMCScripts/CalculatePolymorphism/input/HLA_A_IN4.fasta" \
 --output="/home/ben/MUMCScripts/CalculatePolymorphism/results/HLA_A_IN4_Polymorphism_Results.csv"
 
 python CalculatePolymorphism.py \
 --input="/home/ben/MUMCScripts/CalculatePolymorphism/input/HLA_B_IN4.fasta" \
 --output="/home/ben/MUMCScripts/CalculatePolymorphism/results/HLA_B_IN4_Polymorphism_Results.csv"
 
 python CalculatePolymorphism.py \
 --input="/home/ben/MUMCScripts/CalculatePolymorphism/input/HLA_C_IN4.fasta" \
 --output="/home/ben/MUMCScripts/CalculatePolymorphism/results/HLA_C_IN4_Polymorphism_Results.csv"

python CalculatePolymorphism.py \
 --input="/home/ben/MUMCScripts/CalculatePolymorphism/input/HLA_A_EX3.fasta" \
 --output="/home/ben/MUMCScripts/CalculatePolymorphism/results/HLA_A_IEX3_Polymorphism_Results.csv"
 
 python CalculatePolymorphism.py \
 --input="/home/ben/MUMCScripts/CalculatePolymorphism/input/HLA_B_EX3.fasta" \
 --output="/home/ben/MUMCScripts/CalculatePolymorphism/results/HLA_B_EX3_Polymorphism_Results.csv"
 
 python CalculatePolymorphism.py \
 --input="/home/ben/MUMCScripts/CalculatePolymorphism/input/HLA_C_EX3.fasta" \
 --output="/home/ben/MUMCScripts/CalculatePolymorphism/results/HLA_C_EX3_Polymorphism_Results.csv"

source deactivate
