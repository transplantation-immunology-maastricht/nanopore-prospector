source /home/ben/minionvenv/bin/activate

inputFile="/home/ben/Github/nanopore_prospector/Sample_Data/output/hla_alleles/HLA_Allele_FullLen_Minus_UTRs.fasta"
outputDirectory="/home/ben/Github/nanopore_prospector/Sample_Data/output/homopolymer_results"

python Nanopore_Prospector_Main.py \
 --input=$inputFile \
 --output=$outputDirectory \
 --action="findhomopolymers"


# To use a block comment: https://stackoverflow.com/questions/947897/block-comments-in-a-shell-script
: <<'ENDBLOCK'


inputFile="/home/ben/Github/HLA-allele-analysis/OutputDirectory/5UTR_sequences/HLA_5UTR.fasta"
outputDirectory="/home/ben/Github/find_homopolymers/output"
python find_homopolymers_main.py \
 --input=$inputFile \
 --output=$outputDirectory

inputFile="/home/ben/Github/HLA-allele-analysis/OutputDirectory/EX1_sequences/HLA_EX1.fasta"
outputDirectory="/home/ben/Github/find_homopolymers/output"
python find_homopolymers_main.py \
 --input=$inputFile \
 --output=$outputDirectory

inputFile="/home/ben/Github/HLA-allele-analysis/OutputDirectory/IN1_sequences/HLA_IN1.fasta"
outputDirectory="/home/ben/Github/find_homopolymers/output"
python find_homopolymers_main.py \
 --input=$inputFile \
 --output=$outputDirectory

inputFile="/home/ben/Github/HLA-allele-analysis/OutputDirectory/EX2_sequences/HLA_EX2.fasta"
outputDirectory="/home/ben/Github/find_homopolymers/output"
python find_homopolymers_main.py \
 --input=$inputFile \
 --output=$outputDirectory

inputFile="/home/ben/Github/HLA-allele-analysis/OutputDirectory/IN2_sequences/HLA_IN2.fasta"
outputDirectory="/home/ben/Github/find_homopolymers/output"
python find_homopolymers_main.py \
 --input=$inputFile \
 --output=$outputDirectory

inputFile="/home/ben/Github/HLA-allele-analysis/OutputDirectory/EX3_sequences/HLA_EX3.fasta"
outputDirectory="/home/ben/Github/find_homopolymers/output"
python find_homopolymers_main.py \
 --input=$inputFile \
 --output=$outputDirectory

inputFile="/home/ben/Github/HLA-allele-analysis/OutputDirectory/IN3_sequences/HLA_IN3.fasta"
outputDirectory="/home/ben/Github/find_homopolymers/output"
python find_homopolymers_main.py \
 --input=$inputFile \
 --output=$outputDirectory

inputFile="/home/ben/Github/HLA-allele-analysis/OutputDirectory/EX4_sequences/HLA_EX4.fasta"
outputDirectory="/home/ben/Github/find_homopolymers/output"
python find_homopolymers_main.py \
 --input=$inputFile \
 --output=$outputDirectory

inputFile="/home/ben/Github/HLA-allele-analysis/OutputDirectory/IN4_sequences/HLA_IN4.fasta"
outputDirectory="/home/ben/Github/find_homopolymers/output"
python find_homopolymers_main.py \
 --input=$inputFile \
 --output=$outputDirectory

inputFile="/home/ben/Github/HLA-allele-analysis/OutputDirectory/EX5_sequences/HLA_EX5.fasta"
outputDirectory="/home/ben/Github/find_homopolymers/output"
python find_homopolymers_main.py \
 --input=$inputFile \
 --output=$outputDirectory

inputFile="/home/ben/Github/HLA-allele-analysis/OutputDirectory/IN5_sequences/HLA_IN5.fasta"
outputDirectory="/home/ben/Github/find_homopolymers/output"
python find_homopolymers_main.py \
 --input=$inputFile \
 --output=$outputDirectory

inputFile="/home/ben/Github/HLA-allele-analysis/OutputDirectory/EX6_sequences/HLA_EX6.fasta"
outputDirectory="/home/ben/Github/find_homopolymers/output"
python find_homopolymers_main.py \
 --input=$inputFile \
 --output=$outputDirectory

inputFile="/home/ben/Github/HLA-allele-analysis/OutputDirectory/IN6_sequences/HLA_IN6.fasta"
outputDirectory="/home/ben/Github/find_homopolymers/output"
python find_homopolymers_main.py \
 --input=$inputFile \
 --output=$outputDirectory

inputFile="/home/ben/Github/HLA-allele-analysis/OutputDirectory/EX7_sequences/HLA_EX7.fasta"
outputDirectory="/home/ben/Github/find_homopolymers/output"
python find_homopolymers_main.py \
 --input=$inputFile \
 --output=$outputDirectory

inputFile="/home/ben/Github/HLA-allele-analysis/OutputDirectory/IN7_sequences/HLA_IN7.fasta"
outputDirectory="/home/ben/Github/find_homopolymers/output"
python find_homopolymers_main.py \
 --input=$inputFile \
 --output=$outputDirectory

inputFile="/home/ben/Github/HLA-allele-analysis/OutputDirectory/EX8_sequences/HLA_EX8.fasta"
outputDirectory="/home/ben/Github/find_homopolymers/output"
python find_homopolymers_main.py \
 --input=$inputFile \
 --output=$outputDirectory

inputFile="/home/ben/Github/HLA-allele-analysis/OutputDirectory/3UTR_sequences/HLA_3UTR.fasta"
outputDirectory="/home/ben/Github/find_homopolymers/output"
python find_homopolymers_main.py \
 --input=$inputFile \
 --output=$outputDirectory


ENDBLOCK