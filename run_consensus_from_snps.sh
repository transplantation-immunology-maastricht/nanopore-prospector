source /home/ben/minionvenv/bin/activate

cd /home/ben/Github/nanopore_prospector   

ReferenceSequenceFile="/MinIONData/2019.DRA.2ndRoundBasecalling/GenerateNovelSequences/DRA_01010101.fasta"
SnpFileLocation="/MinIONData/2019.DRA.2ndRoundBasecalling/GenerateNovelSequences/novelSNPTable.csv"
OutputDirectory="/MinIONData/2019.DRA.2ndRoundBasecalling/GenerateNovelSequences/NovelAlleles"

python Nanopore_Prospector_Main.py \
 --inputfile=$SnpFileLocation \
 --outputdirectory=$OutputDirectory \
 --reference=$ReferenceSequenceFile \
 --action="consensusfromsnps"
