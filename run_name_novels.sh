source /home/ben/minionvenv/bin/activate

cd /home/ben/Github/nanopore_prospector   

ReferenceSequenceFileLocation="/MinIONData/2019.DRA.2ndRoundBasecalling/NovelAlleleNaming/DRA_Alleles.fasta"
NovelAlleleFileLocation="/MinIONData/2019.DRA.2ndRoundBasecalling/GenerateNovelSequences/NovelAlleles/GeneratedSequences.fasta"
OutputDirectory="/MinIONData/2019.DRA.2ndRoundBasecalling/NovelAlleleNaming/NovelAlleles"

# Pass the IMGT/HLA Alleles file as the known reference sequence. It contains known alleles. It's the --reference parameter.
# Pass the Novel Sequences as the --inputfile

python Nanopore_Prospector_Main.py \
 --inputfile=$NovelAlleleFileLocation \
 --outputdirectory=$OutputDirectory \
 --reference=$ReferenceSequenceFileLocation \
 --action="namenovels"

