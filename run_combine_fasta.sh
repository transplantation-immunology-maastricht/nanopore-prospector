source /home/ben/minionvenv/bin/activate

InputDirectory="/MinIONData/DP/GenerateDPConsensus/OutputFiles/"
OutputDirectory="/MinIONData/DP/GenerateDPConsensus/CollectedConsensusSequences/"
Pattern="_ConsensusWithAlleles.fasta"

cd /home/ben/Github/nanopore_prospector

python Nanopore_Prospector_Main.py \
 --inputdirectory=$InputDirectory \
 --outputdirectory=$OutputDirectory \
 --pattern=$Pattern \
 --action="combinefastafiles"
