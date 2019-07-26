source /home/ben/minionvenv/bin/activate

#ReadDirectory="/home/ben/Github/nanopore_prospector/33332_analyze_read_quality/A24_Reads"
ReadFile="/MinIONData/DP/GenerateDPConsensus/CollectedConsensusSequences/CombinedSequences.fasta"
ResultsOutputDirectory="/MinIONData/DP/GenerateDPConsensus/ConsensusSNPCalculations"
ReferenceFile="/MinIONData/DP/GenerateDPConsensus/REFERENCE.fasta"

cd /home/ben/Github/nanopore_prospector

python Nanopore_Prospector_Main.py \
 --reads=$ReadFile \
 --outputdirectory=$ResultsOutputDirectory \
 --reference=$ReferenceFile \
  --action="snpanalysis"
