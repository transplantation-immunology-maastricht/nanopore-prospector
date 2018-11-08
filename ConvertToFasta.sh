source /home/ben/minionvenv/bin/activate

inputFile="/home/ben/ben_share/ConvertDPSequences/47470M_DP.fastq"
outputFile="/home/ben/ben_share/ConvertDPSequences/47470M_DP.fasta"
python /home/ben/Github/nanopore_prospector/Nanopore_Prospector_Main.py \
 --inputfile=$inputFile \
 --outputfile=$outputFile \
 --action="fastqtofasta"

