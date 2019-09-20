source /home/ben/minionvenv/bin/activate

inputFile="/home/ben/ben_share/15.10.2018.HaploDP/M180828TJ/8609M_M180828A_BC05.fastq"
outputFile="/home/ben/ben_share/15.10.2018.HaploDP/M180828TJ/8609M_M180828A_BC05.fasta"
python Nanopore_Prospector_Main.py \
 --inputfile=$inputFile \
 --outputfile=$outputFile \
 --action="fastqtofasta"

