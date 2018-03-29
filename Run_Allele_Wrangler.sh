
ReadInputFile="/home/minion/MinIONData/2017.EfiAbstractProject/43785/Analysis/blast_results/SortedReads/HLA-C_15.fastq"
ResultsOutputDirectory="/home/minion/MinIONData/43785_C15_wrangler_test"

source activate minionvironment

cd src

python AlleleWrangler_Main.py \
 --reads=$ReadInputFile \
 --outputdir=$ResultsOutputDirectory \
 --iterations=$NumberIterations \
 --threads=$ThreadCount

source deactivate
