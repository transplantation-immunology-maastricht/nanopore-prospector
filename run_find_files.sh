source /home/ben/minionvenv/bin/activate

InputDirectory="/home/minion/MinIONData/2018.DRA.Basecalling/SNP_Analysis"
OutputDirectory="/home/minion/MinIONData/2018.DRA.Basecalling/Collected_Snp_Arrays"

python find_files_main.py \
 --input=$InputDirectory \
 --output=$OutputDirectory 

