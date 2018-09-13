source /home/ben/minionvenv/bin/activate

InputDirectory="/home/minion/MinIONData/2018.DRA.Basecalling/Collected_Snp_Arrays"
OutputDirectory="/home/minion/MinIONData/2018.DRA.Basecalling/Combined_SNP_Arrays"

python combine_snps_main.py \
 --input=$InputDirectory \
 --output=$OutputDirectory 

#InputDirectory="/home/ben/ben_share/DRA_Analysis/Test_Deletion_Length/Snp_Reports_for_combining"
#OutputDirectory="/home/ben/ben_share/DRA_Analysis/Test_Deletion_Length/combined_snp_reports"

#python combine_snps_main.py \
# --input=$InputDirectory \
# --output=$OutputDirectory 

