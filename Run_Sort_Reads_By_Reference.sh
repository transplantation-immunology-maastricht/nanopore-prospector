# TODO: GPL License
# TODO: Use the environment, instead of conda.

source /home/ben/anaconda2/bin/activate minionvironment



#python Blast_Minion_Reads.py \
#    "/home/ben/MUMCScripts/BlastMinIONReads/inputData/HLA_ClassI_APD.fasta"\
#    "/home/ben/MUMCScripts/BlastMinIONReads/inputData/33332R9_short.fasta"\
#    "/home/ben/MUMCScripts/BlastMinIONReads/output33332_short" 


python Blast_Minion_Reads.py \
    "/home/minion/MinIONData/DianeClassII/HLA_Alleles_APD_ClassII.fasta" \
    "/home/minion/MinIONData/DianeClassII/extracts_min_2000/BC02_TwoDirReads.fasta" \
    "/home/minion/MinIONData/DianeClassII/BC02_Sorted"
 



source /home/ben/anaconda2/bin/deactivate
