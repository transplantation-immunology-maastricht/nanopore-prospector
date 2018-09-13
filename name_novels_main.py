from sys import argv, exc_info
from getopt import getopt, GetoptError

from os.path import isdir, isfile

from combine_snps.name_novels import nameNovels


if __name__=='__main__':

    try:    
        if(True):
            print('Commandline arguments look fine.\nThe hour is at hand.')
            
            IMGTFIleLocation = "/home/minion/MinIONData/2018.DRA.Basecalling/NovelAlleleNaming/IMGT_Alleles.fasta"
            DRAFileLocation = "/home/minion/MinIONData/2018.DRA.Basecalling/NovelAlleleNaming/Ben_Alleles.fasta"
            output = "/home/minion/MinIONData/2018.DRA.Basecalling/NovelAlleleNaming/NovelAlleles"
            
            nameNovels(IMGTFIleLocation, DRAFileLocation, output)
            
            print ('I am done for now, have a nice day.')    
        else:
            print('\nI\'m giving up because I was not satisfied with your commandline arguments.')  
            
    except Exception:
        # Top Level exception handling like a pro.
        # This is not really doing anything.
        print ('Fatal problem during sequence extraction:')
        print (exc_info())
        raise


