import pandas as pd
import os

drop_file="pedigree_data/dum.drop.txt" 
sample_id_file='pedigree_data/mouseIDs_one.txt' 
chrom_num=20 
length=0.5 
spacing= 0.001 
allele_freq=0.2

def genome_dropping(seed):
	map_df = pd.DataFrame(data={'name':[str(i) for i in range(chrom_num)], 
					 'length(cM)':[str(length)]*chrom_num, 
					 'spacing(cM)':[str(spacing)]*chrom_num, 
					 'allele_freq':[str(allele_freq)]*chrom_num})
	map_df.to_csv('pedigree_data/mouse.map.txt', index=False, header=False, sep="\t")

	os.system('g++ -o gdrop -std=c++11 -O4 main.cpp genedrop.cpp -lm')
	os.system('./gdrop -p ' + drop_file + ' -s '+str(seed))

if __name__ == '__main__':
	for i in range(11):
            genome_dropping(i)
