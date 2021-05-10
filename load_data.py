import pandas as pd
import numpy as np

import limix 
from limix.stats import linear_kinship

def load_from_true_geno(seed):
    sample_id_df = pd.read_csv('genome_drop/pedigree_data/mouseIDs_one.txt',
                            index_col=None,
                            sep='\t',
                            dtype="string",
                            header=None)[0].tolist()
    raw_data_file = "genome_drop/pedigree_data/dum.drop_"+str(seed)+".geno_true"
    raw_data = pd.read_csv(raw_data_file,index_col=None,
                                         sep='\t',
                                         dtype="string")
    chromsome_list = raw_data['#chrom']

    geno2score =  {'1': 1, '2': 1, '3': 0, '4': 0}
    M = raw_data.shape[0]
    N = len(sample_id_df)
    G = np.zeros((N,M))
    for i in range(N):
        for j in range(M):
            G[i,j] = int(geno2score[raw_data[sample_id_df[i]][j][0]]+geno2score[raw_data[sample_id_df[i]][j][-1]])

    chromsome_list = np.array([int(chrom) for chrom in chromsome_list])

    return G, chromsome_list

def find_low_maf_allele(G):
    f = np.mean(G, axis=0)/2
    col_index = [(0.1<ff<0.9) for ff in f]
    return col_index
