import numpy as np

from scipy.linalg import block_diag, sqrtm
from numpy_sugar.linalg import economic_qs

import limix 
from limix.stats import linear_kinship
from limix.model.lmm import LMM

from statsmodels.stats.moment_helpers import cov2corr

from r_helpers import find_inv_via_R

def run_limix(Y, X, G, K, etaMax=0.99):
    N, M = G.shape
    _, D = Y.shape
    QS = economic_qs(K)

    Z = np.zeros((M, D))

    for d in range(D):
        lmm = LMM(Y[:,d], X, QS, restricted=False)
        lmm.fit(verbose=False)
        delta=lmm.v0/(lmm.v0+lmm.v1)
        if delta>etaMax:
            lmm = LMM(Y[:,d], X, QS, restricted=True)
            lmm.fit(verbose=False)

        ret=lmm.get_fast_scanner().fast_scan(G,False) 
        wald = ret['effsizes1']/ret['effsizes1_se']  
        Z[:,d]= wald

    return Z.T

# def estimate_C1_via_R(G, chromsome_list, geno_mat = None, Lambda=0.005):
#     N, M = G.shape
#     C1, L = None, None

#     if geno_mat is not None:
#         print(geno_mat.shape)
#         decorr_mat = find_inv_via_R(sqrtm(geno_mat))

#     chromsome_unique, _ = np.unique(chromsome_list, return_inverse=True)
#     centering = np.eye(N) - (1/N)*np.ones((N,N))
#     counter = 1

#     for chrom in chromsome_unique:
#         subset = np.where(chromsome_list==chrom)[0]
#         G_m = G[:,subset]
#         if geno_mat is None:
#             V1_block = G_m.T@centering@G_m/(N-1)
#         else:
#             G_m_decorr = decorr_mat @ G_m
#             V1_block = G_m_decorr.T@centering@G_m_decorr/(N-1)

#         C1_block = cov2corr(V1_block)
#         #dC1_block = np.diag(np.diag(C1_block))
#         #C1_block_shrink = (1-Lambda)*C1_block + Lambda*np.eye(C1_block.shape[0])
#         #L_inv_block = find_inv_via_R(find_chol_via_R(C1_block))
#         if counter == 1:
#             C1 = C1_block
#             #C1_shrink = C1_block_shrink
#             #L_inv = L_inv_block
#             counter += 1
#         else:
#             C1 = block_diag(C1, C1_block)
#             #C1_shrink = block_diag(C1_shrink, C1_block_shrink)
#             #L_inv = block_diag(L_inv, L_inv_block)

#     C1_shrink = (1-Lambda)*C1 + Lambda*np.eye(M)
#     np.savetxt("res_RvsPython_C1/C1_mat.txt",C1)
#     np.savetxt("res_RvsPython_C1/C1_mat_shrink.txt",C1_shrink)

#     return C1,C1_shrink

    #return C1, L_inv

# def estimate_C1_via_py(G, chromsome_list, geno_mat = None, Lambda=0.005):
#     N, M = G.shape
#     C1, L = None, None

#     if geno_mat is not None:
#         print(geno_mat.shape)
#         decorr_mat = np.linalg.inv(sqrtm(geno_mat))

#     chromsome_unique, _ = np.unique(chromsome_list, return_inverse=True)
#     centering = np.eye(N) - (1/N)*np.ones((N,N))
#     counter = 1

#     for chrom in chromsome_unique:
#         subset = np.where(chromsome_list==chrom)[0]
#         G_m = G[:,subset]
#         if geno_mat is None:
#             V1_block = G_m.T@centering@G_m/(N-1)
#         else:
#             G_m_decorr = decorr_mat @ G_m
#             V1_block = G_m_decorr.T@centering@G_m_decorr/(N-1)

#         C1_block = cov2corr(V1_block)
#         #dC1_block = np.diag(np.diag(C1_block))
#         C1_block = (1-Lambda)*C1_block + Lambda*np.eye(C1_block.shape[0])
#         L_inv_block = np.linalg.inv(np.linalg.cholesky(C1_block))
#         if counter == 1:
#             C1 = C1_block
#             L_inv = L_inv_block
#             counter += 1
#         else:
#             C1 = block_diag(C1, C1_block)
#             L_inv = block_diag(L_inv, L_inv_block)

#     return C1, L_inv

# def estimate_C1(G, chromsome_list, Lambda=0.005):
#     N, M = G.shape
#     centering = np.eye(N) - (1/N)*np.ones((N,N))
#     C1, L = None, None

#     chromsome_unique, _ = np.unique(chromsome_list, return_inverse=True)
#     counter = 1

#     for chrom in chromsome_unique:
#         subset = np.where(chromsome_list==chrom)[0]
#         G_m = G[:,subset]
#         V1_block = G_m.T@centering@G_m/(N-1)
#         C1_block = cov2corr(V1_block)

#         C1_block = (1-Lambda)*C1_block + Lambda*np.eye(C1_block.shape[0])
#         L_inv_block = np.linalg.inv(np.linalg.cholesky(C1_block))
#         if counter == 1:
#             C1 = C1_block
#             L_inv = L_inv_block
#             counter += 1
#         else:
#             C1 = block_diag(C1, C1_block)
#             L_inv = block_diag(L_inv, L_inv_block)

#     return C1, L_inv

def estimate_C1(G, chromsome_list, Vmouse=None, Lambda=0.005, UseR=False):
    N, M = G.shape
    C1, C1_inv = None, None

    chromsome_unique, _ = np.unique(chromsome_list, return_inverse=True)
    counter = 1

    if Vmouse is None:
    	Weight_mat = np.eye(N) - (1/N)*np.ones((N,N))
    else:
    	Weight_mat = find_inv_via_R(Vmouse) if UseR else np.linalg.inv(Vmouse)

    for chrom in chromsome_unique:
        subset = np.where(chromsome_list==chrom)[0]
        G_m = G[:,subset]
        V1_block = (G_m-1).T @ Weight_mat @ (G_m-1)
        C1_block = cov2corr(V1_block)

        C1_block = (1-Lambda)*C1_block + Lambda*np.eye(C1_block.shape[0])
        C1_inv_block = find_inv_via_R(C1_block) if UseR else np.linalg.inv(C1_block)
        
        if counter == 1:
            C1 = C1_block
            C1_inv = C1_inv_block
            counter += 1
        else:
            C1 = block_diag(C1, C1_block)
            C1_inv = block_diag(C1_inv, C1_inv_block)
    return C1, C1_inv


def estimate_C2(Z, C1_inv=None):
    D, M = Z.shape
    if C1_inv is None:
        Weight_mat = np.eye(M) - (1/M)*np.ones((M,M))
        return cov2corr(Z@Weight_mat@Z.T)
    return cov2corr(Z@(C1_inv-(C1_inv@np.ones((M,M))@C1_inv)/np.matrix(C1_inv).sum())@Z.T)
