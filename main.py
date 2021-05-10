import random
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from scipy import stats

import limix 
from limix.stats import linear_kinship

from corr_estimation import run_limix, estimate_C1, estimate_C2
from load_data import load_from_true_geno, find_low_maf_allele
from generate_data import generate_rho_u, generate_Trait_Null
from optimal_shrinkage import opt_shrink, matrix_loss

from r_helpers import find_eigen_via_R
from statsmodels.stats.moment_helpers import cov2corr


def compute_performance(D=5000, rho_u_type="spiked", mat_loss="operator", useSample = 5, num_points = 10, h = None, s=None, beta=None):

    G0, _ = load_from_true_geno(seed=0)
    G0 = G0[:,find_low_maf_allele(G0)]
    K0 = linear_kinship(G0, verbose=True)
    N, _ = G0.shape
    X = np.ones((N, 1))

    Id = np.eye(D)

    if rho_u_type=="spiked":
        scale_list = np.linspace(1, 5, num=num_points, endpoint=True)
    elif rho_u_type=="geometric":
        scale_list = np.linspace(0.1, 0.5, num=num_points, endpoint=True)

    ave_loss_S2 = []
    ave_loss_S2_old = []
    ave_loss_C2_naive = []

    for scale in scale_list:
        rho_u = generate_rho_u(rho_u_type, D, scale=scale)
        loss_S2 = []
        loss_S2_old = []
        loss_C2_naive = []

        for seed in range(1,useSample+1):
            G, chromsome_list = load_from_true_geno(seed=seed)
            col_index = find_low_maf_allele(G)
            G = G[:,find_low_maf_allele(G)]
            chromsome_list = chromsome_list[col_index]

            Y = generate_Trait_Null(D, rho_u, rho_u, G, K0, h =h, s = s, beta=beta)

            Z = run_limix(Y, X, G, K0)
            D, M = Z.shape
            Weight_mat = np.eye(M) - (1/M)*np.ones((M,M))
            C2_naive = cov2corr(Z@Weight_mat@Z.T)
            loss_C2_naive.append(matrix_loss(rho_u, C2_naive, mat_loss))

            C1, C1_inv = estimate_C1(G, chromsome_list, UseR=True)
            C2 = estimate_C2(Z, C1_inv)
            S2 = opt_shrink(C2, D/207)
            loss_S2.append(matrix_loss(rho_u, S2, mat_loss))

            C2_old = estimate_C2(Z)
            S2_old = opt_shrink(C2_old, D/207)
            loss_S2_old.append(matrix_loss(rho_u, S2_old, mat_loss))


        ave_loss_S2.append(np.mean(loss_S2))
        ave_loss_S2_old.append(np.mean(loss_S2_old))
        ave_loss_C2_naive.append(np.mean(loss_C2_naive))

        # fname_S2 = "results/"+rho_u_type+"_scale_"+str(scale)+"_"+mat_loss+"loss_S2.txt"
        # fname_S2_old = "results/"+rho_u_type+"_scale_"+str(scale)+"_"+mat_loss+"_loss_S2_old.txt"
        # fname_C2_naive = "results/"+rho_u_type+"_scale_"+str(scale)+"_"+mat_loss+"_loss_C2_naive.txt"

        # np.savetxt(fname_S2, loss_S2)
        # np.savetxt(fname_S2_old, loss_S2_old)
        # np.savetxt(fname_C2_naive, loss_C2_naive)

    fig = plt.figure()
    fig.suptitle(mat_loss)

    plt.plot(scale_list, np.log10(ave_loss_S2), label="optShrink+decorr",marker='o')
    plt.plot(scale_list, np.log10(ave_loss_S2_old), label="optShrink_only",marker='o')
    plt.plot(scale_list, np.log10(ave_loss_C2_naive), label="naive",marker='o')

    plt.xlabel(r"$\rho$")
    plt.ylabel(r'$\log_{10}(Average\,\,Matrix\,\,Loss)$')

    plt.legend(loc='lower right')
    fig_name = "results/"+rho_u_type+"_"+mat_loss+".jpeg"
    plt.savefig(fig_name, dpi=300)

if __name__ == '__main__':
    D = 5000
    h = np.sqrt(np.random.uniform(0,1,D))
    s = np.sqrt(np.random.gamma(1,1,D))
    beta = np.random.normal(0,100,D)
    rho_u_type_list = ["spiked", "geometric"]
    mat_loss_list = ["operator", "frobenius"]
    for rho_u_type in rho_u_type_list:
        for mat_loss in mat_loss_list:
            compute_performance(rho_u_type=rho_u_type, mat_loss=mat_loss)


