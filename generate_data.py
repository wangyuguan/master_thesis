import numpy as np
from scipy.linalg import block_diag, toeplitz
from statsmodels.stats.moment_helpers import cov2corr

def generate_rho_u(corr_type, D, num_spike = 10, scale=0.1, num_block=10):

    if corr_type == "identity":
        return np.eye(D)

    elif corr_type == "spiked":
        spiked = 0
        for i in range(num_spike):
            v = np.random.rand(D, 1)
            v = v/np.linalg.norm(v,2)
            spiked += (2**(-i+1))*(v @ v.T)
        return cov2corr(np.eye(D)+scale*spiked)

    elif corr_type == "geometric":
        rho_u = None
        block_size = int(D/num_block)
        rho_u_block = toeplitz([(scale**i) for i in range(block_size)])
        for i in range(num_block):
            if i == 0:
                rho_u = rho_u_block 
            else:
                rho_u = block_diag(rho_u, rho_u_block)
        return rho_u

def generate_Trait_Null(D, rho_u, rho_e, G, K0, h=None, s=None, beta=None, X=None):
    N, M = G.shape

    delta = np.zeros((M, D))
    Id = np.eye(D)

    if beta is None:
        beta=np.zeros((1, D))
    else:
        beta = beta.reshape((1,D))

    if X is None:
        X=np.ones((N, 1))

    if h is None:
        h=np.sqrt(0.5)
        xi = Id*(h**2)
    else:
        xi = np.diag(h**2)
    
    if s is None:
        s=1
        eta = Id*(s**2)
    else:
        eta = np.diag(s**2)

    #print(eta)
    
    Vu = np.sqrt(eta) @ np.sqrt(xi) @ rho_u @ np.sqrt(xi) @ np.sqrt(eta)
    Ve = np.sqrt(eta) @ np.sqrt(Id-xi) @ rho_e @ np.sqrt(Id-xi) @ np.sqrt(eta)

    mu = np.zeros(D)
    U, Lambda, _ = np.linalg.svd(K0)
    Lambda[Lambda<1e-13] = 1e-13
    _vec_e = np.zeros([D,N])
    for n in range(N):
        cov = Lambda[n]*Vu + Ve 
        _vec_e[:,n] = np.random.multivariate_normal(mu, cov, size=1)
    vec_e = U @ _vec_e.T
    return (X @ beta + G @ delta + vec_e)
    