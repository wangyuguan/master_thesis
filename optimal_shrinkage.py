import numpy as np
from statsmodels.stats.moment_helpers import cov2corr
from scipy.linalg import sqrtm

def cosine_function(l, gamma):

    if l<=(1+np.sqrt(gamma)):
        return 0
    else:
        return np.sqrt(((1-gamma/(l-1)**2)/(1+gamma/(l-1)))) 

def inv_biasing_function(ev, gamma):

    return 0.5*(ev+1-gamma+np.sqrt((ev+1-gamma)**2-4*ev))

def ev_shrink(ev, gamma, loss = "frobenius"):

    if ev>(1+np.sqrt(gamma))**2:
        l = inv_biasing_function(ev, gamma)
        c_square = cosine_function(l, gamma)**2
        s_square = 1 - c_square
        if loss=="operator":
            return l
        elif loss=="frobenius" or loss=="entropy":
            return l*c_square+s_square
        elif loss=="frobenius_precision" or loss=="stein":
            return l/(c_square+l*s_square)
        elif loss == "frechet":
            return (s_square+np.sqrt(l)*c_square)**2
    else:
        return 1

def opt_shrink(C, gamma, loss = "frobenius"):

    U, D, VT = np.linalg.svd(C)

    return cov2corr(U@np.diag(np.array([ev_shrink(ev, gamma, loss) for ev in D]))@VT)

def matrix_loss(A, B, loss="operator", A_inv=None, B_inv=None):

    I = np.eye(A.shape[0])
    if loss == "operator":
        return np.linalg.norm(A-B, ord=2)
    elif loss == "frobenius":
        return np.linalg.norm(A-B, ord='fro')**2
    elif loss == "entropy":
        if B_inv == None:
            B_inv = np.linalg.inv(B)
        return 0.5*(np.trace(B_inv@A - I) - np.log(np.linalg.det(A)/np.linalg.det(B)))
    elif loss == "frobenius_precision":
        if A_inv == None:
            A_inv = np.linalg.inv(A)
        if B_inv == None:
            B_inv = np.linalg.inv(B)
        return np.linalg.norm(A_inv-B_inv, ord='fro')**2
    elif loss == "stein":
        if A_inv == None:
            A_inv = np.linalg.inv(A)
        return 0.5*(np.trace(A_inv@B - I) - np.log(np.linalg.det(B)/np.linalg.det(A)))
    elif loss == "frechet":
        A_sqrt = sqrtm(A)
        B_sqrt = sqrtm(B)
        return np.trace(A+B-2*A_sqrt@B_sqrt)
        