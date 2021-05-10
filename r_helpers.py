import pandas as pd
import numpy as np
import subprocess
import os

def find_eigen_via_R(C):

    np.savetxt("corr_matrix.txt",C)
    command = 'Rscript'
    path2script = 'calculate_ev.R'
    cmd = [command, path2script]
    subprocess.check_output(cmd, universal_newlines=True)

    ret = pd.read_csv("ev_res.txt",delimiter=" ").to_numpy()
    os.system("rm corr_matrix.txt ev_res.txt")

    return ret

def find_chol_via_R(C):

    np.savetxt("mat4chol.txt",C)
    command = 'Rscript'
    path2script = 'calculate_chol.R'
    cmd = [command, path2script]
    subprocess.check_output(cmd, universal_newlines=True)

    ret = np.loadtxt("chol_res.txt")
    os.system("rm mat4chol.txt chol_res.txt")

    return ret

def find_inv_via_R(C):

    np.savetxt("mat4inv.txt",C)
    command = 'Rscript'
    path2script = 'calculate_inv.R'
    cmd = [command, path2script]
    subprocess.check_output(cmd, universal_newlines=True)

    ret = np.loadtxt("inv_res.txt")
    os.system("rm mat4inv.txt inv_res.txt")

    return ret
