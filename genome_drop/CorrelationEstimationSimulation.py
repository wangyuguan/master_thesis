import pandas as pd
import numpy as np
import os
from scipy.stats import norm
from scipy.stats import matrix_normal
import limix 
from limix.stats import linear_kinship
from limix.model.lmm import LMM
from numpy_sugar.linalg import economic_qs
from statsmodels.stats.moment_helpers import cov2corr
from scipy.linalg import block_diag
from datetime import datetime
from random import randrange
import random
from matplotlib import pyplot as plt
from scipy import stats
import seaborn as sns
import statsmodels.api as sm

import rpy2.robjects as robj
import rpy2.robjects.numpy2ri as rpyn

random.seed(20201016)

def getEigenvalues_byR(Z_score):
	nr, nc = Z_score.shape
	Z_score_vec = robj.FloatVector(Z_score.transpose().reshape((Z_score.size)))
	Z_score_r = robj.r.matrix(Z_score_vec, nrow=nr, ncol=nc)
	centering = np.eye(nc) - (1/nc)*np.ones((nc,nc))
	centering_vec = robj.FloatVector(centering.transpose().reshape((centering.size)))
	centering_r = robj.r.matrix(centering_vec, nrow=nc, ncol=nc)
	#eig_r = robj.r['eigen']
	#eigen_val = np.array(eig_r(matr, only.values = True))
	robj.r('''
		f <- function(Z_score_r,centering_r){
			print(dim(Z_score_r))
			print(dim(centering_r))
			eigen(Z_score_r%*%centering_r%*%t(Z_score_r), symmetric = False, only.values = TRUE)
		}
		''')
	r_f = robj.globalenv['f']
	eigen_val = np.asarray(r_f(Z_score_r,centering_r))
	return np.array(list(eigen_val[0]))


#def getEigenvalues_byR(mat):
#	nr, nc = mat.shape
#	matvec = robj.FloatVector(mat.transpose().reshape((mat.size)))
#	matr = robj.r.matrix(matvec, nrow=nr, ncol=nc)
#	#eig_r = robj.r['eigen']
#	#eigen_val = np.array(eig_r(matr, only.values = True))
#	robj.r('''
#		f <- function(matr){
#			eigen(matr, only.values = TRUE)
#		}
#		''')
#	r_f = robj.globalenv['f']
#	eigen_val = np.asarray(r_f(matr))
#	return np.array(list(eigen_val[0]))


class Mouse_eQTLdata_Simulator(object):
	def __init__(self, drop_file, sample_id_file, 
		chrom_num, length, spacing, allele_freq,
		raw_data_file,
		num_trait,  
		H_0, 
		c = 1, 
		U_num_block=None, U_off_diagonal=None, U_indep=True,
		e_num_block=None, e_off_diagonal=None, e_indep=True,
		xi = 0.5, eta=1,
		zero_non_genetic_effect = True):
		self.drop_file = drop_file 
		self.sample_id_file = sample_id_file
		self.raw_data_file = raw_data_file
		self.chrom_num = chrom_num 
		self.length = length
		self.spacing = spacing
		self.allele_freq = allele_freq
		self.num_trait = num_trait
		self.raw_data_file = raw_data_file
		self.H_0 = H_0
		self.c = c
		self.U_num_block = U_num_block
		self.U_off_diagonal = U_off_diagonal
		self.U_indep = U_indep
		self.e_num_block = e_num_block
		self.e_off_diagonal = e_off_diagonal
		self.e_indep = e_indep
		self.xi = xi
		self.eta = eta
		self.zero_non_genetic_effect = zero_non_genetic_effect
		self.compiled = False

	def genome_dropping(self, seed=1):
		map_df = pd.DataFrame(data={'name':[str(i) for i in range(self.chrom_num)], 
					 'length(cM)':[str(self.length)]*self.chrom_num, 
					 'spacing(cM)':[str(self.spacing)]*self.chrom_num, 
					 'allele_freq':[str(self.allele_freq)]*self.chrom_num})
		map_df.to_csv('mouse.map.txt', index=False, header=False, sep="\t")

		if not self.compiled:
			os.system('g++ -o gdrop -std=c++11 -O4 main.cpp genedrop.cpp -lm')
			self.compiled = True

		os.system('./gdrop -p ' + self.drop_file + ' -s '+str(seed))

	def raw2genotype(self):
		true_geno_raw_data = pd.read_csv(self.raw_data_file,
										 index_col=None,
										 sep='\t',
										 dtype="string")
		self.chromsome_list = true_geno_raw_data['#chrom']

		self.sample_id_df = pd.read_csv(self.sample_id_file,
							index_col=None,
							sep='\t',
							dtype="string",
							header=None)[0].tolist()
		
		geno_score_map =  {'1': 1, '2': 1, '3': 0, '4': 0}
		M = true_geno_raw_data.shape[0]
		N = len(self.sample_id_df)
		self.G = np.zeros((N,M))
		for i in range(N):
			for j in range(M):
				self.G[i,j] = geno_score_map[true_geno_raw_data[self.sample_id_df[i]][j][0]]+geno_score_map[true_geno_raw_data[self.sample_id_df[i]][j][-1]] 
		
		self.chromsome_list = np.array([int(chrom) for chrom in self.chromsome_list])
		f = np.mean(self.G, axis=0)/2
		index = [(0.1<ff<0.9) for ff in f]
		self.f = f[index]
		self.G = np.matrix(self.G[:,index])
		self.num_sample, self.num_snps = self.G.shape
		self.chromsome_list = self.chromsome_list[index]

	def check_LD(self):
		geno_corr = np.corrcoef(self.G, rowvar=False)
		one_way_corr_same = []
		one_way_corr_diff = []
		for i in range(self.G.shape[1]-1):
			if self.chromsome_list[i] == self.chromsome_list[i+1]:
				one_way_corr_same.append(geno_corr[i,i+1])
			else:
				one_way_corr_diff.append(geno_corr[i,i+1])

		ave_k_way_corr = []
		for k in range(1,51):
			k_way_corr = []
			for i in range(self.G.shape[1]):
				if (i+k)<self.G.shape[1] and self.chromsome_list[i] == self.chromsome_list[i+k]:
					k_way_corr.append(abs(geno_corr[i,i+k]))
			ave_k_way_corr.append(np.mean(k_way_corr))

		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax.hist(one_way_corr_same, bins=50)
		plt.savefig("check_LD/one_way_corr_same.jpg")
		plt.close('all')

		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax.hist(one_way_corr_diff, bins=10)
		plt.savefig("check_LD/one_way_corr_diff.jpg")
		plt.close('all')

		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax.plot(np.array(range(1,51)),ave_k_way_corr,'.-')
		plt.savefig("ave_abs_k_way_corr.jpg")
		plt.close('all')

	def compute_GRM(self):
		
		self.K = linear_kinship(self.G, verbose=True)
		self.K_corr = cov2corr(self.K)
		Zm = (self.G-2*self.f)/np.sqrt(2*self.f*(1-self.f))
		self.K_man = (1/self.num_snps) * Zm @ Zm.T
		self.K_man_corr = cov2corr(self.K_man)

	def getG0_and_K0(self):
		self.K0 = self.K
		self.G0 = self.G
		self.QS = economic_qs(self.K0)
		
	def generate_rho_u(self):
		if self.U_indep:
			self.rho_u  = np.eye(self.num_trait)

		else:
			block_n = int(self.num_trait/self.U_num_block)
			self.rho_u = np.zeros((self.num_trait, self.num_trait))
			for i in range(self.U_num_block):
				x = np.ones((block_n, block_n))*self.U_off_diagonal
				self.rho_u[i*block_n : (i+1)*block_n, i*block_n : (i+1)*block_n] = x
			for i in range(self.num_trait):
				self.rho_u[i,i] = 1

		return self.rho_u 

	def generate_rho_e(self):
		if self.e_indep:
			self.rho_e  = np.eye(self.num_trait)

		else:
			block_n = int(self.num_trait/self.e_num_block)
			self.rho_u = np.zeros((self.num_trait, self.num_trait))
			for i in range(self.e_num_block):
				x = np.ones((block_n, block_n))*self.e_off_diagonal
				self.rho_e[i*block_n : (i+1)*block_n, i*block_n : (i+1)*block_n] = x
			for i in range(self.num_trait):
				self.rho_e[i,i] = 1

		return self.rho_e

	def generate_heritability(self):
		return np.eye(self.num_trait)*self.xi

	def generate_std_matrix(self):
		return np.eye(self.num_trait)*self.eta

	def generate_covariate(self):
		return np.ones((self.num_sample,1))
		#return np.random.normal(size=self.num_sample*self.c).reshape((self.num_sample,self.c))

	def generate_nongenetic_coefficient(self):
		if self.zero_non_genetic_effect:
			return np.zeros((self.c,self.num_trait))
		else:
			return np.random.normal(self.c*self.num_trait).reshape((self.c,self.num_trait))

	def generate_genetic_coefficient(self):
		if self.H_0:
			return np.zeros((self.num_snps,self.num_trait))
		else:
			return np.random.normal(size=self.num_snps*self.num_trait).reshape((self.num_snps,self.num_trait))

	def generate_parameters(self):
		self.beta = self.generate_nongenetic_coefficient()
		self.delta = self.generate_genetic_coefficient()
		self.xi = self.generate_heritability()
		self.eta = self.generate_std_matrix()
		self.X = self.generate_covariate()
		self.rho_e = self.generate_rho_e()
		self.rho_u = self.generate_rho_u()

	def generate_expressionTrait(self, MatrixVariate=False):
		self.generate_parameters()
		term1 = self.X @ self.beta
		term2 = self.G @ self.delta

		if MatrixVariate:
			try:
				U = matrix_normal(rowcov=self.K, colcov=self.rho_u).rvs(1)
			except:
				K = self.K+1e-8*np.eye(self.num_sample)
				U = matrix_normal(rowcov=K, colcov=self.rho_u).rvs(1)
			epsilon = matrix_normal(rowcov=np.diag(np.ones(self.num_sample)), colcov=self.rho_e).rvs(1)
			term3 = U @ np.sqrt(self.xi)
			term4 = epsilon @ np.sqrt(np.diag(np.ones(self.num_trait))-self.xi)
			term5 = (term3 + term4) @ np.sqrt(self.eta)
			self.Y = term1 + term2 + term5

		else:
			Vu = np.sqrt(self.eta) @ np.sqrt(self.xi) @ self.rho_u @ np.sqrt(self.xi) @ np.sqrt(self.eta)
			Ve = np.sqrt(self.eta) @ np.sqrt(np.eye(self.num_trait)-self.xi) @ self.rho_u @ np.sqrt(np.eye(self.num_trait)-self.xi)@np.sqrt(self.eta)
			mean = np.zeros(self.num_trait)
			U, Lambda, _ = np.linalg.svd(self.K)
			Lambda[Lambda<1e-13] = 1e-13
			_vec_e = np.zeros([self.num_trait,self.num_sample])
			now = datetime.now().time()
			print('******************************* Generating Y:', " now the time is ", now)
			print(self.num_sample)
			for n in range(self.num_sample):
				if (n+1)%20==0:
					print('sampling : {} of {}'.format(n+1,self.num_sample),flush=True)
				cov = Lambda[n]*Vu + Ve 
				_vec_e[:,n] = np.random.multivariate_normal(mean=mean, cov=cov, size=1)
			vec_e = U @ _vec_e.T
			self.Y = term1 + term2 + vec_e

	def check_colCov(self):
		index = np.arange(0,self.num_trait,25)
		self.colCov = np.cov(self.Y,rowvar=False)
		#print(self.rowCov.shape, index.shape)
		self.colCov_subset = self.colCov[index,:][:,index]
		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax = sns.heatmap(self.colCov_subset, center=0)
		plt.savefig("colCov_heatmap.jpg")

		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax = sns.heatmap(np.eye(len(index)), center=0)
		plt.savefig("I_D.jpg")

		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax = sns.heatmap((self.colCov_subset-np.eye(len(index))) , center=0)
		plt.savefig("colCov_Difference_heatmap.jpg")

	def check_colCorr(self):
		index = np.arange(0,self.num_trait,25)
		self.colCorr = cov2corr(self.colCov)
		self.colCorr_subset = self.colCorr[index,:][:,index]
		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax = sns.heatmap(self.colCorr_subset, center=0)
		plt.savefig("colCorr_heatmap.pdf")

	def check_rowCov(self):
		self.rowCov = np.cov(self.Y,rowvar=True)
		print(self.Y.shape, self.rowCov.shape)
		self.true_rowCov = (0.5*self.K0+0.5*np.eye(self.num_sample))
		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax = sns.heatmap(self.rowCov , center=0)
		plt.savefig("rowCov_heatmap.jpg")

		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax = sns.heatmap(self.true_rowCov , center=0)
		plt.savefig("true_rowCov_heatmap.jpg")

		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax = sns.heatmap((self.rowCov-self.true_rowCov) , center=0)
		plt.savefig("rowCov_Difference_heatmap.jpg")
			
	def generateAll(self):
		self.genome_dropping()
		self.raw2genotype()
		self.compute_GRM()
		self.getG0_and_K0()
		self.generate_parameters()
		self.generate_expressionTrait()

class Correlation_Obtainer(object):
	def __init__(self, mouse_eQTLdata, perturbation=1e-10):
		self.perturbation = perturbation
		self.Y = mouse_eQTLdata.Y
		self.X = mouse_eQTLdata.X
		self.G = mouse_eQTLdata.G
		self.K = mouse_eQTLdata.K0
		self.QS = mouse_eQTLdata.QS
		self.num_sample = mouse_eQTLdata.num_sample
		self.num_trait = mouse_eQTLdata.num_trait
		self.num_snps = mouse_eQTLdata.num_snps
		self.chromsome_list = mouse_eQTLdata.chromsome_list
		self.spacing = mouse_eQTLdata.spacing 
		now = datetime.now().time()
		print('sizes of the data: Nsnps=%d, Ntrait=%d, Nsample=%d.'%(self.num_snps,self.num_trait,self.num_sample))

	def run_limix(self, etaMax=0.99, etamin=0.01):

		self.Z = np.zeros((self.num_snps, self.num_trait))
		self.eta = np.ones(self.num_trait)
		self.delta_list = []
		self.var_list = []

		for d in range(self.num_trait):
			if (d+1)%500==0:
				print('runLimix : {} of {}'.format(d+1,self.num_trait),flush=True)
			lmm = LMM(self.Y[:,d], self.X, self.QS, restricted=False)
			lmm.fit(verbose=False)
			delta=lmm.v0/(lmm.v0+lmm.v1)
			if delta>etaMax:
				lmm = LMM(self.Y[:,d], self.X, self.QS, restricted=True)
				lmm.fit(verbose=False)

			delta=lmm.v0/(lmm.v0+lmm.v1)
			variance = lmm.v0+lmm.v1
			self.delta_list.append(delta)
			self.var_list.append(variance)

			ret=lmm.get_fast_scanner().fast_scan(self.G,False) 
			wald = ret['effsizes1']/ret['effsizes1_se']  
			self.Z[:,d]= wald

		self.Z = self.Z.T
		centering = np.eye(res.num_snps) - (1/res.num_snps)*np.ones((res.num_snps,res.num_snps))
		#C1_naive = self.Z @ centering @ self.Z.T
		C1_naive = (self.Z-np.sum(self.Z,axis=1).reshape(-1,1)) @ (self.Z-np.sum(self.Z,axis=1).reshape(-1,1)).T
		print(C1_naive.shape)
		#ev, _ = np.absolute(np.linalg.eig(C1_naive))
		#np.savetxt('C1.csv', C1_naive)

		#ev = getEigenvalues_byR(self.Z)
		#np.savetxt('Spacing%.4fNsnps%degv_naive.txt'%(self.spacing,self.num_snps),ev)
		#fig = plt.figure()
		#plt.hist(self.delta_list, bins=50)
		#plt.savefig("Nsnps%dNtrait%d"%(self.num_snps,self.num_trait)+"_heritability_est.jpg")

		#fig = plt.figure()
		#plt.hist(self.var_list, bins=50)
		#plt.savefig("Nsnps%dNtrait%d"%(self.num_snps,self.num_trait)+"_variance_est.jpg")

	def fit_marginal_model(self):
		self.Z = np.zeros((self.num_snps, self.num_trait))
		
		for d in range(self.num_trait):
			Y_d = self.Y[:,d]
			
			xi_d= self.delta_list[d]
			sigma_d = np.sqrt(self.var_list[d])
			Sigma_d = xi_d*self.K + (1-xi_d)*np.eye(self.num_sample)
			#Sigma_d = .5*self.K + .5*np.eye(self.num_sample)
			Sigma_d_inv = np.linalg.inv(Sigma_d)
			inner = np.linalg.inv(self.X.T @ Sigma_d_inv @ self.X)
			P_d = Sigma_d_inv-Sigma_d_inv @ self.X @ inner @ self.X.T @ Sigma_d_inv
			
			for m in range(self.num_snps):
				G_m = self.G[:,m]
				#print(G_m.shape, P_d.shape, Y_d.shape)
				self.Z[m,d] = (G_m.T @ P_d @ Y_d)/np.sqrt(sigma_d**2 * G_m.T @ P_d @ G_m)
				#self.Z[m,d] = (G_m.T @ P_d @ Y_d)/np.sqrt( G_m.T @ P_d @ G_m)
				 
		now = datetime.now().time()
		print('******************************* Finished computing the Z-score.', " now the time is ", now)
		self.Z = self.Z.T

	def create_diagonastics(self, name=None):
		now = datetime.now().time()
		print('******************************* Creating diagonastic QQplots:', " now the time is ", now)

		fig = plt.figure()
		ax1 = fig.add_subplot(2, 2, 1)
		self.ave_col = np.mean(self.Z, axis=0)
		self.ave_col_scaled = self.ave_col * np.sqrt(self.num_trait)
		ax1.hist(self.ave_col_scaled, bins=100)

		ax2 = fig.add_subplot(2, 2, 2)
		sm.graphics.qqplot(self.ave_col_scaled, ax=ax2, line="45")

		ax3 = fig.add_subplot(2, 2, 3)
		self.ave_row = np.mean(self.Z, axis=1)
		self.ave_row_scaled = self.ave_row * np.sqrt(self.num_snps)
		ax3.hist(self.ave_row_scaled, bins=100)

		ax4 = fig.add_subplot(2, 2, 4)
		sm.graphics.qqplot(self.ave_row_scaled, ax=ax4, line="45")
		plt.savefig("Nsnps%dNtrait%d"%(self.num_snps,self.num_trait)+name+".jpg")

		#fig = plt.figure()
		#ax1 = fig.add_subplot(1, 2, 1)
		#ax1.hist(self.Z.flatten(), bins=100)
		#ax2 = fig.add_subplot(1, 2, 2)
		#sm.qqplot(self.Z.flatten(), ax=ax2, line="45")
		#plt.savefig("Nsnps%dNtrait%d"%(self.num_snps,self.num_trait)+name+"_Zscoresall.pdf")

	def step0_whole(self):
		self.C1_raw = None
		chromsome_unique, _ = np.unique(self.chromsome_list, return_inverse=True)
		centering = np.eye(self.num_sample) - (1/self.num_sample)*np.ones((self.num_sample,self.num_sample))
		now = datetime.now().time()
		self.C1_raw = self.G.T @ centering @ self.G
		self.C1 = cov2corr(self.C1_raw)
		self.C1_ = cov2corr(self.C1_raw+np.eye(self.num_snps)*self.perturbation)

	def step0(self):
		self.C1_raw = None
		chromsome_unique, _ = np.unique(self.chromsome_list, return_inverse=True)
		centering = np.eye(self.num_sample) - (1/self.num_sample)*np.ones((self.num_sample,self.num_sample))
		now = datetime.now().time()
		print('******************************* Start to obtain C1:', " now the time is ", now)
		counter = 1 
		total_block = len(chromsome_unique)

		print('There are %d blocks in total'%total_block)

		largest_evs = []
		smallest_evs = []

		for chrom in chromsome_unique:
			now = datetime.now().time()
			print('computing the %dth block...'%counter, " now the time is ", now)
			subset = np.where(self.chromsome_list==chrom)[0]
			G_m = self.G[:,subset]
			term = G_m.T @ centering @ G_m
			if counter == 1:
				self.C1_raw = term
			else:
				self.C1_raw = block_diag(self.C1_raw, term)
			#ev, _ = np.linalg.eig(term)
			#ev = getEigenvalues_byR(term)
			#log_ev = np.log(ev)
			#log_ev = np.sort(log_ev)[::-1].tolist()
			#fig = plt.figure()

			#ax1 = fig.add_subplot(1, 2, 1)
			#largest_ev = log_ev[:207]
			#ax1.hist(largest_ev, bins=20)
			#largest_evs += largest_ev
			#ax1.set_title('the 207 largest log-eigenvalues')

			#ax2 = fig.add_subplot(1, 2, 2)
			#smallest_ev = log_ev[207:]
			#ax2.hist(smallest_ev, bins=100)
			#smallest_evs += smallest_ev
			#ax2.set_title('the '+str(len(subset)-207)+' smallest log-eigenvalues')

			#plt.savefig("histogram_by_chromosome/chromosome"+str(counter)+".jpg")
			#plt.close('all')
			counter += 1

		#fig = plt.figure()
		#ax1 = fig.add_subplot(1, 2, 1)
		#ax1.hist(largest_evs, bins=20)
		#ax1.set_title('the largest log-eigenvalues')

		#ax2 = fig.add_subplot(1, 2, 2)
		#ax2.hist(smallest_evs, bins=20)
		#ax2.set_title('the smallest log-eigenvalues')
		#plt.savefig("histogram_by_chromosome/all_chromosome.jpg")
		#plt.close('all')

		self.C1 = cov2corr(self.C1_raw)
		self.V1 = cov2corr(self.C1_raw+np.eye(self.num_snps)*self.perturbation)
		#ev, _ = np.linalg.eig(self.V1)
		ev = getEigenvalues_byR(self.V1)
		log_ev = np.log(ev)
		log_ev = np.sort(log_ev)[::-1].tolist()
		np.savetxt('log_ev.txt', log_ev)
		#fig = plt.figure()
		#ax1 = fig.add_subplot(1, 2, 1)
		#ax1.hist(log_ev[:207*20], bins=20)
		#ax1.set_title('the largest ' +str(207*20)+ ' log-eigenvalues')
		
		#ax2 = fig.add_subplot(1, 2, 2)
		#ax2.hist(log_ev[207*20:], bins=20)
		#ax2.set_title('the smallest ' + str(self.num_snps-207*20) + ' log-eigenvalues')
		#plt.savefig("histogram_by_chromosome/C1_logev.jpg")
		#plt.close('all')

		now = datetime.now().time()
		print('******************************* Finished obtainning C1.', " now the time is ", now)

	def step1(self):
		now = datetime.now().time()
		print('******************************* Start to obtain C2:', " now the time is ", now)
		
		#u, s, _ = np.linalg.svd(self.V1)
		#L = u @ np.sqrt(np.diag(s)) @ u.T
		L = np.linalg.cholesky(self.V1)
		Z_update = self.Z @ np.linalg.inv(L).T
		centering = np.eye(self.num_snps) - (1/self.num_snps)*np.ones((self.num_snps,self.num_snps))
		self.C2_raw = Z_update @ centering @ Z_update.T
		self.C2 = cov2corr(self.C2_raw)
		now = datetime.now().time()
		print('******************************* Finished obtainning C2.', "now the time is ", now)

	def fit(self):
		self.fit_marginal_model()
		self.step0()
		self.step1()


if __name__ == '__main__':

	os.system('rm C1.txt C2.txt Z.txt G.txt')
	mouse_eQTLdata = Mouse_eQTLdata_Simulator(drop_file="dum.drop.txt", sample_id_file='mouseIDs_one.txt', 
		chrom_num=20, length=80, spacing= 0.08, allele_freq=0.2, raw_data_file='dum.drop_2.geno_true', num_trait=5000,
		H_0=True)

	# generate the expression traits
	mouse_eQTLdata.genome_dropping(seed=2)
	mouse_eQTLdata.raw2genotype()
	#mouse_eQTLdata.check_LD()
	#from sklearn.datasets import make_spd_matrix
	#mouse_eQTLdata.K = make_spd_matrix(208)+np.diag(range(1,209))*0.01
	#mouse_eQTLdata.QS = economic_qs(mouse_eQTLdata.K)

	mouse_eQTLdata.compute_GRM()
	mouse_eQTLdata.getG0_and_K0()
	mouse_eQTLdata.generate_parameters()
	mouse_eQTLdata.generate_expressionTrait()
	#mouse_eQTLdata.check_rowCov()
	#mouse_eQTLdata.check_colCov()
	#mouse_eQTLdata.check_colCorr()

	#generate the data under null hypothesis
	mouse_eQTLdata.genome_dropping(seed=2020)
	mouse_eQTLdata.raw_data_file='dum.drop_2020.geno_true'
	mouse_eQTLdata.raw2genotype()

	#find Z score and C2
	res = Correlation_Obtainer(mouse_eQTLdata=mouse_eQTLdata)
	
	res.run_limix()
	#res.create_diagonastics(name="by_Limix")
	#res.fit_marginal_model()
	#res.create_diagonastics(name="by_hand")
	#res.step0_whole()
	#res.step0()
	#res.step1()

	#centering = np.eye(res.num_snps) - (1/res.num_snps)*np.ones((res.num_snps,res.num_snps))
	#ev1,_ = np.linalg.eig(res.C2)
	#np.savetxt('ev1.txt', ev1)
	#ev2,_ = np.linalg.eig(res.Z@centering@res.Z.T)
	#np.savetxt('ev2.txt', ev2)

	#ev1,_ = np.linalg.eig(res.C2_raw)
	#np.savetxt('ev_C2raw.txt', ev1)
	
	#ev2,_ = np.linalg.eig(res.Z @ centering @ res.Z.T/res.num_snps)
	#fig = plt.figure()
	#ax1 = fig.add_subplot(1, 2, 1)
	#ax1.hist(ev1, bins=100)
	#ax2 = fig.add_subplot(1, 2, 2)
	#ax2.hist(ev2, bins=100)
	#plt.savefig("eigenvalues.pdf")
	#np.savetxt('ev_ZZT.txt', ev2)

	#ev3,_ = np.linalg.eig(res.C2)
	#np.savetxt('ev_C2.txt', ev3)
	#np.savetxt('G.txt', res.G)
	#np.savetxt('C1.txt', res.C1)
	'''
	np.savetxt('C2.txt', res.C2)
	np.savetxt('Z.txt', res.Z)
	np.savetxt('h.txt', res.hsquare_list)
	np.savetxt('sigma.txt', res.variance_list)
	#np.savetxt('Z.txt', mouse_eQTLdata.Z)

	now = datetime.now().time()
	print('******************************* Computing eigenvalues', "now the time is ", now)
	ev1,_ = np.linalg.eig(res.C2)
	np.savetxt('ev1.txt', res.ev1)
	ev2,_ = np.linalg.eig(res.Z*res.Z.T)
	np.savetxt('ev2.txt', res.ev2)


	now = datetime.now().time()
	print('C1 has shape ', res.C1.shape, ', C2 has shape ', res.C2.shape, ', Z has shape ', res.Z.shape, ', G has shape ', res.G.shape, "now the time is ", now)
	'''

	'''
	mouse_eQTLdata = Mouse_eQTLdata_Simulator(drop_file="dum.drop.txt", sample_id_file='mouseIDs_one.txt', 
		chrom_num=1, length=20, spacing=2, allele_freq=0.2, raw_data_file='dum.drop_1.geno_true', num_trait=1,
		H_0=True)

	mouse_eQTLdata.genome_dropping()
	mouse_eQTLdata.raw2genotype()
	print(mouse_eQTLdata.G)

	mouse_eQTLdata2 = Mouse_eQTLdata_Simulator(drop_file="dum.drop.txt", sample_id_file='mouseIDs_one.txt', 
		chrom_num=1, length=20, spacing=2, allele_freq=0.2, raw_data_file='dum.drop_1.geno_true', num_trait=1,
		H_0=True)

	mouse_eQTLdata2.genome_dropping()
	mouse_eQTLdata2.raw2genotype()
	print(mouse_eQTLdata2.G)
	'''





'''
		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax = sns.heatmap(self.K, center=0)
		ax.set_title('Cov Limix')
		plt.savefig("Cov_Limix.jpg")

		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax = sns.heatmap(self.K_man, center=0)
		ax.set_title('Cov Manual')
		plt.savefig("Cov_Manual.jpg")

		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax = sns.heatmap((self.K-self.K_man), center=0)
		ax.set_title('Cov Difference')
		plt.savefig("Cov_Difference.jpg")

		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax = sns.heatmap(self.K_corr, center=0)
		ax.set_title('Corr Limix')
		plt.savefig("Corr_Limix.jpg")

		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax = sns.heatmap(self.K_man_corr, center=0)
		ax.set_title('Corr Manual')
		plt.savefig("Corr_Manual.jpg")

		fig = plt.figure()
		ax = fig.add_subplot(1,1,1)
		ax = sns.heatmap((self.K_corr-self.K_man_corr), center=0)
		ax.set_title('Corr Difference')
		plt.savefig("Corr_Difference.jpg")
'''


'''
		self.Z = self.Z.T
		fig = plt.figure()
		plt.hist(delta_list, bins=50)
		plt.savefig("heritability_est.jpg")

		fig = plt.figure()
		plt.hist(var_list, bins=50)
		plt.savefig("variance_est.jpg")
'''
