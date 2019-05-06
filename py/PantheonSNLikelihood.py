#
# This module calculates likelihood for the compressed SN.
#

from BaseLikelihood import *
import numpy as np
import scipy.linalg as la


class PantheonSNLikelihood(BaseLikelihood):
    def __init__ (self):
        ## first read data file
        self.name_="PantheonSN"
        da=[x.split() for x in open('data/pantheon_lcparam_full_long_zhel.txt').readlines()[1:]]
        self.zcmb = np.array([float(line[1]) for line in da])
        self.zhelio = np.array([float(line[2]) for line in da])
        self.mag = np.array([float(line[4]) for line in da])
        self.dmag = np.array([float(line[5]) for line in da])
        N=len(self.mag)
        self.syscov=np.loadtxt('data/pantheon_sys_full_long.txt',skiprows=1).reshape((N,N))
        self.cov=np.copy(self.syscov)
        self.cov[np.diag_indices_from(self.cov)]+=self.dmag**2
        self.xdiag=1/self.cov.diagonal() ## diagonal before marginalising constant
        ## add marginalising over a constant
        self.cov+=3**2
        
        self.icov=la.inv(self.cov)
        
    def loglike(self):
        tvec = self.mag-np.array([self.theory_.distance_modulus(z) for z in self.zcmb])

        #print (tvec[:10])
        ## first subtract a rought constant to stabilize marginaliztion of
        ## intrinsic mag.
        tvec-= (tvec*self.xdiag).sum() / (self.xdiag.sum())
        #print (tvec[:10])
        chi2 = np.einsum('i,ij,j',tvec,self.icov,tvec)
        #print ("chi2=",chi2)
        return -chi2/2
    
    
        
        
        
    
