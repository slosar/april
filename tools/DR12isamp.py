#!/usr/bin/env python
from RunBase import *

class DR12Isotropic(CompositeLikelihood):
    ## See 1607.03155
    def __init__(self, alpha_lowz, alphaerr_lowz, alpha_highz, alphaerr_highz, infostr):
        obh2=0.022; Om=0.31; h=0.676; mnu=0.0
        fidtheory=LCDMCosmology(obh2,Om,h,mnu)
        lowz=0.38
        highz=0.61
        rd=fidtheory.rd
        DVL=rd*fidtheory.DVOverrd(lowz)*alpha_lowz
        DVLE=rd*fidtheory.DVOverrd(lowz)*alphaerr_lowz
        DVH=rd*fidtheory.DVOverrd(highz)*alpha_highz
        DVHE=rd*fidtheory.DVOverrd(highz)*alphaerr_highz

        CompositeLikelihood.__init__ (self,[GaussBAODVLikelihood(infostr+"_low",lowz,DVL,DVLE,fidtheory),
                                       GaussBAODVLikelihood(infostr+"_high",highz,DVH,DVHE,fidtheory)])


def getAlpha(z,i):
    alpha,_,_,p2,m2,_=map(float,open('DRMock/result_2_60_BAOfit2D_patchy_win_z%i_postrecon_iso_binsize2_%i.dat'%(z,i)).readlines()[15].split())
    return alpha, (p2+m2)/2/2  ## return 1 sigma average

def dochain(num):
        
    T=ParseModel("NeffLCDM")
    n=num%9+1
    a1,ae1=getAlpha(1,n)
    a2,ae2=getAlpha(3,n)
    print a1,ae1,a2,ae2
    L=DR12Isotropic(a1,ae1,a2,ae2,"TEST%i"%num)
    L.setTheory(T)
    if num/9==0:
        inroot="plikHM_TT_lowTEB"
    else:
        inroot="plikHM_TTTEEE_lowTEB"
    s=CosmoMCImportanceSampler("/data/anze/Planck/base_nnu/%s/base_nnu_%s"%(inroot,inroot),
                           "isamp/%s_patchy%i"%(inroot,n),L)
    s.run()


from multiprocessing import Pool
pool = Pool(processes=4)              
pool.map(dochain, range(18))



