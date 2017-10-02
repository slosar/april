#!/usr/bin/env python
from RunBase import *
import numpy.lib.recfunctions as recfunctions

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


class TwoNeffImportanceSampler(CosmoMCImportanceSampler):
    def __init__ (self, inroot, outroot, like):
        #CosmoMCImportanceSampler.__init__(self,inroot,outroot,like)
        self.dtype=[]
        self.dtype.append(('weight','f8'))
        self.dtype.append(('nloglike','f8'))
        op=open(outroot+".paramnames",'w')
        for line in open(inroot+".paramnames").readlines():
            op.write(line)
            name=line.split()[0]
            self.dtype.append((name,'f8'))
        op.write("NnuLSS\n")
        op.close()
        self.inroot=inroot
        self.outroot=outroot
        self.ofset=None
        self.names=[d[0] for d in self.dtype]
        self.like=like


    def thetaOfh(self,h):
        self.like.updateParams([Parameter("h",h)])#thetap/thetao)])
        return self.like.theory().CMBSimpleVec()[2]
    
    def importanceSample(self,fname,ofname):
        chain=np.loadtxt(fname,dtype=self.dtype)
        ## generate Neff 
        ##
        N=len(chain)
        chain=recfunctions.append_fields(chain,'nnuLSS',np.random.uniform(1,5,N),
                                                 usemask=False)

        llike=[]
        for i,line in enumerate(chain):
            if i%1000==0:
                print "%i/%i..."%(i,N)
            if True:
                self.like.updateParams (self.cosmomc2april(line,True))
            else:
                self.like.updateParams (self.cosmomc2april(line,False))
                h=self.like.theory().h
                thetao=self.thetaOfh(h)
                self.like.updateParams (self.cosmomc2april(line,True))
                thetap=self.thetaOfh(h)
                if thetap>thetao:
                    hlow=h
                    hhigh=h*1.4
                else:
                    hlow=h*0.6
                    hhigh=h

                while True:
                    hmid=(hlow+hhigh)/2
                    thetam=self.thetaOfh(hmid)
                    if abs(thetam/thetao-1)<1e-4:
                        break
                    if thetam>thetao:
                        hlow=hmid
                    else:
                        hhigh=hmid
            # found my h
            llike.append(self.like.loglike())
        llike=np.array(llike)
        if self.ofset is None:
            self.ofset=llike.max()
        rewe=np.exp(llike-self.ofset)
        chain['weight']*=rewe
        print "min/mean/max weight=",rewe.min(), rewe.mean(), rewe.max()
        np.savetxt(ofname,chain)

        
    def cosmomc2april(self,line,NeffLSS=False):
        plist=[Parameter("Obh2",line['omegabh2']),
               Parameter("Om",line['omegam*']),
               Parameter("h",line['H0*']/100.)]
        if "nnu" in self.names:
            if NeffLSS:
                plist.append(Parameter("Nnu",line['nnuLSS']))
            else:
                plist.append(Parameter("Nnu",line['nnu']))
        if "w" in self.names:
            plist.append(Parameter("w",line['w']))
        if "wa" in self.names:
            plist.append(Parameter("wa",line['wa']))
        if "omegak" in self.names:
            plist.append(Parameter("Ok",line['omegak']))
        if "mnu" in self.names:
            plist.append(Parameter("mnu",line['mnu']))
        return plist




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
        droot=inroot
    elif num/9==1:
        inroot="plikHM_TTTEEE_lowTEB"
        droot=inroot
    elif num/9==2:
        inroot="plikHM_TTTEEE_lowTEB_post_lensing"
        droot="plikHM_TTTEEE_lowTEB"
    #s=CosmoMCImportanceSampler("/data/anze/Planck/base_nnu/%s/base_nnu_%s"%(droot,inroot),
    #                       "isamp/%s_patchy%i"%(inroot,n),L)
    s=TwoNeffImportanceSampler("/data/anze/Planck/base_nnu/%s/base_nnu_%s"%(droot,inroot),
                           "isamp_2N_d/%s_patchy%i"%(inroot,n),L)

    s.run()


#dochain(0)
from multiprocessing import Pool
pool = Pool(processes=4)              
pool.map(dochain, range(27)) #just do the last one



