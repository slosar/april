## This is a CDM cosmology with w

from LCDMCosmology import *
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad


class BinnedRhoCosmology(LCDMCosmology):
    def __init__(self, dz=0.2, zmax=1.0,zbinpow=1.4):
        ## two parameters: Om and h
        self.srzbins=np.arange(dz**(1./zbinpow),zmax**(1./zbinpow),dz**(1./zbinpow))
        self.zbinpow=zbinpow
        self.Nb=len(self.srzbins)
        self.rhovals=np.zeros(self.Nb)
        self.pnames=["r%i"%(i+1) for i in range(self.Nb)]
        LCDMCosmology.__init__(self)
        self.integrateOmega()
        
    ## my free parameters. 
    def freeParameters(self):
        rpars=[Parameter(name,self.rhovals[i]) for i,name in enumerate(self.pnames)]
        return rpars+LCDMCosmology.freeParameters(self)

    def updateParams(self,pars):
        ok=LCDMCosmology.updateParams(self,pars)
        if not ok:
            return False
        gotone=False
        for p in pars:
            if p.name in self.pnames:
                i=self.pnames.index(p.name)
                self.rhovals[i]=p.value
                gotone=True
        if gotone:
            self.integrateOmega()
        return True

    def integrateOmega(self):
        abins=np.hstack(([1.0],1./(1+self.srzbins**self.zbinpow),[1e-4]))
        r=np.hstack(([0.0],self.rhovals,[self.rhovals[-1]]))
        self.acut=0.5*(abins[:-1]+abins[1:])
        self.zcut=1./self.acut-1
        self.DEomega=interp1d(abins,r,kind='nearest')
                        

    ##for speed
    def Da_z(self,z):
        zl=0
        r=0
        for zt in self.zcut:
            if z>zt:
                r+=quad(self.DistIntegrand_a,1./(1+zt),1./(1+zl))[0]
                zl=zt
            else:
                r+=quad(self.DistIntegrand_a,1./(1+z),1./(1+zl))[0]
                break
        if self.Curv==0:
            return r
        elif (self.Curv>0):
            q=sqrt(self.Curv)
            ## someone check this eq
            ## Pure ADD has a 1+z fact, but have
            ## comoving one
            return sinh(r*q)/(q)
        else:
            q=sqrt(-self.Curv)
            return sin(r*q)/(q)
                

    
    ## this is relative hsquared as a function of a
    ## i.e. H(z)^2/H(z=0)^2
    def RHSquared_a(self,a):
        NuContrib=self.NuDensity.rho(a)/self.h**2
        return (self.Ocb/a**3+self.Omrad/a**4+NuContrib+(1.0-self.Om)*exp(self.DEomega(a)))


