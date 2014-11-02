#!/usr/bin/env python

from scipy import *
from cosmich import *
import pylab, sys
import string
import math as N

#params1 = {'backend': 'pdf',
#               'axes.labelsize': 20,
#               'text.fontsize': 18,
#               'xtick.labelsize': 20,
#               'ytick.labelsize': 20,
#               'legend.draw_frame': False,
#               'legend.fontsize': 16,
#               'lines.markersize': 6,
#               'font.size': 20,
#               'text.usetex': True}#
#pylab.rcParams.update(params1)


#-----------
#dire = 'final_chains/'
dire='/astro/u/jvazquez/work/SimpleMC/trunk/chains/EDE/'
#dire = '/gpfs01/astro/workarea/jvazquez/chains/'
#dire = '/astro/u/jvazquez/BOSS/cosmomc_july_14/chains/'

#At the moment, 2D-plot is valid for only one couple of parameters, and may use several models/datasets
#At the moment, 1D-plot is valid for only one model, and may use several parameters/datasets

model_1D = 'EarlyDE'
model_2D  = ['EarlyDE','EarlyDE_rd_DE']

extra    = 'phy'
datasetl  = ['phy_BBAO+SN+Planck','rd_DE_phy_BBAO+SN+Planck']
#datasetl  = ['phy_BBAO','BBAO']
#-----------

Plot_2D   = 'False'
param_x   = 'Ode'
param_y   = 'Om'

NBins_2D = 20

xrange_2D = 'False'
xmin, xmax = 0.5, 1
ymin, ymax = 0, 0.8

#-----------

Plot_1D   = 'True'
params_1D      =['Ode']
#params_1D     = ['omegam*','H0*','omegabh2']
#params_1D    = ['6DF*','DR11LOWZ*','MGS*','DR11CMASS*','Lya_Auto_Busca*','Lya_Cross_Andreu*']

NBins_1D = 25

xrange   = 'False'
xmin_1, xmax_1 = 0, 3
xmin_2, xmax_2 = 0, 5
xmin_3, xmax_3 = 0, 5
xmin_4, xmax_4 = 3, 5.5
xmin_5, xmax_5 = 3.5, 5.5

name_fig  = 'Plot_2'

#-------------

def colour(x):
    if x==1: return 'black'
    if x==2: return 'blue'
    if x==3: return 'red'
    if x==4: return 'magenta'
    if x==5: return 'cyan'
    if x==6: return 'orange'
    if x==7: return 'green'
    if x==8: return 'yellow'
    if x>8:  print("Increased colouring") 


def color_legend(leg):
        # """Color legend texts based on color of corresponding lines"""
        for line, txt in zip(leg.get_lines(), leg.get_texts()):
                txt.set_color(line.get_color())

def label_m(name):
    with open('/gpfs01/astro/workarea/jvazquez/chains/LCDM_phy_BBAO+CMBP.paramnames') as inf:
      for line in inf:
        parts = line.split()
        if len(parts) > 1:   
            if name in parts:
                return '$'+str(parts[1])+'$'   

def label_n(name):
    if 'omegabh2' in name: return   '$\Omega_b h^2$'	
    if 'omegam*'  in name: return   '$\Omega_m$'
    if 'H0*'      in name: return   '$H_0$'
    if 'DR11LOWZ*' in name: return 'LOWZ'
    if 'DR11CMASS*' in name: return 'DR11CMASS'
    if 'Lya_Auto_Busca*' in name: return 'DR11LyaAuto'
    if 'Lya_Cross_Andreu*' in name: return 'DR11LyaCross'
    if 'MGS*' in name: return 'MGS'
    if '6DF*' in name: return '6DF'
 
    if 'Obh2' in name: return '$\\Omega_bh^2$'
    if 'Om'   in name: return '$\\Omega_m$'
    if 'h'    in name: return '$h$'
    if 'Ok'   in name: return '$\\Omega_k$'
    if 'w'    in name: return '$w$'
    if 'wa'   in name: return '$w_a$'
    if 'Pr'   in name: return '$c/(H_0 r_d)$'
    if 'mnu'  in name: return '$\\sum m_\\nu$'
    if 'beta' in name: return '$\\beta$'
    if 'Sp1'  in name: return '$w(z_1)$'
    if 'Sp2'  in name: return '$w(z_2)$' 
    if 'Sp3'  in name: return '$w(z_3)$'
    if 'Sp4'  in name: return '$w(z_4)$'
    if 'Lambda' in name: return '$\\lambda$'
    if 'Omr'  in name:  return '$\\Omega_{r_x}$'
    if 'Ode'  in name: return '$\\Omega^e_{\\rm de}$'
    if 'DR11LOWZ*' in name: return 'LOWZ'
   

def cosmodata2(datsets):
	cosmodata=''
	if 'phy' in datasets:
		cosmodata= 'SimpleMC'
 	else:
		cosmodata= 'CosmoMC'
	return cosmodata

def cosmodata(datasets):
        cosmodata=''
        if 'BBAO' in datasets:
                if '+BBAO' in datasets:
                   cosmodata=cosmodata+'+'
                cosmodata=cosmodata+'BAO'
        if 'GBAO' in datasets:
                if '+GBAO' in datasets:
                   cosmodata=cosmodata+'+'
                cosmodata=cosmodata+'Gal BAO'
        if 'LBAO' in datasets:
                if '+LBAO' in datasets:
                   cosmodata=cosmodata+'+'
                cosmodata=cosmodata+'Lya BAO'
        if 'SN' in datasets:
                if '+SN' in datasets:
                   cosmodata=cosmodata+'+'
                cosmodata=cosmodata+'SNe'
        if "Planck"  in datasets:
                if '+Planck' in datasets:
                   cosmodata=cosmodata+'+'
                cosmodata=cosmodata+'Planck'           
        if '6dFGS' in datasets:
                if '+6dFGS' in datasets:
                   cosmodata=cosmodata+'+'
                cosmodata=cosmodata+'6dFGS'
        return cosmodata

a=0
if 'True' in Plot_2D:
  for model in model_2D:    
     for datasets in datasetl:
        a+=1             

	C=cosmochain(dire+ model+'_'+extra+'_'+datasets)
        C.Plot2D(param_x, param_y, filled=colour(a), label=cosmodata(datasets), N=NBins_2D)

  if 'True' in xrange_2D:   
     pylab.xlim(xmin,xmax)
     pylab.ylim(ymin,ymax)

  if 'EarlyDE' in model_2D:
      Ode = arange(0,0.5,0.01)
      Om=[0.3/(1+Odes) for Odes in Ode]
      Om2=[0.3*(1-Odes) for Odes in Ode]

      pylab.plot(Ode,Om,'r--')
      pylab.plot(Ode,Om2,'r-')

  leg=pylab.legend()
  leg.draw_frame(False)		# No box & colour legend
  color_legend(leg)

  pylab.xlabel(label_n(param_x))
  pylab.ylabel(label_n(param_y))
  pylab.savefig(name_fig+'_2D.pdf')
  pylab.show()





a=0
if 'True' in Plot_1D:
   for datasets in datasetl:
      a+=1
      b=0
      for params in params_1D:
        b+=1

	pylab.subplot(1,len(params_1D),b)
	C=cosmochain(dire+ model_1D+'_'+datasets,'auto')	
        #C=cosmochain(dire+ model_1D+'_'+extra+'_'+datasets,'auto')

	if 'phy' in datasets:
	    if 'H0' in params:
		C[params]= 100*C[params]	
#	    if '*' in params:
#		C[params]=-2*C[params]	
        xx,yy=C.GetHisto(params,NormPeak=True,nbins=NBins_1D)	
        pylab.plot(xx, yy, colour(a), label=cosmodata2(datasets))

        if 'True' in xrange:
            xmin="xmin_"+str(b)
            xmax="xmax_"+str(b) 
            pylab.xlim(eval(xmin),eval(xmax))

        pylab.xlabel(label_n(str(params)))
        pylab.ylabel('prob.')

   pylab.tight_layout() 
   pylab.legend(loc='upper right')
   pylab.savefig(name_fig+'_1D.pdf')
   pylab.show()
   
else:
   print 'Nothing else to do'


