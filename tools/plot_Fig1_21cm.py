#!/usr/bin/env python

#Error bars mainly got from  table 3,
#http://arxiv.org/pdf/1108.2635v1.pdf
#and MGS paper

from RunBase import *
import matplotlib.pyplot as plt
import matplotlib.ticker
import pylab 
import math as N



do21cm=True
plaw=0.25
#for division by log(1+z) use this
#plaw=-1


params1 = {'backend': 'pdf',
               'axes.labelsize': 20,
#               'text.fontsize': 18,
               'xtick.labelsize': 20,
               'ytick.labelsize': 20,
#               'legend.draw_frame': False,
               'legend.fontsize': 16,
               'lines.markersize': 6,
               'font.size': 20,
               'text.usetex': True}#
pylab.rcParams.update(params1)


#Planck best fit cosmology
T=LCDMCosmology(Obh2=0.022,Om=0.31,h=0.676)
#PLK-15
#T=LCDMCosmology(Obh2=0.02225,Om=0.3156,h=0.6727)

zLOWZ  = 0.32 
zCMASS = 0.57
zLyaA  = 2.33
zLyaC  = 2.40

z6dFGS   = 0.106
zMGS     = 0.15
zSDSS1   = 0.2
zSDSS2   = 0.35
zWiggleZ1=0.44
zWiggleZ2= 0.6
zWiggleZ3= 0.73

z_CMB = 1090.43

zCombBAO1 = 0.38
zCombBAO2 = 0.51
zCombBAO3 = 0.61

zEBQSO=1.52


rd_EHtoCAMB =153.19/149.28
rd_fid_DR12 = 147.78
rd_fid_DR7  = 151.84  

zl=arange(0,8,0.01)

def fixer(z):
    if plaw>0:
        return z**plaw
    else:
        return log(1.+z)


y1=[T.DaOverrd(z)/fixer(z)   for z in zl]
y2=[T.HIOverrd(z)*z/fixer(z) for z in zl]
y3=[T.DVOverrd(z)/fixer(z)   for z in zl]

fig = plt.figure(1,figsize=(9,6))
sfig= plt.figure(2,figsize=(9,6))
ax = fig.add_subplot(1,1,1)
sax= sfig.add_subplot(1,1,1)
plt.figure(1)

l1,=plt.plot(zl,y1,'r-',lw=2)
#l2,=plt.plot(zl,y3,'b-',lw=2)
l3,=plt.plot(zl,y2,'g-',lw=2)

if plaw>0:
    legend1=pylab.legend([l1,l3],["$%sD_%s(z)/r_d\\sqrt{z}$"%st for st in [('','M'),('z','H')]],loc="lower center")
else:
    legend1=pylab.legend([l1,l3],["$D_%s(z)/r_d\log(1+z)$"%st for st in ['A','H']],loc="lower center")


pylab.gca().add_artist(legend1)
def color_legend(leg):
        # """Color legend texts based on color of corresponding lines"""
        for line, txt in zip(leg.get_lines(), leg.get_texts()):
                txt.set_color(line.get_color())

legend1.draw_frame(False)
color_legend(legend1)

def ersys(x, y):
 	return sqrt(x**2 + y**2)


### Plotting -  Error bars 
def plot_errorbar(z,val, yerr=0, color='blue', fmt='x', markersize=0,label=None, empty=True, alpha=1):
    if empty:
        mfc='white'
        lw=1
    else:
        mfc=color
        lw=2
    pylab.errorbar(z,val/fixer(z), yerr=yerr/fixer(z), color=color, fmt=fmt, markersize=markersize, lw=lw, capthick=lw,capsize=2+2*lw,markerfacecolor=mfc, alpha=alpha)
    if label>0:
        if (mfc=='white'):
            pylab.plot ([],[],fmt,color=color,label=label,markersize=markersize,markerfacecolor=mfc)
        else:
            pylab.plot ([],[],fmt,color=color,label=label,markersize=markersize)

pylab.legend(loc="lower right")

fmt1 = '^'
fmt2 = 's'
empty1= True
empty2= False
alpha= 1.0

######## CURRRENT ERRORS
#Errorbars from DR12 Full-shape
fact = (300000./rd_fid_DR12)

#666/148.651,    yerr=25/148.651
plot_errorbar(z6dFGS,    2.97*rd_EHtoCAMB,   yerr=rd_EHtoCAMB*0.015/0.336**2,  color ='blue', fmt='o', markersize=6, label="current\ generation", empty=True,alpha=alpha)
plot_errorbar(zMGS,      4.464,    yerr=0.168,               color ='blue', fmt='o', markersize=6, alpha=alpha)
plot_errorbar(zSDSS1,    5.2493*rd_EHtoCAMB, yerr=rd_EHtoCAMB*0.0061/0.1905**2,color ='blue', fmt='o', markersize=6, alpha=alpha)
plot_errorbar(zSDSS2,    1348./rd_fid_DR12, yerr=26./rd_fid_DR12 ,color ='blue', fmt='o', markersize=6, alpha=alpha)
plot_errorbar(zWiggleZ1, 1695./rd_fid_DR12 ,yerr=82./rd_fid_DR12 ,color ='blue', fmt='o', markersize=6,  alpha=alpha)
plot_errorbar(zWiggleZ2, 2194./rd_fid_DR12 ,yerr=100./rd_fid_DR12 ,color ='blue', fmt='o', markersize=6,  alpha=alpha)
plot_errorbar(zWiggleZ3, 2486./rd_fid_DR12 ,yerr=85./rd_fid_DR12 ,color ='blue', fmt='o', markersize=6,  alpha=alpha)
plot_errorbar(zCombBAO1,  1512.4/rd_fid_DR12,     yerr=ersys(22.5, 11.0)/rd_fid_DR12,       #
       	color ='blue', fmt='o', markersize=6)
plot_errorbar(zCombBAO2,  1975.2/rd_fid_DR12,     yerr=ersys(26.6, 14.1)/rd_fid_DR12,
		color ='blue', fmt='o', markersize=6)
plot_errorbar(zCombBAO3,  2306.7/rd_fid_DR12,  	 yerr=ersys(33.2, 16.7)/rd_fid_DR12, #
       	color ='blue', fmt='o', markersize=6)
plot_errorbar(zCombBAO1,  fact*zCombBAO1/81.21,       yerr=fact*zCombBAO1*ersys(2.17, 0.97)/(81.21)**2,  #
		color ='blue', fmt='o', markersize=6)
plot_errorbar(zCombBAO2,  fact*zCombBAO2/90.90,       yerr=fact*zCombBAO2*ersys(2.07, 1.08)/(90.90)**2,    
		color ='blue', fmt='o', markersize=6)
plot_errorbar(zCombBAO3,  fact*zCombBAO3/98.96,       yerr=fact*zCombBAO3*ersys(2.21, 1.18)/(98.96)**2, 
		color ='blue', fmt='o', markersize=6)
#plot_errorbar(zEBQSO,  3855/rd_fid_DR12, yerr=170/rd_fid_DR12,  color ='blue', fmt='o', markersize=6, alpha=alpha)
plot_errorbar(zLyaA,  37.77,  yerr=2.13,  color ='blue', fmt='o', markersize=6)
plot_errorbar(zLyaA,  9.07*zLyaA,       yerr=0.31*zLyaA,       color ='blue', fmt='o', markersize=6)
plot_errorbar(zLyaC,  35.7,   yerr=1.5,    color ='blue', fmt='o', markersize=6)
plot_errorbar(zLyaC,  9.01*zLyaC,        yerr=0.32*zLyaC,        color ='blue', fmt='o', markersize=6)



def plotFutureErrorbar (z,dae,he,label,fmt,color,markersize=6):
    da=T.DaOverrd(z)
    hi=T.HIOverrd(z)*z
    plot_errorbar(z,  da,   yerr=da*dae,   color=color, fmt=fmt, markersize=markersize,label=label, empty=False)
    plot_errorbar(z,  hi,        yerr=hi*he, color=color, fmt=fmt, markersize=markersize, empty=False)


desi="""0.15 & 2.78 & 5.34 & 1.95 & 5.24 & 13.79 & 0.23 & 376 &  50 &   8 & 7.51 & 3.60 \\ 
0.25 & 1.87 & 3.51 & 1.30 & 3.24 & 8.19 & 0.56 & 347 & 125 &  23 & 5.24 & 2.55 \\ 
0.35 & 1.45 & 2.69 & 1.00 & 2.58 & 6.35 & 0.99 & 291 & 222 &  31 & 4.44 & 2.17 \\ 
0.45 & 1.19 & 2.20 & 0.82 & 2.36 & 5.74 & 1.46 & 285 & 332 &  31 & 3.92 & 1.91 \\ 
0.55 & 1.01 & 1.85 & 0.70 & 2.42 & 5.90 & 1.94 & 431 & 448 &  32 & 3.31 & 1.60 \\ 
0.65 & 0.87 & 1.60 & 0.60 & 2.58 & 6.34 & 2.42 & 722 & 563 &  34 & 2.80 & 1.34 \\ 
0.75 & 0.77 & 1.41 & 0.53 & 2.77 & 6.85 & 2.87 & 1112 & 675 &  37 & 2.47 & 1.18 \\ 
0.85 & 0.76 & 1.35 & 0.52 & 2.05 & 5.17 & 3.29 & 1333 & 471 &  44 & 2.34 & 1.11 \\ 
0.95 & 0.88 & 1.42 & 0.58 & 1.03 & 2.76 & 3.67 & 1401 &  91 &  50 & 2.34 & 1.13 \\ 
1.05 & 0.91 & 1.41 & 0.59 & 0.82 & 2.24 & 4.01 & 1469 &  11 &  56 & 2.32 & 1.12 \\ 
1.15 & 0.91 & 1.38 & 0.58 & 0.75 & 2.05 & 4.31 & 1483 &   0 &  62 & 2.30 & 1.12 \\ 
1.25 & 0.91 & 1.36 & 0.58 & 0.69 & 1.86 & 4.57 & 1421 &   0 &  69 & 2.32 & 1.14 \\ 
1.35 & 1.00 & 1.46 & 0.64 & 0.53 & 1.42 & 4.80 & 1120 &   0 &  75 & 2.45 & 1.26 \\ 
1.45 & 1.17 & 1.66 & 0.74 & 0.38 & 1.00 & 4.99 & 775 &   0 &  81 & 2.71 & 1.47 \\ 
1.55 & 1.50 & 2.04 & 0.93 & 0.25 & 0.63 & 5.15 & 460 &   0 &  83 & 3.22 & 1.89 \\ 
1.65 & 2.36 & 3.15 & 1.45 & 0.13 & 0.33 & 5.29 & 179 &   0 &  80 & 4.63 & 3.06 \\ 
1.75 & 3.62 & 4.87 & 2.23 & 0.08 & 0.19 & 5.40 &  49 &   0 &  77 & 7.17 & 5.14 \\ 
1.85 & 4.79 & 6.55 & 2.98 & 0.06 & 0.13 & 5.49 &   0 &   0 &  74 & 10.26 & 7.66 \\ """
desi=np.array(([[float(x) for x in s.split('&')[:3]] for s in desi.replace('\\','').split('\n')]))
lab='DESI'
for z,dae,he in desi:
    plotFutureErrorbar(z,dae/100,he/100,lab,'o','green')
    lab=None
plt.figure(2)
plt.plot(desi[:,0],desi[:,1],'g-', lw=2,label='DESI')
plt.plot(desi[:,0],desi[:,2],'g--',lw=2)
plt.figure(1)

euclid="""0.65 & 1.23 & 1.89 & 0.79 & 0.75 & 2.24 & 2.59 & 1100 \\ 
0.75 & 0.83 & 1.42 & 0.56 & 1.69 & 5.03 & 3.07 & 2950 \\ 
0.85 & 0.74 & 1.27 & 0.50 & 1.90 & 5.60 & 3.52 & 3800 \\ 
0.95 & 0.71 & 1.19 & 0.48 & 1.75 & 5.11 & 3.93 & 3900 \\ 
1.05 & 0.70 & 1.14 & 0.46 & 1.55 & 4.48 & 4.29 & 3775 \\ 
1.15 & 0.70 & 1.12 & 0.46 & 1.35 & 3.85 & 4.62 & 3525 \\ 
1.25 & 0.70 & 1.10 & 0.46 & 1.17 & 3.31 & 4.90 & 3250 \\ 
1.35 & 0.73 & 1.11 & 0.47 & 0.98 & 2.74 & 5.14 & 2850 \\ 
1.45 & 0.78 & 1.16 & 0.50 & 0.78 & 2.15 & 5.35 & 2350 \\ 
1.55 & 0.87 & 1.24 & 0.55 & 0.59 & 1.62 & 5.52 & 1850 \\ 
1.65 & 1.01 & 1.40 & 0.63 & 0.43 & 1.16 & 5.66 & 1375 \\ 
1.75 & 1.23 & 1.64 & 0.75 & 0.30 & 0.80 & 5.78 & 975 \\ 
1.85 & 1.61 & 2.07 & 0.97 & 0.20 & 0.52 & 5.88 & 650 \\ 
1.95 & 2.32 & 2.90 & 1.38 & 0.12 & 0.31 & 5.95 & 400 \\ 
2.05 & 5.32 & 6.39 & 3.11 & 0.04 & 0.12 & 6.01 & 150 \\ """
euclid=np.array(([[float(x) for x in s.split('&')[:3]] for s in euclid.replace('\\','').split('\n')]))
lab='EUCLID'
for z,dae,he in euclid:
    plotFutureErrorbar(z,dae/100,he/100,lab,'^','cyan')
    lab=None
plt.figure(2)
plt.plot(euclid[:,0],euclid[:,1],'-',lw=2,color='cyan', label='EUCLID')
plt.plot(euclid[:,0],euclid[:,2],'--',lw=2,color='cyan')
plt.figure(1)


wfirst="""1.05 & 1.51 & 2.72 & 1.03 & 4.37 & 12.60 & 0.57 & 10623 \\ 
1.15 & 1.43 & 2.56 & 0.98 & 4.50 & 12.85 & 0.62 & 11776 \\ 
1.25 & 1.35 & 2.42 & 0.92 & 5.00 & 14.13 & 0.65 & 13877 \\ 
1.35 & 1.29 & 2.30 & 0.88 & 5.33 & 14.90 & 0.69 & 15527 \\ 
1.45 & 1.24 & 2.21 & 0.85 & 5.58 & 15.42 & 0.71 & 16890 \\ 
1.55 & 1.23 & 2.16 & 0.84 & 5.04 & 13.79 & 0.74 & 15759 \\ 
1.65 & 1.25 & 2.15 & 0.84 & 4.15 & 11.23 & 0.76 & 13305 \\ 
1.75 & 1.28 & 2.16 & 0.86 & 3.33 & 8.94 & 0.77 & 10918 \\ 
1.85 & 1.33 & 2.19 & 0.88 & 2.61 & 6.94 & 0.78 & 8697 \\ 
1.95 & 1.41 & 2.27 & 0.93 & 1.99 & 5.25 & 0.79 & 6718 \\ 
2.05 & 2.51 & 3.52 & 1.57 & 0.47 & 1.23 & 0.80 & 1610 \\ 
2.15 & 2.60 & 3.62 & 1.62 & 0.44 & 1.14 & 0.81 & 1509 \\ 
2.25 & 2.74 & 3.78 & 1.70 & 0.40 & 1.02 & 0.81 & 1368 \\ 
2.35 & 3.02 & 4.09 & 1.86 & 0.33 & 0.85 & 0.81 & 1156 \\ 
2.45 & 3.38 & 4.52 & 2.08 & 0.28 & 0.70 & 0.81 & 960 \\ 
2.55 & 3.87 & 5.11 & 2.36 & 0.23 & 0.57 & 0.81 & 781 \\ 
2.65 & 4.52 & 5.90 & 2.75 & 0.18 & 0.45 & 0.81 & 626 \\ 
2.75 & 5.41 & 6.99 & 3.27 & 0.14 & 0.35 & 0.81 & 490 \\ """
wfirst=np.array(([[float(x) for x in s.split('&')[:3]] for s in wfirst.replace('\\','').split('\n')]))
lab='WFIRST'
for z,dae,he in wfirst:
    plotFutureErrorbar(z,dae/100,he/100,lab,'>','orange')
    lab=None

plt.figure(2)
plt.plot(wfirst[:,0],wfirst[:,1],'-',lw=2,color='orange', label='WFIRST')
plt.plot(wfirst[:,0],wfirst[:,2],'--',lw=2,color='orange')
plt.figure(1)




if do21cm:
## run printbaoerrors.py in RSDCalib
    obuljen=np.loadtxt('petite_errors.txt')
    lab='PUMA'
    for z,dae1,he1 in obuljen:
        plotFutureErrorbar(z,dae1*0.01,he1*0.01,lab,'*','red',markersize=8)
        lab=None
    
plt.figure(2)
plt.plot(obuljen[:,0],obuljen[:,1],'-',lw=2,color='red', label='21cm Stage II')
plt.plot(obuljen[:,0],obuljen[:,2],'--',lw=2,color='red')
plt.figure(1)


ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())


plt.yticks(list(range(0, 50, 10)))



if plaw>0:
    pylab.ylabel("${\\rm distance}/r_d\\sqrt{z}$")
    pylab.xlabel("$z$")
    pylab.legend(loc=(-0.02,0.65), numpoints=1, frameon=False)
    pylab.ylim(0,45)
    pylab.xlim(0.00,6.5)
    pylab.tight_layout()

    if do21cm:
        pylab.savefig("Fig1_21cm.pdf")
    else:
        pylab.savefig("Fig1_21cmSans.pdf")
else:
    pylab.legend(loc='lower left', numpoints=1)
    pylab.ylabel("$D(z)/r_d \log(1+z)$")
    pylab.ylim(15,35)
    pylab.xlim(0.08,8.0)
    pylab.xlabel("$z$")
    plt.semilogx()
    pylab.tight_layout()
    pylab.savefig("Fig1_21cm_v2.pdf")


plt.figure(2)
plt.xlim(0,6.5)
plt.plot([],[],'k-',label='$D_A$')
plt.plot([],[],'k--',label='$D_H$')
plt.legend(frameon=False)



plt.xlabel("$z$")
plt.ylabel("${\\rm precision}\ [\%]$")
plt.tight_layout()
plt.savefig("baoprecision.pdf")
    
#pylab.show()

