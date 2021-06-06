#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2020-12-07 10:42:03 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

import uproot as up
import scipy as sci
import scipy.optimize as opt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import ticker
from collections import namedtuple
from sys import path
import time,math,sys,itertools

sys.path.insert(0,'/home/trottar/bin/python/')
import root2py as r2p

kinematics = sys.argv[1]

# kinematics="pi_n_18on275"
# kinematics="pi_n_10on100"
# kinematics="pi_n_5on100"
# kinematics="pi_n_5on41"
# kinematics="k_lambda_5on100"
# kinematics="k_lambda_18on275"

rootName="/home/trottar/ResearchNP/JLEIC/USERS/trottar/OUTPUTS/%s.root" % kinematics

tree = up.open(rootName)["Evnts"]
branch = r2p.pyBranch(tree)

# Define phyisics data
s_e = branch.findBranch("invts","s_e")
s_q = branch.findBranch("invts","s_q")
Q2 = branch.findBranch("invts","Q2")
xBj = branch.findBranch("invts","xBj")
t = -branch.findBranch("invts","tSpectator")
tPrime = -branch.findBranch("invts","tPrime")
y_D = branch.findBranch("invts","y_D")
nu = branch.findBranch("invts","nu")
TwoPdotk = branch.findBranch("invts","TwoPdotk")
# xBj = tree.array("TDIS_xbj")
TDIS_xbj = tree.array("TDIS_xbj")
sigma_dis = tree.array("sigma_dis")*(1e-5) # used to compare to HERA cross-section
TDIS_y = tree.array("TDIS_y")
ppix_Lab = tree.array("ppix_Lab") # Long pion momentum
ppiy_Lab = tree.array("ppiy_Lab")
ppiz_Lab = tree.array("ppiz_Lab")
EpiE_Lab = tree.array("EpiE_Lab")
pprx_inc = tree.array("pprx_inc") # Long proton momentum
ppry_inc = tree.array("ppry_inc")
pprz_inc = tree.array("pprz_inc")
EprE_inc = tree.array("EprE_inc")
EnE_Lab = tree.array("EnE_Lab")
# pprx_inc = tree.array("pnx_Lab") # Long proton momentum
# ppry_inc = tree.array("pny_Lab")
# pprz_inc = tree.array("pnz_Lab")
# EprE_inc = tree.array("EnE_Lab")
y = tree.array("TDIS_y")
fpi = tree.array("fpi")
f2N = tree.array("f2N")
xpi = tree.array("xpi")
# xpi = TDIS_xbj
ypi = tree.array("ypi")
tpi = tree.array("tpi")
escat = tree.array("EScatRest")
pprz_inc = tree.array("pprz_inc")
ppiz_Lab = tree.array("ppiz_Lab")

# binx  = 0.001
binx  = 0.01
binQ2 = 10.0

print("\nx binning:",binx)
print("Q^2 binning:",binQ2,"\n")

# Create cut dictionary
cutDict = {}

# More coarse, set bin size
xpiarray = np.arange(binx/2,1.0,binx).tolist()

# Less coarse
for i,x in enumerate(xpiarray):
    xpitmp = '{"xpicut%i" : ((%0.5f <= xpi) & (xpi <= %0.5f))}' % (i,xpiarray[i]-binx/2,xpiarray[i]+binx/2) # for proper binning
    # (i,xpiarray[i]-binx/2,xpiarray[i]+binx/2) # no binning
    # (i,xpiarray[i]-binx/20,xpiarray[i]+binx/20) # for proper binning
    print('{"xpicut%i" : ((%0.5f <= xpi) & (xpi <= %0.5f))}' % (i,xpiarray[i]-binx/2,xpiarray[i]+binx/2))
    cutDict.update(eval(xpitmp))
c = r2p.pyPlot(cutDict)

xpicut = []
for i,evt in enumerate(xpiarray):
    xpicut.append("xpicut%s" % i)
    tmp  = "xpicut_%s = [\"xpicut%s\"]" % (i,i)
    exec(tmp)
    
s_e = c.applyCuts(s_e,xpicut,either=True)
s_q = c.applyCuts(s_q,xpicut,either=True)
Q2 = c.applyCuts(Q2,xpicut,either=True)
xBj = c.applyCuts(xBj,xpicut,either=True)
t = c.applyCuts(t,xpicut,either=True)
tPrime = c.applyCuts(tPrime,xpicut,either=True)
y_D = c.applyCuts(y_D,xpicut,either=True)
nu = c.applyCuts(nu,xpicut,either=True)
TwoPdotk = c.applyCuts(TwoPdotk,xpicut,either=True)
TDIS_xbj = c.applyCuts(TDIS_xbj,xpicut,either=True)
sigma_dis = c.applyCuts(sigma_dis,xpicut,either=True)
TDIS_y = c.applyCuts(TDIS_y,xpicut,either=True)
ppix_Lab = c.applyCuts(ppix_Lab,xpicut,either=True)
ppiy_Lab = c.applyCuts(ppiy_Lab,xpicut,either=True)
ppiz_Lab = c.applyCuts(ppiz_Lab,xpicut,either=True)
EpiE_Lab = c.applyCuts(EpiE_Lab,xpicut,either=True)
pprx_inc = c.applyCuts(pprx_inc,xpicut,either=True)
ppry_inc = c.applyCuts(ppry_inc,xpicut,either=True)
pprz_inc = c.applyCuts(pprz_inc,xpicut,either=True)
EprE_inc = c.applyCuts(EprE_inc,xpicut,either=True)
EnE_Lab = c.applyCuts(EnE_Lab,xpicut,either=True)
y = c.applyCuts(y,xpicut,either=True)
fpi = c.applyCuts(fpi,xpicut,either=True)
f2N = c.applyCuts(f2N,xpicut,either=True)
ypi = c.applyCuts(ypi,xpicut,either=True)
tpi = c.applyCuts(tpi,xpicut,either=True)
escat = c.applyCuts(escat,xpicut,either=True)
# pprz_inc = c.applyCuts(pprz_inc,xpicut,either=True)
# ppiz_Lab = c.applyCuts(ppiz_Lab,xpicut,either=True)

xpi = c.applyCuts(xpi,xpicut,either=True)

Q2array = np.arange(binQ2/2,1000.0,binQ2).tolist()
for i,x in enumerate(Q2array) :
    Q2tmp = '{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,Q2array[i]-binQ2/2,Q2array[i]+binQ2/2)
    # (i,Q2array[i]-binQ2/20,Q2array[i]+binQ2/20) # for proper binning
    print('{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,Q2array[i]-binQ2/20,Q2array[i]+binQ2/20))
    cutDict.update(eval(Q2tmp))
c = r2p.pyPlot(cutDict)

Q2cut = []
for i,evt in enumerate(Q2array):
    Q2cut.append("Q2cut%s" % i)

s_e = c.applyCuts(s_e,Q2cut,either=True)
s_q = c.applyCuts(s_q,Q2cut,either=True)
xBj = c.applyCuts(xBj,Q2cut,either=True)
t = c.applyCuts(t,Q2cut,either=True)
tPrime = c.applyCuts(tPrime,Q2cut,either=True)
y_D = c.applyCuts(y_D,Q2cut,either=True)
nu = c.applyCuts(nu,Q2cut,either=True)
TwoPdotk = c.applyCuts(TwoPdotk,Q2cut,either=True)
TDIS_xbj = c.applyCuts(TDIS_xbj,Q2cut,either=True)
sigma_dis = c.applyCuts(sigma_dis,Q2cut,either=True)
TDIS_y = c.applyCuts(TDIS_y,Q2cut,either=True)
ppix_Lab = c.applyCuts(ppix_Lab,Q2cut,either=True)
ppiy_Lab = c.applyCuts(ppiy_Lab,Q2cut,either=True)
ppiz_Lab = c.applyCuts(ppiz_Lab,Q2cut,either=True)
EpiE_Lab = c.applyCuts(EpiE_Lab,Q2cut,either=True)
pprx_inc = c.applyCuts(pprx_inc,Q2cut,either=True)
ppry_inc = c.applyCuts(ppry_inc,Q2cut,either=True)
pprz_inc = c.applyCuts(pprz_inc,Q2cut,either=True)
EprE_inc = c.applyCuts(EprE_inc,Q2cut,either=True)
EnE_Lab = c.applyCuts(EnE_Lab,Q2cut,either=True)
y = c.applyCuts(y,Q2cut,either=True)
fpi = c.applyCuts(fpi,Q2cut,either=True)
f2N = c.applyCuts(f2N,Q2cut,either=True)
xpi = c.applyCuts(xpi,Q2cut,either=True)
ypi = c.applyCuts(ypi,Q2cut,either=True)
tpi = c.applyCuts(tpi,Q2cut,either=True)
escat = c.applyCuts(escat,Q2cut,either=True)
# pprz_inc = c.applyCuts(pprz_inc,Q2cut,either=True)
# ppiz_Lab = c.applyCuts(ppiz_Lab,Q2cut,either=True)

Q2 = c.applyCuts(Q2,Q2cut,either=True)

sigma_tdis = sigma_dis*(fpi/f2N)
pNz_Lab = pprz_inc-ppiz_Lab # Long neutron momentum
xL = EnE_Lab/EprE_inc # Frac of proton momentum
# xL = 1-(xBj/xpi) # Frac of proton momentum
Q2_new = xBj/TwoPdotk
piQ2 = s_q/(xpi*ypi)
pNz_Lab = pprz_inc-ppiz_Lab # Long neutron momentum
ppiz_frac = ppiz_Lab/pprz_inc # Frac of proton momentum
pNz_frac = pNz_Lab/pprz_inc # Frac of proton momentum
yplus = 1+((1-y)*(1-y))
tot_sigma = (sigma_tdis)*((TDIS_xbj*(Q2*Q2)*(137)*(137))/(2*math.pi*yplus))
# red_sig  = np.array([1.397,1.258,1.137,1.055,0.944,0.836,0.706,0.519])
# red_y = np.array([0.657,0.416,0.263,0.163,0.103,0.066,0.033,0.008])
# red_x = np.array([2.53e-4,4.0e-4,6.32e-4,1.02e-3,1.61e-3,2.52e-3,5.0e-3,2.10e-2])
# red_Q2 = np.array([15,15,15,15,15,15,15,15])
# red_yplus = 1+(1-red_y)*(1-red_y)
# my_sigma = (red_sig)*((2*math.pi*red_yplus)/(red_x*(red_Q2*red_Q2)*(137)*(137)))

Q2binarray = [7,15,30,60,120,240,480,1000]
for i,x in enumerate(Q2binarray) :
    Q2tmp = '{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,Q2binarray[i]-binQ2/20,Q2binarray[i]+binQ2/20)
    print('{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,Q2binarray[i]-binQ2/20,Q2binarray[i]+binQ2/20))
    cutDict.update(eval(Q2tmp))

# xarray = np.arange(0.05,1.0,0.1).tolist() # table, large x (usual)
xarray = np.arange(0.005,1.0,0.001).tolist() # table, small x 
for i,x in enumerate(xarray):
    xtmp = '{"xcut%i" : ((%0.4f <= xpi) & (xpi <= %0.4f))}' % (i,xarray[i]-0.005,xarray[i]+0.005)
    print('{"xcut%i" : ((%0.4f <= xpi) & (xpi <= %0.4f))}' % (i,xarray[i]-0.005,xarray[i]+0.005))
    cutDict.update(eval(xtmp))

xcut = []
for i,evt in enumerate(xarray):
    xcut.append("xcut%s" % i)
    tmp  = "xcut_%s = [\"xcut%s\"]" % (i,i)
    exec(tmp)

    
ytmp = '{"ycut" : ((0.01 <= y) & (y <= 0.95))}'
ttmp = '{"tcut" : ((-1.00 <= t) & (t <= 0.00))}'
cutDict.update(eval(ttmp))
cutDict.update(eval(ytmp))
c = r2p.pyPlot(cutDict)

ycut1 = ["ycut"]
tcut1 = ["tcut"]
cut7 = ["Q2cut0","tcut","ycut"]
cut15 = ["Q2cut1","tcut","ycut"]
cut30 = ["Q2cut2","tcut","ycut"]
cut60 = ["Q2cut3","tcut","ycut"]
cut120 = ["Q2cut4","tcut","ycut"]
cut240 = ["Q2cut5","tcut","ycut"]
cut480 = ["Q2cut6","tcut","ycut"]
cut1000 = ["Q2cut7","tcut","ycut"]

cutx0_15 = ["Q2cut1","tcut","ycut","xcut0"]
cutx1_15 = ["Q2cut1","tcut","ycut","xcut1"]
cutx2_15 = ["Q2cut1","tcut","ycut","xcut2"]
cutx3_15 = ["Q2cut1","tcut","ycut","xcut3"]
cutx4_15 = ["Q2cut1","tcut","ycut","xcut4"]
cutx5_15 = ["Q2cut1","tcut","ycut","xcut5"]
cutx6_15 = ["Q2cut1","tcut","ycut","xcut6"]
cutx7_15 = ["Q2cut1","tcut","ycut","xcut7"]
cutx8_15 = ["Q2cut1","tcut","ycut","xcut8"]
cutx9_15 = ["Q2cut1","tcut","ycut","xcut9"]

cutx0_30 = ["Q2cut2","tcut","ycut","xcut0"]
cutx1_30 = ["Q2cut2","tcut","ycut","xcut1"]
cutx2_30 = ["Q2cut2","tcut","ycut","xcut2"]
cutx3_30 = ["Q2cut2","tcut","ycut","xcut3"]
cutx4_30 = ["Q2cut2","tcut","ycut","xcut4"]
cutx5_30 = ["Q2cut2","tcut","ycut","xcut5"]
cutx6_30 = ["Q2cut2","tcut","ycut","xcut6"]
cutx7_30 = ["Q2cut2","tcut","ycut","xcut7"]
cutx8_30 = ["Q2cut2","tcut","ycut","xcut8"]
cutx9_30 = ["Q2cut2","tcut","ycut","xcut9"]

cutx0_60 = ["Q2cut3","tcut","ycut","xcut0"]
cutx1_60 = ["Q2cut3","tcut","ycut","xcut1"]
cutx2_60 = ["Q2cut3","tcut","ycut","xcut2"]
cutx3_60 = ["Q2cut3","tcut","ycut","xcut3"]
cutx4_60 = ["Q2cut3","tcut","ycut","xcut4"]
cutx5_60 = ["Q2cut3","tcut","ycut","xcut5"]
cutx6_60 = ["Q2cut3","tcut","ycut","xcut6"]
cutx7_60 = ["Q2cut3","tcut","ycut","xcut7"]
cutx8_60 = ["Q2cut3","tcut","ycut","xcut8"]
cutx9_60 = ["Q2cut3","tcut","ycut","xcut9"]

cutx0_120 = ["Q2cut4","tcut","ycut","xcut0"]
cutx1_120 = ["Q2cut4","tcut","ycut","xcut1"]
cutx2_120 = ["Q2cut4","tcut","ycut","xcut2"]
cutx3_120 = ["Q2cut4","tcut","ycut","xcut3"]
cutx4_120 = ["Q2cut4","tcut","ycut","xcut4"]
cutx5_120 = ["Q2cut4","tcut","ycut","xcut5"]
cutx6_120 = ["Q2cut4","tcut","ycut","xcut6"]
cutx7_120 = ["Q2cut4","tcut","ycut","xcut7"]
cutx8_120 = ["Q2cut4","tcut","ycut","xcut8"]
cutx9_120 = ["Q2cut4","tcut","ycut","xcut9"]


## For table, small x

cutxpi001_15 = ["Q2cut1","tcut","ycut","xcut1"]
cutxpi01_15 = ["Q2cut1","tcut","ycut","xcut10"]

cutxpi01_30 = ["Q2cut2","tcut","ycut","xcut10"]

cutxpi01_60 = ["Q2cut3","tcut","ycut","xcut10"]

## For table, large x

cutxpi1_15 = ["Q2cut1","tcut","ycut","xcut0"]
cutxpi2_15 = ["Q2cut1","tcut","ycut","xcut1"]
cutxpi3_15 = ["Q2cut1","tcut","ycut","xcut2"]
cutxpi4_15 = ["Q2cut1","tcut","ycut","xcut3"]
cutxpi5_15 = ["Q2cut1","tcut","ycut","xcut4"]
cutxpi6_15 = ["Q2cut1","tcut","ycut","xcut5"]
cutxpi7_15 = ["Q2cut1","tcut","ycut","xcut6"]

cutxpi1_30 = ["Q2cut2","tcut","ycut","xcut0"]
cutxpi2_30 = ["Q2cut2","tcut","ycut","xcut1"]
cutxpi3_30 = ["Q2cut2","tcut","ycut","xcut2"]
cutxpi4_30 = ["Q2cut2","tcut","ycut","xcut3"]
cutxpi5_30 = ["Q2cut2","tcut","ycut","xcut4"]
cutxpi6_30 = ["Q2cut2","tcut","ycut","xcut5"]
cutxpi7_30 = ["Q2cut2","tcut","ycut","xcut6"]

cutxpi1_60 = ["Q2cut3","tcut","ycut","xcut0"]
cutxpi2_60 = ["Q2cut3","tcut","ycut","xcut1"]
cutxpi3_60 = ["Q2cut3","tcut","ycut","xcut2"]
cutxpi4_60 = ["Q2cut3","tcut","ycut","xcut3"]
cutxpi5_60 = ["Q2cut3","tcut","ycut","xcut4"]
cutxpi6_60 = ["Q2cut3","tcut","ycut","xcut5"]
cutxpi7_60 = ["Q2cut3","tcut","ycut","xcut6"]

cutxpi1_120 = ["Q2cut4","tcut","ycut","xcut0"]
cutxpi2_120 = ["Q2cut4","tcut","ycut","xcut1"]
cutxpi3_120 = ["Q2cut4","tcut","ycut","xcut2"]
cutxpi4_120 = ["Q2cut4","tcut","ycut","xcut3"]
cutxpi5_120 = ["Q2cut4","tcut","ycut","xcut4"]
cutxpi6_120 = ["Q2cut4","tcut","ycut","xcut5"]
cutxpi7_120 = ["Q2cut4","tcut","ycut","xcut6"]

cutxpi1_240 = ["Q2cut5","tcut","ycut","xcut0"]
cutxpi2_240 = ["Q2cut5","tcut","ycut","xcut1"]
cutxpi3_240 = ["Q2cut5","tcut","ycut","xcut2"]
cutxpi4_240 = ["Q2cut5","tcut","ycut","xcut3"]
cutxpi5_240 = ["Q2cut5","tcut","ycut","xcut4"]
cutxpi6_240 = ["Q2cut5","tcut","ycut","xcut5"]
cutxpi7_240 = ["Q2cut5","tcut","ycut","xcut6"]

cutxpi1_480 = ["Q2cut6","tcut","ycut","xcut0"]
cutxpi2_480 = ["Q2cut6","tcut","ycut","xcut1"]
cutxpi3_480 = ["Q2cut6","tcut","ycut","xcut2"]
cutxpi4_480 = ["Q2cut6","tcut","ycut","xcut3"]
cutxpi5_480 = ["Q2cut6","tcut","ycut","xcut4"]
cutxpi6_480 = ["Q2cut6","tcut","ycut","xcut5"]
cutxpi7_480 = ["Q2cut6","tcut","ycut","xcut6"]


def phase_space():

    phaseSpace = c.densityPlot(xpi, Q2, '$Q^2$ vs $x_{\pi}$','$x_{\pi}$','$Q^{2}$', 200, 200,  c, 0., 1.0, 0., 1000.)
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(0.,1.)
    plt.ylim(0.,1000.)

    phaseSpace[1].savefig('OUTPUTS/phase_space.png')

def sigmavxpi_Plot():
    
    f = plt.figure(figsize=(11.69,8.27))
    plt.style.use('classic')
    
    ax = f.add_subplot(421)
    xpiscat1 = ax.errorbar(c.applyCuts(xpi,cut7),c.applyCuts(tot_sigma,cut7),yerr=np.sqrt(c.applyCuts(tot_sigma,cut7))/200,fmt='o',label='$Q^2$=7 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=10 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    plt.ylabel('$\sigma_{TDIS}$ ($nb/GeV^{2}$)')

    ax = f.add_subplot(422)
    xpiscat2 = ax.errorbar(c.applyCuts(xpi,cut15),c.applyCuts(tot_sigma,cut15),yerr=np.sqrt(c.applyCuts(tot_sigma,cut15))/150,fmt='o',label='$Q^2$=15 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=20 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    plt.title('$\sigma_{TDIS}$ vs $x_\pi$', fontsize =20)
    
    ax = f.add_subplot(423)
    xpiscat3 = ax.errorbar(c.applyCuts(xpi,cut30),c.applyCuts(tot_sigma,cut30),yerr=np.sqrt(c.applyCuts(tot_sigma,cut30))/200,fmt='o',label='$Q^2$=30 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    
    ax = f.add_subplot(424)
    xpiscat4 = ax.errorbar(c.applyCuts(xpi,cut60),c.applyCuts(tot_sigma,cut60),yerr=np.sqrt(c.applyCuts(tot_sigma,cut60))/200,fmt='o',label='$Q^2$=60 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=40 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    ax = f.add_subplot(425)
    xpiscat5 = ax.errorbar(c.applyCuts(xpi,cut120),c.applyCuts(tot_sigma,cut120),yerr=np.sqrt(c.applyCuts(tot_sigma,cut120))/200,fmt='o',label='$Q^2$=120 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=50 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    ax = f.add_subplot(426)
    xpiscat6 = ax.errorbar(c.applyCuts(xpi,cut240),c.applyCuts(tot_sigma,cut240),yerr=np.sqrt(c.applyCuts(tot_sigma,cut240))/200,fmt='o',label='$Q^2$=240 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=70 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(427)
    xpiscat7 = ax.errorbar(c.applyCuts(xpi,cut480),c.applyCuts(tot_sigma,cut480),yerr=np.sqrt(c.applyCuts(tot_sigma,cut480))/200,fmt='o',label='$Q^2$=480 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=80 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))    
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    ax = f.add_subplot(428)
    xpiscat8 = ax.errorbar(c.applyCuts(xpi,cut1000),c.applyCuts(tot_sigma,cut1000),yerr=np.sqrt(c.applyCuts(tot_sigma,cut1000))/200,fmt='o',label='$Q^2$=1000 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=90 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))    
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    plt.xlabel('$x_\pi$')
    plt.tight_layout()

    plt.style.use('default')


def fpivxpi_Plot():
    
    f = plt.figure(figsize=(11.69,8.27))
    plt.rcParams.update({'font.size': 15})
    plt.style.use('classic')
    
    ax = f.add_subplot(221)
    xpiscat4 = ax.errorbar(c.applyCuts(xpi,cut60),c.applyCuts(fpi,cut60),yerr=np.sqrt(c.applyCuts(lumi,cut60))/c.applyCuts(lumi,cut60),fmt='.',label='$Q^2$=60 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    plt.plot([0.001,0.01,0.1],[1.2,0.45,0.25], label="GRV fit",color="y")
    plt.xscale('log')
    plt.ylim(-0.1,0.3)
    plt.xlim(1e-3,1.)
    ax.text(0.25, 0.65, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.set_yticks([-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3])
    ax.set_xticks([1e-2,1e-1])
    
    plt.ylabel('$F^{\pi}_{2}$', fontsize=20)
    
    ax = f.add_subplot(222)
    xpiscat5 = ax.errorbar(c.applyCuts(xpi,cut120),c.applyCuts(fpi,cut120),yerr=np.sqrt(c.applyCuts(lumi,cut120))/c.applyCuts(lumi,cut120),fmt='.',label='$Q^2$=120 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.5,0.25,0.15], label="GRV fit",color="y")
    plt.xscale('log')
    plt.ylim(-0.1,0.3)
    plt.xlim(1e-3,1.)
    ax.text(0.25, 0.65, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.set_yticks([-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3])
    ax.set_xticks([1e-2,1e-1])
    
    ax = f.add_subplot(223)
    xpiscat6 = ax.errorbar(c.applyCuts(xpi,cut240),c.applyCuts(fpi,cut240),yerr=np.sqrt(c.applyCuts(lumi,cut240))/c.applyCuts(lumi,cut240),fmt='.',label='$Q^2$=240 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit",color="y")
    plt.xscale('log')
    plt.ylim(-0.1,0.3)
    plt.xlim(1e-3,1.)
    ax.text(0.25, 0.65, '$Q^2$=240 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.set_yticks([-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3])
    ax.set_xticks([1e-2,1e-1])

    ax = f.add_subplot(224)
    xpiscat7 = ax.errorbar(c.applyCuts(xpi,cut480),c.applyCuts(fpi,cut480),yerr=np.sqrt(c.applyCuts(lumi,cut480))/c.applyCuts(lumi,cut480),fmt='.',label='$Q^2$=480 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit",color="y")
    plt.xscale('log')
    plt.ylim(-0.1,0.3)
    plt.xlim(1e-3,1.)
    ax.text(0.25, 0.65, '$Q^2$=480 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.set_yticks([-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3])
    ax.set_xticks([1e-2,1e-1])

    plt.xlabel('$x_\pi$', fontsize=20)    
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.0,wspace=0.0)

    plt.style.use('default')

    f.savefig('OUTPUTS/fpivxpi.png')

def fpivxpi_Plot_nolog():
    
    f = plt.figure(figsize=(11.69,8.27))
    
    ax = f.add_subplot(221)
    xpiscat4 = ax.errorbar(c.applyCuts(xpi,cut60),c.applyCuts(fpi,cut60),yerr=np.sqrt(c.applyCuts(lumi,cut60))/c.applyCuts(lumi,cut60),fmt='.',label='$Q^2$=60 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    plt.plot([0.001,0.01,0.1],[1.2,0.45,0.25], label="GRV fit",color="y")
    plt.xlim(0.,1.)
    plt.ylim(-0.1,0.3)
    ax.text(0.60, 0.85, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    plt.ylabel('$F^{\pi}_{2}$', fontsize=20)

    ax = f.add_subplot(222)
    xpiscat5 = ax.errorbar(c.applyCuts(xpi,cut120),c.applyCuts(fpi,cut120),yerr=np.sqrt(c.applyCuts(lumi,cut120))/c.applyCuts(lumi,cut120),fmt='.',label='$Q^2$=120 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.5,0.25,0.15], label="GRV fit",color="y")
    plt.xlim(0.,1.)
    plt.ylim(-0.1,0.3)
    ax.text(0.60, 0.85, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    plt.title('$F^{\pi}_{2}$ vs $x_\pi$', fontsize =20)
    
    ax = f.add_subplot(223)
    xpiscat6 = ax.errorbar(c.applyCuts(xpi,cut240),c.applyCuts(fpi,cut240),yerr=np.sqrt(c.applyCuts(lumi,cut240))/c.applyCuts(lumi,cut240),fmt='.',label='$Q^2$=240 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit",color="y")
    plt.xlim(0.,1.)
    plt.ylim(-0.1,0.3)
    ax.text(0.60, 0.85, '$Q^2$=240 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    # ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))

    ax = f.add_subplot(224)
    xpiscat7 = ax.errorbar(c.applyCuts(xpi,cut480),c.applyCuts(fpi,cut480),yerr=np.sqrt(c.applyCuts(lumi,cut480))/c.applyCuts(lumi,cut480),fmt='.',label='$Q^2$=480 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit",color="y")
    plt.xlim(1e-4,1.)
    plt.ylim(-0.1,0.3)
    ax.text(0.60, 0.85, '$Q^2$=480 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    
    plt.xlabel('$x_\pi$', fontsize=20)
    plt.tight_layout()

    plt.subplots_adjust(hspace=0.0,wspace=0.0)

    f.savefig('OUTPUTS/fpivxpi_nolog.png')

    
def fpivt_Plot():

    f = plt.figure(figsize=(11.69,8.27))
    plt.style.use('classic')
    
    ax = f.add_subplot(221)
    tscat4 = ax.errorbar(c.applyCuts(-t,cut60),c.applyCuts(fpi,cut60),yerr=np.sqrt(c.applyCuts(lumi,cut60))/c.applyCuts(lumi,cut60),fmt='.',label='$Q^2$=60 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    # plt.xscale('log')
    # plt.yscale('log')
    plt.ylim(-0.1,0.3)
    plt.xlim(0.0,1.0)
    ax.text(0.65, 0.85, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.set_yticks([-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3])
    ax.set_xticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])

    plt.ylabel('$F^{\pi}_{2}$', fontsize=20)
    
    ax = f.add_subplot(222)
    tscat5 = ax.errorbar(c.applyCuts(-t,cut120),c.applyCuts(fpi,cut120),yerr=np.sqrt(c.applyCuts(lumi,cut120))/c.applyCuts(lumi,cut120),fmt='.',label='$Q^2$=120 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    # plt.xscale('log')
    # plt.yscale('log')
    plt.ylim(-0.1,0.3)
    plt.xlim(0.0,1.0)
    ax.text(0.65, 0.85, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.set_yticks([-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3])
    ax.set_xticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])

    ax = f.add_subplot(223)
    tscat6 = ax.errorbar(c.applyCuts(-t,cut240),c.applyCuts(fpi,cut240),yerr=np.sqrt(c.applyCuts(lumi,cut240))/c.applyCuts(lumi,cut240),fmt='.',label='$Q^2$=240 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    # plt.xscale('log')
    # plt.yscale('log')
    plt.ylim(-0.1,0.3)
    plt.xlim(0.0,1.0)
    ax.text(0.65, 0.85, '$Q^2$=240 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.set_yticks([-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3])
    ax.set_xticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])

    ax = f.add_subplot(224)
    tscat7 = ax.errorbar(c.applyCuts(-t,cut480),c.applyCuts(fpi,cut480),yerr=np.sqrt(c.applyCuts(lumi,cut480))/c.applyCuts(lumi,cut480),fmt='.',label='$Q^2$=480 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    # plt.xscale('log')
    # plt.yscale('log')
    plt.ylim(-0.1,0.3)
    plt.xlim(0.0,1.0)
    ax.text(0.65, 0.85, '$Q^2$=480 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.set_yticks([-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3])
    ax.set_xticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])
    
    plt.xlabel('-t', fontsize=20)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.0,wspace=0.0)

    f.savefig('OUTPUTS/fpivt.png')

def fpivt_xbin_Plot():

    f = plt.figure(figsize=(11.69,8.27))
    
    ax = f.add_subplot(221)
    tscat0 = ax.errorbar(c.applyCuts(t,cutx0_15),c.applyCuts(fpi,cutx0_15),yerr=np.sqrt(c.applyCuts(lumi,cutx0_15))/c.applyCuts(lumi,cutx0_15),fmt='.',label='$x_{\pi}$=(0.0,0.1)',ecolor='cyan',capsize=2, capthick=2,color='b')
    tscat1 = ax.errorbar(c.applyCuts(t,cutx1_15),c.applyCuts(fpi,cutx1_15),yerr=np.sqrt(c.applyCuts(lumi,cutx1_15))/c.applyCuts(lumi,cutx1_15),fmt='.',label='$x_{\pi}$=(0.1,0.2)',ecolor='cyan',capsize=2, capthick=2,color='g')
    tscat2 = ax.errorbar(c.applyCuts(t,cutx2_15),c.applyCuts(fpi,cutx2_15),yerr=np.sqrt(c.applyCuts(lumi,cutx2_15))/c.applyCuts(lumi,cutx2_15),fmt='.',label='$x_{\pi}$=(0.2,0.3)',ecolor='cyan',capsize=2, capthick=2,color='r')
    tscat3 = ax.errorbar(c.applyCuts(t,cutx3_15),c.applyCuts(fpi,cutx3_15),yerr=np.sqrt(c.applyCuts(lumi,cutx3_15))/c.applyCuts(lumi,cutx3_15),fmt='.',label='$x_{\pi}$=(0.3,0.4)',ecolor='cyan',capsize=2, capthick=2,color='c')
    tscat4 = ax.errorbar(c.applyCuts(t,cutx4_15),c.applyCuts(fpi,cutx4_15),yerr=np.sqrt(c.applyCuts(lumi,cutx4_15))/c.applyCuts(lumi,cutx4_15),fmt='.',label='$x_{\pi}$=(0.4,0.5)',ecolor='cyan',capsize=2, capthick=2,color='m')
    tscat5 = ax.errorbar(c.applyCuts(t,cutx5_15),c.applyCuts(fpi,cutx5_15),yerr=np.sqrt(c.applyCuts(lumi,cutx5_15))/c.applyCuts(lumi,cutx5_15),fmt='.',label='$x_{\pi}$=(0.5,0.6)',ecolor='cyan',capsize=2, capthick=2,color='y')
    tscat6 = ax.errorbar(c.applyCuts(t,cutx6_15),c.applyCuts(fpi,cutx6_15),yerr=np.sqrt(c.applyCuts(lumi,cutx6_15))/c.applyCuts(lumi,cutx6_15),fmt='.',label='$x_{\pi}$=(0.6,0.7)',ecolor='cyan',capsize=2, capthick=2,color='purple')
    tscat7 = ax.errorbar(c.applyCuts(t,cutx7_15),c.applyCuts(fpi,cutx7_15),yerr=np.sqrt(c.applyCuts(lumi,cutx7_15))/c.applyCuts(lumi,cutx7_15),fmt='.',label='$x_{\pi}$=(0.7,0.8)',ecolor='cyan',capsize=2, capthick=2,color='aqua')
    tscat8 = ax.errorbar(c.applyCuts(t,cutx8_15),c.applyCuts(fpi,cutx8_15),yerr=np.sqrt(c.applyCuts(lumi,cutx8_15))/c.applyCuts(lumi,cutx8_15),fmt='.',label='$x_{\pi}$=(0.8,0.9)',ecolor='cyan',capsize=2, capthick=2,color='gray')
    tscat9 = ax.errorbar(c.applyCuts(t,cutx9_15),c.applyCuts(fpi,cutx9_15),yerr=np.sqrt(c.applyCuts(lumi,cutx9_15))/c.applyCuts(lumi,cutx9_15),fmt='.',label='$x_{\pi}$=(0.9,1.0)',ecolor='cyan',capsize=2, capthick=2,color='pink')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(-1.,0.)
    plt.ylim(-0.1,0.3)
    ax.text(0.85, 0.85, '$Q^2$=15 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='lower'))

    plt.ylabel('$F^{\pi}_{2}$', fontsize=20)
    
    ax = f.add_subplot(222)
    tscat0 = ax.errorbar(c.applyCuts(t,cutx0_30),c.applyCuts(fpi,cutx0_30),yerr=np.sqrt(c.applyCuts(lumi,cutx0_30))/c.applyCuts(lumi,cutx0_30),fmt='.',label='$x_{\pi}$=(0.0,0.1)',ecolor='cyan',capsize=2, capthick=2,color='b')
    tscat1 = ax.errorbar(c.applyCuts(t,cutx1_30),c.applyCuts(fpi,cutx1_30),yerr=np.sqrt(c.applyCuts(lumi,cutx1_30))/c.applyCuts(lumi,cutx1_30),fmt='.',label='$x_{\pi}$=(0.1,0.2)',ecolor='cyan',capsize=2, capthick=2,color='g')
    tscat2 = ax.errorbar(c.applyCuts(t,cutx2_30),c.applyCuts(fpi,cutx2_30),yerr=np.sqrt(c.applyCuts(lumi,cutx2_30))/c.applyCuts(lumi,cutx2_30),fmt='.',label='$x_{\pi}$=(0.2,0.3)',ecolor='cyan',capsize=2, capthick=2,color='r')
    tscat3 = ax.errorbar(c.applyCuts(t,cutx3_30),c.applyCuts(fpi,cutx3_30),yerr=np.sqrt(c.applyCuts(lumi,cutx3_30))/c.applyCuts(lumi,cutx3_30),fmt='.',label='$x_{\pi}$=(0.3,0.4)',ecolor='cyan',capsize=2, capthick=2,color='c')
    tscat4 = ax.errorbar(c.applyCuts(t,cutx4_30),c.applyCuts(fpi,cutx4_30),yerr=np.sqrt(c.applyCuts(lumi,cutx4_30))/c.applyCuts(lumi,cutx4_30),fmt='.',label='$x_{\pi}$=(0.4,0.5)',ecolor='cyan',capsize=2, capthick=2,color='m')
    tscat5 = ax.errorbar(c.applyCuts(t,cutx5_30),c.applyCuts(fpi,cutx5_30),yerr=np.sqrt(c.applyCuts(lumi,cutx5_30))/c.applyCuts(lumi,cutx5_30),fmt='.',label='$x_{\pi}$=(0.5,0.6)',ecolor='cyan',capsize=2, capthick=2,color='y')
    tscat6 = ax.errorbar(c.applyCuts(t,cutx6_30),c.applyCuts(fpi,cutx6_30),yerr=np.sqrt(c.applyCuts(lumi,cutx6_30))/c.applyCuts(lumi,cutx6_30),fmt='.',label='$x_{\pi}$=(0.6,0.7)',ecolor='cyan',capsize=2, capthick=2,color='purple')
    tscat7 = ax.errorbar(c.applyCuts(t,cutx7_30),c.applyCuts(fpi,cutx7_30),yerr=np.sqrt(c.applyCuts(lumi,cutx7_30))/c.applyCuts(lumi,cutx7_30),fmt='.',label='$x_{\pi}$=(0.7,0.8)',ecolor='cyan',capsize=2, capthick=2,color='aqua')
    tscat8 = ax.errorbar(c.applyCuts(t,cutx8_30),c.applyCuts(fpi,cutx8_30),yerr=np.sqrt(c.applyCuts(lumi,cutx8_30))/c.applyCuts(lumi,cutx8_30),fmt='.',label='$x_{\pi}$=(0.8,0.9)',ecolor='cyan',capsize=2, capthick=2,color='gray')
    tscat9 = ax.errorbar(c.applyCuts(t,cutx9_30),c.applyCuts(fpi,cutx9_30),yerr=np.sqrt(c.applyCuts(lumi,cutx9_30))/c.applyCuts(lumi,cutx9_30),fmt='.',label='$x_{\pi}$=(0.9,1.0)',ecolor='cyan',capsize=2, capthick=2,color='pink')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(-1.,0.)
    plt.ylim(-0.1,0.3)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.text(0.85, 0.85, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    plt.title('$F^{\pi}_{2}$ vs -t', fontsize =20)
    
    ax = f.add_subplot(223)
    tscat0 = ax.errorbar(c.applyCuts(t,cutx0_60),c.applyCuts(fpi,cutx0_60),yerr=np.sqrt(c.applyCuts(lumi,cutx0_60))/c.applyCuts(lumi,cutx0_60),fmt='.',label='$x_{\pi}$=(0.0,0.1)',ecolor='cyan',capsize=2, capthick=2,color='b')
    tscat1 = ax.errorbar(c.applyCuts(t,cutx1_60),c.applyCuts(fpi,cutx1_60),yerr=np.sqrt(c.applyCuts(lumi,cutx1_60))/c.applyCuts(lumi,cutx1_60),fmt='.',label='$x_{\pi}$=(0.1,0.2)',ecolor='cyan',capsize=2, capthick=2,color='g')
    tscat2 = ax.errorbar(c.applyCuts(t,cutx2_60),c.applyCuts(fpi,cutx2_60),yerr=np.sqrt(c.applyCuts(lumi,cutx2_60))/c.applyCuts(lumi,cutx2_60),fmt='.',label='$x_{\pi}$=(0.2,0.3)',ecolor='cyan',capsize=2, capthick=2,color='r')
    tscat3 = ax.errorbar(c.applyCuts(t,cutx3_60),c.applyCuts(fpi,cutx3_60),yerr=np.sqrt(c.applyCuts(lumi,cutx3_60))/c.applyCuts(lumi,cutx3_60),fmt='.',label='$x_{\pi}$=(0.3,0.4)',ecolor='cyan',capsize=2, capthick=2,color='c')
    tscat4 = ax.errorbar(c.applyCuts(t,cutx4_60),c.applyCuts(fpi,cutx4_60),yerr=np.sqrt(c.applyCuts(lumi,cutx4_60))/c.applyCuts(lumi,cutx4_60),fmt='.',label='$x_{\pi}$=(0.4,0.5)',ecolor='cyan',capsize=2, capthick=2,color='m')
    tscat5 = ax.errorbar(c.applyCuts(t,cutx5_60),c.applyCuts(fpi,cutx5_60),yerr=np.sqrt(c.applyCuts(lumi,cutx5_60))/c.applyCuts(lumi,cutx5_60),fmt='.',label='$x_{\pi}$=(0.5,0.6)',ecolor='cyan',capsize=2, capthick=2,color='y')
    tscat6 = ax.errorbar(c.applyCuts(t,cutx6_60),c.applyCuts(fpi,cutx6_60),yerr=np.sqrt(c.applyCuts(lumi,cutx6_60))/c.applyCuts(lumi,cutx6_60),fmt='.',label='$x_{\pi}$=(0.6,0.7)',ecolor='cyan',capsize=2, capthick=2,color='purple')
    tscat7 = ax.errorbar(c.applyCuts(t,cutx7_60),c.applyCuts(fpi,cutx7_60),yerr=np.sqrt(c.applyCuts(lumi,cutx7_60))/c.applyCuts(lumi,cutx7_60),fmt='.',label='$x_{\pi}$=(0.7,0.8)',ecolor='cyan',capsize=2, capthick=2,color='aqua')
    tscat8 = ax.errorbar(c.applyCuts(t,cutx8_60),c.applyCuts(fpi,cutx8_60),yerr=np.sqrt(c.applyCuts(lumi,cutx8_60))/c.applyCuts(lumi,cutx8_60),fmt='.',label='$x_{\pi}$=(0.8,0.9)',ecolor='cyan',capsize=2, capthick=2,color='gray')
    tscat9 = ax.errorbar(c.applyCuts(t,cutx9_60),c.applyCuts(fpi,cutx9_60),yerr=np.sqrt(c.applyCuts(lumi,cutx9_60))/c.applyCuts(lumi,cutx9_60),fmt='.',label='$x_{\pi}$=(0.9,1.0)',ecolor='cyan',capsize=2, capthick=2,color='pink')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(-1.,0.)
    plt.ylim(-0.1,0.3)
    ax.text(0.85, 0.85, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    
    ax = f.add_subplot(224)
    tscat0 = ax.errorbar(c.applyCuts(t,cutx0_120),c.applyCuts(fpi,cutx0_120),yerr=np.sqrt(c.applyCuts(lumi,cutx0_120))/c.applyCuts(lumi,cutx0_120),fmt='.',label='$x_{\pi}$=(0.0,0.1)',ecolor='cyan',capsize=2, capthick=2,color='b')
    tscat1 = ax.errorbar(c.applyCuts(t,cutx1_120),c.applyCuts(fpi,cutx1_120),yerr=np.sqrt(c.applyCuts(lumi,cutx1_120))/c.applyCuts(lumi,cutx1_120),fmt='.',label='$x_{\pi}$=(0.1,0.2)',ecolor='cyan',capsize=2, capthick=2,color='g')
    tscat2 = ax.errorbar(c.applyCuts(t,cutx2_120),c.applyCuts(fpi,cutx2_120),yerr=np.sqrt(c.applyCuts(lumi,cutx2_120))/c.applyCuts(lumi,cutx2_120),fmt='.',label='$x_{\pi}$=(0.2,0.3)',ecolor='cyan',capsize=2, capthick=2,color='r')
    tscat3 = ax.errorbar(c.applyCuts(t,cutx3_120),c.applyCuts(fpi,cutx3_120),yerr=np.sqrt(c.applyCuts(lumi,cutx3_120))/c.applyCuts(lumi,cutx3_120),fmt='.',label='$x_{\pi}$=(0.3,0.4)',ecolor='cyan',capsize=2, capthick=2,color='c')
    tscat4 = ax.errorbar(c.applyCuts(t,cutx4_120),c.applyCuts(fpi,cutx4_120),yerr=np.sqrt(c.applyCuts(lumi,cutx4_120))/c.applyCuts(lumi,cutx4_120),fmt='.',label='$x_{\pi}$=(0.4,0.5)',ecolor='cyan',capsize=2, capthick=2,color='m')
    tscat5 = ax.errorbar(c.applyCuts(t,cutx5_120),c.applyCuts(fpi,cutx5_120),yerr=np.sqrt(c.applyCuts(lumi,cutx5_120))/c.applyCuts(lumi,cutx5_120),fmt='.',label='$x_{\pi}$=(0.5,0.6)',ecolor='cyan',capsize=2, capthick=2,color='y')
    tscat6 = ax.errorbar(c.applyCuts(t,cutx6_120),c.applyCuts(fpi,cutx6_120),yerr=np.sqrt(c.applyCuts(lumi,cutx6_120))/c.applyCuts(lumi,cutx6_120),fmt='.',label='$x_{\pi}$=(0.6,0.7)',ecolor='cyan',capsize=2, capthick=2,color='purple')
    tscat7 = ax.errorbar(c.applyCuts(t,cutx7_120),c.applyCuts(fpi,cutx7_120),yerr=np.sqrt(c.applyCuts(lumi,cutx7_120))/c.applyCuts(lumi,cutx7_120),fmt='.',label='$x_{\pi}$=(0.7,0.8)',ecolor='cyan',capsize=2, capthick=2,color='aqua')
    tscat8 = ax.errorbar(c.applyCuts(t,cutx8_120),c.applyCuts(fpi,cutx8_120),yerr=np.sqrt(c.applyCuts(lumi,cutx8_120))/c.applyCuts(lumi,cutx8_120),fmt='.',label='$x_{\pi}$=(0.8,0.9)',ecolor='cyan',capsize=2, capthick=2,color='gray')
    tscat9 = ax.errorbar(c.applyCuts(t,cutx9_120),c.applyCuts(fpi,cutx9_120),yerr=np.sqrt(c.applyCuts(lumi,cutx9_120))/c.applyCuts(lumi,cutx9_120),fmt='.',label='$x_{\pi}$=(0.9,1.0)',ecolor='cyan',capsize=2, capthick=2,color='pink')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(-1.,0.)
    plt.ylim(-0.1,0.3)
    ax.text(0.85, 0.85, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    plt.xlabel('-t', fontsize=20)
    plt.tight_layout()

    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    
    f.savefig('OUTPUTS/fpivt_xbin.png')

def fpivtPrime_Plot():

    f = plt.figure(figsize=(11.69,8.27))
    
    ax = f.add_subplot(421)
    tscat1 = ax.errorbar(c.applyCuts(tPrime,cut7),c.applyCuts(fpi,cut7),yerr=np.sqrt(c.applyCuts(lumi,cut7))/c.applyCuts(lumi,cut7),fmt='.',label='$Q^2$=7 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    # plt.yscale('log')
    plt.xlim(-1.0,0.)
    plt.ylim(-0.1,0.3)
    ax.text(0.85, 0.45, '$Q^2$=7 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    plt.ylabel('$F^{\pi}_{2}$', fontsize=20)

    ax = f.add_subplot(422)
    tscat2 = ax.errorbar(c.applyCuts(tPrime,cut15),c.applyCuts(fpi,cut15),yerr=np.sqrt(c.applyCuts(lumi,cut15))/c.applyCuts(lumi,cut15),fmt='.',label='$Q^2$=15 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    # plt.yscale('log')
    plt.xlim(-1.0,0.)
    plt.ylim(-0.1,0.3)
    ax.text(0.85, 0.45, '$Q^2$=15 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_formatter(plt.NullFormatter())

    plt.title('$F^{\pi}_{2}$ vs -t\'', fontsize =20)
    
    ax = f.add_subplot(423)
    tscat3 = ax.errorbar(c.applyCuts(tPrime,cut30),c.applyCuts(fpi,cut30),yerr=np.sqrt(c.applyCuts(lumi,cut30))/c.applyCuts(lumi,cut30),fmt='.',label='$Q^2$=30 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    # plt.yscale('log')
    plt.xlim(-1.0,0.)
    plt.ylim(-0.1,0.3)
    ax.text(0.85, 0.45, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    
    ax = f.add_subplot(424)
    tscat4 = ax.errorbar(c.applyCuts(tPrime,cut60),c.applyCuts(fpi,cut60),yerr=np.sqrt(c.applyCuts(lumi,cut60))/c.applyCuts(lumi,cut60),fmt='.',label='$Q^2$=60 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    # plt.yscale('log')
    plt.xlim(-1.0,0.)
    plt.ylim(-0.1,0.3)
    ax.text(0.85, 0.45, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    
    ax = f.add_subplot(425)
    tscat5 = ax.errorbar(c.applyCuts(tPrime,cut120),c.applyCuts(fpi,cut120),yerr=np.sqrt(c.applyCuts(lumi,cut120))/c.applyCuts(lumi,cut120),fmt='.',label='$Q^2$=120 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    # plt.yscale('log')
    plt.xlim(-1.0,0.)
    plt.ylim(-0.1,0.3)
    ax.text(0.85, 0.45, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    ax = f.add_subplot(426)
    tscat6 = ax.errorbar(c.applyCuts(tPrime,cut240),c.applyCuts(fpi,cut240),yerr=np.sqrt(c.applyCuts(lumi,cut240))/c.applyCuts(lumi,cut240),fmt='.',label='$Q^2$=240 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    # plt.yscale('log')
    plt.xlim(-1.0,0.)
    plt.ylim(-0.1,0.3)
    ax.text(0.85, 0.45, '$Q^2$=240 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(427)
    tscat7 = ax.errorbar(c.applyCuts(tPrime,cut480),c.applyCuts(fpi,cut480),yerr=np.sqrt(c.applyCuts(lumi,cut480))/c.applyCuts(lumi,cut480),fmt='.',label='$Q^2$=480 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    # plt.yscale('log')
    plt.xlim(-1.0,0.)
    plt.ylim(-0.1,0.3)
    ax.text(0.85, 0.45, '$Q^2$=480 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))

    ax = f.add_subplot(428)
    tscat8 = ax.errorbar(c.applyCuts(tPrime,cut1000),c.applyCuts(fpi,cut1000),yerr=np.sqrt(c.applyCuts(lumi,cut1000))/c.applyCuts(lumi,cut1000),fmt='.',label='$Q^2$=1000 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    # plt.yscale('log')
    plt.xlim(-1.0,0.)
    plt.ylim(-0.1,0.3)
    ax.text(0.85, 0.45, '$Q^2$=1000 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    plt.xlabel('-t\'', fontsize=20)
    plt.tight_layout()

    plt.subplots_adjust(hspace=0.0,wspace=0.0)

    f.savefig('OUTPUTS/fpivtPrime.png')

def fpivtPrime_xbin_Plot():

    f = plt.figure(figsize=(11.69,8.27))
    
    ax = f.add_subplot(221)
    tscat0 = ax.errorbar(c.applyCuts(tPrime,cutx0_15),c.applyCuts(fpi,cutx0_15),yerr=np.sqrt(c.applyCuts(lumi,cutx0_15))/c.applyCuts(lumi,cutx0_15),fmt='.',label='$x_{\pi}$=(0.0,0.1)',ecolor='cyan',capsize=2, capthick=2,color='b')
    tscat1 = ax.errorbar(c.applyCuts(tPrime,cutx1_15),c.applyCuts(fpi,cutx1_15),yerr=np.sqrt(c.applyCuts(lumi,cutx1_15))/c.applyCuts(lumi,cutx1_15),fmt='.',label='$x_{\pi}$=(0.1,0.2)',ecolor='cyan',capsize=2, capthick=2,color='g')
    tscat2 = ax.errorbar(c.applyCuts(tPrime,cutx2_15),c.applyCuts(fpi,cutx2_15),yerr=np.sqrt(c.applyCuts(lumi,cutx2_15))/c.applyCuts(lumi,cutx2_15),fmt='.',label='$x_{\pi}$=(0.2,0.3)',ecolor='cyan',capsize=2, capthick=2,color='r')
    tscat3 = ax.errorbar(c.applyCuts(tPrime,cutx3_15),c.applyCuts(fpi,cutx3_15),yerr=np.sqrt(c.applyCuts(lumi,cutx3_15))/c.applyCuts(lumi,cutx3_15),fmt='.',label='$x_{\pi}$=(0.3,0.4)',ecolor='cyan',capsize=2, capthick=2,color='c')
    tscat4 = ax.errorbar(c.applyCuts(tPrime,cutx4_15),c.applyCuts(fpi,cutx4_15),yerr=np.sqrt(c.applyCuts(lumi,cutx4_15))/c.applyCuts(lumi,cutx4_15),fmt='.',label='$x_{\pi}$=(0.4,0.5)',ecolor='cyan',capsize=2, capthick=2,color='m')
    tscat5 = ax.errorbar(c.applyCuts(tPrime,cutx5_15),c.applyCuts(fpi,cutx5_15),yerr=np.sqrt(c.applyCuts(lumi,cutx5_15))/c.applyCuts(lumi,cutx5_15),fmt='.',label='$x_{\pi}$=(0.5,0.6)',ecolor='cyan',capsize=2, capthick=2,color='y')
    tscat6 = ax.errorbar(c.applyCuts(tPrime,cutx6_15),c.applyCuts(fpi,cutx6_15),yerr=np.sqrt(c.applyCuts(lumi,cutx6_15))/c.applyCuts(lumi,cutx6_15),fmt='.',label='$x_{\pi}$=(0.6,0.7)',ecolor='cyan',capsize=2, capthick=2,color='purple')
    tscat7 = ax.errorbar(c.applyCuts(tPrime,cutx7_15),c.applyCuts(fpi,cutx7_15),yerr=np.sqrt(c.applyCuts(lumi,cutx7_15))/c.applyCuts(lumi,cutx7_15),fmt='.',label='$x_{\pi}$=(0.7,0.8)',ecolor='cyan',capsize=2, capthick=2,color='aqua')
    tscat8 = ax.errorbar(c.applyCuts(tPrime,cutx8_15),c.applyCuts(fpi,cutx8_15),yerr=np.sqrt(c.applyCuts(lumi,cutx8_15))/c.applyCuts(lumi,cutx8_15),fmt='.',label='$x_{\pi}$=(0.8,0.9)',ecolor='cyan',capsize=2, capthick=2,color='gray')
    tscat9 = ax.errorbar(c.applyCuts(tPrime,cutx9_15),c.applyCuts(fpi,cutx9_15),yerr=np.sqrt(c.applyCuts(lumi,cutx9_15))/c.applyCuts(lumi,cutx9_15),fmt='.',label='$x_{\pi}$=(0.9,1.0)',ecolor='cyan',capsize=2, capthick=2,color='pink')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(-1.,0.)
    plt.ylim(-0.1,0.3)
    ax.text(0.65, 0.65, '$Q^2$=15 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='lower'))

    plt.ylabel('$F^{\pi}_{2}$', fontsize=20)
    
    ax = f.add_subplot(222)
    tscat0 = ax.errorbar(c.applyCuts(tPrime,cutx0_30),c.applyCuts(fpi,cutx0_30),yerr=np.sqrt(c.applyCuts(lumi,cutx0_30))/c.applyCuts(lumi,cutx0_30),fmt='.',label='$x_{\pi}$=(0.0,0.1)',ecolor='cyan',capsize=2, capthick=2,color='b')
    tscat1 = ax.errorbar(c.applyCuts(tPrime,cutx1_30),c.applyCuts(fpi,cutx1_30),yerr=np.sqrt(c.applyCuts(lumi,cutx1_30))/c.applyCuts(lumi,cutx1_30),fmt='.',label='$x_{\pi}$=(0.1,0.2)',ecolor='cyan',capsize=2, capthick=2,color='g')
    tscat2 = ax.errorbar(c.applyCuts(tPrime,cutx2_30),c.applyCuts(fpi,cutx2_30),yerr=np.sqrt(c.applyCuts(lumi,cutx2_30))/c.applyCuts(lumi,cutx2_30),fmt='.',label='$x_{\pi}$=(0.2,0.3)',ecolor='cyan',capsize=2, capthick=2,color='r')
    tscat3 = ax.errorbar(c.applyCuts(tPrime,cutx3_30),c.applyCuts(fpi,cutx3_30),yerr=np.sqrt(c.applyCuts(lumi,cutx3_30))/c.applyCuts(lumi,cutx3_30),fmt='.',label='$x_{\pi}$=(0.3,0.4)',ecolor='cyan',capsize=2, capthick=2,color='c')
    tscat4 = ax.errorbar(c.applyCuts(tPrime,cutx4_30),c.applyCuts(fpi,cutx4_30),yerr=np.sqrt(c.applyCuts(lumi,cutx4_30))/c.applyCuts(lumi,cutx4_30),fmt='.',label='$x_{\pi}$=(0.4,0.5)',ecolor='cyan',capsize=2, capthick=2,color='m')
    tscat5 = ax.errorbar(c.applyCuts(tPrime,cutx5_30),c.applyCuts(fpi,cutx5_30),yerr=np.sqrt(c.applyCuts(lumi,cutx5_30))/c.applyCuts(lumi,cutx5_30),fmt='.',label='$x_{\pi}$=(0.5,0.6)',ecolor='cyan',capsize=2, capthick=2,color='y')
    tscat6 = ax.errorbar(c.applyCuts(tPrime,cutx6_30),c.applyCuts(fpi,cutx6_30),yerr=np.sqrt(c.applyCuts(lumi,cutx6_30))/c.applyCuts(lumi,cutx6_30),fmt='.',label='$x_{\pi}$=(0.6,0.7)',ecolor='cyan',capsize=2, capthick=2,color='purple')
    tscat7 = ax.errorbar(c.applyCuts(tPrime,cutx7_30),c.applyCuts(fpi,cutx7_30),yerr=np.sqrt(c.applyCuts(lumi,cutx7_30))/c.applyCuts(lumi,cutx7_30),fmt='.',label='$x_{\pi}$=(0.7,0.8)',ecolor='cyan',capsize=2, capthick=2,color='aqua')
    tscat8 = ax.errorbar(c.applyCuts(tPrime,cutx8_30),c.applyCuts(fpi,cutx8_30),yerr=np.sqrt(c.applyCuts(lumi,cutx8_30))/c.applyCuts(lumi,cutx8_30),fmt='.',label='$x_{\pi}$=(0.8,0.9)',ecolor='cyan',capsize=2, capthick=2,color='gray')
    tscat9 = ax.errorbar(c.applyCuts(tPrime,cutx9_30),c.applyCuts(fpi,cutx9_30),yerr=np.sqrt(c.applyCuts(lumi,cutx9_30))/c.applyCuts(lumi,cutx9_30),fmt='.',label='$x_{\pi}$=(0.9,1.0)',ecolor='cyan',capsize=2, capthick=2,color='pink')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(-1.,0.)
    plt.ylim(-0.1,0.3)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.text(0.65, 0.65, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    plt.title('$F^{\pi}_{2}$ vs -t\'', fontsize =20)
    
    ax = f.add_subplot(223)
    tscat0 = ax.errorbar(c.applyCuts(tPrime,cutx0_60),c.applyCuts(fpi,cutx0_60),yerr=np.sqrt(c.applyCuts(lumi,cutx0_60))/c.applyCuts(lumi,cutx0_60),fmt='.',label='$x_{\pi}$=(0.0,0.1)',ecolor='cyan',capsize=2, capthick=2,color='b')
    tscat1 = ax.errorbar(c.applyCuts(tPrime,cutx1_60),c.applyCuts(fpi,cutx1_60),yerr=np.sqrt(c.applyCuts(lumi,cutx1_60))/c.applyCuts(lumi,cutx1_60),fmt='.',label='$x_{\pi}$=(0.1,0.2)',ecolor='cyan',capsize=2, capthick=2,color='g')
    tscat2 = ax.errorbar(c.applyCuts(tPrime,cutx2_60),c.applyCuts(fpi,cutx2_60),yerr=np.sqrt(c.applyCuts(lumi,cutx2_60))/c.applyCuts(lumi,cutx2_60),fmt='.',label='$x_{\pi}$=(0.2,0.3)',ecolor='cyan',capsize=2, capthick=2,color='r')
    tscat3 = ax.errorbar(c.applyCuts(tPrime,cutx3_60),c.applyCuts(fpi,cutx3_60),yerr=np.sqrt(c.applyCuts(lumi,cutx3_60))/c.applyCuts(lumi,cutx3_60),fmt='.',label='$x_{\pi}$=(0.3,0.4)',ecolor='cyan',capsize=2, capthick=2,color='c')
    tscat4 = ax.errorbar(c.applyCuts(tPrime,cutx4_60),c.applyCuts(fpi,cutx4_60),yerr=np.sqrt(c.applyCuts(lumi,cutx4_60))/c.applyCuts(lumi,cutx4_60),fmt='.',label='$x_{\pi}$=(0.4,0.5)',ecolor='cyan',capsize=2, capthick=2,color='m')
    tscat5 = ax.errorbar(c.applyCuts(tPrime,cutx5_60),c.applyCuts(fpi,cutx5_60),yerr=np.sqrt(c.applyCuts(lumi,cutx5_60))/c.applyCuts(lumi,cutx5_60),fmt='.',label='$x_{\pi}$=(0.5,0.6)',ecolor='cyan',capsize=2, capthick=2,color='y')
    tscat6 = ax.errorbar(c.applyCuts(tPrime,cutx6_60),c.applyCuts(fpi,cutx6_60),yerr=np.sqrt(c.applyCuts(lumi,cutx6_60))/c.applyCuts(lumi,cutx6_60),fmt='.',label='$x_{\pi}$=(0.6,0.7)',ecolor='cyan',capsize=2, capthick=2,color='purple')
    tscat7 = ax.errorbar(c.applyCuts(tPrime,cutx7_60),c.applyCuts(fpi,cutx7_60),yerr=np.sqrt(c.applyCuts(lumi,cutx7_60))/c.applyCuts(lumi,cutx7_60),fmt='.',label='$x_{\pi}$=(0.7,0.8)',ecolor='cyan',capsize=2, capthick=2,color='aqua')
    tscat8 = ax.errorbar(c.applyCuts(tPrime,cutx8_60),c.applyCuts(fpi,cutx8_60),yerr=np.sqrt(c.applyCuts(lumi,cutx8_60))/c.applyCuts(lumi,cutx8_60),fmt='.',label='$x_{\pi}$=(0.8,0.9)',ecolor='cyan',capsize=2, capthick=2,color='gray')
    tscat9 = ax.errorbar(c.applyCuts(tPrime,cutx9_60),c.applyCuts(fpi,cutx9_60),yerr=np.sqrt(c.applyCuts(lumi,cutx9_60))/c.applyCuts(lumi,cutx9_60),fmt='.',label='$x_{\pi}$=(0.9,1.0)',ecolor='cyan',capsize=2, capthick=2,color='pink')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(-1.,0.)
    plt.ylim(-0.1,0.3)
    ax.text(0.65, 0.65, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    
    ax = f.add_subplot(224)
    tscat0 = ax.errorbar(c.applyCuts(tPrime,cutx0_120),c.applyCuts(fpi,cutx0_120),yerr=np.sqrt(c.applyCuts(lumi,cutx0_120))/c.applyCuts(lumi,cutx0_120),fmt='.',label='$x_{\pi}$=(0.0,0.1)',ecolor='cyan',capsize=2, capthick=2,color='b')
    tscat1 = ax.errorbar(c.applyCuts(tPrime,cutx1_120),c.applyCuts(fpi,cutx1_120),yerr=np.sqrt(c.applyCuts(lumi,cutx1_120))/c.applyCuts(lumi,cutx1_120),fmt='.',label='$x_{\pi}$=(0.1,0.2)',ecolor='cyan',capsize=2, capthick=2,color='g')
    tscat2 = ax.errorbar(c.applyCuts(tPrime,cutx2_120),c.applyCuts(fpi,cutx2_120),yerr=np.sqrt(c.applyCuts(lumi,cutx2_120))/c.applyCuts(lumi,cutx2_120),fmt='.',label='$x_{\pi}$=(0.2,0.3)',ecolor='cyan',capsize=2, capthick=2,color='r')
    tscat3 = ax.errorbar(c.applyCuts(tPrime,cutx3_120),c.applyCuts(fpi,cutx3_120),yerr=np.sqrt(c.applyCuts(lumi,cutx3_120))/c.applyCuts(lumi,cutx3_120),fmt='.',label='$x_{\pi}$=(0.3,0.4)',ecolor='cyan',capsize=2, capthick=2,color='c')
    tscat4 = ax.errorbar(c.applyCuts(tPrime,cutx4_120),c.applyCuts(fpi,cutx4_120),yerr=np.sqrt(c.applyCuts(lumi,cutx4_120))/c.applyCuts(lumi,cutx4_120),fmt='.',label='$x_{\pi}$=(0.4,0.5)',ecolor='cyan',capsize=2, capthick=2,color='m')
    tscat5 = ax.errorbar(c.applyCuts(tPrime,cutx5_120),c.applyCuts(fpi,cutx5_120),yerr=np.sqrt(c.applyCuts(lumi,cutx5_120))/c.applyCuts(lumi,cutx5_120),fmt='.',label='$x_{\pi}$=(0.5,0.6)',ecolor='cyan',capsize=2, capthick=2,color='y')
    tscat6 = ax.errorbar(c.applyCuts(tPrime,cutx6_120),c.applyCuts(fpi,cutx6_120),yerr=np.sqrt(c.applyCuts(lumi,cutx6_120))/c.applyCuts(lumi,cutx6_120),fmt='.',label='$x_{\pi}$=(0.6,0.7)',ecolor='cyan',capsize=2, capthick=2,color='purple')
    tscat7 = ax.errorbar(c.applyCuts(tPrime,cutx7_120),c.applyCuts(fpi,cutx7_120),yerr=np.sqrt(c.applyCuts(lumi,cutx7_120))/c.applyCuts(lumi,cutx7_120),fmt='.',label='$x_{\pi}$=(0.7,0.8)',ecolor='cyan',capsize=2, capthick=2,color='aqua')
    tscat8 = ax.errorbar(c.applyCuts(tPrime,cutx8_120),c.applyCuts(fpi,cutx8_120),yerr=np.sqrt(c.applyCuts(lumi,cutx8_120))/c.applyCuts(lumi,cutx8_120),fmt='.',label='$x_{\pi}$=(0.8,0.9)',ecolor='cyan',capsize=2, capthick=2,color='gray')
    tscat9 = ax.errorbar(c.applyCuts(tPrime,cutx9_120),c.applyCuts(fpi,cutx9_120),yerr=np.sqrt(c.applyCuts(lumi,cutx9_120))/c.applyCuts(lumi,cutx9_120),fmt='.',label='$x_{\pi}$=(0.9,1.0)',ecolor='cyan',capsize=2, capthick=2,color='pink')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(-1.,0.)
    plt.ylim(-0.1,0.3)
    ax.text(0.65, 0.65, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    plt.xlabel('-t\'', fontsize=20)
    plt.tight_layout()

    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    
    f.savefig('OUTPUTS/fpivtPrime_xbin.png')    

def plot3D():
    f = plt.figure(figsize=(11.69,8.27))
    ax = f.gca(projection='3d')
    
    xQfpi = ax.scatter3D(Q2, xpi, fpi, c=fpi)
    ax.set_xlabel('$Q^2$')
    ax.set_ylabel('$x_{\pi}$')
    ax.set_zlabel('$F^{\pi}_{2}$')
    ax.invert_xaxis()
    f.colorbar(xQfpi)
    
def Lumi():

    # Luminosity
    
    # all binned events
    sig_all = tot_sigma*1e5 # the 1e5 corrects the scaling up top
    evts_all = len(tot_sigma)*[len(tot_sigma)]

    lumi = [(e)/((s)*(binQ2)*(binx)) for e,s in zip(evts_all,sig_all)]
    print("---------------------------------\n")
    # print("\nLuminosity: ", lumi)

    nevt = [100/(l*1e-6) for l in lumi] # The 1e-6 converts properly, integrated luminosiy: 100 fb^-1
    nevt  = np.asarray(nevt)
    print("\nEvents expected running at 100 $fb^{-1}$: ", nevt)

    # f = plt.figure(figsize=(11.69,8.27))
    
    # ax = f.add_subplot(331)
    # ax.hist(nevt[0],bins=c.setbin(nevt[0],200),label='$Q^2$=7 $GeV^2$',histtype='step', alpha=0.5, stacked=True, fill=True)
    # plt.xscale('log')
    # ax.text(0.60, 0.85, '$Q^2$=7 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')

    # ax = f.add_subplot(332)
    # ax.hist(nevt[1],bins=c.setbin(nevt[1],200),label='$Q^2$=15 $GeV^2$',histtype='step', alpha=0.5, stacked=True, fill=True)
    # plt.xscale('log')
    # ax.text(0.60, 0.85, '$Q^2$=15 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    
    # ax = f.add_subplot(333)
    # ax.hist(nevt[2],bins=c.setbin(nevt[2],200),label='$Q^2$=30 $GeV^2$',histtype='step', alpha=0.5, stacked=True, fill=True)
    # plt.xscale('log')
    # ax.text(0.60, 0.85, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')

    # plt.title('Events expected running at 10 $fb^{-1}$', fontsize =20)
    
    # ax = f.add_subplot(334)
    # ax.hist(nevt[3],bins=c.setbin(nevt[3],200),label='$Q^2$=60 $GeV^2$',histtype='step', alpha=0.5, stacked=True, fill=True)
    # plt.xscale('log')
    # ax.text(0.60, 0.85, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    
    # ax = f.add_subplot(335)
    # ax.hist(nevt[4],bins=c.setbin(nevt[4],200),label='$Q^2$=120 $GeV^2$',histtype='step', alpha=0.5, stacked=True, fill=True)
    # plt.xscale('log')
    # ax.text(0.60, 0.85, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    
    # ax = f.add_subplot(336)
    # ax.hist(nevt[5],bins=c.setbin(nevt[5],200),label='$Q^2$=240 $GeV^2$',histtype='step', alpha=0.5, stacked=True, fill=True)
    # plt.xscale('log')
    # ax.text(0.60, 0.85, '$Q^2$=240 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')

    # ax = f.add_subplot(337)
    # ax.hist(nevt[6],bins=c.setbin(nevt[6],200),label='$Q^2$=480 $GeV^2$',histtype='step', alpha=0.5, stacked=True, fill=True)
    # plt.xscale('log')
    # ax.text(0.60, 0.85, '$Q^2$=480 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')

    # ax = f.add_subplot(338)
    # ax.hist(nevt[7],bins=c.setbin(nevt[7],200),label='$Q^2$=1000 $GeV^2$',histtype='step', alpha=0.5, stacked=True, fill=True)
    # plt.xscale('log')
    # ax.text(0.60, 0.85, '$Q^2$=1000 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')

    # plt.tight_layout()

    # return [nevt_7,nevt_15,nevt_30,nevt_60,nevt_120,nevt_240,nevt_480,nevt_1000,nevt]
    return nevt

lumi = Lumi()

def theory_table():

    #################################
    ymin_cutxpi001_15 = min(c.applyCuts(y,cutxpi001_15))
    ymin_cutxpi01_15 = min(c.applyCuts(y,cutxpi01_15))
    
    ymax_cutxpi001_15 = max(c.applyCuts(y,cutxpi001_15))
    ymax_cutxpi01_15 = max(c.applyCuts(y,cutxpi01_15))

    yieldmin_cutxpi001_15 = float(c.applyCuts(lumi,cutxpi001_15)[
        np.where(c.applyCuts(y,cutxpi001_15)==(ymin_cutxpi001_15))])
    yieldmin_cutxpi01_15 = float(c.applyCuts(lumi,cutxpi01_15)[
        np.where(c.applyCuts(y,cutxpi01_15)==(ymin_cutxpi01_15))])

    yieldmax_cutxpi001_15 = float(c.applyCuts(lumi,cutxpi001_15)[
        np.where(c.applyCuts(y,cutxpi001_15)==(ymax_cutxpi001_15))])
    yieldmax_cutxpi01_15 = float(c.applyCuts(lumi,cutxpi01_15)[
        np.where(c.applyCuts(y,cutxpi01_15)==(ymax_cutxpi01_15))])

    xLmin_cutxpi001_15 = float(c.applyCuts(xL,cutxpi001_15)[
        np.where(c.applyCuts(y,cutxpi001_15)==(ymin_cutxpi001_15))])
    xLmin_cutxpi01_15 = float(c.applyCuts(xL,cutxpi01_15)[
        np.where(c.applyCuts(y,cutxpi01_15)==(ymin_cutxpi01_15))])
    
    xLmax_cutxpi001_15 = float(c.applyCuts(xL,cutxpi001_15)[
        np.where(c.applyCuts(y,cutxpi001_15)==(ymax_cutxpi001_15))])
    xLmax_cutxpi01_15 = float(c.applyCuts(xL,cutxpi01_15)[
        np.where(c.applyCuts(y,cutxpi01_15)==(ymax_cutxpi01_15))])

    xpimin_cutxpi001_15 = float(c.applyCuts(xpi,cutxpi001_15)[
        np.where(c.applyCuts(y,cutxpi001_15)==(ymin_cutxpi001_15))])
    xpimin_cutxpi01_15 = float(c.applyCuts(xpi,cutxpi01_15)[
        np.where(c.applyCuts(y,cutxpi01_15)==(ymin_cutxpi01_15))])
    
    xpimax_cutxpi001_15 = float(c.applyCuts(xpi,cutxpi001_15)[
        np.where(c.applyCuts(y,cutxpi001_15)==(ymax_cutxpi001_15))])
    xpimax_cutxpi01_15 = float(c.applyCuts(xpi,cutxpi01_15)[
        np.where(c.applyCuts(y,cutxpi01_15)==(ymax_cutxpi01_15))])

    Q2min_cutxpi001_15 = float(c.applyCuts(Q2,cutxpi001_15)[
        np.where(c.applyCuts(y,cutxpi001_15)==(ymin_cutxpi001_15))])
    Q2min_cutxpi01_15 = float(c.applyCuts(Q2,cutxpi01_15)[
        np.where(c.applyCuts(y,cutxpi01_15)==(ymin_cutxpi01_15))])
    
    Q2max_cutxpi001_15 = float(c.applyCuts(Q2,cutxpi001_15)[
        np.where(c.applyCuts(y,cutxpi001_15)==(ymax_cutxpi001_15))])
    Q2max_cutxpi01_15 = float(c.applyCuts(Q2,cutxpi01_15)[
        np.where(c.applyCuts(y,cutxpi01_15)==(ymax_cutxpi01_15))])
    #################################    
    
    #################################
    ymin_cutxpi01_30 = min(c.applyCuts(y,cutxpi01_30))
    
    ymax_cutxpi01_30 = max(c.applyCuts(y,cutxpi01_30))
    
    yieldmin_cutxpi01_30 = float(c.applyCuts(lumi,cutxpi01_30)[
        np.where(c.applyCuts(y,cutxpi01_30)==(ymin_cutxpi01_30))])

    yieldmax_cutxpi01_30 = float(c.applyCuts(lumi,cutxpi01_30)[
        np.where(c.applyCuts(y,cutxpi01_30)==(ymax_cutxpi01_30))])

    xLmin_cutxpi01_30 = float(c.applyCuts(xL,cutxpi01_30)[
        np.where(c.applyCuts(y,cutxpi01_30)==(ymin_cutxpi01_30))])

    xLmax_cutxpi01_30 = float(c.applyCuts(xL,cutxpi01_30)[
        np.where(c.applyCuts(y,cutxpi01_30)==(ymax_cutxpi01_30))])

    xpimin_cutxpi01_30 = float(c.applyCuts(xpi,cutxpi01_30)[
        np.where(c.applyCuts(y,cutxpi01_30)==(ymin_cutxpi01_30))])

    xpimax_cutxpi01_30 = float(c.applyCuts(xpi,cutxpi01_30)[
        np.where(c.applyCuts(y,cutxpi01_30)==(ymax_cutxpi01_30))])

    Q2min_cutxpi01_30 = float(c.applyCuts(Q2,cutxpi01_30)[
        np.where(c.applyCuts(y,cutxpi01_30)==(ymin_cutxpi01_30))])

    Q2max_cutxpi01_30 = float(c.applyCuts(Q2,cutxpi01_30)[
        np.where(c.applyCuts(y,cutxpi01_30)==(ymax_cutxpi01_30))])

    #################################
    
    #################################
    ymin_cutxpi01_60 = min(c.applyCuts(y,cutxpi01_60))
    
    ymax_cutxpi01_60 = max(c.applyCuts(y,cutxpi01_60))

    yieldmin_cutxpi01_60 = float(c.applyCuts(lumi,cutxpi01_60)[
        np.where(c.applyCuts(y,cutxpi01_60)==(ymin_cutxpi01_60))])

    yieldmax_cutxpi01_60 = float(c.applyCuts(lumi,cutxpi01_60)[
        np.where(c.applyCuts(y,cutxpi01_60)==(ymax_cutxpi01_60))])

    xLmin_cutxpi01_60 = float(c.applyCuts(xL,cutxpi01_60)[
        np.where(c.applyCuts(y,cutxpi01_60)==(ymin_cutxpi01_60))])
    
    xLmax_cutxpi01_60 = float(c.applyCuts(xL,cutxpi01_60)[
        np.where(c.applyCuts(y,cutxpi01_60)==(ymax_cutxpi01_60))])

    xpimin_cutxpi01_60 = float(c.applyCuts(xpi,cutxpi01_60)[
        np.where(c.applyCuts(y,cutxpi01_60)==(ymin_cutxpi01_60))])
    
    xpimax_cutxpi01_60 = float(c.applyCuts(xpi,cutxpi01_60)[
        np.where(c.applyCuts(y,cutxpi01_60)==(ymax_cutxpi01_60))])

    Q2min_cutxpi01_60 = float(c.applyCuts(Q2,cutxpi01_60)[
        np.where(c.applyCuts(y,cutxpi01_60)==(ymin_cutxpi01_60))])
    
    Q2max_cutxpi01_60 = float(c.applyCuts(Q2,cutxpi01_60)[
        np.where(c.applyCuts(y,cutxpi01_60)==(ymax_cutxpi01_60))])
    #################################
    
    table_dict = {
    
            "15" : 
             {"0.001" :{
                 "Q2" : [Q2min_cutxpi001_15,Q2max_cutxpi001_15],
                 "y bin" : [ymin_cutxpi001_15,ymax_cutxpi001_15],
                 "Yield (100 fb^-1)" : [yieldmin_cutxpi001_15,yieldmax_cutxpi001_15],
                 "xL" : [xLmin_cutxpi001_15,xLmax_cutxpi001_15],
                 "xpi": [xpimin_cutxpi001_15,xpimax_cutxpi001_15]
             },
              "0.01":{
                  "Q2" : [Q2min_cutxpi01_15,Q2max_cutxpi01_15],
                  "y bin" : [ymin_cutxpi01_15,ymax_cutxpi01_15],
                  "Yield (100 fb^-1)" : [yieldmin_cutxpi01_15,yieldmax_cutxpi01_15],
                  "xL" : [xLmin_cutxpi01_15,xLmax_cutxpi01_15],
                  "xpi": [xpimin_cutxpi01_15,xpimax_cutxpi01_15]
             }
            },
            "30" : 
             {"0.01":{
                 "Q2" : [Q2min_cutxpi01_30,Q2max_cutxpi01_30],
                 "y bin" : [ymin_cutxpi01_30,ymax_cutxpi01_30],
                 "Yield (100 fb^-1)" : [yieldmin_cutxpi01_30,yieldmax_cutxpi01_30],
                 "xL" : [xLmin_cutxpi01_30,xLmax_cutxpi01_30],
                 "xpi": [xpimin_cutxpi01_30,xpimax_cutxpi01_30]
             }
            },
            "60" :
             {"0.01":{
                 "Q2" : [Q2min_cutxpi01_60,Q2max_cutxpi01_60],
                 "y bin" : [ymin_cutxpi01_60,ymax_cutxpi01_60],
                 "Yield (100 fb^-1)" : [yieldmin_cutxpi01_60,yieldmax_cutxpi01_60],
                 "xL" : [xLmin_cutxpi01_60,xLmax_cutxpi01_60],
                 "xpi": [xpimin_cutxpi01_60,xpimax_cutxpi01_60]
             }
            }
    }

    #################################
        
    
    # ymin_cutxpi1_15 = min(c.applyCuts(y,cutxpi1_15))
    # ymin_cutxpi2_15 = min(c.applyCuts(y,cutxpi2_15))
    # ymin_cutxpi3_15 = min(c.applyCuts(y,cutxpi3_15))
    # # ymin_cutxpi4_15 = min(c.applyCuts(y,cutxpi4_15))
    # # ymin_cutxpi5_15 = min(c.applyCuts(y,cutxpi5_15))
    # # ymin_cutxpi6_15 = min(c.applyCuts(y,cutxpi6_15))
    # # ymin_cutxpi7_15 = min(c.applyCuts(y,cutxpi7_15))
    
    # ymax_cutxpi1_15 = max(c.applyCuts(y,cutxpi1_15))
    # ymax_cutxpi2_15 = max(c.applyCuts(y,cutxpi2_15))
    # ymax_cutxpi3_15 = max(c.applyCuts(y,cutxpi3_15))
    # # ymax_cutxpi4_15 = max(c.applyCuts(y,cutxpi4_15))
    # # ymax_cutxpi5_15 = max(c.applyCuts(y,cutxpi5_15))
    # # ymax_cutxpi6_15 = max(c.applyCuts(y,cutxpi6_15))
    # # ymax_cutxpi7_15 = max(c.applyCuts(y,cutxpi7_15))

    # yieldmin_cutxpi1_15 = float(c.applyCuts(lumi,cutxpi1_15)[
    #     np.where(c.applyCuts(y,cutxpi1_15)==(ymin_cutxpi1_15))])
    # yieldmin_cutxpi2_15 = float(c.applyCuts(lumi,cutxpi2_15)[
    #     np.where(c.applyCuts(y,cutxpi2_15)==(ymin_cutxpi2_15))])
    # yieldmin_cutxpi3_15 = float(c.applyCuts(lumi,cutxpi3_15)[
    #     np.where(c.applyCuts(y,cutxpi3_15)==(ymin_cutxpi3_15))])
    # # yieldmin_cutxpi4_15 = float(c.applyCuts(lumi,cutxpi4_15)[
    # #     np.where(c.applyCuts(y,cutxpi4_15)==(ymin_cutxpi4_15))])
    # # yieldmin_cutxpi5_15 = float(c.applyCuts(lumi,cutxpi5_15)[
    # #     np.where(c.applyCuts(y,cutxpi5_15)==(ymin_cutxpi5_15))])
    # # yieldmin_cutxpi6_15 = float(c.applyCuts(lumi,cutxpi6_15)[
    # #     np.where(c.applyCuts(y,cutxpi6_15)==(ymin_cutxpi6_15))])
    # # yieldmin_cutxpi7_15 = float(c.applyCuts(lumi,cutxpi7_15)[
    # #     np.where(c.applyCuts(y,cutxpi7_15)==(ymin_cutxpi7_15))])
    
    # yieldmax_cutxpi1_15 = float(c.applyCuts(lumi,cutxpi1_15)[
    #     np.where(c.applyCuts(y,cutxpi1_15)==(ymax_cutxpi1_15))])
    # yieldmax_cutxpi2_15 = float(c.applyCuts(lumi,cutxpi2_15)[
    #     np.where(c.applyCuts(y,cutxpi2_15)==(ymax_cutxpi2_15))])
    # yieldmax_cutxpi3_15 = float(c.applyCuts(lumi,cutxpi3_15)[
    #     np.where(c.applyCuts(y,cutxpi3_15)==(ymax_cutxpi3_15))])
    # # yieldmax_cutxpi4_15 = float(c.applyCuts(lumi,cutxpi4_15)[
    # #     np.where(c.applyCuts(y,cutxpi4_15)==(ymax_cutxpi4_15))])
    # # yieldmax_cutxpi5_15 = float(c.applyCuts(lumi,cutxpi5_15)[
    # #     np.where(c.applyCuts(y,cutxpi5_15)==(ymax_cutxpi5_15))])
    # # yieldmax_cutxpi6_15 = float(c.applyCuts(lumi,cutxpi6_15)[
    # #     np.where(c.applyCuts(y,cutxpi6_15)==(ymax_cutxpi6_15))])
    # # yieldmax_cutxpi7_15 = float(c.applyCuts(lumi,cutxpi7_15)[
    # #     np.where(c.applyCuts(y,cutxpi7_15)==(ymax_cutxpi7_15))])

    # xLmin_cutxpi1_15 = float(c.applyCuts(xL,cutxpi1_15)[
    #     np.where(c.applyCuts(y,cutxpi1_15)==(ymin_cutxpi1_15))])
    # xLmin_cutxpi2_15 = float(c.applyCuts(xL,cutxpi2_15)[
    #     np.where(c.applyCuts(y,cutxpi2_15)==(ymin_cutxpi2_15))])
    # xLmin_cutxpi3_15 = float(c.applyCuts(xL,cutxpi3_15)[
    #     np.where(c.applyCuts(y,cutxpi3_15)==(ymin_cutxpi3_15))])
    # # xLmin_cutxpi4_15 = float(c.applyCuts(xL,cutxpi4_15)[
    # #     np.where(c.applyCuts(y,cutxpi4_15)==(ymin_cutxpi4_15))])
    # # xLmin_cutxpi5_15 = float(c.applyCuts(xL,cutxpi5_15)[
    # #     np.where(c.applyCuts(y,cutxpi5_15)==(ymin_cutxpi5_15))])
    # # xLmin_cutxpi6_15 = float(c.applyCuts(xL,cutxpi6_15)[
    # #     np.where(c.applyCuts(y,cutxpi6_15)==(ymin_cutxpi6_15))])
    # # xLmin_cutxpi7_15 = float(c.applyCuts(xL,cutxpi7_15)[
    # #     np.where(c.applyCuts(y,cutxpi7_15)==(ymin_cutxpi7_15))])

    # xLmax_cutxpi1_15 = float(c.applyCuts(xL,cutxpi1_15)[
    #     np.where(c.applyCuts(y,cutxpi1_15)==(ymax_cutxpi1_15))])
    # xLmax_cutxpi2_15 = float(c.applyCuts(xL,cutxpi2_15)[
    #     np.where(c.applyCuts(y,cutxpi2_15)==(ymax_cutxpi2_15))])
    # xLmax_cutxpi3_15 = float(c.applyCuts(xL,cutxpi3_15)[
    #     np.where(c.applyCuts(y,cutxpi3_15)==(ymax_cutxpi3_15))])
    # # xLmax_cutxpi4_15 = float(c.applyCuts(xL,cutxpi4_15)[
    # #     np.where(c.applyCuts(y,cutxpi4_15)==(ymax_cutxpi4_15))])
    # # xLmax_cutxpi5_15 = float(c.applyCuts(xL,cutxpi5_15)[
    # #     np.where(c.applyCuts(y,cutxpi5_15)==(ymax_cutxpi5_15))])
    # # xLmax_cutxpi6_15 = float(c.applyCuts(xL,cutxpi6_15)[
    # #     np.where(c.applyCuts(y,cutxpi6_15)==(ymax_cutxpi6_15))])
    # # xLmax_cutxpi7_15 = float(c.applyCuts(xL,cutxpi7_15)[
    # #     np.where(c.applyCuts(y,cutxpi7_15)==(ymax_cutxpi7_15))])

    # xpimin_cutxpi1_15 = float(c.applyCuts(xpi,cutxpi1_15)[
    #     np.where(c.applyCuts(y,cutxpi1_15)==(ymin_cutxpi1_15))])
    # xpimin_cutxpi2_15 = float(c.applyCuts(xpi,cutxpi2_15)[
    #     np.where(c.applyCuts(y,cutxpi2_15)==(ymin_cutxpi2_15))])
    # xpimin_cutxpi3_15 = float(c.applyCuts(xpi,cutxpi3_15)[
    #     np.where(c.applyCuts(y,cutxpi3_15)==(ymin_cutxpi3_15))])
    # # xpimin_cutxpi4_15 = float(c.applyCuts(xpi,cutxpi4_15)[
    # #     np.where(c.applyCuts(y,cutxpi4_15)==(ymin_cutxpi4_15))])
    # # xpimin_cutxpi5_15 = float(c.applyCuts(xpi,cutxpi5_15)[
    # #     np.where(c.applyCuts(y,cutxpi5_15)==(ymin_cutxpi5_15))])
    # # xpimin_cutxpi6_15 = float(c.applyCuts(xpi,cutxpi6_15)[
    # #     np.where(c.applyCuts(y,cutxpi6_15)==(ymin_cutxpi6_15))])
    # # xpimin_cutxpi7_15 = float(c.applyCuts(xpi,cutxpi7_15)[
    # #     np.where(c.applyCuts(y,cutxpi7_15)==(ymin_cutxpi7_15))])

    # xpimax_cutxpi1_15 = float(c.applyCuts(xpi,cutxpi1_15)[
    #     np.where(c.applyCuts(y,cutxpi1_15)==(ymax_cutxpi1_15))])
    # xpimax_cutxpi2_15 = float(c.applyCuts(xpi,cutxpi2_15)[
    #     np.where(c.applyCuts(y,cutxpi2_15)==(ymax_cutxpi2_15))])
    # xpimax_cutxpi3_15 = float(c.applyCuts(xpi,cutxpi3_15)[
    #     np.where(c.applyCuts(y,cutxpi3_15)==(ymax_cutxpi3_15))])
    # # xpimax_cutxpi4_15 = float(c.applyCuts(xpi,cutxpi4_15)[
    # #     np.where(c.applyCuts(y,cutxpi4_15)==(ymax_cutxpi4_15))])
    # # xpimax_cutxpi5_15 = float(c.applyCuts(xpi,cutxpi5_15)[
    # #     np.where(c.applyCuts(y,cutxpi5_15)==(ymax_cutxpi5_15))])
    # # xpimax_cutxpi6_15 = float(c.applyCuts(xpi,cutxpi6_15)[
    # #     np.where(c.applyCuts(y,cutxpi6_15)==(ymax_cutxpi6_15))])
    # # xpimax_cutxpi7_15 = float(c.applyCuts(xpi,cutxpi7_15)[
    # #     np.where(c.applyCuts(y,cutxpi7_15)==(ymax_cutxpi7_15))])

    # Q2min_cutxpi1_15 = float(c.applyCuts(Q2,cutxpi1_15)[
    #     np.where(c.applyCuts(y,cutxpi1_15)==(ymin_cutxpi1_15))])
    # Q2min_cutxpi2_15 = float(c.applyCuts(Q2,cutxpi2_15)[
    #     np.where(c.applyCuts(y,cutxpi2_15)==(ymin_cutxpi2_15))])
    # Q2min_cutxpi3_15 = float(c.applyCuts(Q2,cutxpi3_15)[
    #     np.where(c.applyCuts(y,cutxpi3_15)==(ymin_cutxpi3_15))])
    # # Q2min_cutxpi4_15 = float(c.applyCuts(Q2,cutxpi4_15)[
    # #     np.where(c.applyCuts(y,cutxpi4_15)==(ymin_cutxpi4_15))])
    # # Q2min_cutxpi5_15 = float(c.applyCuts(Q2,cutxpi5_15)[
    # #     np.where(c.applyCuts(y,cutxpi5_15)==(ymin_cutxpi5_15))])
    # # Q2min_cutxpi6_15 = float(c.applyCuts(Q2,cutxpi6_15)[
    # #     np.where(c.applyCuts(y,cutxpi6_15)==(ymin_cutxpi6_15))])
    # # Q2min_cutxpi7_15 = float(c.applyCuts(Q2,cutxpi7_15)[
    # #     np.where(c.applyCuts(y,cutxpi7_15)==(ymin_cutxpi7_15))])

    # Q2max_cutxpi1_15 = float(c.applyCuts(Q2,cutxpi1_15)[
    #     np.where(c.applyCuts(y,cutxpi1_15)==(ymax_cutxpi1_15))])
    # Q2max_cutxpi2_15 = float(c.applyCuts(Q2,cutxpi2_15)[
    #     np.where(c.applyCuts(y,cutxpi2_15)==(ymax_cutxpi2_15))])
    # Q2max_cutxpi3_15 = float(c.applyCuts(Q2,cutxpi3_15)[
    #     np.where(c.applyCuts(y,cutxpi3_15)==(ymax_cutxpi3_15))])
    # # Q2max_cutxpi4_15 = float(c.applyCuts(Q2,cutxpi4_15)[
    # #     np.where(c.applyCuts(y,cutxpi4_15)==(ymax_cutxpi4_15))])
    # # Q2max_cutxpi5_15 = float(c.applyCuts(Q2,cutxpi5_15)[
    # #     np.where(c.applyCuts(y,cutxpi5_15)==(ymax_cutxpi5_15))])
    # # Q2max_cutxpi6_15 = float(c.applyCuts(Q2,cutxpi6_15)[
    # #     np.where(c.applyCuts(y,cutxpi6_15)==(ymax_cutxpi6_15))])
    # # Q2max_cutxpi7_15 = float(c.applyCuts(Q2,cutxpi7_15)[
    # #     np.where(c.applyCuts(y,cutxpi7_15)==(ymax_cutxpi7_15))])
    
    # #################################

    # #################################
        
    
    # ymin_cutxpi1_30 = min(c.applyCuts(y,cutxpi1_30))
    # ymin_cutxpi2_30 = min(c.applyCuts(y,cutxpi2_30))
    # ymin_cutxpi3_30 = min(c.applyCuts(y,cutxpi3_30))
    # ymin_cutxpi4_30 = min(c.applyCuts(y,cutxpi4_30))
    # ymin_cutxpi5_30 = min(c.applyCuts(y,cutxpi5_30))
    # ymin_cutxpi6_30 = min(c.applyCuts(y,cutxpi6_30))
    # # ymin_cutxpi7_30 = min(c.applyCuts(y,cutxpi7_30))
    
    # ymax_cutxpi1_30 = max(c.applyCuts(y,cutxpi1_30))
    # ymax_cutxpi2_30 = max(c.applyCuts(y,cutxpi2_30))
    # ymax_cutxpi3_30 = max(c.applyCuts(y,cutxpi3_30))
    # ymax_cutxpi4_30 = max(c.applyCuts(y,cutxpi4_30))
    # ymax_cutxpi5_30 = max(c.applyCuts(y,cutxpi5_30))
    # ymax_cutxpi6_30 = max(c.applyCuts(y,cutxpi6_30))
    # # ymax_cutxpi7_30 = max(c.applyCuts(y,cutxpi7_30))

    # yieldmin_cutxpi1_30 = float(c.applyCuts(lumi,cutxpi1_30)[
    #     np.where(c.applyCuts(y,cutxpi1_30)==(ymin_cutxpi1_30))])
    # yieldmin_cutxpi2_30 = float(c.applyCuts(lumi,cutxpi2_30)[
    #     np.where(c.applyCuts(y,cutxpi2_30)==(ymin_cutxpi2_30))])
    # yieldmin_cutxpi3_30 = float(c.applyCuts(lumi,cutxpi3_30)[
    #     np.where(c.applyCuts(y,cutxpi3_30)==(ymin_cutxpi3_30))])
    # yieldmin_cutxpi4_30 = float(c.applyCuts(lumi,cutxpi4_30)[
    #     np.where(c.applyCuts(y,cutxpi4_30)==(ymin_cutxpi4_30))])
    # yieldmin_cutxpi5_30 = float(c.applyCuts(lumi,cutxpi5_30)[
    #     np.where(c.applyCuts(y,cutxpi5_30)==(ymin_cutxpi5_30))])
    # yieldmin_cutxpi6_30 = float(c.applyCuts(lumi,cutxpi6_30)[
    #     np.where(c.applyCuts(y,cutxpi6_30)==(ymin_cutxpi6_30))])
    # # yieldmin_cutxpi7_30 = float(c.applyCuts(lumi,cutxpi7_30)[
    # #     np.where(c.applyCuts(y,cutxpi7_30)==(ymin_cutxpi7_30))])
    
    # yieldmax_cutxpi1_30 = float(c.applyCuts(lumi,cutxpi1_30)[
    #     np.where(c.applyCuts(y,cutxpi1_30)==(ymax_cutxpi1_30))])
    # yieldmax_cutxpi2_30 = float(c.applyCuts(lumi,cutxpi2_30)[
    #     np.where(c.applyCuts(y,cutxpi2_30)==(ymax_cutxpi2_30))])
    # yieldmax_cutxpi3_30 = float(c.applyCuts(lumi,cutxpi3_30)[
    #     np.where(c.applyCuts(y,cutxpi3_30)==(ymax_cutxpi3_30))])
    # yieldmax_cutxpi4_30 = float(c.applyCuts(lumi,cutxpi4_30)[
    #     np.where(c.applyCuts(y,cutxpi4_30)==(ymax_cutxpi4_30))])
    # yieldmax_cutxpi5_30 = float(c.applyCuts(lumi,cutxpi5_30)[
    #     np.where(c.applyCuts(y,cutxpi5_30)==(ymax_cutxpi5_30))])
    # yieldmax_cutxpi6_30 = float(c.applyCuts(lumi,cutxpi6_30)[
    #     np.where(c.applyCuts(y,cutxpi6_30)==(ymax_cutxpi6_30))])
    # # yieldmax_cutxpi7_30 = float(c.applyCuts(lumi,cutxpi7_30)[
    # #     np.where(c.applyCuts(y,cutxpi7_30)==(ymax_cutxpi7_30))])

    # xLmin_cutxpi1_30 = float(c.applyCuts(xL,cutxpi1_30)[
    #     np.where(c.applyCuts(y,cutxpi1_30)==(ymin_cutxpi1_30))])
    # xLmin_cutxpi2_30 = float(c.applyCuts(xL,cutxpi2_30)[
    #     np.where(c.applyCuts(y,cutxpi2_30)==(ymin_cutxpi2_30))])
    # xLmin_cutxpi3_30 = float(c.applyCuts(xL,cutxpi3_30)[
    #     np.where(c.applyCuts(y,cutxpi3_30)==(ymin_cutxpi3_30))])
    # xLmin_cutxpi4_30 = float(c.applyCuts(xL,cutxpi4_30)[
    #     np.where(c.applyCuts(y,cutxpi4_30)==(ymin_cutxpi4_30))])
    # xLmin_cutxpi5_30 = float(c.applyCuts(xL,cutxpi5_30)[
    #     np.where(c.applyCuts(y,cutxpi5_30)==(ymin_cutxpi5_30))])
    # xLmin_cutxpi6_30 = float(c.applyCuts(xL,cutxpi6_30)[
    #     np.where(c.applyCuts(y,cutxpi6_30)==(ymin_cutxpi6_30))])
    # # xLmin_cutxpi7_30 = float(c.applyCuts(xL,cutxpi7_30)[
    # #     np.where(c.applyCuts(y,cutxpi7_30)==(ymin_cutxpi7_30))])

    # xLmax_cutxpi1_30 = float(c.applyCuts(xL,cutxpi1_30)[
    #     np.where(c.applyCuts(y,cutxpi1_30)==(ymax_cutxpi1_30))])
    # xLmax_cutxpi2_30 = float(c.applyCuts(xL,cutxpi2_30)[
    #     np.where(c.applyCuts(y,cutxpi2_30)==(ymax_cutxpi2_30))])
    # xLmax_cutxpi3_30 = float(c.applyCuts(xL,cutxpi3_30)[
    #     np.where(c.applyCuts(y,cutxpi3_30)==(ymax_cutxpi3_30))])
    # xLmax_cutxpi4_30 = float(c.applyCuts(xL,cutxpi4_30)[
    #     np.where(c.applyCuts(y,cutxpi4_30)==(ymax_cutxpi4_30))])
    # xLmax_cutxpi5_30 = float(c.applyCuts(xL,cutxpi5_30)[
    #     np.where(c.applyCuts(y,cutxpi5_30)==(ymax_cutxpi5_30))])
    # xLmax_cutxpi6_30 = float(c.applyCuts(xL,cutxpi6_30)[
    #     np.where(c.applyCuts(y,cutxpi6_30)==(ymax_cutxpi6_30))])
    # # xLmax_cutxpi7_30 = float(c.applyCuts(xL,cutxpi7_30)[
    # #     np.where(c.applyCuts(y,cutxpi7_30)==(ymax_cutxpi7_30))])

    # xpimin_cutxpi1_30 = float(c.applyCuts(xpi,cutxpi1_30)[
    #     np.where(c.applyCuts(y,cutxpi1_30)==(ymin_cutxpi1_30))])
    # xpimin_cutxpi2_30 = float(c.applyCuts(xpi,cutxpi2_30)[
    #     np.where(c.applyCuts(y,cutxpi2_30)==(ymin_cutxpi2_30))])
    # xpimin_cutxpi3_30 = float(c.applyCuts(xpi,cutxpi3_30)[
    #     np.where(c.applyCuts(y,cutxpi3_30)==(ymin_cutxpi3_30))])
    # xpimin_cutxpi4_30 = float(c.applyCuts(xpi,cutxpi4_30)[
    #     np.where(c.applyCuts(y,cutxpi4_30)==(ymin_cutxpi4_30))])
    # xpimin_cutxpi5_30 = float(c.applyCuts(xpi,cutxpi5_30)[
    #     np.where(c.applyCuts(y,cutxpi5_30)==(ymin_cutxpi5_30))])
    # xpimin_cutxpi6_30 = float(c.applyCuts(xpi,cutxpi6_30)[
    #     np.where(c.applyCuts(y,cutxpi6_30)==(ymin_cutxpi6_30))])
    # # xpimin_cutxpi7_30 = float(c.applyCuts(xpi,cutxpi7_30)[
    # #     np.where(c.applyCuts(y,cutxpi7_30)==(ymin_cutxpi7_30))])

    # xpimax_cutxpi1_30 = float(c.applyCuts(xpi,cutxpi1_30)[
    #     np.where(c.applyCuts(y,cutxpi1_30)==(ymax_cutxpi1_30))])
    # xpimax_cutxpi2_30 = float(c.applyCuts(xpi,cutxpi2_30)[
    #     np.where(c.applyCuts(y,cutxpi2_30)==(ymax_cutxpi2_30))])
    # xpimax_cutxpi3_30 = float(c.applyCuts(xpi,cutxpi3_30)[
    #     np.where(c.applyCuts(y,cutxpi3_30)==(ymax_cutxpi3_30))])
    # xpimax_cutxpi4_30 = float(c.applyCuts(xpi,cutxpi4_30)[
    #     np.where(c.applyCuts(y,cutxpi4_30)==(ymax_cutxpi4_30))])
    # xpimax_cutxpi5_30 = float(c.applyCuts(xpi,cutxpi5_30)[
    #     np.where(c.applyCuts(y,cutxpi5_30)==(ymax_cutxpi5_30))])
    # xpimax_cutxpi6_30 = float(c.applyCuts(xpi,cutxpi6_30)[
    #     np.where(c.applyCuts(y,cutxpi6_30)==(ymax_cutxpi6_30))])
    # # xpimax_cutxpi7_30 = float(c.applyCuts(xpi,cutxpi7_30)[
    # #     np.where(c.applyCuts(y,cutxpi7_30)==(ymax_cutxpi7_30))])

    # Q2min_cutxpi1_30 = float(c.applyCuts(Q2,cutxpi1_30)[
    #     np.where(c.applyCuts(y,cutxpi1_30)==(ymin_cutxpi1_30))])
    # Q2min_cutxpi2_30 = float(c.applyCuts(Q2,cutxpi2_30)[
    #     np.where(c.applyCuts(y,cutxpi2_30)==(ymin_cutxpi2_30))])
    # Q2min_cutxpi3_30 = float(c.applyCuts(Q2,cutxpi3_30)[
    #     np.where(c.applyCuts(y,cutxpi3_30)==(ymin_cutxpi3_30))])
    # Q2min_cutxpi4_30 = float(c.applyCuts(Q2,cutxpi4_30)[
    #     np.where(c.applyCuts(y,cutxpi4_30)==(ymin_cutxpi4_30))])
    # Q2min_cutxpi5_30 = float(c.applyCuts(Q2,cutxpi5_30)[
    #     np.where(c.applyCuts(y,cutxpi5_30)==(ymin_cutxpi5_30))])
    # Q2min_cutxpi6_30 = float(c.applyCuts(Q2,cutxpi6_30)[
    #     np.where(c.applyCuts(y,cutxpi6_30)==(ymin_cutxpi6_30))])
    # # Q2min_cutxpi7_30 = float(c.applyCuts(Q2,cutxpi7_30)[
    # #     np.where(c.applyCuts(y,cutxpi7_30)==(ymin_cutxpi7_30))])

    # Q2max_cutxpi1_30 = float(c.applyCuts(Q2,cutxpi1_30)[
    #     np.where(c.applyCuts(y,cutxpi1_30)==(ymax_cutxpi1_30))])
    # Q2max_cutxpi2_30 = float(c.applyCuts(Q2,cutxpi2_30)[
    #     np.where(c.applyCuts(y,cutxpi2_30)==(ymax_cutxpi2_30))])
    # Q2max_cutxpi3_30 = float(c.applyCuts(Q2,cutxpi3_30)[
    #     np.where(c.applyCuts(y,cutxpi3_30)==(ymax_cutxpi3_30))])
    # Q2max_cutxpi4_30 = float(c.applyCuts(Q2,cutxpi4_30)[
    #     np.where(c.applyCuts(y,cutxpi4_30)==(ymax_cutxpi4_30))])
    # Q2max_cutxpi5_30 = float(c.applyCuts(Q2,cutxpi5_30)[
    #     np.where(c.applyCuts(y,cutxpi5_30)==(ymax_cutxpi5_30))])
    # Q2max_cutxpi6_30 = float(c.applyCuts(Q2,cutxpi6_30)[
    #     np.where(c.applyCuts(y,cutxpi6_30)==(ymax_cutxpi6_30))])
    # # Q2max_cutxpi7_30 = float(c.applyCuts(Q2,cutxpi7_30)[
    # #     np.where(c.applyCuts(y,cutxpi7_30)==(ymax_cutxpi7_30))])
    
    # #################################

    # #################################
        
    
    # ymin_cutxpi1_60 = min(c.applyCuts(y,cutxpi1_60))
    # ymin_cutxpi2_60 = min(c.applyCuts(y,cutxpi2_60))
    # ymin_cutxpi3_60 = min(c.applyCuts(y,cutxpi3_60))
    # ymin_cutxpi4_60 = min(c.applyCuts(y,cutxpi4_60))
    # ymin_cutxpi5_60 = min(c.applyCuts(y,cutxpi5_60))
    # ymin_cutxpi6_60 = min(c.applyCuts(y,cutxpi6_60))
    # ymin_cutxpi7_60 = min(c.applyCuts(y,cutxpi7_60))
    
    # ymax_cutxpi1_60 = max(c.applyCuts(y,cutxpi1_60))
    # ymax_cutxpi2_60 = max(c.applyCuts(y,cutxpi2_60))
    # ymax_cutxpi3_60 = max(c.applyCuts(y,cutxpi3_60))
    # ymax_cutxpi4_60 = max(c.applyCuts(y,cutxpi4_60))
    # ymax_cutxpi5_60 = max(c.applyCuts(y,cutxpi5_60))
    # ymax_cutxpi6_60 = max(c.applyCuts(y,cutxpi6_60))
    # ymax_cutxpi7_60 = max(c.applyCuts(y,cutxpi7_60))

    # yieldmin_cutxpi1_60 = float(c.applyCuts(lumi,cutxpi1_60)[
    #     np.where(c.applyCuts(y,cutxpi1_60)==(ymin_cutxpi1_60))])
    # yieldmin_cutxpi2_60 = float(c.applyCuts(lumi,cutxpi2_60)[
    #     np.where(c.applyCuts(y,cutxpi2_60)==(ymin_cutxpi2_60))])
    # yieldmin_cutxpi3_60 = float(c.applyCuts(lumi,cutxpi3_60)[
    #     np.where(c.applyCuts(y,cutxpi3_60)==(ymin_cutxpi3_60))])
    # yieldmin_cutxpi4_60 = float(c.applyCuts(lumi,cutxpi4_60)[
    #     np.where(c.applyCuts(y,cutxpi4_60)==(ymin_cutxpi4_60))])
    # yieldmin_cutxpi5_60 = float(c.applyCuts(lumi,cutxpi5_60)[
    #     np.where(c.applyCuts(y,cutxpi5_60)==(ymin_cutxpi5_60))])
    # yieldmin_cutxpi6_60 = float(c.applyCuts(lumi,cutxpi6_60)[
    #     np.where(c.applyCuts(y,cutxpi6_60)==(ymin_cutxpi6_60))])
    # yieldmin_cutxpi7_60 = float(c.applyCuts(lumi,cutxpi7_60)[
    #     np.where(c.applyCuts(y,cutxpi7_60)==(ymin_cutxpi7_60))])
    
    # yieldmax_cutxpi1_60 = float(c.applyCuts(lumi,cutxpi1_60)[
    #     np.where(c.applyCuts(y,cutxpi1_60)==(ymax_cutxpi1_60))])
    # yieldmax_cutxpi2_60 = float(c.applyCuts(lumi,cutxpi2_60)[
    #     np.where(c.applyCuts(y,cutxpi2_60)==(ymax_cutxpi2_60))])
    # yieldmax_cutxpi3_60 = float(c.applyCuts(lumi,cutxpi3_60)[
    #     np.where(c.applyCuts(y,cutxpi3_60)==(ymax_cutxpi3_60))])
    # yieldmax_cutxpi4_60 = float(c.applyCuts(lumi,cutxpi4_60)[
    #     np.where(c.applyCuts(y,cutxpi4_60)==(ymax_cutxpi4_60))])
    # yieldmax_cutxpi5_60 = float(c.applyCuts(lumi,cutxpi5_60)[
    #     np.where(c.applyCuts(y,cutxpi5_60)==(ymax_cutxpi5_60))])
    # yieldmax_cutxpi6_60 = float(c.applyCuts(lumi,cutxpi6_60)[
    #     np.where(c.applyCuts(y,cutxpi6_60)==(ymax_cutxpi6_60))])
    # yieldmax_cutxpi7_60 = float(c.applyCuts(lumi,cutxpi7_60)[
    #     np.where(c.applyCuts(y,cutxpi7_60)==(ymax_cutxpi7_60))])

    # xLmin_cutxpi1_60 = float(c.applyCuts(xL,cutxpi1_60)[
    #     np.where(c.applyCuts(y,cutxpi1_60)==(ymin_cutxpi1_60))])
    # xLmin_cutxpi2_60 = float(c.applyCuts(xL,cutxpi2_60)[
    #     np.where(c.applyCuts(y,cutxpi2_60)==(ymin_cutxpi2_60))])
    # xLmin_cutxpi3_60 = float(c.applyCuts(xL,cutxpi3_60)[
    #     np.where(c.applyCuts(y,cutxpi3_60)==(ymin_cutxpi3_60))])
    # xLmin_cutxpi4_60 = float(c.applyCuts(xL,cutxpi4_60)[
    #     np.where(c.applyCuts(y,cutxpi4_60)==(ymin_cutxpi4_60))])
    # xLmin_cutxpi5_60 = float(c.applyCuts(xL,cutxpi5_60)[
    #     np.where(c.applyCuts(y,cutxpi5_60)==(ymin_cutxpi5_60))])
    # xLmin_cutxpi6_60 = float(c.applyCuts(xL,cutxpi6_60)[
    #     np.where(c.applyCuts(y,cutxpi6_60)==(ymin_cutxpi6_60))])
    # xLmin_cutxpi7_60 = float(c.applyCuts(xL,cutxpi7_60)[
    #     np.where(c.applyCuts(y,cutxpi7_60)==(ymin_cutxpi7_60))])

    # xLmax_cutxpi1_60 = float(c.applyCuts(xL,cutxpi1_60)[
    #     np.where(c.applyCuts(y,cutxpi1_60)==(ymax_cutxpi1_60))])
    # xLmax_cutxpi2_60 = float(c.applyCuts(xL,cutxpi2_60)[
    #     np.where(c.applyCuts(y,cutxpi2_60)==(ymax_cutxpi2_60))])
    # xLmax_cutxpi3_60 = float(c.applyCuts(xL,cutxpi3_60)[
    #     np.where(c.applyCuts(y,cutxpi3_60)==(ymax_cutxpi3_60))])
    # xLmax_cutxpi4_60 = float(c.applyCuts(xL,cutxpi4_60)[
    #     np.where(c.applyCuts(y,cutxpi4_60)==(ymax_cutxpi4_60))])
    # xLmax_cutxpi5_60 = float(c.applyCuts(xL,cutxpi5_60)[
    #     np.where(c.applyCuts(y,cutxpi5_60)==(ymax_cutxpi5_60))])
    # xLmax_cutxpi6_60 = float(c.applyCuts(xL,cutxpi6_60)[
    #     np.where(c.applyCuts(y,cutxpi6_60)==(ymax_cutxpi6_60))])
    # xLmax_cutxpi7_60 = float(c.applyCuts(xL,cutxpi7_60)[
    #     np.where(c.applyCuts(y,cutxpi7_60)==(ymax_cutxpi7_60))])

    # xpimin_cutxpi1_60 = float(c.applyCuts(xpi,cutxpi1_60)[
    #     np.where(c.applyCuts(y,cutxpi1_60)==(ymin_cutxpi1_60))])
    # xpimin_cutxpi2_60 = float(c.applyCuts(xpi,cutxpi2_60)[
    #     np.where(c.applyCuts(y,cutxpi2_60)==(ymin_cutxpi2_60))])
    # xpimin_cutxpi3_60 = float(c.applyCuts(xpi,cutxpi3_60)[
    #     np.where(c.applyCuts(y,cutxpi3_60)==(ymin_cutxpi3_60))])
    # xpimin_cutxpi4_60 = float(c.applyCuts(xpi,cutxpi4_60)[
    #     np.where(c.applyCuts(y,cutxpi4_60)==(ymin_cutxpi4_60))])
    # xpimin_cutxpi5_60 = float(c.applyCuts(xpi,cutxpi5_60)[
    #     np.where(c.applyCuts(y,cutxpi5_60)==(ymin_cutxpi5_60))])
    # xpimin_cutxpi6_60 = float(c.applyCuts(xpi,cutxpi6_60)[
    #     np.where(c.applyCuts(y,cutxpi6_60)==(ymin_cutxpi6_60))])
    # xpimin_cutxpi7_60 = float(c.applyCuts(xpi,cutxpi7_60)[
    #     np.where(c.applyCuts(y,cutxpi7_60)==(ymin_cutxpi7_60))])

    # xpimax_cutxpi1_60 = float(c.applyCuts(xpi,cutxpi1_60)[
    #     np.where(c.applyCuts(y,cutxpi1_60)==(ymax_cutxpi1_60))])
    # xpimax_cutxpi2_60 = float(c.applyCuts(xpi,cutxpi2_60)[
    #     np.where(c.applyCuts(y,cutxpi2_60)==(ymax_cutxpi2_60))])
    # xpimax_cutxpi3_60 = float(c.applyCuts(xpi,cutxpi3_60)[
    #     np.where(c.applyCuts(y,cutxpi3_60)==(ymax_cutxpi3_60))])
    # xpimax_cutxpi4_60 = float(c.applyCuts(xpi,cutxpi4_60)[
    #     np.where(c.applyCuts(y,cutxpi4_60)==(ymax_cutxpi4_60))])
    # xpimax_cutxpi5_60 = float(c.applyCuts(xpi,cutxpi5_60)[
    #     np.where(c.applyCuts(y,cutxpi5_60)==(ymax_cutxpi5_60))])
    # xpimax_cutxpi6_60 = float(c.applyCuts(xpi,cutxpi6_60)[
    #     np.where(c.applyCuts(y,cutxpi6_60)==(ymax_cutxpi6_60))])
    # xpimax_cutxpi7_60 = float(c.applyCuts(xpi,cutxpi7_60)[
    #     np.where(c.applyCuts(y,cutxpi7_60)==(ymax_cutxpi7_60))])

    # Q2min_cutxpi1_60 = float(c.applyCuts(Q2,cutxpi1_60)[
    #     np.where(c.applyCuts(y,cutxpi1_60)==(ymin_cutxpi1_60))])
    # Q2min_cutxpi2_60 = float(c.applyCuts(Q2,cutxpi2_60)[
    #     np.where(c.applyCuts(y,cutxpi2_60)==(ymin_cutxpi2_60))])
    # Q2min_cutxpi3_60 = float(c.applyCuts(Q2,cutxpi3_60)[
    #     np.where(c.applyCuts(y,cutxpi3_60)==(ymin_cutxpi3_60))])
    # Q2min_cutxpi4_60 = float(c.applyCuts(Q2,cutxpi4_60)[
    #     np.where(c.applyCuts(y,cutxpi4_60)==(ymin_cutxpi4_60))])
    # Q2min_cutxpi5_60 = float(c.applyCuts(Q2,cutxpi5_60)[
    #     np.where(c.applyCuts(y,cutxpi5_60)==(ymin_cutxpi5_60))])
    # Q2min_cutxpi6_60 = float(c.applyCuts(Q2,cutxpi6_60)[
    #     np.where(c.applyCuts(y,cutxpi6_60)==(ymin_cutxpi6_60))])
    # Q2min_cutxpi7_60 = float(c.applyCuts(Q2,cutxpi7_60)[
    #     np.where(c.applyCuts(y,cutxpi7_60)==(ymin_cutxpi7_60))])

    # Q2max_cutxpi1_60 = float(c.applyCuts(Q2,cutxpi1_60)[
    #     np.where(c.applyCuts(y,cutxpi1_60)==(ymax_cutxpi1_60))])
    # Q2max_cutxpi2_60 = float(c.applyCuts(Q2,cutxpi2_60)[
    #     np.where(c.applyCuts(y,cutxpi2_60)==(ymax_cutxpi2_60))])
    # Q2max_cutxpi3_60 = float(c.applyCuts(Q2,cutxpi3_60)[
    #     np.where(c.applyCuts(y,cutxpi3_60)==(ymax_cutxpi3_60))])
    # Q2max_cutxpi4_60 = float(c.applyCuts(Q2,cutxpi4_60)[
    #     np.where(c.applyCuts(y,cutxpi4_60)==(ymax_cutxpi4_60))])
    # Q2max_cutxpi5_60 = float(c.applyCuts(Q2,cutxpi5_60)[
    #     np.where(c.applyCuts(y,cutxpi5_60)==(ymax_cutxpi5_60))])
    # Q2max_cutxpi6_60 = float(c.applyCuts(Q2,cutxpi6_60)[
    #     np.where(c.applyCuts(y,cutxpi6_60)==(ymax_cutxpi6_60))])
    # Q2max_cutxpi7_60 = float(c.applyCuts(Q2,cutxpi7_60)[
    #     np.where(c.applyCuts(y,cutxpi7_60)==(ymax_cutxpi7_60))])
    
    # #################################

    # #################################
        
    
    # ymin_cutxpi1_120 = min(c.applyCuts(y,cutxpi1_120))
    # ymin_cutxpi2_120 = min(c.applyCuts(y,cutxpi2_120))
    # ymin_cutxpi3_120 = min(c.applyCuts(y,cutxpi3_120))
    # ymin_cutxpi4_120 = min(c.applyCuts(y,cutxpi4_120))
    # ymin_cutxpi5_120 = min(c.applyCuts(y,cutxpi5_120))
    # ymin_cutxpi6_120 = min(c.applyCuts(y,cutxpi6_120))
    # ymin_cutxpi7_120 = min(c.applyCuts(y,cutxpi7_120))
    
    # ymax_cutxpi1_120 = max(c.applyCuts(y,cutxpi1_120))
    # ymax_cutxpi2_120 = max(c.applyCuts(y,cutxpi2_120))
    # ymax_cutxpi3_120 = max(c.applyCuts(y,cutxpi3_120))
    # ymax_cutxpi4_120 = max(c.applyCuts(y,cutxpi4_120))
    # ymax_cutxpi5_120 = max(c.applyCuts(y,cutxpi5_120))
    # ymax_cutxpi6_120 = max(c.applyCuts(y,cutxpi6_120))
    # ymax_cutxpi7_120 = max(c.applyCuts(y,cutxpi7_120))

    # yieldmin_cutxpi1_120 = float(c.applyCuts(lumi,cutxpi1_120)[
    #     np.where(c.applyCuts(y,cutxpi1_120)==(ymin_cutxpi1_120))])
    # yieldmin_cutxpi2_120 = float(c.applyCuts(lumi,cutxpi2_120)[
    #     np.where(c.applyCuts(y,cutxpi2_120)==(ymin_cutxpi2_120))])
    # yieldmin_cutxpi3_120 = float(c.applyCuts(lumi,cutxpi3_120)[
    #     np.where(c.applyCuts(y,cutxpi3_120)==(ymin_cutxpi3_120))])
    # yieldmin_cutxpi4_120 = float(c.applyCuts(lumi,cutxpi4_120)[
    #     np.where(c.applyCuts(y,cutxpi4_120)==(ymin_cutxpi4_120))])
    # yieldmin_cutxpi5_120 = float(c.applyCuts(lumi,cutxpi5_120)[
    #     np.where(c.applyCuts(y,cutxpi5_120)==(ymin_cutxpi5_120))])
    # yieldmin_cutxpi6_120 = float(c.applyCuts(lumi,cutxpi6_120)[
    #     np.where(c.applyCuts(y,cutxpi6_120)==(ymin_cutxpi6_120))])
    # yieldmin_cutxpi7_120 = float(c.applyCuts(lumi,cutxpi7_120)[
    #     np.where(c.applyCuts(y,cutxpi7_120)==(ymin_cutxpi7_120))])
    
    # yieldmax_cutxpi1_120 = float(c.applyCuts(lumi,cutxpi1_120)[
    #     np.where(c.applyCuts(y,cutxpi1_120)==(ymax_cutxpi1_120))])
    # yieldmax_cutxpi2_120 = float(c.applyCuts(lumi,cutxpi2_120)[
    #     np.where(c.applyCuts(y,cutxpi2_120)==(ymax_cutxpi2_120))])
    # yieldmax_cutxpi3_120 = float(c.applyCuts(lumi,cutxpi3_120)[
    #     np.where(c.applyCuts(y,cutxpi3_120)==(ymax_cutxpi3_120))])
    # yieldmax_cutxpi4_120 = float(c.applyCuts(lumi,cutxpi4_120)[
    #     np.where(c.applyCuts(y,cutxpi4_120)==(ymax_cutxpi4_120))])
    # yieldmax_cutxpi5_120 = float(c.applyCuts(lumi,cutxpi5_120)[
    #     np.where(c.applyCuts(y,cutxpi5_120)==(ymax_cutxpi5_120))])
    # yieldmax_cutxpi6_120 = float(c.applyCuts(lumi,cutxpi6_120)[
    #     np.where(c.applyCuts(y,cutxpi6_120)==(ymax_cutxpi6_120))])
    # yieldmax_cutxpi7_120 = float(c.applyCuts(lumi,cutxpi7_120)[
    #     np.where(c.applyCuts(y,cutxpi7_120)==(ymax_cutxpi7_120))])

    # xLmin_cutxpi1_120 = float(c.applyCuts(xL,cutxpi1_120)[
    #     np.where(c.applyCuts(y,cutxpi1_120)==(ymin_cutxpi1_120))])
    # xLmin_cutxpi2_120 = float(c.applyCuts(xL,cutxpi2_120)[
    #     np.where(c.applyCuts(y,cutxpi2_120)==(ymin_cutxpi2_120))])
    # xLmin_cutxpi3_120 = float(c.applyCuts(xL,cutxpi3_120)[
    #     np.where(c.applyCuts(y,cutxpi3_120)==(ymin_cutxpi3_120))])
    # xLmin_cutxpi4_120 = float(c.applyCuts(xL,cutxpi4_120)[
    #     np.where(c.applyCuts(y,cutxpi4_120)==(ymin_cutxpi4_120))])
    # xLmin_cutxpi5_120 = float(c.applyCuts(xL,cutxpi5_120)[
    #     np.where(c.applyCuts(y,cutxpi5_120)==(ymin_cutxpi5_120))])
    # xLmin_cutxpi6_120 = float(c.applyCuts(xL,cutxpi6_120)[
    #     np.where(c.applyCuts(y,cutxpi6_120)==(ymin_cutxpi6_120))])
    # xLmin_cutxpi7_120 = float(c.applyCuts(xL,cutxpi7_120)[
    #     np.where(c.applyCuts(y,cutxpi7_120)==(ymin_cutxpi7_120))])

    # xLmax_cutxpi1_120 = float(c.applyCuts(xL,cutxpi1_120)[
    #     np.where(c.applyCuts(y,cutxpi1_120)==(ymax_cutxpi1_120))])
    # xLmax_cutxpi2_120 = float(c.applyCuts(xL,cutxpi2_120)[
    #     np.where(c.applyCuts(y,cutxpi2_120)==(ymax_cutxpi2_120))])
    # xLmax_cutxpi3_120 = float(c.applyCuts(xL,cutxpi3_120)[
    #     np.where(c.applyCuts(y,cutxpi3_120)==(ymax_cutxpi3_120))])
    # xLmax_cutxpi4_120 = float(c.applyCuts(xL,cutxpi4_120)[
    #     np.where(c.applyCuts(y,cutxpi4_120)==(ymax_cutxpi4_120))])
    # xLmax_cutxpi5_120 = float(c.applyCuts(xL,cutxpi5_120)[
    #     np.where(c.applyCuts(y,cutxpi5_120)==(ymax_cutxpi5_120))])
    # xLmax_cutxpi6_120 = float(c.applyCuts(xL,cutxpi6_120)[
    #     np.where(c.applyCuts(y,cutxpi6_120)==(ymax_cutxpi6_120))])
    # xLmax_cutxpi7_120 = float(c.applyCuts(xL,cutxpi7_120)[
    #     np.where(c.applyCuts(y,cutxpi7_120)==(ymax_cutxpi7_120))])

    # xpimin_cutxpi1_120 = float(c.applyCuts(xpi,cutxpi1_120)[
    #     np.where(c.applyCuts(y,cutxpi1_120)==(ymin_cutxpi1_120))])
    # xpimin_cutxpi2_120 = float(c.applyCuts(xpi,cutxpi2_120)[
    #     np.where(c.applyCuts(y,cutxpi2_120)==(ymin_cutxpi2_120))])
    # xpimin_cutxpi3_120 = float(c.applyCuts(xpi,cutxpi3_120)[
    #     np.where(c.applyCuts(y,cutxpi3_120)==(ymin_cutxpi3_120))])
    # xpimin_cutxpi4_120 = float(c.applyCuts(xpi,cutxpi4_120)[
    #     np.where(c.applyCuts(y,cutxpi4_120)==(ymin_cutxpi4_120))])
    # xpimin_cutxpi5_120 = float(c.applyCuts(xpi,cutxpi5_120)[
    #     np.where(c.applyCuts(y,cutxpi5_120)==(ymin_cutxpi5_120))])
    # xpimin_cutxpi6_120 = float(c.applyCuts(xpi,cutxpi6_120)[
    #     np.where(c.applyCuts(y,cutxpi6_120)==(ymin_cutxpi6_120))])
    # xpimin_cutxpi7_120 = float(c.applyCuts(xpi,cutxpi7_120)[
    #     np.where(c.applyCuts(y,cutxpi7_120)==(ymin_cutxpi7_120))])

    # xpimax_cutxpi1_120 = float(c.applyCuts(xpi,cutxpi1_120)[
    #     np.where(c.applyCuts(y,cutxpi1_120)==(ymax_cutxpi1_120))])
    # xpimax_cutxpi2_120 = float(c.applyCuts(xpi,cutxpi2_120)[
    #     np.where(c.applyCuts(y,cutxpi2_120)==(ymax_cutxpi2_120))])
    # xpimax_cutxpi3_120 = float(c.applyCuts(xpi,cutxpi3_120)[
    #     np.where(c.applyCuts(y,cutxpi3_120)==(ymax_cutxpi3_120))])
    # xpimax_cutxpi4_120 = float(c.applyCuts(xpi,cutxpi4_120)[
    #     np.where(c.applyCuts(y,cutxpi4_120)==(ymax_cutxpi4_120))])
    # xpimax_cutxpi5_120 = float(c.applyCuts(xpi,cutxpi5_120)[
    #     np.where(c.applyCuts(y,cutxpi5_120)==(ymax_cutxpi5_120))])
    # xpimax_cutxpi6_120 = float(c.applyCuts(xpi,cutxpi6_120)[
    #     np.where(c.applyCuts(y,cutxpi6_120)==(ymax_cutxpi6_120))])
    # xpimax_cutxpi7_120 = float(c.applyCuts(xpi,cutxpi7_120)[
    #     np.where(c.applyCuts(y,cutxpi7_120)==(ymax_cutxpi7_120))])

    # Q2min_cutxpi1_120 = float(c.applyCuts(Q2,cutxpi1_120)[
    #     np.where(c.applyCuts(y,cutxpi1_120)==(ymin_cutxpi1_120))])
    # Q2min_cutxpi2_120 = float(c.applyCuts(Q2,cutxpi2_120)[
    #     np.where(c.applyCuts(y,cutxpi2_120)==(ymin_cutxpi2_120))])
    # Q2min_cutxpi3_120 = float(c.applyCuts(Q2,cutxpi3_120)[
    #     np.where(c.applyCuts(y,cutxpi3_120)==(ymin_cutxpi3_120))])
    # Q2min_cutxpi4_120 = float(c.applyCuts(Q2,cutxpi4_120)[
    #     np.where(c.applyCuts(y,cutxpi4_120)==(ymin_cutxpi4_120))])
    # Q2min_cutxpi5_120 = float(c.applyCuts(Q2,cutxpi5_120)[
    #     np.where(c.applyCuts(y,cutxpi5_120)==(ymin_cutxpi5_120))])
    # Q2min_cutxpi6_120 = float(c.applyCuts(Q2,cutxpi6_120)[
    #     np.where(c.applyCuts(y,cutxpi6_120)==(ymin_cutxpi6_120))])
    # Q2min_cutxpi7_120 = float(c.applyCuts(Q2,cutxpi7_120)[
    #     np.where(c.applyCuts(y,cutxpi7_120)==(ymin_cutxpi7_120))])

    # Q2max_cutxpi1_120 = float(c.applyCuts(Q2,cutxpi1_120)[
    #     np.where(c.applyCuts(y,cutxpi1_120)==(ymax_cutxpi1_120))])
    # Q2max_cutxpi2_120 = float(c.applyCuts(Q2,cutxpi2_120)[
    #     np.where(c.applyCuts(y,cutxpi2_120)==(ymax_cutxpi2_120))])
    # Q2max_cutxpi3_120 = float(c.applyCuts(Q2,cutxpi3_120)[
    #     np.where(c.applyCuts(y,cutxpi3_120)==(ymax_cutxpi3_120))])
    # Q2max_cutxpi4_120 = float(c.applyCuts(Q2,cutxpi4_120)[
    #     np.where(c.applyCuts(y,cutxpi4_120)==(ymax_cutxpi4_120))])
    # Q2max_cutxpi5_120 = float(c.applyCuts(Q2,cutxpi5_120)[
    #     np.where(c.applyCuts(y,cutxpi5_120)==(ymax_cutxpi5_120))])
    # Q2max_cutxpi6_120 = float(c.applyCuts(Q2,cutxpi6_120)[
    #     np.where(c.applyCuts(y,cutxpi6_120)==(ymax_cutxpi6_120))])
    # Q2max_cutxpi7_120 = float(c.applyCuts(Q2,cutxpi7_120)[
    #     np.where(c.applyCuts(y,cutxpi7_120)==(ymax_cutxpi7_120))])
    
    # #################################

    # #################################
        
    
    # ymin_cutxpi1_240 = min(c.applyCuts(y,cutxpi1_240))
    # ymin_cutxpi2_240 = min(c.applyCuts(y,cutxpi2_240))
    # ymin_cutxpi3_240 = min(c.applyCuts(y,cutxpi3_240))
    # ymin_cutxpi4_240 = min(c.applyCuts(y,cutxpi4_240))
    # ymin_cutxpi5_240 = min(c.applyCuts(y,cutxpi5_240))
    # ymin_cutxpi6_240 = min(c.applyCuts(y,cutxpi6_240))
    # ymin_cutxpi7_240 = min(c.applyCuts(y,cutxpi7_240))
    
    # ymax_cutxpi1_240 = max(c.applyCuts(y,cutxpi1_240))
    # ymax_cutxpi2_240 = max(c.applyCuts(y,cutxpi2_240))
    # ymax_cutxpi3_240 = max(c.applyCuts(y,cutxpi3_240))
    # ymax_cutxpi4_240 = max(c.applyCuts(y,cutxpi4_240))
    # ymax_cutxpi5_240 = max(c.applyCuts(y,cutxpi5_240))
    # ymax_cutxpi6_240 = max(c.applyCuts(y,cutxpi6_240))
    # ymax_cutxpi7_240 = max(c.applyCuts(y,cutxpi7_240))

    # yieldmin_cutxpi1_240 = float(c.applyCuts(lumi,cutxpi1_240)[
    #     np.where(c.applyCuts(y,cutxpi1_240)==(ymin_cutxpi1_240))])
    # yieldmin_cutxpi2_240 = float(c.applyCuts(lumi,cutxpi2_240)[
    #     np.where(c.applyCuts(y,cutxpi2_240)==(ymin_cutxpi2_240))])
    # yieldmin_cutxpi3_240 = float(c.applyCuts(lumi,cutxpi3_240)[
    #     np.where(c.applyCuts(y,cutxpi3_240)==(ymin_cutxpi3_240))])
    # yieldmin_cutxpi4_240 = float(c.applyCuts(lumi,cutxpi4_240)[
    #     np.where(c.applyCuts(y,cutxpi4_240)==(ymin_cutxpi4_240))])
    # yieldmin_cutxpi5_240 = float(c.applyCuts(lumi,cutxpi5_240)[
    #     np.where(c.applyCuts(y,cutxpi5_240)==(ymin_cutxpi5_240))])
    # yieldmin_cutxpi6_240 = float(c.applyCuts(lumi,cutxpi6_240)[
    #     np.where(c.applyCuts(y,cutxpi6_240)==(ymin_cutxpi6_240))])
    # yieldmin_cutxpi7_240 = float(c.applyCuts(lumi,cutxpi7_240)[
    #     np.where(c.applyCuts(y,cutxpi7_240)==(ymin_cutxpi7_240))])
    
    # yieldmax_cutxpi1_240 = float(c.applyCuts(lumi,cutxpi1_240)[
    #     np.where(c.applyCuts(y,cutxpi1_240)==(ymax_cutxpi1_240))])
    # yieldmax_cutxpi2_240 = float(c.applyCuts(lumi,cutxpi2_240)[
    #     np.where(c.applyCuts(y,cutxpi2_240)==(ymax_cutxpi2_240))])
    # yieldmax_cutxpi3_240 = float(c.applyCuts(lumi,cutxpi3_240)[
    #     np.where(c.applyCuts(y,cutxpi3_240)==(ymax_cutxpi3_240))])
    # yieldmax_cutxpi4_240 = float(c.applyCuts(lumi,cutxpi4_240)[
    #     np.where(c.applyCuts(y,cutxpi4_240)==(ymax_cutxpi4_240))])
    # yieldmax_cutxpi5_240 = float(c.applyCuts(lumi,cutxpi5_240)[
    #     np.where(c.applyCuts(y,cutxpi5_240)==(ymax_cutxpi5_240))])
    # yieldmax_cutxpi6_240 = float(c.applyCuts(lumi,cutxpi6_240)[
    #     np.where(c.applyCuts(y,cutxpi6_240)==(ymax_cutxpi6_240))])
    # yieldmax_cutxpi7_240 = float(c.applyCuts(lumi,cutxpi7_240)[
    #     np.where(c.applyCuts(y,cutxpi7_240)==(ymax_cutxpi7_240))])

    # xLmin_cutxpi1_240 = float(c.applyCuts(xL,cutxpi1_240)[
    #     np.where(c.applyCuts(y,cutxpi1_240)==(ymin_cutxpi1_240))])
    # xLmin_cutxpi2_240 = float(c.applyCuts(xL,cutxpi2_240)[
    #     np.where(c.applyCuts(y,cutxpi2_240)==(ymin_cutxpi2_240))])
    # xLmin_cutxpi3_240 = float(c.applyCuts(xL,cutxpi3_240)[
    #     np.where(c.applyCuts(y,cutxpi3_240)==(ymin_cutxpi3_240))])
    # xLmin_cutxpi4_240 = float(c.applyCuts(xL,cutxpi4_240)[
    #     np.where(c.applyCuts(y,cutxpi4_240)==(ymin_cutxpi4_240))])
    # xLmin_cutxpi5_240 = float(c.applyCuts(xL,cutxpi5_240)[
    #     np.where(c.applyCuts(y,cutxpi5_240)==(ymin_cutxpi5_240))])
    # xLmin_cutxpi6_240 = float(c.applyCuts(xL,cutxpi6_240)[
    #     np.where(c.applyCuts(y,cutxpi6_240)==(ymin_cutxpi6_240))])
    # xLmin_cutxpi7_240 = float(c.applyCuts(xL,cutxpi7_240)[
    #     np.where(c.applyCuts(y,cutxpi7_240)==(ymin_cutxpi7_240))])

    # xLmax_cutxpi1_240 = float(c.applyCuts(xL,cutxpi1_240)[
    #     np.where(c.applyCuts(y,cutxpi1_240)==(ymax_cutxpi1_240))])
    # xLmax_cutxpi2_240 = float(c.applyCuts(xL,cutxpi2_240)[
    #     np.where(c.applyCuts(y,cutxpi2_240)==(ymax_cutxpi2_240))])
    # xLmax_cutxpi3_240 = float(c.applyCuts(xL,cutxpi3_240)[
    #     np.where(c.applyCuts(y,cutxpi3_240)==(ymax_cutxpi3_240))])
    # xLmax_cutxpi4_240 = float(c.applyCuts(xL,cutxpi4_240)[
    #     np.where(c.applyCuts(y,cutxpi4_240)==(ymax_cutxpi4_240))])
    # xLmax_cutxpi5_240 = float(c.applyCuts(xL,cutxpi5_240)[
    #     np.where(c.applyCuts(y,cutxpi5_240)==(ymax_cutxpi5_240))])
    # xLmax_cutxpi6_240 = float(c.applyCuts(xL,cutxpi6_240)[
    #     np.where(c.applyCuts(y,cutxpi6_240)==(ymax_cutxpi6_240))])
    # xLmax_cutxpi7_240 = float(c.applyCuts(xL,cutxpi7_240)[
    #     np.where(c.applyCuts(y,cutxpi7_240)==(ymax_cutxpi7_240))])

    # xpimin_cutxpi1_240 = float(c.applyCuts(xpi,cutxpi1_240)[
    #     np.where(c.applyCuts(y,cutxpi1_240)==(ymin_cutxpi1_240))])
    # xpimin_cutxpi2_240 = float(c.applyCuts(xpi,cutxpi2_240)[
    #     np.where(c.applyCuts(y,cutxpi2_240)==(ymin_cutxpi2_240))])
    # xpimin_cutxpi3_240 = float(c.applyCuts(xpi,cutxpi3_240)[
    #     np.where(c.applyCuts(y,cutxpi3_240)==(ymin_cutxpi3_240))])
    # xpimin_cutxpi4_240 = float(c.applyCuts(xpi,cutxpi4_240)[
    #     np.where(c.applyCuts(y,cutxpi4_240)==(ymin_cutxpi4_240))])
    # xpimin_cutxpi5_240 = float(c.applyCuts(xpi,cutxpi5_240)[
    #     np.where(c.applyCuts(y,cutxpi5_240)==(ymin_cutxpi5_240))])
    # xpimin_cutxpi6_240 = float(c.applyCuts(xpi,cutxpi6_240)[
    #     np.where(c.applyCuts(y,cutxpi6_240)==(ymin_cutxpi6_240))])
    # xpimin_cutxpi7_240 = float(c.applyCuts(xpi,cutxpi7_240)[
    #     np.where(c.applyCuts(y,cutxpi7_240)==(ymin_cutxpi7_240))])

    # xpimax_cutxpi1_240 = float(c.applyCuts(xpi,cutxpi1_240)[
    #     np.where(c.applyCuts(y,cutxpi1_240)==(ymax_cutxpi1_240))])
    # xpimax_cutxpi2_240 = float(c.applyCuts(xpi,cutxpi2_240)[
    #     np.where(c.applyCuts(y,cutxpi2_240)==(ymax_cutxpi2_240))])
    # xpimax_cutxpi3_240 = float(c.applyCuts(xpi,cutxpi3_240)[
    #     np.where(c.applyCuts(y,cutxpi3_240)==(ymax_cutxpi3_240))])
    # xpimax_cutxpi4_240 = float(c.applyCuts(xpi,cutxpi4_240)[
    #     np.where(c.applyCuts(y,cutxpi4_240)==(ymax_cutxpi4_240))])
    # xpimax_cutxpi5_240 = float(c.applyCuts(xpi,cutxpi5_240)[
    #     np.where(c.applyCuts(y,cutxpi5_240)==(ymax_cutxpi5_240))])
    # xpimax_cutxpi6_240 = float(c.applyCuts(xpi,cutxpi6_240)[
    #     np.where(c.applyCuts(y,cutxpi6_240)==(ymax_cutxpi6_240))])
    # xpimax_cutxpi7_240 = float(c.applyCuts(xpi,cutxpi7_240)[
    #     np.where(c.applyCuts(y,cutxpi7_240)==(ymax_cutxpi7_240))])

    # Q2min_cutxpi1_240 = float(c.applyCuts(Q2,cutxpi1_240)[
    #     np.where(c.applyCuts(y,cutxpi1_240)==(ymin_cutxpi1_240))])
    # Q2min_cutxpi2_240 = float(c.applyCuts(Q2,cutxpi2_240)[
    #     np.where(c.applyCuts(y,cutxpi2_240)==(ymin_cutxpi2_240))])
    # Q2min_cutxpi3_240 = float(c.applyCuts(Q2,cutxpi3_240)[
    #     np.where(c.applyCuts(y,cutxpi3_240)==(ymin_cutxpi3_240))])
    # Q2min_cutxpi4_240 = float(c.applyCuts(Q2,cutxpi4_240)[
    #     np.where(c.applyCuts(y,cutxpi4_240)==(ymin_cutxpi4_240))])
    # Q2min_cutxpi5_240 = float(c.applyCuts(Q2,cutxpi5_240)[
    #     np.where(c.applyCuts(y,cutxpi5_240)==(ymin_cutxpi5_240))])
    # Q2min_cutxpi6_240 = float(c.applyCuts(Q2,cutxpi6_240)[
    #     np.where(c.applyCuts(y,cutxpi6_240)==(ymin_cutxpi6_240))])
    # Q2min_cutxpi7_240 = float(c.applyCuts(Q2,cutxpi7_240)[
    #     np.where(c.applyCuts(y,cutxpi7_240)==(ymin_cutxpi7_240))])

    # Q2max_cutxpi1_240 = float(c.applyCuts(Q2,cutxpi1_240)[
    #     np.where(c.applyCuts(y,cutxpi1_240)==(ymax_cutxpi1_240))])
    # Q2max_cutxpi2_240 = float(c.applyCuts(Q2,cutxpi2_240)[
    #     np.where(c.applyCuts(y,cutxpi2_240)==(ymax_cutxpi2_240))])
    # Q2max_cutxpi3_240 = float(c.applyCuts(Q2,cutxpi3_240)[
    #     np.where(c.applyCuts(y,cutxpi3_240)==(ymax_cutxpi3_240))])
    # Q2max_cutxpi4_240 = float(c.applyCuts(Q2,cutxpi4_240)[
    #     np.where(c.applyCuts(y,cutxpi4_240)==(ymax_cutxpi4_240))])
    # Q2max_cutxpi5_240 = float(c.applyCuts(Q2,cutxpi5_240)[
    #     np.where(c.applyCuts(y,cutxpi5_240)==(ymax_cutxpi5_240))])
    # Q2max_cutxpi6_240 = float(c.applyCuts(Q2,cutxpi6_240)[
    #     np.where(c.applyCuts(y,cutxpi6_240)==(ymax_cutxpi6_240))])
    # Q2max_cutxpi7_240 = float(c.applyCuts(Q2,cutxpi7_240)[
    #     np.where(c.applyCuts(y,cutxpi7_240)==(ymax_cutxpi7_240))])
    
    # #################################

    # #################################
        
    
    # # ymin_cutxpi1_480 = min(c.applyCuts(y,cutxpi1_480))
    # ymin_cutxpi2_480 = min(c.applyCuts(y,cutxpi2_480))
    # ymin_cutxpi3_480 = min(c.applyCuts(y,cutxpi3_480))
    # ymin_cutxpi4_480 = min(c.applyCuts(y,cutxpi4_480))
    # ymin_cutxpi5_480 = min(c.applyCuts(y,cutxpi5_480))
    # ymin_cutxpi6_480 = min(c.applyCuts(y,cutxpi6_480))
    # ymin_cutxpi7_480 = min(c.applyCuts(y,cutxpi7_480))
    
    # # ymax_cutxpi1_480 = max(c.applyCuts(y,cutxpi1_480))
    # ymax_cutxpi2_480 = max(c.applyCuts(y,cutxpi2_480))
    # ymax_cutxpi3_480 = max(c.applyCuts(y,cutxpi3_480))
    # ymax_cutxpi4_480 = max(c.applyCuts(y,cutxpi4_480))
    # ymax_cutxpi5_480 = max(c.applyCuts(y,cutxpi5_480))
    # ymax_cutxpi6_480 = max(c.applyCuts(y,cutxpi6_480))
    # ymax_cutxpi7_480 = max(c.applyCuts(y,cutxpi7_480))

    # # yieldmin_cutxpi1_480 = float(c.applyCuts(lumi,cutxpi1_480)[
    # #     np.where(c.applyCuts(y,cutxpi1_480)==(ymin_cutxpi1_480))])
    # yieldmin_cutxpi2_480 = float(c.applyCuts(lumi,cutxpi2_480)[
    #     np.where(c.applyCuts(y,cutxpi2_480)==(ymin_cutxpi2_480))])
    # yieldmin_cutxpi3_480 = float(c.applyCuts(lumi,cutxpi3_480)[
    #     np.where(c.applyCuts(y,cutxpi3_480)==(ymin_cutxpi3_480))])
    # yieldmin_cutxpi4_480 = float(c.applyCuts(lumi,cutxpi4_480)[
    #     np.where(c.applyCuts(y,cutxpi4_480)==(ymin_cutxpi4_480))])
    # yieldmin_cutxpi5_480 = float(c.applyCuts(lumi,cutxpi5_480)[
    #     np.where(c.applyCuts(y,cutxpi5_480)==(ymin_cutxpi5_480))])
    # yieldmin_cutxpi6_480 = float(c.applyCuts(lumi,cutxpi6_480)[
    #     np.where(c.applyCuts(y,cutxpi6_480)==(ymin_cutxpi6_480))])
    # yieldmin_cutxpi7_480 = float(c.applyCuts(lumi,cutxpi7_480)[
    #     np.where(c.applyCuts(y,cutxpi7_480)==(ymin_cutxpi7_480))])
    
    # # yieldmax_cutxpi1_480 = float(c.applyCuts(lumi,cutxpi1_480)[
    # #     np.where(c.applyCuts(y,cutxpi1_480)==(ymax_cutxpi1_480))])
    # yieldmax_cutxpi2_480 = float(c.applyCuts(lumi,cutxpi2_480)[
    #     np.where(c.applyCuts(y,cutxpi2_480)==(ymax_cutxpi2_480))])
    # yieldmax_cutxpi3_480 = float(c.applyCuts(lumi,cutxpi3_480)[
    #     np.where(c.applyCuts(y,cutxpi3_480)==(ymax_cutxpi3_480))])
    # yieldmax_cutxpi4_480 = float(c.applyCuts(lumi,cutxpi4_480)[
    #     np.where(c.applyCuts(y,cutxpi4_480)==(ymax_cutxpi4_480))])
    # yieldmax_cutxpi5_480 = float(c.applyCuts(lumi,cutxpi5_480)[
    #     np.where(c.applyCuts(y,cutxpi5_480)==(ymax_cutxpi5_480))])
    # yieldmax_cutxpi6_480 = float(c.applyCuts(lumi,cutxpi6_480)[
    #     np.where(c.applyCuts(y,cutxpi6_480)==(ymax_cutxpi6_480))])
    # yieldmax_cutxpi7_480 = float(c.applyCuts(lumi,cutxpi7_480)[
    #     np.where(c.applyCuts(y,cutxpi7_480)==(ymax_cutxpi7_480))])

    # # xLmin_cutxpi1_480 = float(c.applyCuts(xL,cutxpi1_480)[
    # #     np.where(c.applyCuts(y,cutxpi1_480)==(ymin_cutxpi1_480))])
    # xLmin_cutxpi2_480 = float(c.applyCuts(xL,cutxpi2_480)[
    #     np.where(c.applyCuts(y,cutxpi2_480)==(ymin_cutxpi2_480))])
    # xLmin_cutxpi3_480 = float(c.applyCuts(xL,cutxpi3_480)[
    #     np.where(c.applyCuts(y,cutxpi3_480)==(ymin_cutxpi3_480))])
    # xLmin_cutxpi4_480 = float(c.applyCuts(xL,cutxpi4_480)[
    #     np.where(c.applyCuts(y,cutxpi4_480)==(ymin_cutxpi4_480))])
    # xLmin_cutxpi5_480 = float(c.applyCuts(xL,cutxpi5_480)[
    #     np.where(c.applyCuts(y,cutxpi5_480)==(ymin_cutxpi5_480))])
    # xLmin_cutxpi6_480 = float(c.applyCuts(xL,cutxpi6_480)[
    #     np.where(c.applyCuts(y,cutxpi6_480)==(ymin_cutxpi6_480))])
    # xLmin_cutxpi7_480 = float(c.applyCuts(xL,cutxpi7_480)[
    #     np.where(c.applyCuts(y,cutxpi7_480)==(ymin_cutxpi7_480))])

    # # xLmax_cutxpi1_480 = float(c.applyCuts(xL,cutxpi1_480)[
    # #     np.where(c.applyCuts(y,cutxpi1_480)==(ymax_cutxpi1_480))])
    # xLmax_cutxpi2_480 = float(c.applyCuts(xL,cutxpi2_480)[
    #     np.where(c.applyCuts(y,cutxpi2_480)==(ymax_cutxpi2_480))])
    # xLmax_cutxpi3_480 = float(c.applyCuts(xL,cutxpi3_480)[
    #     np.where(c.applyCuts(y,cutxpi3_480)==(ymax_cutxpi3_480))])
    # xLmax_cutxpi4_480 = float(c.applyCuts(xL,cutxpi4_480)[
    #     np.where(c.applyCuts(y,cutxpi4_480)==(ymax_cutxpi4_480))])
    # xLmax_cutxpi5_480 = float(c.applyCuts(xL,cutxpi5_480)[
    #     np.where(c.applyCuts(y,cutxpi5_480)==(ymax_cutxpi5_480))])
    # xLmax_cutxpi6_480 = float(c.applyCuts(xL,cutxpi6_480)[
    #     np.where(c.applyCuts(y,cutxpi6_480)==(ymax_cutxpi6_480))])
    # xLmax_cutxpi7_480 = float(c.applyCuts(xL,cutxpi7_480)[
    #     np.where(c.applyCuts(y,cutxpi7_480)==(ymax_cutxpi7_480))])

    # # xpimin_cutxpi1_480 = float(c.applyCuts(xpi,cutxpi1_480)[
    # #     np.where(c.applyCuts(y,cutxpi1_480)==(ymin_cutxpi1_480))])
    # xpimin_cutxpi2_480 = float(c.applyCuts(xpi,cutxpi2_480)[
    #     np.where(c.applyCuts(y,cutxpi2_480)==(ymin_cutxpi2_480))])
    # xpimin_cutxpi3_480 = float(c.applyCuts(xpi,cutxpi3_480)[
    #     np.where(c.applyCuts(y,cutxpi3_480)==(ymin_cutxpi3_480))])
    # xpimin_cutxpi4_480 = float(c.applyCuts(xpi,cutxpi4_480)[
    #     np.where(c.applyCuts(y,cutxpi4_480)==(ymin_cutxpi4_480))])
    # xpimin_cutxpi5_480 = float(c.applyCuts(xpi,cutxpi5_480)[
    #     np.where(c.applyCuts(y,cutxpi5_480)==(ymin_cutxpi5_480))])
    # xpimin_cutxpi6_480 = float(c.applyCuts(xpi,cutxpi6_480)[
    #     np.where(c.applyCuts(y,cutxpi6_480)==(ymin_cutxpi6_480))])
    # xpimin_cutxpi7_480 = float(c.applyCuts(xpi,cutxpi7_480)[
    #     np.where(c.applyCuts(y,cutxpi7_480)==(ymin_cutxpi7_480))])

    # # xpimax_cutxpi1_480 = float(c.applyCuts(xpi,cutxpi1_480)[
    # #     np.where(c.applyCuts(y,cutxpi1_480)==(ymax_cutxpi1_480))])
    # xpimax_cutxpi2_480 = float(c.applyCuts(xpi,cutxpi2_480)[
    #     np.where(c.applyCuts(y,cutxpi2_480)==(ymax_cutxpi2_480))])
    # xpimax_cutxpi3_480 = float(c.applyCuts(xpi,cutxpi3_480)[
    #     np.where(c.applyCuts(y,cutxpi3_480)==(ymax_cutxpi3_480))])
    # xpimax_cutxpi4_480 = float(c.applyCuts(xpi,cutxpi4_480)[
    #     np.where(c.applyCuts(y,cutxpi4_480)==(ymax_cutxpi4_480))])
    # xpimax_cutxpi5_480 = float(c.applyCuts(xpi,cutxpi5_480)[
    #     np.where(c.applyCuts(y,cutxpi5_480)==(ymax_cutxpi5_480))])
    # xpimax_cutxpi6_480 = float(c.applyCuts(xpi,cutxpi6_480)[
    #     np.where(c.applyCuts(y,cutxpi6_480)==(ymax_cutxpi6_480))])
    # xpimax_cutxpi7_480 = float(c.applyCuts(xpi,cutxpi7_480)[
    #     np.where(c.applyCuts(y,cutxpi7_480)==(ymax_cutxpi7_480))])

    # # Q2min_cutxpi1_480 = float(c.applyCuts(Q2,cutxpi1_480)[
    # #     np.where(c.applyCuts(y,cutxpi1_480)==(ymin_cutxpi1_480))])
    # Q2min_cutxpi2_480 = float(c.applyCuts(Q2,cutxpi2_480)[
    #     np.where(c.applyCuts(y,cutxpi2_480)==(ymin_cutxpi2_480))])
    # Q2min_cutxpi3_480 = float(c.applyCuts(Q2,cutxpi3_480)[
    #     np.where(c.applyCuts(y,cutxpi3_480)==(ymin_cutxpi3_480))])
    # Q2min_cutxpi4_480 = float(c.applyCuts(Q2,cutxpi4_480)[
    #     np.where(c.applyCuts(y,cutxpi4_480)==(ymin_cutxpi4_480))])
    # Q2min_cutxpi5_480 = float(c.applyCuts(Q2,cutxpi5_480)[
    #     np.where(c.applyCuts(y,cutxpi5_480)==(ymin_cutxpi5_480))])
    # Q2min_cutxpi6_480 = float(c.applyCuts(Q2,cutxpi6_480)[
    #     np.where(c.applyCuts(y,cutxpi6_480)==(ymin_cutxpi6_480))])
    # Q2min_cutxpi7_480 = float(c.applyCuts(Q2,cutxpi7_480)[
    #     np.where(c.applyCuts(y,cutxpi7_480)==(ymin_cutxpi7_480))])

    # # Q2max_cutxpi1_480 = float(c.applyCuts(Q2,cutxpi1_480)[
    # #     np.where(c.applyCuts(y,cutxpi1_480)==(ymax_cutxpi1_480))])
    # Q2max_cutxpi2_480 = float(c.applyCuts(Q2,cutxpi2_480)[
    #     np.where(c.applyCuts(y,cutxpi2_480)==(ymax_cutxpi2_480))])
    # Q2max_cutxpi3_480 = float(c.applyCuts(Q2,cutxpi3_480)[
    #     np.where(c.applyCuts(y,cutxpi3_480)==(ymax_cutxpi3_480))])
    # Q2max_cutxpi4_480 = float(c.applyCuts(Q2,cutxpi4_480)[
    #     np.where(c.applyCuts(y,cutxpi4_480)==(ymax_cutxpi4_480))])
    # Q2max_cutxpi5_480 = float(c.applyCuts(Q2,cutxpi5_480)[
    #     np.where(c.applyCuts(y,cutxpi5_480)==(ymax_cutxpi5_480))])
    # Q2max_cutxpi6_480 = float(c.applyCuts(Q2,cutxpi6_480)[
    #     np.where(c.applyCuts(y,cutxpi6_480)==(ymax_cutxpi6_480))])
    # Q2max_cutxpi7_480 = float(c.applyCuts(Q2,cutxpi7_480)[
    #     np.where(c.applyCuts(y,cutxpi7_480)==(ymax_cutxpi7_480))])
    
    # #################################

    
    
    
    # table_dict = {
    
    #         "15" : 
    #          {"0.1" :{
    #              "Q2" : [Q2min_cutxpi1_15,Q2max_cutxpi1_15],
    #              "y bin" : [ymin_cutxpi1_15,ymax_cutxpi1_15],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi1_15,yieldmax_cutxpi1_15],
    #              "xL" : [xLmin_cutxpi1_15,xLmax_cutxpi1_15],
    #              "xpi": [xpimin_cutxpi1_15,xpimax_cutxpi1_15]
    #          },
    #          "0.2" :{
    #              "Q2" : [Q2min_cutxpi2_15,Q2max_cutxpi2_15],
    #              "y bin" : [ymin_cutxpi2_15,ymax_cutxpi2_15],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi2_15,yieldmax_cutxpi2_15],
    #              "xL" : [xLmin_cutxpi2_15,xLmax_cutxpi2_15],
    #              "xpi": [xpimin_cutxpi2_15,xpimax_cutxpi2_15]
    #          },
    #          "0.3" :{
    #              "Q2" : [Q2min_cutxpi3_15,Q2max_cutxpi3_15],
    #              "y bin" : [ymin_cutxpi3_15,ymax_cutxpi3_15],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi3_15,yieldmax_cutxpi3_15],
    #              "xL" : [xLmin_cutxpi3_15,xLmax_cutxpi3_15],
    #              "xpi": [xpimin_cutxpi3_15,xpimax_cutxpi3_15]
    #          }
    #         },
    #         "30" : 
    #          {"0.1" :{
    #              "Q2" : [Q2min_cutxpi1_30,Q2max_cutxpi1_30],
    #              "y bin" : [ymin_cutxpi1_30,ymax_cutxpi1_30],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi1_30,yieldmax_cutxpi1_30],
    #              "xL" : [xLmin_cutxpi1_30,xLmax_cutxpi1_30],
    #              "xpi": [xpimin_cutxpi1_30,xpimax_cutxpi1_30]
    #          },
    #          "0.2" :{
    #              "Q2" : [Q2min_cutxpi2_30,Q2max_cutxpi2_30],
    #              "y bin" : [ymin_cutxpi2_30,ymax_cutxpi2_30],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi2_30,yieldmax_cutxpi2_30],
    #              "xL" : [xLmin_cutxpi2_30,xLmax_cutxpi2_30],
    #              "xpi": [xpimin_cutxpi2_30,xpimax_cutxpi2_30]
    #          },
    #          "0.3" :{
    #              "Q2" : [Q2min_cutxpi3_30,Q2max_cutxpi3_30],
    #              "y bin" : [ymin_cutxpi3_30,ymax_cutxpi3_30],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi3_30,yieldmax_cutxpi3_30],
    #              "xL" : [xLmin_cutxpi3_30,xLmax_cutxpi3_30],
    #              "xpi": [xpimin_cutxpi3_30,xpimax_cutxpi3_30]
    #          },
    #          "0.4" :{
    #              "Q2" : [Q2min_cutxpi4_30,Q2max_cutxpi4_30],
    #              "y bin" : [ymin_cutxpi4_30,ymax_cutxpi4_30],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi4_30,yieldmax_cutxpi4_30],
    #              "xL" : [xLmin_cutxpi4_30,xLmax_cutxpi4_30],
    #              "xpi": [xpimin_cutxpi4_30,xpimax_cutxpi4_30]
    #          },
    #          "0.5" :{
    #              "Q2" : [Q2min_cutxpi5_30,Q2max_cutxpi5_30],
    #              "y bin" : [ymin_cutxpi5_30,ymax_cutxpi5_30],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi5_30,yieldmax_cutxpi5_30],
    #              "xL" : [xLmin_cutxpi5_30,xLmax_cutxpi5_30],
    #              "xpi": [xpimin_cutxpi5_30,xpimax_cutxpi5_30]
    #          },
    #          "0.6" :{
    #              "Q2" : [Q2min_cutxpi6_30,Q2max_cutxpi6_30],
    #              "y bin" : [ymin_cutxpi6_30,ymax_cutxpi6_30],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi6_30,yieldmax_cutxpi6_30],
    #              "xL" : [xLmin_cutxpi6_30,xLmax_cutxpi6_30],
    #              "xpi": [xpimin_cutxpi6_30,xpimax_cutxpi6_30]
    #          }
    #         },
    #         "60" : 
    #          {"0.1" :{
    #              "Q2" : [Q2min_cutxpi1_60,Q2max_cutxpi1_60],
    #              "y bin" : [ymin_cutxpi1_60,ymax_cutxpi1_60],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi1_60,yieldmax_cutxpi1_60],
    #              "xL" : [xLmin_cutxpi1_60,xLmax_cutxpi1_60],
    #              "xpi": [xpimin_cutxpi1_60,xpimax_cutxpi1_60]
    #          },
    #          "0.2" :{
    #              "Q2" : [Q2min_cutxpi2_60,Q2max_cutxpi2_60],
    #              "y bin" : [ymin_cutxpi2_60,ymax_cutxpi2_60],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi2_60,yieldmax_cutxpi2_60],
    #              "xL" : [xLmin_cutxpi2_60,xLmax_cutxpi2_60],
    #              "xpi": [xpimin_cutxpi2_60,xpimax_cutxpi2_60]
    #          },
    #          "0.3" :{
    #              "Q2" : [Q2min_cutxpi3_60,Q2max_cutxpi3_60],
    #              "y bin" : [ymin_cutxpi3_60,ymax_cutxpi3_60],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi3_60,yieldmax_cutxpi3_60],
    #              "xL" : [xLmin_cutxpi3_60,xLmax_cutxpi3_60],
    #              "xpi": [xpimin_cutxpi3_60,xpimax_cutxpi3_60]
    #          },
    #          "0.4" :{
    #              "Q2" : [Q2min_cutxpi4_60,Q2max_cutxpi4_60],
    #              "y bin" : [ymin_cutxpi4_60,ymax_cutxpi4_60],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi4_60,yieldmax_cutxpi4_60],
    #              "xL" : [xLmin_cutxpi4_60,xLmax_cutxpi4_60],
    #              "xpi": [xpimin_cutxpi4_60,xpimax_cutxpi4_60]
    #          },
    #          "0.5" :{
    #              "Q2" : [Q2min_cutxpi5_60,Q2max_cutxpi5_60],
    #              "y bin" : [ymin_cutxpi5_60,ymax_cutxpi5_60],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi5_60,yieldmax_cutxpi5_60],
    #              "xL" : [xLmin_cutxpi5_60,xLmax_cutxpi5_60],
    #              "xpi": [xpimin_cutxpi5_60,xpimax_cutxpi5_60]
    #          },
    #          "0.6" :{
    #              "Q2" : [Q2min_cutxpi6_60,Q2max_cutxpi6_60],
    #              "y bin" : [ymin_cutxpi6_60,ymax_cutxpi6_60],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi6_60,yieldmax_cutxpi6_60],
    #              "xL" : [xLmin_cutxpi6_60,xLmax_cutxpi6_60],
    #              "xpi": [xpimin_cutxpi6_60,xpimax_cutxpi6_60]
    #          },
    #          "0.7" :{
    #              "Q2" : [Q2min_cutxpi7_60,Q2max_cutxpi7_60],
    #              "y bin" : [ymin_cutxpi7_60,ymax_cutxpi7_60],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi7_60,yieldmax_cutxpi7_60],
    #              "xL" : [xLmin_cutxpi7_60,xLmax_cutxpi7_60],
    #              "xpi": [xpimin_cutxpi7_60,xpimax_cutxpi7_60]
    #          }
    #         },
    #         "120" : 
    #          {"0.1" :{
    #              "Q2" : [Q2min_cutxpi1_120,Q2max_cutxpi1_120],
    #              "y bin" : [ymin_cutxpi1_120,ymax_cutxpi1_120],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi1_120,yieldmax_cutxpi1_120],
    #              "xL" : [xLmin_cutxpi1_120,xLmax_cutxpi1_120],
    #              "xpi": [xpimin_cutxpi1_120,xpimax_cutxpi1_120]
    #          },
    #          "0.2" :{
    #              "Q2" : [Q2min_cutxpi2_120,Q2max_cutxpi2_120],
    #              "y bin" : [ymin_cutxpi2_120,ymax_cutxpi2_120],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi2_120,yieldmax_cutxpi2_120],
    #              "xL" : [xLmin_cutxpi2_120,xLmax_cutxpi2_120],
    #              "xpi": [xpimin_cutxpi2_120,xpimax_cutxpi2_120]
    #          },
    #          "0.3" :{
    #              "Q2" : [Q2min_cutxpi3_120,Q2max_cutxpi3_120],
    #              "y bin" : [ymin_cutxpi3_120,ymax_cutxpi3_120],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi3_120,yieldmax_cutxpi3_120],
    #              "xL" : [xLmin_cutxpi3_120,xLmax_cutxpi3_120],
    #              "xpi": [xpimin_cutxpi3_120,xpimax_cutxpi3_120]
    #          },
    #          "0.4" :{
    #              "Q2" : [Q2min_cutxpi4_120,Q2max_cutxpi4_120],
    #              "y bin" : [ymin_cutxpi4_120,ymax_cutxpi4_120],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi4_120,yieldmax_cutxpi4_120],
    #              "xL" : [xLmin_cutxpi4_120,xLmax_cutxpi4_120],
    #              "xpi": [xpimin_cutxpi4_120,xpimax_cutxpi4_120]
    #          },
    #          "0.5" :{
    #              "Q2" : [Q2min_cutxpi5_120,Q2max_cutxpi5_120],
    #              "y bin" : [ymin_cutxpi5_120,ymax_cutxpi5_120],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi5_120,yieldmax_cutxpi5_120],
    #              "xL" : [xLmin_cutxpi5_120,xLmax_cutxpi5_120],
    #              "xpi": [xpimin_cutxpi5_120,xpimax_cutxpi5_120]
    #          },
    #          "0.6" :{
    #              "Q2" : [Q2min_cutxpi6_120,Q2max_cutxpi6_120],
    #              "y bin" : [ymin_cutxpi6_120,ymax_cutxpi6_120],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi6_120,yieldmax_cutxpi6_120],
    #              "xL" : [xLmin_cutxpi6_120,xLmax_cutxpi6_120],
    #              "xpi": [xpimin_cutxpi6_120,xpimax_cutxpi6_120]
    #          },
    #          "0.7" :{
    #              "Q2" : [Q2min_cutxpi7_120,Q2max_cutxpi7_120],
    #              "y bin" : [ymin_cutxpi7_120,ymax_cutxpi7_120],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi7_120,yieldmax_cutxpi7_120],
    #              "xL" : [xLmin_cutxpi7_120,xLmax_cutxpi7_120],
    #              "xpi": [xpimin_cutxpi7_120,xpimax_cutxpi7_120]
    #          }
    #         },
    #         "240" : 
    #          {"0.1" :{
    #              "Q2" : [Q2min_cutxpi1_240,Q2max_cutxpi1_240],
    #              "y bin" : [ymin_cutxpi1_240,ymax_cutxpi1_240],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi1_240,yieldmax_cutxpi1_240],
    #              "xL" : [xLmin_cutxpi1_240,xLmax_cutxpi1_240],
    #              "xpi": [xpimin_cutxpi1_240,xpimax_cutxpi1_240]
    #          },
    #          "0.2" :{
    #              "Q2" : [Q2min_cutxpi2_240,Q2max_cutxpi2_240],
    #              "y bin" : [ymin_cutxpi2_240,ymax_cutxpi2_240],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi2_240,yieldmax_cutxpi2_240],
    #              "xL" : [xLmin_cutxpi2_240,xLmax_cutxpi2_240],
    #              "xpi": [xpimin_cutxpi2_240,xpimax_cutxpi2_240]
    #          },
    #          "0.3" :{
    #              "Q2" : [Q2min_cutxpi3_240,Q2max_cutxpi3_240],
    #              "y bin" : [ymin_cutxpi3_240,ymax_cutxpi3_240],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi3_240,yieldmax_cutxpi3_240],
    #              "xL" : [xLmin_cutxpi3_240,xLmax_cutxpi3_240],
    #              "xpi": [xpimin_cutxpi3_240,xpimax_cutxpi3_240]
    #          },
    #          "0.4" :{
    #              "Q2" : [Q2min_cutxpi4_240,Q2max_cutxpi4_240],
    #              "y bin" : [ymin_cutxpi4_240,ymax_cutxpi4_240],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi4_240,yieldmax_cutxpi4_240],
    #              "xL" : [xLmin_cutxpi4_240,xLmax_cutxpi4_240],
    #              "xpi": [xpimin_cutxpi4_240,xpimax_cutxpi4_240]
    #          },
    #          "0.5" :{
    #              "Q2" : [Q2min_cutxpi5_240,Q2max_cutxpi5_240],
    #              "y bin" : [ymin_cutxpi5_240,ymax_cutxpi5_240],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi5_240,yieldmax_cutxpi5_240],
    #              "xL" : [xLmin_cutxpi5_240,xLmax_cutxpi5_240],
    #              "xpi": [xpimin_cutxpi5_240,xpimax_cutxpi5_240]
    #          },
    #          "0.6" :{
    #              "Q2" : [Q2min_cutxpi6_240,Q2max_cutxpi6_240],
    #              "y bin" : [ymin_cutxpi6_240,ymax_cutxpi6_240],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi6_240,yieldmax_cutxpi6_240],
    #              "xL" : [xLmin_cutxpi6_240,xLmax_cutxpi6_240],
    #              "xpi": [xpimin_cutxpi6_240,xpimax_cutxpi6_240]
    #          },
    #          "0.7" :{
    #              "Q2" : [Q2min_cutxpi7_240,Q2max_cutxpi7_240],
    #              "y bin" : [ymin_cutxpi7_240,ymax_cutxpi7_240],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi7_240,yieldmax_cutxpi7_240],
    #              "xL" : [xLmin_cutxpi7_240,xLmax_cutxpi7_240],
    #              "xpi": [xpimin_cutxpi7_240,xpimax_cutxpi7_240]
    #          }
    #         },
    #         "480" : 
    #          {"0.2" :{
    #              "Q2" : [Q2min_cutxpi2_480,Q2max_cutxpi2_480],
    #              "y bin" : [ymin_cutxpi2_480,ymax_cutxpi2_480],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi2_480,yieldmax_cutxpi2_480],
    #              "xL" : [xLmin_cutxpi2_480,xLmax_cutxpi2_480],
    #              "xpi": [xpimin_cutxpi2_480,xpimax_cutxpi2_480]
    #          },
    #          "0.3" :{
    #              "Q2" : [Q2min_cutxpi3_480,Q2max_cutxpi3_480],
    #              "y bin" : [ymin_cutxpi3_480,ymax_cutxpi3_480],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi3_480,yieldmax_cutxpi3_480],
    #              "xL" : [xLmin_cutxpi3_480,xLmax_cutxpi3_480],
    #              "xpi": [xpimin_cutxpi3_480,xpimax_cutxpi3_480]
    #          },
    #          "0.4" :{
    #              "Q2" : [Q2min_cutxpi4_480,Q2max_cutxpi4_480],
    #              "y bin" : [ymin_cutxpi4_480,ymax_cutxpi4_480],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi4_480,yieldmax_cutxpi4_480],
    #              "xL" : [xLmin_cutxpi4_480,xLmax_cutxpi4_480],
    #              "xpi": [xpimin_cutxpi4_480,xpimax_cutxpi4_480]
    #          },
    #          "0.5" :{
    #              "Q2" : [Q2min_cutxpi5_480,Q2max_cutxpi5_480],
    #              "y bin" : [ymin_cutxpi5_480,ymax_cutxpi5_480],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi5_480,yieldmax_cutxpi5_480],
    #              "xL" : [xLmin_cutxpi5_480,xLmax_cutxpi5_480],
    #              "xpi": [xpimin_cutxpi5_480,xpimax_cutxpi5_480]
    #          },
    #          "0.6" :{
    #              "Q2" : [Q2min_cutxpi6_480,Q2max_cutxpi6_480],
    #              "y bin" : [ymin_cutxpi6_480,ymax_cutxpi6_480],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi6_480,yieldmax_cutxpi6_480],
    #              "xL" : [xLmin_cutxpi6_480,xLmax_cutxpi6_480],
    #              "xpi": [xpimin_cutxpi6_480,xpimax_cutxpi6_480]
    #          },
    #          "0.7" :{
    #              "Q2" : [Q2min_cutxpi7_480,Q2max_cutxpi7_480],
    #              "y bin" : [ymin_cutxpi7_480,ymax_cutxpi7_480],
    #              "Yield (100 fb^-1)" : [yieldmin_cutxpi7_480,yieldmax_cutxpi7_480],
    #              "xL" : [xLmin_cutxpi7_480,xLmax_cutxpi7_480],
    #              "xpi": [xpimin_cutxpi7_480,xpimax_cutxpi7_480]
    #          }
    #     }
    # }
    
    f = plt.figure(figsize=(11.69,8.27))
    plt.scatter(c.applyCuts(y,cutxpi001_15),c.applyCuts(lumi,cutxpi001_15),label='0.001')
    # plt.scatter(c.applyCuts(y,cutxpi01_15),c.applyCuts(lumi,cutxpi01_15),label='0.01')
    plt.scatter(c.applyCuts(Q2/(5400*xpi),cutxpi001_15),c.applyCuts(lumi,cutxpi001_15),label='0.001, Q2/sx')
    plt.title('15', fontsize =20)
    plt.xlabel('y')
    plt.ylabel('yield')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    
    f = plt.figure(figsize=(11.69,8.27))
    plt.scatter(c.applyCuts(y,cutxpi1_15),c.applyCuts(lumi,cutxpi1_15),label='0.1')
    plt.scatter(c.applyCuts(Q2/(5400*xpi),cutxpi1_15),c.applyCuts(lumi,cutxpi1_15),label='0.1, Q2/sx')
    # plt.scatter(c.applyCuts(y,cutxpi2_15),c.applyCuts(lumi,cutxpi2_15),label='0.2')
    # plt.scatter(c.applyCuts(y,cutxpi3_15),c.applyCuts(lumi,cutxpi3_15),label='0.3')
    # plt.scatter(c.applyCuts(y,cutxpi4_15),c.applyCuts(lumi,cutxpi4_15),label='0.4')
    # plt.scatter(c.applyCuts(y,cutxpi5_15),c.applyCuts(lumi,cutxpi5_15),label='0.5')
    # plt.scatter(c.applyCuts(y,cutxpi6_15),c.applyCuts(lumi,cutxpi6_15),label='0.6')
    # plt.scatter(c.applyCuts(y,cutxpi7_15),c.applyCuts(lumi,cutxpi7_15),label='0.7')
    plt.title('15', fontsize =20)
    plt.xlabel('y')
    plt.ylabel('yield')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    print("\n\n",np.average(c.applyCuts(xL,cutxpi1_15)))
    print(np.average(c.applyCuts(xL,cutxpi2_15)))
    print(np.average(c.applyCuts(xL,cutxpi3_15)))
    print(np.average(c.applyCuts(xL,cutxpi4_15)))
    print(np.average(c.applyCuts(xL,cutxpi5_15)))
    print(np.average(c.applyCuts(xL,cutxpi6_15)))
    print(np.average(c.applyCuts(xL,cutxpi7_15)),"\n\n")
    
    f = plt.figure(figsize=(11.69,8.27))
    plt.scatter(c.applyCuts(y,cutxpi1_30),c.applyCuts(lumi,cutxpi1_30),label='0.1')
    plt.scatter(c.applyCuts(y,cutxpi2_30),c.applyCuts(lumi,cutxpi2_30),label='0.2')
    plt.scatter(c.applyCuts(y,cutxpi3_30),c.applyCuts(lumi,cutxpi3_30),label='0.3')
    plt.scatter(c.applyCuts(y,cutxpi4_30),c.applyCuts(lumi,cutxpi4_30),label='0.4')
    plt.scatter(c.applyCuts(y,cutxpi5_30),c.applyCuts(lumi,cutxpi5_30),label='0.5')
    plt.scatter(c.applyCuts(y,cutxpi6_30),c.applyCuts(lumi,cutxpi6_30),label='0.6')
    plt.scatter(c.applyCuts(y,cutxpi7_30),c.applyCuts(lumi,cutxpi7_30),label='0.7')
    plt.title('30', fontsize =20)
    plt.xlabel('y')
    plt.ylabel('yield')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    print("\n\n",np.average(c.applyCuts(xL,cutxpi1_30)))
    print(np.average(c.applyCuts(xL,cutxpi2_30)))
    print(np.average(c.applyCuts(xL,cutxpi3_30)))
    print(np.average(c.applyCuts(xL,cutxpi4_30)))
    print(np.average(c.applyCuts(xL,cutxpi5_30)))
    print(np.average(c.applyCuts(xL,cutxpi6_30)))
    print(np.average(c.applyCuts(xL,cutxpi7_30)),"\n\n")
    
    f = plt.figure(figsize=(11.69,8.27))
    plt.scatter(c.applyCuts(y,cutxpi1_60),c.applyCuts(lumi,cutxpi1_60),label='0.1# ')
    plt.scatter(c.applyCuts(y,cutxpi2_60),c.applyCuts(lumi,cutxpi2_60),label='0.2')
    plt.scatter(c.applyCuts(y,cutxpi3_60),c.applyCuts(lumi,cutxpi3_60),label='0.3')
    plt.scatter(c.applyCuts(y,cutxpi4_60),c.applyCuts(lumi,cutxpi4_60),label='0.4')
    plt.scatter(c.applyCuts(y,cutxpi5_60),c.applyCuts(lumi,cutxpi5_60),label='0.5')
    plt.scatter(c.applyCuts(y,cutxpi6_60),c.applyCuts(lumi,cutxpi6_60),label='0.6')
    plt.scatter(c.applyCuts(y,cutxpi7_60),c.applyCuts(lumi,cutxpi7_60),label='0.7')
    plt.title('60', fontsize =20)
    plt.xlabel('y')
    plt.ylabel('yield')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    print("\n\n",np.average(c.applyCuts(xL,cutxpi1_60)))
    print(np.average(c.applyCuts(xL,cutxpi2_60)))
    print(np.average(c.applyCuts(xL,cutxpi3_60)))
    print(np.average(c.applyCuts(xL,cutxpi4_60)))
    print(np.average(c.applyCuts(xL,cutxpi5_60)))
    print(np.average(c.applyCuts(xL,cutxpi6_60)))
    print(np.average(c.applyCuts(xL,cutxpi7_60)),"\n\n")
    
    f = plt.figure(figsize=(11.69,8.27))
    plt.scatter(c.applyCuts(y,cutxpi1_120),c.applyCuts(lumi,cutxpi1_120),label='0.1')
    plt.scatter(c.applyCuts(y,cutxpi2_120),c.applyCuts(lumi,cutxpi2_120),label='0.2')
    plt.scatter(c.applyCuts(y,cutxpi3_120),c.applyCuts(lumi,cutxpi3_120),label='0.3')
    plt.scatter(c.applyCuts(y,cutxpi4_120),c.applyCuts(lumi,cutxpi4_120),label='0.4')
    plt.scatter(c.applyCuts(y,cutxpi5_120),c.applyCuts(lumi,cutxpi5_120),label='0.5')
    plt.scatter(c.applyCuts(y,cutxpi6_120),c.applyCuts(lumi,cutxpi6_120),label='0.6')
    plt.scatter(c.applyCuts(y,cutxpi7_120),c.applyCuts(lumi,cutxpi7_120),label='0.7')
    plt.title('120', fontsize =20)
    plt.xlabel('y')
    plt.ylabel('yield')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    print("\n\n",np.average(c.applyCuts(xL,cutxpi1_120)))
    print(np.average(c.applyCuts(xL,cutxpi2_120)))
    print(np.average(c.applyCuts(xL,cutxpi3_120)))
    print(np.average(c.applyCuts(xL,cutxpi4_120)))
    print(np.average(c.applyCuts(xL,cutxpi5_120)))
    print(np.average(c.applyCuts(xL,cutxpi6_120)))
    print(np.average(c.applyCuts(xL,cutxpi7_120)),"\n\n")
    
    f = plt.figure(figsize=(11.69,8.27))
    plt.scatter(c.applyCuts(y,cutxpi1_240),c.applyCuts(lumi,cutxpi1_240),label='0.1')
    plt.scatter(c.applyCuts(y,cutxpi2_240),c.applyCuts(lumi,cutxpi2_240),label='0.2')
    plt.scatter(c.applyCuts(y,cutxpi3_240),c.applyCuts(lumi,cutxpi3_240),label='0.3')
    plt.scatter(c.applyCuts(y,cutxpi4_240),c.applyCuts(lumi,cutxpi4_240),label='0.4')
    plt.scatter(c.applyCuts(y,cutxpi5_240),c.applyCuts(lumi,cutxpi5_240),label='0.5')
    plt.scatter(c.applyCuts(y,cutxpi6_240),c.applyCuts(lumi,cutxpi6_240),label='0.6')
    plt.scatter(c.applyCuts(y,cutxpi7_240),c.applyCuts(lumi,cutxpi7_240),label='0.7')
    plt.title('240', fontsize =20)
    plt.xlabel('y')
    plt.ylabel('yield')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    print("\n\n",np.average(c.applyCuts(xL,cutxpi1_240)))
    print(np.average(c.applyCuts(xL,cutxpi2_240)))
    print(np.average(c.applyCuts(xL,cutxpi3_240)))
    print(np.average(c.applyCuts(xL,cutxpi4_240)))
    print(np.average(c.applyCuts(xL,cutxpi5_240)))
    print(np.average(c.applyCuts(xL,cutxpi6_240)))
    print(np.average(c.applyCuts(xL,cutxpi7_240)),"\n\n")

    f = plt.figure(figsize=(11.69,8.27))
    plt.scatter(c.applyCuts(y,cutxpi1_480),c.applyCuts(lumi,cutxpi1_480),label='0.1')
    plt.scatter(c.applyCuts(y,cutxpi2_480),c.applyCuts(lumi,cutxpi2_480),label='0.2')
    plt.scatter(c.applyCuts(y,cutxpi3_480),c.applyCuts(lumi,cutxpi3_480),label='0.3')
    plt.scatter(c.applyCuts(y,cutxpi4_480),c.applyCuts(lumi,cutxpi4_480),label='0.4')
    plt.scatter(c.applyCuts(y,cutxpi5_480),c.applyCuts(lumi,cutxpi5_480),label='0.5')
    plt.scatter(c.applyCuts(y,cutxpi6_480),c.applyCuts(lumi,cutxpi6_480),label='0.6')
    plt.scatter(c.applyCuts(y,cutxpi7_480),c.applyCuts(lumi,cutxpi7_480),label='0.7')
    plt.title('480', fontsize =20)
    plt.xlabel('y')
    plt.ylabel('yield')
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    print("\n\n",np.average(c.applyCuts(xL,cutxpi1_480)))
    print(np.average(c.applyCuts(xL,cutxpi2_480)))
    print(np.average(c.applyCuts(xL,cutxpi3_480)))
    print(np.average(c.applyCuts(xL,cutxpi4_480)))
    print(np.average(c.applyCuts(xL,cutxpi5_480)))
    print(np.average(c.applyCuts(xL,cutxpi6_480)))
    print(np.average(c.applyCuts(xL,cutxpi7_480)),"\n\n")

    ################################################
    # print(np.average(c.applyCuts(Q2/(5400*xBj),cutxpi001_15)))
    # print(np.average(c.applyCuts(Q2/(5400*xBj),cutxpi01_15)))
    # print(np.average(c.applyCuts(Q2/(5400*xBj),cutxpi01_30)))
    # print(np.average(c.applyCuts(Q2/(5400*xBj),cutxpi01_60)))

    # f = plt.figure(figsize=(11.69,8.27))
    # plt.scatter(c.applyCuts(Q2/(5400*xBj),cutxpi1_15),c.applyCuts(y,cutxpi1_15),label='0.1')
    # plt.scatter(c.applyCuts(Q2/(5400*xBj),cutxpi2_15),c.applyCuts(y,cutxpi2_15),label='0.2')
    # plt.scatter(c.applyCuts(Q2/(5400*xBj),cutxpi3_15),c.applyCuts(y,cutxpi3_15),label='0.3')
    # plt.xlabel('y=Q2/sx')
    # plt.ylabel('y')
    # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # f = plt.figure(figsize=(11.69,8.27))
    # plt.scatter(c.applyCuts(1-(xBj/xpi),cutxpi1_15),c.applyCuts(xL,cutxpi1_15),label='0.1')
    # plt.scatter(c.applyCuts(1-(xBj/xpi),cutxpi2_15),c.applyCuts(xL,cutxpi2_15),label='0.2')
    # plt.scatter(c.applyCuts(1-(xBj/xpi),cutxpi3_15),c.applyCuts(xL,cutxpi3_15),label='0.3')
    # plt.xlabel('xL=1-(xBj/xpi)')
    # plt.ylabel('xL')
    # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # f = plt.figure(figsize=(11.69,8.27))
    # plt.scatter(c.applyCuts(Q2/(5400*y),cutxpi1_15),c.applyCuts(xBj,cutxpi1_15),label='0.1')
    # plt.scatter(c.applyCuts(Q2/(5400*y),cutxpi2_15),c.applyCuts(xBj,cutxpi2_15),label='0.2')
    # plt.scatter(c.applyCuts(Q2/(5400*y),cutxpi3_15),c.applyCuts(xBj,cutxpi3_15),label='0.3')
    # plt.xlabel('xBj=Q2/(5400*y)')
    # plt.ylabel('xBj')
    # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    # f = plt.figure(figsize=(11.69,8.27))
    # plt.scatter(c.applyCuts(xBj/(1-xL),cutxpi1_15),c.applyCuts(xpi,cutxpi1_15),label='0.1')
    # plt.scatter(c.applyCuts(xBj/(1-xL),cutxpi2_15),c.applyCuts(xpi,cutxpi2_15),label='0.2')
    # plt.scatter(c.applyCuts(xBj/(1-xL),cutxpi3_15),c.applyCuts(xpi,cutxpi3_15),label='0.3')
    # plt.xlabel('xpi=xBj/(1-xL)')
    # plt.ylabel('xpi')
    # plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    
    tmp = []
    frames = []

    for key,val in table_dict.items():
        tmp.append(key)
        frames.append(pd.DataFrame.from_dict(val,orient='index'))    
    theory_df_table = pd.concat(frames,keys=tmp)        

    # theory_df_table = theory_df_table.T
    
    print(theory_df_table)
    return theory_df_table
    

    
def main() :

    # phase_space()
    # fpivxpi_Plot()
    # # fpivxpi_Plot_nolog()
    # fpivt_Plot()
    # # fpivt_xbin_Plot()
    # # fpivtPrime_Plot()
    # # fpivtPrime_xbin_Plot()
    # # sigmavxpi_Plot()
    # # plot3D()
    # plt.show()

    df_table = theory_table()

    plt.show()
    #render dataframe as html
    html_table = df_table.to_html()
    csv_table = df_table.to_csv()

    #write html to file
    text_file = open("theory_table.html", "w")
    text_file.write(html_table)
    text_file.close()

    text_file = open("theory_table.csv", "w")
    text_file.write(csv_table)
    text_file.close()

    
if __name__=='__main__': main()
