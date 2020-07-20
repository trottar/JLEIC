#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2020-07-19 23:54:00 trottar"
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
y = tree.array("y")
fpi = tree.array("fpi")
f2N = tree.array("f2N")
xpi = tree.array("xpi")
ypi = tree.array("ypi")
tpi = tree.array("tpi")
escat = tree.array("EScatRest")
pprz_inc = tree.array("pprz_inc")
ppiz_Lab = tree.array("ppiz_Lab")

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
    xpitmp = '{"xpicut%i" : ((%0.5f <= xpi) & (xpi <= %0.5f))}' % (i,xpiarray[i]-binx/20,xpiarray[i]+binx/20)
    print('{"xpicut%i" : ((%0.5f <= xpi) & (xpi <= %0.5f))}' % (i,xpiarray[i]-binx/20,xpiarray[i]+binx/20))
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
    Q2tmp = '{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,Q2array[i]-binQ2/20,Q2array[i]+binQ2/20)
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
xL = pNz_Lab/pprz_inc # Frac of proton momentum
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
    Q2tmp = '{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,Q2array[i]-binQ2/2,Q2array[i]+binQ2/2)
    print('{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,Q2array[i]-binQ2/2,Q2array[i]+binQ2/2))
    cutDict.update(eval(Q2tmp))

xarray = np.arange(0.05,1.0,0.1).tolist()
for i,x in enumerate(xarray):
    xtmp = '{"xcut%i" : ((%0.2f <= xpi) & (xpi <= %0.2f))}' % (i,xarray[i]-0.05,xarray[i]+0.05)
    print('{"xcut%i" : ((%0.2f <= xpi) & (xpi <= %0.2f))}' % (i,xarray[i]-0.05,xarray[i]+0.05))
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
cut120 = ["Q2cut3","tcut","ycut"]
cut240 = ["Q2cut4","tcut","ycut"]
cut480 = ["Q2cut5","tcut","ycut"]
cut1000 = ["Q2cut6","tcut","ycut"]

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

    plt.ylabel('d$\sigma_{TDIS}$ ($nb/GeV^{2}$)')

    ax = f.add_subplot(422)
    xpiscat2 = ax.errorbar(c.applyCuts(xpi,cut15),c.applyCuts(tot_sigma,cut15),yerr=np.sqrt(c.applyCuts(tot_sigma,cut15))/150,fmt='o',label='$Q^2$=15 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=20 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    plt.title('d$\sigma_{TDIS}$ vs $x_\pi$', fontsize =20)
    
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

    print(type(lumi))
    print(type(xpi))
    
    ax = f.add_subplot(421)
    xpiscat1 = ax.errorbar(c.applyCuts(xpi,cut7),c.applyCuts(fpi,cut7),yerr=np.sqrt(c.applyCuts(lumi,cut7))/c.applyCuts(lumi,cut7),fmt='.',label='$Q^2$=7 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.plot([0.001,0.01,0.1],[0.55,0.25,0.20], label="GRV fit",color="y")
    plt.xscale('log')
    plt.xlim(1e-4,1.)
    plt.ylim(0.,0.5)
    ax.text(0.60, 0.85, '$Q^2$=7 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    plt.ylabel('$F^{\pi}_{2}$')

    ax = f.add_subplot(422)
    error1 = np.sqrt(c.applyCuts(lumi,cut15))/c.applyCuts(lumi,cut15)
    xpiscat2 = ax.errorbar(c.applyCuts(xpi,cut15),c.applyCuts(fpi,cut15),yerr=np.sqrt(c.applyCuts(lumi,cut15))/c.applyCuts(lumi,cut15),fmt='.',label='$Q^2$=15 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.plot([0.001,0.01,0.1],[0.80,0.30,0.20], label="GRV fit",color="y")
    plt.xscale('log')
    plt.xlim(1e-4,1.)
    plt.ylim(0.,0.5)
    ax.text(0.60, 0.85, '$Q^2$=15 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    plt.title('$F^{\pi}_{2}$ vs $x_\pi$', fontsize =20)
    
    ax = f.add_subplot(423)
    xpiscat3 = ax.errorbar(c.applyCuts(xpi,cut30),c.applyCuts(fpi,cut30),yerr=np.sqrt(c.applyCuts(lumi,cut30))/c.applyCuts(lumi,cut30),fmt='.',label='$Q^2$=30 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.plot([0.001,0.01,0.1],[0.85,0.45,0.20], label="GRV fit",color="y")
    plt.xscale('log')
    plt.xlim(1e-4,1.)
    plt.ylim(0.,0.5)
    ax.text(0.60, 0.85, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    
    ax = f.add_subplot(424)
    xpiscat4 = ax.errorbar(c.applyCuts(xpi,cut60),c.applyCuts(fpi,cut60),yerr=np.sqrt(c.applyCuts(lumi,cut60))/c.applyCuts(lumi,cut60),fmt='.',label='$Q^2$=60 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.plot([0.001,0.01,0.1],[1.2,0.45,0.25], label="GRV fit",color="y")
    plt.xscale('log')
    plt.xlim(1e-4,1.)
    plt.ylim(0.,0.5)
    ax.text(0.60, 0.85, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    ax = f.add_subplot(425)
    xpiscat5 = ax.errorbar(c.applyCuts(xpi,cut120),c.applyCuts(fpi,cut120),yerr=np.sqrt(c.applyCuts(lumi,cut120))/c.applyCuts(lumi,cut120),fmt='.',label='$Q^2$=120 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.5,0.25,0.15], label="GRV fit",color="y")
    plt.xscale('log')
    plt.xlim(1e-4,1.)
    plt.ylim(0.,0.5)
    ax.text(0.60, 0.85, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))

    ax = f.add_subplot(426)
    xpiscat6 = ax.errorbar(c.applyCuts(xpi,cut240),c.applyCuts(fpi,cut240),yerr=np.sqrt(c.applyCuts(lumi,cut240))/c.applyCuts(lumi,cut240),fmt='.',label='$Q^2$=240 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit",color="y")
    plt.xscale('log')
    plt.xlim(1e-4,1.)
    plt.ylim(0.,0.5)
    ax.text(0.60, 0.85, '$Q^2$=240 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(427)
    xpiscat7 = ax.errorbar(c.applyCuts(xpi,cut480),c.applyCuts(fpi,cut480),yerr=np.sqrt(c.applyCuts(lumi,cut480))/c.applyCuts(lumi,cut480),fmt='.',label='$Q^2$=480 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit",color="y")
    plt.xscale('log')
    plt.xlim(1e-4,1.)
    plt.ylim(0.,0.5)
    ax.text(0.60, 0.85, '$Q^2$=480 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    # ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))

    ax = f.add_subplot(428)
    xpiscat8 = ax.errorbar(c.applyCuts(xpi,cut1000),c.applyCuts(fpi,cut1000),yerr=np.sqrt(c.applyCuts(lumi,cut1000))/c.applyCuts(lumi,cut1000),fmt='.',label='$Q^2$=1000 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit",color="y")
    plt.xscale('log')
    plt.xlim(1e-4,1.)
    plt.ylim(0.,0.5)
    ax.text(0.60, 0.85, '$Q^2$=1000 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    
    plt.xlabel('log($x_\pi$)')
    plt.tight_layout()

    plt.subplots_adjust(hspace=0.0,wspace=0.0)

    f.savefig('OUTPUTS/fpivxpi.png')

def fpivxpi_Plot_nolog():
    
    f = plt.figure(figsize=(11.69,8.27))
    
    ax = f.add_subplot(421)
    xpiscat1 = ax.errorbar(c.applyCuts(xpi,cut7),c.applyCuts(fpi,cut7),yerr=np.sqrt(c.applyCuts(lumi,cut7))/c.applyCuts(lumi,cut7),fmt='.',label='$Q^2$=7 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.plot([0.001,0.01,0.1],[0.55,0.25,0.20], label="GRV fit",color="y")
    plt.xlim(0.,1.)
    plt.ylim(0.,0.5)
    ax.text(0.60, 0.85, '$Q^2$=7 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    plt.ylabel('$F^{\pi}_{2}$')

    ax = f.add_subplot(422)
    xpiscat2 = ax.errorbar(c.applyCuts(xpi,cut15),c.applyCuts(fpi,cut15),yerr=np.sqrt(c.applyCuts(lumi,cut15))/c.applyCuts(lumi,cut15),fmt='.',label='$Q^2$=15 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.plot([0.001,0.01,0.1],[0.80,0.30,0.20], label="GRV fit",color="y")
    plt.xlim(0.,1.)
    plt.ylim(0.,0.5)
    ax.text(0.60, 0.85, '$Q^2$=15 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    plt.title('$F^{\pi}_{2}$ vs $x_\pi$', fontsize =20)
    
    ax = f.add_subplot(423)
    xpiscat3 = ax.errorbar(c.applyCuts(xpi,cut30),c.applyCuts(fpi,cut30),yerr=np.sqrt(c.applyCuts(lumi,cut30))/c.applyCuts(lumi,cut30),fmt='.',label='$Q^2$=30 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.plot([0.001,0.01,0.1],[0.85,0.45,0.20], label="GRV fit",color="y")
    plt.xlim(0.,1.)
    plt.ylim(0.,0.5)
    ax.text(0.60, 0.85, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    
    ax = f.add_subplot(424)
    xpiscat4 = ax.errorbar(c.applyCuts(xpi,cut60),c.applyCuts(fpi,cut60),yerr=np.sqrt(c.applyCuts(lumi,cut60))/c.applyCuts(lumi,cut60),fmt='.',label='$Q^2$=60 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.plot([0.001,0.01,0.1],[1.2,0.45,0.25], label="GRV fit",color="y")
    plt.xlim(0.,1.)
    plt.ylim(0.,0.5)
    ax.text(0.60, 0.85, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    ax = f.add_subplot(425)
    xpiscat5 = ax.errorbar(c.applyCuts(xpi,cut120),c.applyCuts(fpi,cut120),yerr=np.sqrt(c.applyCuts(lumi,cut120))/c.applyCuts(lumi,cut120),fmt='.',label='$Q^2$=120 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.5,0.25,0.15], label="GRV fit",color="y")
    plt.xlim(0.,1.)
    plt.ylim(0.,0.5)
    ax.text(0.60, 0.85, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))

    ax = f.add_subplot(426)
    xpiscat6 = ax.errorbar(c.applyCuts(xpi,cut240),c.applyCuts(fpi,cut240),yerr=np.sqrt(c.applyCuts(lumi,cut240))/c.applyCuts(lumi,cut240),fmt='.',label='$Q^2$=240 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit",color="y")
    plt.xlim(0.,1.)
    plt.ylim(0.,0.5)
    ax.text(0.60, 0.85, '$Q^2$=240 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(427)
    xpiscat7 = ax.errorbar(c.applyCuts(xpi,cut480),c.applyCuts(fpi,cut480),yerr=np.sqrt(c.applyCuts(lumi,cut480))/c.applyCuts(lumi,cut480),fmt='.',label='$Q^2$=480 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit",color="y")
    plt.xlim(0.,1.)
    plt.ylim(0.,0.5)
    ax.text(0.60, 0.85, '$Q^2$=480 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower',nbins=5))

    ax = f.add_subplot(428)
    xpiscat8 = ax.errorbar(c.applyCuts(xpi,cut1000),c.applyCuts(fpi,cut1000),yerr=np.sqrt(c.applyCuts(lumi,cut1000))/c.applyCuts(lumi,cut1000),fmt='.',label='$Q^2$=1000 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit",color="y")
    plt.xlim(0.,1.)
    plt.ylim(0.,0.5)
    ax.text(0.60, 0.85, '$Q^2$=1000 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower',nbins=5))
    
    plt.xlabel('$x_\pi$')
    plt.tight_layout()

    plt.subplots_adjust(hspace=0.0,wspace=0.0)

    f.savefig('OUTPUTS/fpivxpi_nolog.png')

    
def fpivt_Plot():

    f = plt.figure(figsize=(11.69,8.27))
    
    ax = f.add_subplot(421)
    tscat1 = ax.errorbar(c.applyCuts(t,cut7),c.applyCuts(fpi,cut7),yerr=np.sqrt(c.applyCuts(lumi,cut7))/c.applyCuts(lumi,cut7),fmt='.',label='$Q^2$=7 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=7 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    plt.ylabel('log($F^{\pi}_{2}$)')

    ax = f.add_subplot(422)
    tscat2 = ax.errorbar(c.applyCuts(t,cut15),c.applyCuts(fpi,cut15),yerr=np.sqrt(c.applyCuts(lumi,cut15))/c.applyCuts(lumi,cut15),fmt='.',label='$Q^2$=15 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=15 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_formatter(plt.NullFormatter())

    plt.title('log($F^{\pi}_{2}$) vs -t', fontsize =20)
    
    ax = f.add_subplot(423)
    tscat3 = ax.errorbar(c.applyCuts(t,cut30),c.applyCuts(fpi,cut30),yerr=np.sqrt(c.applyCuts(lumi,cut30))/c.applyCuts(lumi,cut30),fmt='.',label='$Q^2$=30 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    
    ax = f.add_subplot(424)
    tscat4 = ax.errorbar(c.applyCuts(t,cut60),c.applyCuts(fpi,cut60),yerr=np.sqrt(c.applyCuts(lumi,cut60))/c.applyCuts(lumi,cut60),fmt='.',label='$Q^2$=60 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    
    ax = f.add_subplot(425)
    tscat5 = ax.errorbar(c.applyCuts(t,cut120),c.applyCuts(fpi,cut120),yerr=np.sqrt(c.applyCuts(lumi,cut120))/c.applyCuts(lumi,cut120),fmt='.',label='$Q^2$=120 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    ax = f.add_subplot(426)
    tscat6 = ax.errorbar(c.applyCuts(t,cut240),c.applyCuts(fpi,cut240),yerr=np.sqrt(c.applyCuts(lumi,cut240))/c.applyCuts(lumi,cut240),fmt='.',label='$Q^2$=240 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=240 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(427)
    tscat7 = ax.errorbar(c.applyCuts(t,cut480),c.applyCuts(fpi,cut480),yerr=np.sqrt(c.applyCuts(lumi,cut480))/c.applyCuts(lumi,cut480),fmt='.',label='$Q^2$=480 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=480 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))

    ax = f.add_subplot(428)
    tscat8 = ax.errorbar(c.applyCuts(t,cut1000),c.applyCuts(fpi,cut1000),yerr=np.sqrt(c.applyCuts(lumi,cut1000))/c.applyCuts(lumi,cut1000),fmt='.',label='$Q^2$=1000 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=1000 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    plt.xlabel('-t')
    plt.tight_layout()

    plt.subplots_adjust(hspace=0.0,wspace=0.0)

    f.savefig('OUTPUTS/fpivt.png')

def fpivt_xbin_Plot():

    f = plt.figure(figsize=(11.69,8.27))
    
    ax = f.add_subplot(221)
    tscat0 = ax.errorbar(c.applyCuts(t,cutx0_15),c.applyCuts(fpi,cutx0_15),yerr=np.sqrt(c.applyCuts(lumi,cutx0_15))/c.applyCuts(lumi,cutx0_15),fmt='.',label='$x_{\pi}$=(0.0,0.1)',ecolor='black',capsize=2, capthick=2,color='b')
    tscat1 = ax.errorbar(c.applyCuts(t,cutx1_15),c.applyCuts(fpi,cutx1_15),yerr=np.sqrt(c.applyCuts(lumi,cutx1_15))/c.applyCuts(lumi,cutx1_15),fmt='.',label='$x_{\pi}$=(0.1,0.2)',ecolor='black',capsize=2, capthick=2,color='g')
    tscat2 = ax.errorbar(c.applyCuts(t,cutx2_15),c.applyCuts(fpi,cutx2_15),yerr=np.sqrt(c.applyCuts(lumi,cutx2_15))/c.applyCuts(lumi,cutx2_15),fmt='.',label='$x_{\pi}$=(0.2,0.3)',ecolor='black',capsize=2, capthick=2,color='r')
    tscat3 = ax.errorbar(c.applyCuts(t,cutx3_15),c.applyCuts(fpi,cutx3_15),yerr=np.sqrt(c.applyCuts(lumi,cutx3_15))/c.applyCuts(lumi,cutx3_15),fmt='.',label='$x_{\pi}$=(0.3,0.4)',ecolor='black',capsize=2, capthick=2,color='c')
    tscat4 = ax.errorbar(c.applyCuts(t,cutx4_15),c.applyCuts(fpi,cutx4_15),yerr=np.sqrt(c.applyCuts(lumi,cutx4_15))/c.applyCuts(lumi,cutx4_15),fmt='.',label='$x_{\pi}$=(0.4,0.5)',ecolor='black',capsize=2, capthick=2,color='m')
    tscat5 = ax.errorbar(c.applyCuts(t,cutx5_15),c.applyCuts(fpi,cutx5_15),yerr=np.sqrt(c.applyCuts(lumi,cutx5_15))/c.applyCuts(lumi,cutx5_15),fmt='.',label='$x_{\pi}$=(0.5,0.6)',ecolor='black',capsize=2, capthick=2,color='y')
    tscat6 = ax.errorbar(c.applyCuts(t,cutx6_15),c.applyCuts(fpi,cutx6_15),yerr=np.sqrt(c.applyCuts(lumi,cutx6_15))/c.applyCuts(lumi,cutx6_15),fmt='.',label='$x_{\pi}$=(0.6,0.7)',ecolor='black',capsize=2, capthick=2,color='purple')
    tscat7 = ax.errorbar(c.applyCuts(t,cutx7_15),c.applyCuts(fpi,cutx7_15),yerr=np.sqrt(c.applyCuts(lumi,cutx7_15))/c.applyCuts(lumi,cutx7_15),fmt='.',label='$x_{\pi}$=(0.7,0.8)',ecolor='black',capsize=2, capthick=2,color='aqua')
    tscat8 = ax.errorbar(c.applyCuts(t,cutx8_15),c.applyCuts(fpi,cutx8_15),yerr=np.sqrt(c.applyCuts(lumi,cutx8_15))/c.applyCuts(lumi,cutx8_15),fmt='.',label='$x_{\pi}$=(0.8,0.9)',ecolor='black',capsize=2, capthick=2,color='gray')
    tscat9 = ax.errorbar(c.applyCuts(t,cutx9_15),c.applyCuts(fpi,cutx9_15),yerr=np.sqrt(c.applyCuts(lumi,cutx9_15))/c.applyCuts(lumi,cutx9_15),fmt='.',label='$x_{\pi}$=(0.9,1.0)',ecolor='black',capsize=2, capthick=2,color='pink')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(0.,0.3)
    ax.text(0.65, 0.85, '$Q^2$=15 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='lower'))

    plt.ylabel('$F^{\pi}_{2}$')
    
    ax = f.add_subplot(222)
    tscat0 = ax.errorbar(c.applyCuts(t,cutx0_30),c.applyCuts(fpi,cutx0_30),yerr=np.sqrt(c.applyCuts(lumi,cutx0_30))/c.applyCuts(lumi,cutx0_30),fmt='.',label='$x_{\pi}$=(0.0,0.1)',ecolor='black',capsize=2, capthick=2,color='b')
    tscat1 = ax.errorbar(c.applyCuts(t,cutx1_30),c.applyCuts(fpi,cutx1_30),yerr=np.sqrt(c.applyCuts(lumi,cutx1_30))/c.applyCuts(lumi,cutx1_30),fmt='.',label='$x_{\pi}$=(0.1,0.2)',ecolor='black',capsize=2, capthick=2,color='g')
    tscat2 = ax.errorbar(c.applyCuts(t,cutx2_30),c.applyCuts(fpi,cutx2_30),yerr=np.sqrt(c.applyCuts(lumi,cutx2_30))/c.applyCuts(lumi,cutx2_30),fmt='.',label='$x_{\pi}$=(0.2,0.3)',ecolor='black',capsize=2, capthick=2,color='r')
    tscat3 = ax.errorbar(c.applyCuts(t,cutx3_30),c.applyCuts(fpi,cutx3_30),yerr=np.sqrt(c.applyCuts(lumi,cutx3_30))/c.applyCuts(lumi,cutx3_30),fmt='.',label='$x_{\pi}$=(0.3,0.4)',ecolor='black',capsize=2, capthick=2,color='c')
    tscat4 = ax.errorbar(c.applyCuts(t,cutx4_30),c.applyCuts(fpi,cutx4_30),yerr=np.sqrt(c.applyCuts(lumi,cutx4_30))/c.applyCuts(lumi,cutx4_30),fmt='.',label='$x_{\pi}$=(0.4,0.5)',ecolor='black',capsize=2, capthick=2,color='m')
    tscat5 = ax.errorbar(c.applyCuts(t,cutx5_30),c.applyCuts(fpi,cutx5_30),yerr=np.sqrt(c.applyCuts(lumi,cutx5_30))/c.applyCuts(lumi,cutx5_30),fmt='.',label='$x_{\pi}$=(0.5,0.6)',ecolor='black',capsize=2, capthick=2,color='y')
    tscat6 = ax.errorbar(c.applyCuts(t,cutx6_30),c.applyCuts(fpi,cutx6_30),yerr=np.sqrt(c.applyCuts(lumi,cutx6_30))/c.applyCuts(lumi,cutx6_30),fmt='.',label='$x_{\pi}$=(0.6,0.7)',ecolor='black',capsize=2, capthick=2,color='purple')
    tscat7 = ax.errorbar(c.applyCuts(t,cutx7_30),c.applyCuts(fpi,cutx7_30),yerr=np.sqrt(c.applyCuts(lumi,cutx7_30))/c.applyCuts(lumi,cutx7_30),fmt='.',label='$x_{\pi}$=(0.7,0.8)',ecolor='black',capsize=2, capthick=2,color='aqua')
    tscat8 = ax.errorbar(c.applyCuts(t,cutx8_30),c.applyCuts(fpi,cutx8_30),yerr=np.sqrt(c.applyCuts(lumi,cutx8_30))/c.applyCuts(lumi,cutx8_30),fmt='.',label='$x_{\pi}$=(0.8,0.9)',ecolor='black',capsize=2, capthick=2,color='gray')
    tscat9 = ax.errorbar(c.applyCuts(t,cutx9_30),c.applyCuts(fpi,cutx9_30),yerr=np.sqrt(c.applyCuts(lumi,cutx9_30))/c.applyCuts(lumi,cutx9_30),fmt='.',label='$x_{\pi}$=(0.9,1.0)',ecolor='black',capsize=2, capthick=2,color='pink')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(0.,0.3)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.text(0.65, 0.85, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    plt.title('$F^{\pi}_{2}$ vs -t', fontsize =20)
    
    ax = f.add_subplot(223)
    tscat0 = ax.errorbar(c.applyCuts(t,cutx0_60),c.applyCuts(fpi,cutx0_60),yerr=np.sqrt(c.applyCuts(lumi,cutx0_60))/c.applyCuts(lumi,cutx0_60),fmt='.',label='$x_{\pi}$=(0.0,0.1)',ecolor='black',capsize=2, capthick=2,color='b')
    tscat1 = ax.errorbar(c.applyCuts(t,cutx1_60),c.applyCuts(fpi,cutx1_60),yerr=np.sqrt(c.applyCuts(lumi,cutx1_60))/c.applyCuts(lumi,cutx1_60),fmt='.',label='$x_{\pi}$=(0.1,0.2)',ecolor='black',capsize=2, capthick=2,color='g')
    tscat2 = ax.errorbar(c.applyCuts(t,cutx2_60),c.applyCuts(fpi,cutx2_60),yerr=np.sqrt(c.applyCuts(lumi,cutx2_60))/c.applyCuts(lumi,cutx2_60),fmt='.',label='$x_{\pi}$=(0.2,0.3)',ecolor='black',capsize=2, capthick=2,color='r')
    tscat3 = ax.errorbar(c.applyCuts(t,cutx3_60),c.applyCuts(fpi,cutx3_60),yerr=np.sqrt(c.applyCuts(lumi,cutx3_60))/c.applyCuts(lumi,cutx3_60),fmt='.',label='$x_{\pi}$=(0.3,0.4)',ecolor='black',capsize=2, capthick=2,color='c')
    tscat4 = ax.errorbar(c.applyCuts(t,cutx4_60),c.applyCuts(fpi,cutx4_60),yerr=np.sqrt(c.applyCuts(lumi,cutx4_60))/c.applyCuts(lumi,cutx4_60),fmt='.',label='$x_{\pi}$=(0.4,0.5)',ecolor='black',capsize=2, capthick=2,color='m')
    tscat5 = ax.errorbar(c.applyCuts(t,cutx5_60),c.applyCuts(fpi,cutx5_60),yerr=np.sqrt(c.applyCuts(lumi,cutx5_60))/c.applyCuts(lumi,cutx5_60),fmt='.',label='$x_{\pi}$=(0.5,0.6)',ecolor='black',capsize=2, capthick=2,color='y')
    tscat6 = ax.errorbar(c.applyCuts(t,cutx6_60),c.applyCuts(fpi,cutx6_60),yerr=np.sqrt(c.applyCuts(lumi,cutx6_60))/c.applyCuts(lumi,cutx6_60),fmt='.',label='$x_{\pi}$=(0.6,0.7)',ecolor='black',capsize=2, capthick=2,color='purple')
    tscat7 = ax.errorbar(c.applyCuts(t,cutx7_60),c.applyCuts(fpi,cutx7_60),yerr=np.sqrt(c.applyCuts(lumi,cutx7_60))/c.applyCuts(lumi,cutx7_60),fmt='.',label='$x_{\pi}$=(0.7,0.8)',ecolor='black',capsize=2, capthick=2,color='aqua')
    tscat8 = ax.errorbar(c.applyCuts(t,cutx8_60),c.applyCuts(fpi,cutx8_60),yerr=np.sqrt(c.applyCuts(lumi,cutx8_60))/c.applyCuts(lumi,cutx8_60),fmt='.',label='$x_{\pi}$=(0.8,0.9)',ecolor='black',capsize=2, capthick=2,color='gray')
    tscat9 = ax.errorbar(c.applyCuts(t,cutx9_60),c.applyCuts(fpi,cutx9_60),yerr=np.sqrt(c.applyCuts(lumi,cutx9_60))/c.applyCuts(lumi,cutx9_60),fmt='.',label='$x_{\pi}$=(0.9,1.0)',ecolor='black',capsize=2, capthick=2,color='pink')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(0.,0.3)
    ax.text(0.65, 0.85, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    
    ax = f.add_subplot(224)
    tscat0 = ax.errorbar(c.applyCuts(t,cutx0_120),c.applyCuts(fpi,cutx0_120),yerr=np.sqrt(c.applyCuts(lumi,cutx0_120))/c.applyCuts(lumi,cutx0_120),fmt='.',label='$x_{\pi}$=(0.0,0.1)',ecolor='black',capsize=2, capthick=2,color='b')
    tscat1 = ax.errorbar(c.applyCuts(t,cutx1_120),c.applyCuts(fpi,cutx1_120),yerr=np.sqrt(c.applyCuts(lumi,cutx1_120))/c.applyCuts(lumi,cutx1_120),fmt='.',label='$x_{\pi}$=(0.1,0.2)',ecolor='black',capsize=2, capthick=2,color='g')
    tscat2 = ax.errorbar(c.applyCuts(t,cutx2_120),c.applyCuts(fpi,cutx2_120),yerr=np.sqrt(c.applyCuts(lumi,cutx2_120))/c.applyCuts(lumi,cutx2_120),fmt='.',label='$x_{\pi}$=(0.2,0.3)',ecolor='black',capsize=2, capthick=2,color='r')
    tscat3 = ax.errorbar(c.applyCuts(t,cutx3_120),c.applyCuts(fpi,cutx3_120),yerr=np.sqrt(c.applyCuts(lumi,cutx3_120))/c.applyCuts(lumi,cutx3_120),fmt='.',label='$x_{\pi}$=(0.3,0.4)',ecolor='black',capsize=2, capthick=2,color='c')
    tscat4 = ax.errorbar(c.applyCuts(t,cutx4_120),c.applyCuts(fpi,cutx4_120),yerr=np.sqrt(c.applyCuts(lumi,cutx4_120))/c.applyCuts(lumi,cutx4_120),fmt='.',label='$x_{\pi}$=(0.4,0.5)',ecolor='black',capsize=2, capthick=2,color='m')
    tscat5 = ax.errorbar(c.applyCuts(t,cutx5_120),c.applyCuts(fpi,cutx5_120),yerr=np.sqrt(c.applyCuts(lumi,cutx5_120))/c.applyCuts(lumi,cutx5_120),fmt='.',label='$x_{\pi}$=(0.5,0.6)',ecolor='black',capsize=2, capthick=2,color='y')
    tscat6 = ax.errorbar(c.applyCuts(t,cutx6_120),c.applyCuts(fpi,cutx6_120),yerr=np.sqrt(c.applyCuts(lumi,cutx6_120))/c.applyCuts(lumi,cutx6_120),fmt='.',label='$x_{\pi}$=(0.6,0.7)',ecolor='black',capsize=2, capthick=2,color='purple')
    tscat7 = ax.errorbar(c.applyCuts(t,cutx7_120),c.applyCuts(fpi,cutx7_120),yerr=np.sqrt(c.applyCuts(lumi,cutx7_120))/c.applyCuts(lumi,cutx7_120),fmt='.',label='$x_{\pi}$=(0.7,0.8)',ecolor='black',capsize=2, capthick=2,color='aqua')
    tscat8 = ax.errorbar(c.applyCuts(t,cutx8_120),c.applyCuts(fpi,cutx8_120),yerr=np.sqrt(c.applyCuts(lumi,cutx8_120))/c.applyCuts(lumi,cutx8_120),fmt='.',label='$x_{\pi}$=(0.8,0.9)',ecolor='black',capsize=2, capthick=2,color='gray')
    tscat9 = ax.errorbar(c.applyCuts(t,cutx9_120),c.applyCuts(fpi,cutx9_120),yerr=np.sqrt(c.applyCuts(lumi,cutx9_120))/c.applyCuts(lumi,cutx9_120),fmt='.',label='$x_{\pi}$=(0.9,1.0)',ecolor='black',capsize=2, capthick=2,color='pink')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(0.,0.3)
    ax.text(0.65, 0.85, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    plt.xlabel('-t')
    plt.tight_layout()

    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    
    f.savefig('OUTPUTS/fpivt_xbin.png')

def fpivtPrime_Plot():

    f = plt.figure(figsize=(11.69,8.27))
    
    ax = f.add_subplot(421)
    tscat1 = ax.errorbar(c.applyCuts(tPrime,cut7),c.applyCuts(fpi,cut7),yerr=np.sqrt(c.applyCuts(lumi,cut7))/c.applyCuts(lumi,cut7),fmt='.',label='$Q^2$=7 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=7 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    plt.ylabel('log($F^{\pi}_{2}$)')

    ax = f.add_subplot(422)
    tscat2 = ax.errorbar(c.applyCuts(tPrime,cut15),c.applyCuts(fpi,cut15),yerr=np.sqrt(c.applyCuts(lumi,cut15))/c.applyCuts(lumi,cut15),fmt='.',label='$Q^2$=15 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=15 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_formatter(plt.NullFormatter())

    plt.title('log($F^{\pi}_{2}$) vs -t\'', fontsize =20)
    
    ax = f.add_subplot(423)
    tscat3 = ax.errorbar(c.applyCuts(tPrime,cut30),c.applyCuts(fpi,cut30),yerr=np.sqrt(c.applyCuts(lumi,cut30))/c.applyCuts(lumi,cut30),fmt='.',label='$Q^2$=30 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    
    ax = f.add_subplot(424)
    tscat4 = ax.errorbar(c.applyCuts(tPrime,cut60),c.applyCuts(fpi,cut60),yerr=np.sqrt(c.applyCuts(lumi,cut60))/c.applyCuts(lumi,cut60),fmt='.',label='$Q^2$=60 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    
    ax = f.add_subplot(425)
    tscat5 = ax.errorbar(c.applyCuts(tPrime,cut120),c.applyCuts(fpi,cut120),yerr=np.sqrt(c.applyCuts(lumi,cut120))/c.applyCuts(lumi,cut120),fmt='.',label='$Q^2$=120 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    ax = f.add_subplot(426)
    tscat6 = ax.errorbar(c.applyCuts(tPrime,cut240),c.applyCuts(fpi,cut240),yerr=np.sqrt(c.applyCuts(lumi,cut240))/c.applyCuts(lumi,cut240),fmt='.',label='$Q^2$=240 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=240 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(427)
    tscat7 = ax.errorbar(c.applyCuts(tPrime,cut480),c.applyCuts(fpi,cut480),yerr=np.sqrt(c.applyCuts(lumi,cut480))/c.applyCuts(lumi,cut480),fmt='.',label='$Q^2$=480 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=480 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))

    ax = f.add_subplot(428)
    tscat8 = ax.errorbar(c.applyCuts(tPrime,cut1000),c.applyCuts(fpi,cut1000),yerr=np.sqrt(c.applyCuts(lumi,cut1000))/c.applyCuts(lumi,cut1000),fmt='.',label='$Q^2$=1000 $GeV^2$',ecolor='black',capsize=2, capthick=2)
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=1000 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    plt.xlabel('-t\'')
    plt.tight_layout()

    plt.subplots_adjust(hspace=0.0,wspace=0.0)

    f.savefig('OUTPUTS/fpivtPrime.png')

def fpivtPrime_xbin_Plot():

    f = plt.figure(figsize=(11.69,8.27))
    
    ax = f.add_subplot(221)
    tscat0 = ax.errorbar(c.applyCuts(tPrime,cutx0_15),c.applyCuts(fpi,cutx0_15),yerr=np.sqrt(c.applyCuts(lumi,cutx0_15))/c.applyCuts(lumi,cutx0_15),fmt='.',label='$x_{\pi}$=(0.0,0.1)',ecolor='black',capsize=2, capthick=2,color='b')
    tscat1 = ax.errorbar(c.applyCuts(tPrime,cutx1_15),c.applyCuts(fpi,cutx1_15),yerr=np.sqrt(c.applyCuts(lumi,cutx1_15))/c.applyCuts(lumi,cutx1_15),fmt='.',label='$x_{\pi}$=(0.1,0.2)',ecolor='black',capsize=2, capthick=2,color='g')
    tscat2 = ax.errorbar(c.applyCuts(tPrime,cutx2_15),c.applyCuts(fpi,cutx2_15),yerr=np.sqrt(c.applyCuts(lumi,cutx2_15))/c.applyCuts(lumi,cutx2_15),fmt='.',label='$x_{\pi}$=(0.2,0.3)',ecolor='black',capsize=2, capthick=2,color='r')
    tscat3 = ax.errorbar(c.applyCuts(tPrime,cutx3_15),c.applyCuts(fpi,cutx3_15),yerr=np.sqrt(c.applyCuts(lumi,cutx3_15))/c.applyCuts(lumi,cutx3_15),fmt='.',label='$x_{\pi}$=(0.3,0.4)',ecolor='black',capsize=2, capthick=2,color='c')
    tscat4 = ax.errorbar(c.applyCuts(tPrime,cutx4_15),c.applyCuts(fpi,cutx4_15),yerr=np.sqrt(c.applyCuts(lumi,cutx4_15))/c.applyCuts(lumi,cutx4_15),fmt='.',label='$x_{\pi}$=(0.4,0.5)',ecolor='black',capsize=2, capthick=2,color='m')
    tscat5 = ax.errorbar(c.applyCuts(tPrime,cutx5_15),c.applyCuts(fpi,cutx5_15),yerr=np.sqrt(c.applyCuts(lumi,cutx5_15))/c.applyCuts(lumi,cutx5_15),fmt='.',label='$x_{\pi}$=(0.5,0.6)',ecolor='black',capsize=2, capthick=2,color='y')
    tscat6 = ax.errorbar(c.applyCuts(tPrime,cutx6_15),c.applyCuts(fpi,cutx6_15),yerr=np.sqrt(c.applyCuts(lumi,cutx6_15))/c.applyCuts(lumi,cutx6_15),fmt='.',label='$x_{\pi}$=(0.6,0.7)',ecolor='black',capsize=2, capthick=2,color='purple')
    tscat7 = ax.errorbar(c.applyCuts(tPrime,cutx7_15),c.applyCuts(fpi,cutx7_15),yerr=np.sqrt(c.applyCuts(lumi,cutx7_15))/c.applyCuts(lumi,cutx7_15),fmt='.',label='$x_{\pi}$=(0.7,0.8)',ecolor='black',capsize=2, capthick=2,color='aqua')
    tscat8 = ax.errorbar(c.applyCuts(tPrime,cutx8_15),c.applyCuts(fpi,cutx8_15),yerr=np.sqrt(c.applyCuts(lumi,cutx8_15))/c.applyCuts(lumi,cutx8_15),fmt='.',label='$x_{\pi}$=(0.8,0.9)',ecolor='black',capsize=2, capthick=2,color='gray')
    tscat9 = ax.errorbar(c.applyCuts(tPrime,cutx9_15),c.applyCuts(fpi,cutx9_15),yerr=np.sqrt(c.applyCuts(lumi,cutx9_15))/c.applyCuts(lumi,cutx9_15),fmt='.',label='$x_{\pi}$=(0.9,1.0)',ecolor='black',capsize=2, capthick=2,color='pink')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(0.,0.25)
    ax.text(0.65, 0.85, '$Q^2$=15 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='lower'))

    plt.ylabel('$F^{\pi}_{2}$')
    
    ax = f.add_subplot(222)
    tscat0 = ax.errorbar(c.applyCuts(tPrime,cutx0_30),c.applyCuts(fpi,cutx0_30),yerr=np.sqrt(c.applyCuts(lumi,cutx0_30))/c.applyCuts(lumi,cutx0_30),fmt='.',label='$x_{\pi}$=(0.0,0.1)',ecolor='black',capsize=2, capthick=2,color='b')
    tscat1 = ax.errorbar(c.applyCuts(tPrime,cutx1_30),c.applyCuts(fpi,cutx1_30),yerr=np.sqrt(c.applyCuts(lumi,cutx1_30))/c.applyCuts(lumi,cutx1_30),fmt='.',label='$x_{\pi}$=(0.1,0.2)',ecolor='black',capsize=2, capthick=2,color='g')
    tscat2 = ax.errorbar(c.applyCuts(tPrime,cutx2_30),c.applyCuts(fpi,cutx2_30),yerr=np.sqrt(c.applyCuts(lumi,cutx2_30))/c.applyCuts(lumi,cutx2_30),fmt='.',label='$x_{\pi}$=(0.2,0.3)',ecolor='black',capsize=2, capthick=2,color='r')
    tscat3 = ax.errorbar(c.applyCuts(tPrime,cutx3_30),c.applyCuts(fpi,cutx3_30),yerr=np.sqrt(c.applyCuts(lumi,cutx3_30))/c.applyCuts(lumi,cutx3_30),fmt='.',label='$x_{\pi}$=(0.3,0.4)',ecolor='black',capsize=2, capthick=2,color='c')
    tscat4 = ax.errorbar(c.applyCuts(tPrime,cutx4_30),c.applyCuts(fpi,cutx4_30),yerr=np.sqrt(c.applyCuts(lumi,cutx4_30))/c.applyCuts(lumi,cutx4_30),fmt='.',label='$x_{\pi}$=(0.4,0.5)',ecolor='black',capsize=2, capthick=2,color='m')
    tscat5 = ax.errorbar(c.applyCuts(tPrime,cutx5_30),c.applyCuts(fpi,cutx5_30),yerr=np.sqrt(c.applyCuts(lumi,cutx5_30))/c.applyCuts(lumi,cutx5_30),fmt='.',label='$x_{\pi}$=(0.5,0.6)',ecolor='black',capsize=2, capthick=2,color='y')
    tscat6 = ax.errorbar(c.applyCuts(tPrime,cutx6_30),c.applyCuts(fpi,cutx6_30),yerr=np.sqrt(c.applyCuts(lumi,cutx6_30))/c.applyCuts(lumi,cutx6_30),fmt='.',label='$x_{\pi}$=(0.6,0.7)',ecolor='black',capsize=2, capthick=2,color='purple')
    tscat7 = ax.errorbar(c.applyCuts(tPrime,cutx7_30),c.applyCuts(fpi,cutx7_30),yerr=np.sqrt(c.applyCuts(lumi,cutx7_30))/c.applyCuts(lumi,cutx7_30),fmt='.',label='$x_{\pi}$=(0.7,0.8)',ecolor='black',capsize=2, capthick=2,color='aqua')
    tscat8 = ax.errorbar(c.applyCuts(tPrime,cutx8_30),c.applyCuts(fpi,cutx8_30),yerr=np.sqrt(c.applyCuts(lumi,cutx8_30))/c.applyCuts(lumi,cutx8_30),fmt='.',label='$x_{\pi}$=(0.8,0.9)',ecolor='black',capsize=2, capthick=2,color='gray')
    tscat9 = ax.errorbar(c.applyCuts(tPrime,cutx9_30),c.applyCuts(fpi,cutx9_30),yerr=np.sqrt(c.applyCuts(lumi,cutx9_30))/c.applyCuts(lumi,cutx9_30),fmt='.',label='$x_{\pi}$=(0.9,1.0)',ecolor='black',capsize=2, capthick=2,color='pink')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(0.,0.25)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    ax.text(0.65, 0.85, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    plt.title('$F^{\pi}_{2}$ vs -t\'', fontsize =20)
    
    ax = f.add_subplot(223)
    tscat0 = ax.errorbar(c.applyCuts(tPrime,cutx0_60),c.applyCuts(fpi,cutx0_60),yerr=np.sqrt(c.applyCuts(lumi,cutx0_60))/c.applyCuts(lumi,cutx0_60),fmt='.',label='$x_{\pi}$=(0.0,0.1)',ecolor='black',capsize=2, capthick=2,color='b')
    tscat1 = ax.errorbar(c.applyCuts(tPrime,cutx1_60),c.applyCuts(fpi,cutx1_60),yerr=np.sqrt(c.applyCuts(lumi,cutx1_60))/c.applyCuts(lumi,cutx1_60),fmt='.',label='$x_{\pi}$=(0.1,0.2)',ecolor='black',capsize=2, capthick=2,color='g')
    tscat2 = ax.errorbar(c.applyCuts(tPrime,cutx2_60),c.applyCuts(fpi,cutx2_60),yerr=np.sqrt(c.applyCuts(lumi,cutx2_60))/c.applyCuts(lumi,cutx2_60),fmt='.',label='$x_{\pi}$=(0.2,0.3)',ecolor='black',capsize=2, capthick=2,color='r')
    tscat3 = ax.errorbar(c.applyCuts(tPrime,cutx3_60),c.applyCuts(fpi,cutx3_60),yerr=np.sqrt(c.applyCuts(lumi,cutx3_60))/c.applyCuts(lumi,cutx3_60),fmt='.',label='$x_{\pi}$=(0.3,0.4)',ecolor='black',capsize=2, capthick=2,color='c')
    tscat4 = ax.errorbar(c.applyCuts(tPrime,cutx4_60),c.applyCuts(fpi,cutx4_60),yerr=np.sqrt(c.applyCuts(lumi,cutx4_60))/c.applyCuts(lumi,cutx4_60),fmt='.',label='$x_{\pi}$=(0.4,0.5)',ecolor='black',capsize=2, capthick=2,color='m')
    tscat5 = ax.errorbar(c.applyCuts(tPrime,cutx5_60),c.applyCuts(fpi,cutx5_60),yerr=np.sqrt(c.applyCuts(lumi,cutx5_60))/c.applyCuts(lumi,cutx5_60),fmt='.',label='$x_{\pi}$=(0.5,0.6)',ecolor='black',capsize=2, capthick=2,color='y')
    tscat6 = ax.errorbar(c.applyCuts(tPrime,cutx6_60),c.applyCuts(fpi,cutx6_60),yerr=np.sqrt(c.applyCuts(lumi,cutx6_60))/c.applyCuts(lumi,cutx6_60),fmt='.',label='$x_{\pi}$=(0.6,0.7)',ecolor='black',capsize=2, capthick=2,color='purple')
    tscat7 = ax.errorbar(c.applyCuts(tPrime,cutx7_60),c.applyCuts(fpi,cutx7_60),yerr=np.sqrt(c.applyCuts(lumi,cutx7_60))/c.applyCuts(lumi,cutx7_60),fmt='.',label='$x_{\pi}$=(0.7,0.8)',ecolor='black',capsize=2, capthick=2,color='aqua')
    tscat8 = ax.errorbar(c.applyCuts(tPrime,cutx8_60),c.applyCuts(fpi,cutx8_60),yerr=np.sqrt(c.applyCuts(lumi,cutx8_60))/c.applyCuts(lumi,cutx8_60),fmt='.',label='$x_{\pi}$=(0.8,0.9)',ecolor='black',capsize=2, capthick=2,color='gray')
    tscat9 = ax.errorbar(c.applyCuts(tPrime,cutx9_60),c.applyCuts(fpi,cutx9_60),yerr=np.sqrt(c.applyCuts(lumi,cutx9_60))/c.applyCuts(lumi,cutx9_60),fmt='.',label='$x_{\pi}$=(0.9,1.0)',ecolor='black',capsize=2, capthick=2,color='pink')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(0.,0.25)
    ax.text(0.65, 0.85, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    
    ax = f.add_subplot(224)
    tscat0 = ax.errorbar(c.applyCuts(tPrime,cutx0_120),c.applyCuts(fpi,cutx0_120),yerr=np.sqrt(c.applyCuts(lumi,cutx0_120))/c.applyCuts(lumi,cutx0_120),fmt='.',label='$x_{\pi}$=(0.0,0.1)',ecolor='black',capsize=2, capthick=2,color='b')
    tscat1 = ax.errorbar(c.applyCuts(tPrime,cutx1_120),c.applyCuts(fpi,cutx1_120),yerr=np.sqrt(c.applyCuts(lumi,cutx1_120))/c.applyCuts(lumi,cutx1_120),fmt='.',label='$x_{\pi}$=(0.1,0.2)',ecolor='black',capsize=2, capthick=2,color='g')
    tscat2 = ax.errorbar(c.applyCuts(tPrime,cutx2_120),c.applyCuts(fpi,cutx2_120),yerr=np.sqrt(c.applyCuts(lumi,cutx2_120))/c.applyCuts(lumi,cutx2_120),fmt='.',label='$x_{\pi}$=(0.2,0.3)',ecolor='black',capsize=2, capthick=2,color='r')
    tscat3 = ax.errorbar(c.applyCuts(tPrime,cutx3_120),c.applyCuts(fpi,cutx3_120),yerr=np.sqrt(c.applyCuts(lumi,cutx3_120))/c.applyCuts(lumi,cutx3_120),fmt='.',label='$x_{\pi}$=(0.3,0.4)',ecolor='black',capsize=2, capthick=2,color='c')
    tscat4 = ax.errorbar(c.applyCuts(tPrime,cutx4_120),c.applyCuts(fpi,cutx4_120),yerr=np.sqrt(c.applyCuts(lumi,cutx4_120))/c.applyCuts(lumi,cutx4_120),fmt='.',label='$x_{\pi}$=(0.4,0.5)',ecolor='black',capsize=2, capthick=2,color='m')
    tscat5 = ax.errorbar(c.applyCuts(tPrime,cutx5_120),c.applyCuts(fpi,cutx5_120),yerr=np.sqrt(c.applyCuts(lumi,cutx5_120))/c.applyCuts(lumi,cutx5_120),fmt='.',label='$x_{\pi}$=(0.5,0.6)',ecolor='black',capsize=2, capthick=2,color='y')
    tscat6 = ax.errorbar(c.applyCuts(tPrime,cutx6_120),c.applyCuts(fpi,cutx6_120),yerr=np.sqrt(c.applyCuts(lumi,cutx6_120))/c.applyCuts(lumi,cutx6_120),fmt='.',label='$x_{\pi}$=(0.6,0.7)',ecolor='black',capsize=2, capthick=2,color='purple')
    tscat7 = ax.errorbar(c.applyCuts(tPrime,cutx7_120),c.applyCuts(fpi,cutx7_120),yerr=np.sqrt(c.applyCuts(lumi,cutx7_120))/c.applyCuts(lumi,cutx7_120),fmt='.',label='$x_{\pi}$=(0.7,0.8)',ecolor='black',capsize=2, capthick=2,color='aqua')
    tscat8 = ax.errorbar(c.applyCuts(tPrime,cutx8_120),c.applyCuts(fpi,cutx8_120),yerr=np.sqrt(c.applyCuts(lumi,cutx8_120))/c.applyCuts(lumi,cutx8_120),fmt='.',label='$x_{\pi}$=(0.8,0.9)',ecolor='black',capsize=2, capthick=2,color='gray')
    tscat9 = ax.errorbar(c.applyCuts(tPrime,cutx9_120),c.applyCuts(fpi,cutx9_120),yerr=np.sqrt(c.applyCuts(lumi,cutx9_120))/c.applyCuts(lumi,cutx9_120),fmt='.',label='$x_{\pi}$=(0.9,1.0)',ecolor='black',capsize=2, capthick=2,color='pink')
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(0.,0.25)
    ax.text(0.65, 0.85, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    plt.xlabel('-t\'')
    plt.tight_layout()

    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    
    f.savefig('OUTPUTS/fpivtPrime_xbin.png')    

    
def Lumi():

    # Luminosity
    
    # all binned events
    sig_all = tot_sigma*1e5 # the 1e5 corrects the adjustment up top
    evts_all = len(tot_sigma)*[len(tot_sigma)]

    lumi = [(e)/((s)*(binQ2)*(binx)) for e,s in zip(evts_all,sig_all)]
    print("---------------------------------\n")
    print("\nLuminosity: ", lumi)

    nevt = [10/(l*1e-6) for l in lumi] # The 1e-6 converts properly and the 1e-5 corrects for the adjustment to the total cross-section
    nevt  = np.asarray(nevt)
    print("\nEvents expected running at 10 $fb^{-1}$: ", nevt)

    # f = plt.figure(figsize=(11.69,8.27))
    
    # ax = f.add_subplot(331)
    # ax.hist(nevt[0],bins=c.setbin(nevt[0],200),label='$Q^2$=7 $GeV^2$',histtype='step', alpha=0.5, stacked=True, fill=True)
    # plt.xscale('log')
    # ax.text(0.60, 0.85, '$Q^2$=7 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')

    # ax = f.add_subplot(332)
    # ax.hist(nevt[1],bins=c.setbin(nevt[1],200),label='$Q^2$=15 $GeV^2$',histtype='step', alpha=0.5, stacked=True, fill=True)
    # plt.xscale('log')
    # ax.text(0.60, 0.85, '$Q^2$=15 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    
    # ax = f.add_subplot(333)
    # ax.hist(nevt[2],bins=c.setbin(nevt[2],200),label='$Q^2$=30 $GeV^2$',histtype='step', alpha=0.5, stacked=True, fill=True)
    # plt.xscale('log')
    # ax.text(0.60, 0.85, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')

    # plt.title('Events expected running at 10 $fb^{-1}$', fontsize =20)
    
    # ax = f.add_subplot(334)
    # ax.hist(nevt[3],bins=c.setbin(nevt[3],200),label='$Q^2$=60 $GeV^2$',histtype='step', alpha=0.5, stacked=True, fill=True)
    # plt.xscale('log')
    # ax.text(0.60, 0.85, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    
    # ax = f.add_subplot(335)
    # ax.hist(nevt[4],bins=c.setbin(nevt[4],200),label='$Q^2$=120 $GeV^2$',histtype='step', alpha=0.5, stacked=True, fill=True)
    # plt.xscale('log')
    # ax.text(0.60, 0.85, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    
    # ax = f.add_subplot(336)
    # ax.hist(nevt[5],bins=c.setbin(nevt[5],200),label='$Q^2$=240 $GeV^2$',histtype='step', alpha=0.5, stacked=True, fill=True)
    # plt.xscale('log')
    # ax.text(0.60, 0.85, '$Q^2$=240 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')

    # ax = f.add_subplot(337)
    # ax.hist(nevt[6],bins=c.setbin(nevt[6],200),label='$Q^2$=480 $GeV^2$',histtype='step', alpha=0.5, stacked=True, fill=True)
    # plt.xscale('log')
    # ax.text(0.60, 0.85, '$Q^2$=480 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')

    # ax = f.add_subplot(338)
    # ax.hist(nevt[7],bins=c.setbin(nevt[7],200),label='$Q^2$=1000 $GeV^2$',histtype='step', alpha=0.5, stacked=True, fill=True)
    # plt.xscale('log')
    # ax.text(0.60, 0.85, '$Q^2$=1000 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')

    # plt.tight_layout()

    # return [nevt_7,nevt_15,nevt_30,nevt_60,nevt_120,nevt_240,nevt_480,nevt_1000,nevt]
    return nevt

lumi = Lumi()    
    
def main() :

    phase_space()
    fpivxpi_Plot()
    fpivxpi_Plot_nolog()
    fpivt_Plot()
    fpivt_xbin_Plot()
    fpivtPrime_Plot()
    fpivtPrime_xbin_Plot()
    sigmavxpi_Plot()
    plt.show()
    
if __name__=='__main__': main()
