#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2020-07-15 15:39:54 trottar"
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
import time,math,sys

sys.path.insert(0,'/home/trottar/bin/python/')
import root2py as r2p

kinematics = sys.argv[1]

# kinematics="k_lambda_18on275"
# kinematics="k_lambda_10on100"
# kinematics="k_lambda_5on100"
# kinematics="k_lambda_5on41"

rootName="/home/trottar/ResearchNP/JLEIC/USERS/trottar/OUTPUTS/%s.root" % kinematics

tree = up.open(rootName)["Evnts"]
branch = r2p.pyBranch(tree)

# Define phyisics data
s_e = branch.findBranch("invts","s_e")
s_q = branch.findBranch("invts","s_q")
Q2 = branch.findBranch("invts","Q2")
xBj = branch.findBranch("invts","xBj")
t = -branch.findBranch("invts","tSpectator")
y_D = branch.findBranch("invts","y_D")
nu = branch.findBranch("invts","nu")
TwoPdotk = branch.findBranch("invts","TwoPdotk")
TDIS_xbj = tree.array("TDIS_xbj")
sigma_dis = tree.array("sigma_dis")*(1e-5)
TDIS_y = tree.array("TDIS_y")
pkx_Lab = tree.array("pkx_Lab") # Long kaon momentum
pky_Lab = tree.array("pky_Lab")
pkz_Lab = tree.array("pkz_Lab")
EkE_Lab = tree.array("EkE_Lab")
pprx_inc = tree.array("pprx_inc") # Long proton momentum
ppry_inc = tree.array("ppry_inc")
pprz_inc = tree.array("pprz_inc")
EprE_inc = tree.array("EprE_inc")
y = tree.array("y")
fk = tree.array("fk")
f2N = tree.array("f2N")
xk = tree.array("xk")
yk = tree.array("yk")
tk = tree.array("tk")
escat = tree.array("EScatRest")
pprz_inc = tree.array("pprz_inc")
pkz_Lab = tree.array("pkz_Lab")

# More coarse, set bin size
xkarray = np.arange(0.005,1.0,0.01).tolist()

cutDict = {}

# Less coarse
for i,x in enumerate(xkarray):
    xktmp = '{"xkcut%i" : ((%0.5f <= xk) & (xk <= %0.5f))}' % (i,xkarray[i]-0.0005,xkarray[i]+0.0005)
    print('{"xkcut%i" : ((%0.5f <= xk) & (xk <= %0.5f))}' % (i,xkarray[i]-0.0005,xkarray[i]+0.0005))
    cutDict.update(eval(xktmp))
c = r2p.pyPlot(cutDict)

xkcut = []
for i,evt in enumerate(xkarray):
    xkcut.append("xkcut%s" % i)
    tmp  = "xkcut_%s = [\"xkcut%s\"]" % (i,i)
    exec(tmp)
    
s_e = c.applyCuts(s_e,xkcut,either=True)
s_q = c.applyCuts(s_q,xkcut,either=True)
Q2 = c.applyCuts(Q2,xkcut,either=True)
xBj = c.applyCuts(xBj,xkcut,either=True)
t = c.applyCuts(t,xkcut,either=True)
y_D = c.applyCuts(y_D,xkcut,either=True)
nu = c.applyCuts(nu,xkcut,either=True)
TwoPdotk = c.applyCuts(TwoPdotk,xkcut,either=True)
TDIS_xbj = c.applyCuts(TDIS_xbj,xkcut,either=True)
sigma_dis = c.applyCuts(sigma_dis,xkcut,either=True)
TDIS_y = c.applyCuts(TDIS_y,xkcut,either=True)
pkx_Lab = c.applyCuts(pkx_Lab,xkcut,either=True)
pky_Lab = c.applyCuts(pky_Lab,xkcut,either=True)
pkz_Lab = c.applyCuts(pkz_Lab,xkcut,either=True)
EkE_Lab = c.applyCuts(EkE_Lab,xkcut,either=True)
pprx_inc = c.applyCuts(pprx_inc,xkcut,either=True)
ppry_inc = c.applyCuts(ppry_inc,xkcut,either=True)
pprz_inc = c.applyCuts(pprz_inc,xkcut,either=True)
EprE_inc = c.applyCuts(EprE_inc,xkcut,either=True)
y = c.applyCuts(y,xkcut,either=True)
fk = c.applyCuts(fk,xkcut,either=True)
f2N = c.applyCuts(f2N,xkcut,either=True)
yk = c.applyCuts(yk,xkcut,either=True)
tk = c.applyCuts(tk,xkcut,either=True)
escat = c.applyCuts(escat,xkcut,either=True)
# pprz_inc = c.applyCuts(pprz_inc,xkcut,either=True)
# pkz_Lab = c.applyCuts(pkz_Lab,xkcut,either=True)

xk = c.applyCuts(xk,xkcut,either=True)

# More coarse, set bin size
Q2array = np.arange(5.0,1000.0,10.0).tolist()

# Less coarse
for i,x in enumerate(Q2array) :
    Q2tmp = '{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,Q2array[i]-0.5,Q2array[i]+0.5)
    print('{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,Q2array[i]-0.5,Q2array[i]+0.5))
    cutDict.update(eval(Q2tmp))
c = r2p.pyPlot(cutDict)

Q2cut = []
for i,evt in enumerate(Q2array):
    Q2cut.append("Q2cut%s" % i)

s_e = c.applyCuts(s_e,Q2cut,either=True)
s_q = c.applyCuts(s_q,Q2cut,either=True)
xBj = c.applyCuts(xBj,Q2cut,either=True)
t = c.applyCuts(t,Q2cut,either=True)
y_D = c.applyCuts(y_D,Q2cut,either=True)
nu = c.applyCuts(nu,Q2cut,either=True)
TwoPdotk = c.applyCuts(TwoPdotk,Q2cut,either=True)
TDIS_xbj = c.applyCuts(TDIS_xbj,Q2cut,either=True)
sigma_dis = c.applyCuts(sigma_dis,Q2cut,either=True)
TDIS_y = c.applyCuts(TDIS_y,Q2cut,either=True)
pkx_Lab = c.applyCuts(pkx_Lab,Q2cut,either=True)
pky_Lab = c.applyCuts(pky_Lab,Q2cut,either=True)
pkz_Lab = c.applyCuts(pkz_Lab,Q2cut,either=True)
EkE_Lab = c.applyCuts(EkE_Lab,Q2cut,either=True)
pprx_inc = c.applyCuts(pprx_inc,Q2cut,either=True)
ppry_inc = c.applyCuts(ppry_inc,Q2cut,either=True)
pprz_inc = c.applyCuts(pprz_inc,Q2cut,either=True)
EprE_inc = c.applyCuts(EprE_inc,Q2cut,either=True)
y = c.applyCuts(y,Q2cut,either=True)
fk = c.applyCuts(fk,Q2cut,either=True)
f2N = c.applyCuts(f2N,Q2cut,either=True)
xk = c.applyCuts(xk,Q2cut,either=True)
yk = c.applyCuts(yk,Q2cut,either=True)
tk = c.applyCuts(tk,Q2cut,either=True)
escat = c.applyCuts(escat,Q2cut,either=True)
# pprz_inc = c.applyCuts(pprz_inc,Q2cut,either=True)
# pkz_Lab = c.applyCuts(pkz_Lab,Q2cut,either=True)

Q2 = c.applyCuts(Q2,Q2cut,either=True)

sigma_tdis = sigma_dis*(fk/f2N)
pNz_Lab = pprz_inc-pkz_Lab # Long neutron momentum
xL = pNz_Lab/pprz_inc # Frac of proton momentum
Q2_new = xBj/TwoPdotk
kQ2 = s_q/(xk*yk)
pNz_Lab = pprz_inc-pkz_Lab # Long neutron momentum
pkz_frac = pkz_Lab/pprz_inc # Frac of proton momentum
pNz_frac = pNz_Lab/pprz_inc # Frac of proton momentum
yplus = 1+((1-y)*(1-y))
tot_sigma = (sigma_tdis)*((TDIS_xbj*(Q2*Q2)*(137)*(137))/(2*math.pi*yplus))
# red_sig  = np.array([1.397,1.258,1.137,1.055,0.944,0.836,0.706,0.519])
# red_y = np.array([0.657,0.416,0.263,0.163,0.103,0.066,0.033,0.008])
# red_x = np.array([2.53e-4,4.0e-4,6.32e-4,1.02e-3,1.61e-3,2.52e-3,5.0e-3,2.10e-2])
# red_Q2 = np.array([15,15,15,15,15,15,15,15])
# red_yplus = 1+(1-red_y)*(1-red_y)
# my_sigma = (red_sig)*((2*math.pi*red_yplus)/(red_x*(red_Q2*red_Q2)*(137)*(137)))

Q2array = [7,15,30,60,120,240,480,1000]

for i,x in enumerate(Q2array) :
    Q2tmp = '{"Q2cut%i" : ((%0.1f < Q2) & (Q2 < %0.1f))}' % (i,Q2array[i]-5.0,Q2array[i]+5.0)
    print('{"Q2cut%i" : ((%0.1f < Q2) & (Q2 < %0.1f))}' % (i,Q2array[i]-5.0,Q2array[i]+5.0))
    cutDict.update(eval(Q2tmp))

ytmp = '{"ycut" : ((0.01 < y) & (y < 0.95))}'
ttmp = '{"tcut" : ((-1.00 < t) & (t < 0.00))}'
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

def phase_space():

    phaseSpace = c.densityPlot(xk, Q2, '$Q^2$ vs $x_{\k}$','$x_{\k}$','$Q^{2}$', 200, 200,  c, 0., 1.0, 0., 1000.)
    # plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(0.,1.)
    plt.ylim(0.,1000.)

    phaseSpace[1].savefig('OUTPUTS/phase_space.png')

def sigmavxk_Plot():
    
    f = plt.figure(figsize=(11.69,8.27))
    plt.style.use('classic')
    
    ax = f.add_subplot(421)
    xkscat1 = ax.errorbar(c.applyCuts(xk,cut7),c.applyCuts(tot_sigma,cut7),yerr=np.sqrt(c.applyCuts(tot_sigma,cut7))/200,fmt='o',label='$Q^2$=7 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=10 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    plt.ylabel('d$\sigma_{TDIS}$ ($nb/GeV^{2}$)')

    ax = f.add_subplot(422)
    xkscat2 = ax.errorbar(c.applyCuts(xk,cut15),c.applyCuts(tot_sigma,cut15),yerr=np.sqrt(c.applyCuts(tot_sigma,cut15))/150,fmt='o',label='$Q^2$=15 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=20 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    plt.title('d$\sigma_{TDIS}$ vs $x_\k$', fontsize =20)
    
    ax = f.add_subplot(423)
    xkscat3 = ax.errorbar(c.applyCuts(xk,cut30),c.applyCuts(tot_sigma,cut30),yerr=np.sqrt(c.applyCuts(tot_sigma,cut30))/200,fmt='o',label='$Q^2$=30 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    
    ax = f.add_subplot(424)
    xkscat4 = ax.errorbar(c.applyCuts(xk,cut60),c.applyCuts(tot_sigma,cut60),yerr=np.sqrt(c.applyCuts(tot_sigma,cut60))/200,fmt='o',label='$Q^2$=60 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=40 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    ax = f.add_subplot(425)
    xkscat5 = ax.errorbar(c.applyCuts(xk,cut120),c.applyCuts(tot_sigma,cut120),yerr=np.sqrt(c.applyCuts(tot_sigma,cut120))/200,fmt='o',label='$Q^2$=120 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=50 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    ax = f.add_subplot(426)
    xkscat6 = ax.errorbar(c.applyCuts(xk,cut240),c.applyCuts(tot_sigma,cut240),yerr=np.sqrt(c.applyCuts(tot_sigma,cut240))/200,fmt='o',label='$Q^2$=240 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=70 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(427)
    xkscat7 = ax.errorbar(c.applyCuts(xk,cut480),c.applyCuts(tot_sigma,cut480),yerr=np.sqrt(c.applyCuts(tot_sigma,cut480))/200,fmt='o',label='$Q^2$=480 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=80 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))    
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    ax = f.add_subplot(428)
    xkscat8 = ax.errorbar(c.applyCuts(xk,cut1000),c.applyCuts(tot_sigma,cut1000),yerr=np.sqrt(c.applyCuts(tot_sigma,cut1000))/200,fmt='o',label='$Q^2$=1000 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=90 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))    
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    plt.xlabel('$x_\k$')
    plt.tight_layout()

    plt.style.use('default')

    
def fkvxk_Plot():

    lumi = Lumi()
    
    f = plt.figure(figsize=(11.69,8.27))
    
    ax = f.add_subplot(421)
    xkscat1 = ax.errorbar(c.applyCuts(xk,cut7),c.applyCuts(fk,cut7),yerr=np.sqrt(lumi[1-1])/lumi[1-1],fmt='.',label='$Q^2$=7 $GeV^2$')
    plt.plot([0.001,0.01,0.1],[0.55,0.25,0.20], label="GRV fit")
    plt.xscale('log')
    plt.xlim(1e-4,1.)
    plt.ylim(0,0.5)
    ax.text(0.60, 0.85, '$Q^2$=7 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    plt.ylabel('$F^{k}_{2}$')

    ax = f.add_subplot(422)
    xkscat2 = ax.errorbar(c.applyCuts(xk,cut15),c.applyCuts(fk,cut15),yerr=np.sqrt(lumi[2-1])/lumi[2-1],fmt='.',label='$Q^2$=15 $GeV^2$')
    plt.plot([0.001,0.01,0.1],[0.80,0.30,0.20], label="GRV fit")
    plt.xscale('log')
    plt.xlim(1e-4,1.)
    plt.ylim(0,0.5)
    ax.text(0.60, 0.85, '$Q^2$=15 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    plt.title('$F^{k}_{2}$ vs $x_\k$', fontsize =20)
    
    ax = f.add_subplot(423)
    xkscat3 = ax.errorbar(c.applyCuts(xk,cut30),c.applyCuts(fk,cut30),yerr=np.sqrt(lumi[3-1])/lumi[3-1],fmt='.',label='$Q^2$=30 $GeV^2$')
    plt.plot([0.001,0.01,0.1],[0.85,0.45,0.20], label="GRV fit")
    plt.xscale('log')
    plt.xlim(1e-4,1.)
    plt.ylim(0,0.5)
    ax.text(0.60, 0.85, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    
    ax = f.add_subplot(424)
    xkscat4 = ax.errorbar(c.applyCuts(xk,cut60),c.applyCuts(fk,cut60),yerr=np.sqrt(lumi[4-1])/lumi[4-1],fmt='.',label='$Q^2$=60 $GeV^2$')
    plt.plot([0.001,0.01,0.1],[1.2,0.45,0.25], label="GRV fit")
    plt.xscale('log')
    plt.xlim(1e-4,1.)
    plt.ylim(0,0.5)
    ax.text(0.60, 0.85, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    ax = f.add_subplot(425)
    xkscat5 = ax.errorbar(c.applyCuts(xk,cut120),c.applyCuts(fk,cut120),yerr=np.sqrt(lumi[5-1])/lumi[5-1],fmt='.',label='$Q^2$=120 $GeV^2$')
    plt.plot([0.01,0.1,0.3],[0.5,0.25,0.15], label="GRV fit")
    plt.xscale('log')
    plt.xlim(1e-4,1.)
    plt.ylim(0,0.5)
    ax.text(0.60, 0.85, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))

    ax = f.add_subplot(426)
    xkscat6 = ax.errorbar(c.applyCuts(xk,cut240),c.applyCuts(fk,cut240),yerr=np.sqrt(lumi[6-1])/lumi[6-1],fmt='.',label='$Q^2$=240 $GeV^2$')
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit")
    plt.xscale('log')
    plt.xlim(1e-4,1.)
    plt.ylim(0,0.5)
    ax.text(0.60, 0.85, '$Q^2$=240 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(427)
    xkscat7 = ax.errorbar(c.applyCuts(xk,cut480),c.applyCuts(fk,cut480),yerr=np.sqrt(lumi[7-1])/lumi[7-1],fmt='.',label='$Q^2$=480 $GeV^2$')
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit")
    plt.xscale('log')
    plt.xlim(1e-4,1.)
    plt.ylim(0,0.5)
    ax.text(0.60, 0.85, '$Q^2$=480 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    # ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))

    ax = f.add_subplot(428)
    xkscat8 = ax.errorbar(c.applyCuts(xk,cut1000),c.applyCuts(fk,cut1000),yerr=np.sqrt(lumi[8-1])/lumi[8-1],fmt='.',label='$Q^2$=1000 $GeV^2$')
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit")
    plt.xscale('log')
    plt.xlim(1e-4,1.)
    plt.ylim(0,0.5)
    ax.text(0.60, 0.85, '$Q^2$=1000 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    
    plt.xlabel('log($x_\k$)')
    plt.tight_layout()

    plt.subplots_adjust(hspace=0.0,wspace=0.0)

    f.savefig('OUTPUTS/fkvxk.png')

def fkvxk_Plot_nolog():

    lumi = Lumi()
    
    f = plt.figure(figsize=(11.69,8.27))
    
    ax = f.add_subplot(421)
    xkscat1 = ax.errorbar(c.applyCuts(xk,cut7),c.applyCuts(fk,cut7),yerr=np.sqrt(lumi[1-1])/lumi[1-1],fmt='.',label='$Q^2$=7 $GeV^2$')
    plt.plot([0.001,0.01,0.1],[0.55,0.25,0.20], label="GRV fit")
    plt.xlim(0.,1.)
    plt.ylim(0,0.5)
    ax.text(0.60, 0.85, '$Q^2$=7 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    plt.ylabel('$F^{k}_{2}$')

    ax = f.add_subplot(422)
    xkscat2 = ax.errorbar(c.applyCuts(xk,cut15),c.applyCuts(fk,cut15),yerr=np.sqrt(lumi[2-1])/lumi[2-1],fmt='.',label='$Q^2$=15 $GeV^2$')
    plt.plot([0.001,0.01,0.1],[0.80,0.30,0.20], label="GRV fit")
    plt.xlim(0.,1.)
    plt.ylim(0,0.5)
    ax.text(0.60, 0.85, '$Q^2$=15 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    plt.title('$F^{k}_{2}$ vs $x_\k$', fontsize =20)
    
    ax = f.add_subplot(423)
    xkscat3 = ax.errorbar(c.applyCuts(xk,cut30),c.applyCuts(fk,cut30),yerr=np.sqrt(lumi[3-1])/lumi[3-1],fmt='.',label='$Q^2$=30 $GeV^2$')
    plt.plot([0.001,0.01,0.1],[0.85,0.45,0.20], label="GRV fit")
    plt.xlim(0.,1.)
    plt.ylim(0,0.5)
    ax.text(0.60, 0.85, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    
    ax = f.add_subplot(424)
    xkscat4 = ax.errorbar(c.applyCuts(xk,cut60),c.applyCuts(fk,cut60),yerr=np.sqrt(lumi[4-1])/lumi[4-1],fmt='.',label='$Q^2$=60 $GeV^2$')
    plt.plot([0.001,0.01,0.1],[1.2,0.45,0.25], label="GRV fit")
    plt.xlim(0.,1.)
    plt.ylim(0,0.5)
    ax.text(0.60, 0.85, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    ax = f.add_subplot(425)
    xkscat5 = ax.errorbar(c.applyCuts(xk,cut120),c.applyCuts(fk,cut120),yerr=np.sqrt(lumi[5-1])/lumi[5-1],fmt='.',label='$Q^2$=120 $GeV^2$')
    plt.plot([0.01,0.1,0.3],[0.5,0.25,0.15], label="GRV fit")
    plt.xlim(0.,1.)
    plt.ylim(0,0.5)
    ax.text(0.60, 0.85, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))

    ax = f.add_subplot(426)
    xkscat6 = ax.errorbar(c.applyCuts(xk,cut240),c.applyCuts(fk,cut240),yerr=np.sqrt(lumi[6-1])/lumi[6-1],fmt='.',label='$Q^2$=240 $GeV^2$')
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit")
    plt.xlim(0.,1.)
    plt.ylim(0,0.5)
    ax.text(0.60, 0.85, '$Q^2$=240 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(427)
    xkscat7 = ax.errorbar(c.applyCuts(xk,cut480),c.applyCuts(fk,cut480),yerr=np.sqrt(lumi[7-1])/lumi[7-1],fmt='.',label='$Q^2$=480 $GeV^2$')
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit")
    plt.xlim(0.,1.)
    plt.ylim(0,0.5)
    ax.text(0.60, 0.85, '$Q^2$=480 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower',nbins=5))

    ax = f.add_subplot(428)
    xkscat8 = ax.errorbar(c.applyCuts(xk,cut1000),c.applyCuts(fk,cut1000),yerr=np.sqrt(lumi[8-1])/lumi[8-1],fmt='.',label='$Q^2$=1000 $GeV^2$')
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit")
    plt.xlim(0.,1.)
    plt.ylim(0,0.5)
    ax.text(0.60, 0.85, '$Q^2$=1000 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower',nbins=5))
    
    plt.xlabel('$x_\k$')
    plt.tight_layout()

    plt.subplots_adjust(hspace=0.0,wspace=0.0)

    f.savefig('OUTPUTS/fkvxk_nolog.png')

    
def fkvt_Plot():
    
    lumi = Lumi()

    f = plt.figure(figsize=(11.69,8.27))
    
    ax = f.add_subplot(421)
    tscat1 = ax.errorbar(c.applyCuts(t,cut7),c.applyCuts(fk,cut7),yerr=np.sqrt(lumi[1-1])/lumi[1-1],fmt='.',label='$Q^2$=7 $GeV^2$')
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=7 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    plt.ylabel('log($F^{k}_{2}$)')

    ax = f.add_subplot(422)
    tscat2 = ax.errorbar(c.applyCuts(t,cut15),c.applyCuts(fk,cut15),yerr=np.sqrt(lumi[2-1])/lumi[2-1],fmt='.',label='$Q^2$=15 $GeV^2$')
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=15 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_formatter(plt.NullFormatter())

    plt.title('log($F^{k}_{2}$) vs -t', fontsize =20)
    
    ax = f.add_subplot(423)
    tscat3 = ax.errorbar(c.applyCuts(t,cut30),c.applyCuts(fk,cut30),yerr=np.sqrt(lumi[3-1])/lumi[3-1],fmt='.',label='$Q^2$=30 $GeV^2$')
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    
    ax = f.add_subplot(424)
    tscat4 = ax.errorbar(c.applyCuts(t,cut60),c.applyCuts(fk,cut60),yerr=np.sqrt(lumi[4-1])/lumi[4-1],fmt='.',label='$Q^2$=60 $GeV^2$')
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    
    ax = f.add_subplot(425)
    tscat5 = ax.errorbar(c.applyCuts(t,cut120),c.applyCuts(fk,cut120),yerr=np.sqrt(lumi[5-1])/lumi[5-1],fmt='.',label='$Q^2$=120 $GeV^2$')
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    ax = f.add_subplot(426)
    tscat6 = ax.errorbar(c.applyCuts(t,cut240),c.applyCuts(fk,cut240),yerr=np.sqrt(lumi[6-1])/lumi[6-1],fmt='.',label='$Q^2$=240 $GeV^2$')
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=240 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(427)
    tscat7 = ax.errorbar(c.applyCuts(t,cut480),c.applyCuts(fk,cut480),yerr=np.sqrt(lumi[7-1])/lumi[7-1],fmt='.',label='$Q^2$=480 $GeV^2$')
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=480 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))

    ax = f.add_subplot(428)
    tscat8 = ax.errorbar(c.applyCuts(t,cut1000),c.applyCuts(fk,cut1000),yerr=np.sqrt(lumi[8-1])/lumi[8-1],fmt='.',label='$Q^2$=1000 $GeV^2$')
    plt.yscale('log')
    plt.xlim(-0.4,0.)
    plt.ylim(1e-8,1.0)
    ax.text(0.65, 0.65, '$Q^2$=1000 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    plt.xlabel('-t')
    plt.tight_layout()

    plt.subplots_adjust(hspace=0.0,wspace=0.0)

    f.savefig('OUTPUTS/fkvt.png')

    
def bin_sigma(xbin,cut):
    
    sigma = c.applyCuts(tot_sigma,cut)
    xk_cut = c.applyCuts(xk,cut)

    xsigma = [sig for sig,x in zip(sigma,xk_cut) if xbin-0.0005 <= x <= xbin+0.0005]
    
    return xsigma

def Lumi():

    # Luminosity
    lumi = []

    for i,q in enumerate(Q2array):
        lumi_xk = []
        k = 0
        for j,x in enumerate(xkarray):
            xkcut = x

            tmp  = "cut%i" % (q)
            cut = eval(tmp)
            
            sigval = bin_sigma(xkcut,cut)
            
            evts = len(sigval)
            
            if evts == 0:
                k+=1
                continue
        
            binx = 0.01 # bin x size
            binQ2 = 10.0 # bin Q2 size
            print("---------------------------------")
            print("Events in bin: ",evts)
            print("\nbinx size: ", binx, "\nbinQ2 size: ", binQ2)
            print("\nbinx: ", xkcut, "\nbinQ2: ", cut[0],"\n")

            avgSig = np.average(sigval)
            
            for l,sig in enumerate(sigval):
                print("Sigma: ", sig)
                if evts/((avgSig)*(binQ2)*(binx)) < 1e4:
                    lumi_xk.append(evts/((avgSig)*(binQ2)*(binx)))
                else:
                    lumi_xk.append(np.average(lumi_xk))
                # print("Luminosity: ", lumi_xk)
        
            print("---------------------------------\n")
        lumi.append(np.asarray(lumi_xk))

    print("\nLuminosity: ", lumi)

    # nevt = [10/(lum*1e-6) for i,lum in enumerate(lumi)]
    nevt = [10/(lum*1e-6) for i,lum in enumerate(lumi)]
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

    return nevt

    
def main() :

    phase_space()
    fkvxk_Plot()
    fkvxk_Plot_nolog()
    fkvt_Plot()
    sigmavxk_Plot()
    plt.show()
    
if __name__=='__main__': main()
