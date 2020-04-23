#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2020-04-22 16:11:31 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

import uproot as up
import scipy as sci
from scipy import optimize
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from matplotlib.ticker import MaxNLocator
from collections import namedtuple
from sys import path
import time,math,sys

sys.path.insert(0,'/home/trottar/bin/python/')
import root2py as r2p

# rootName="/home/trottar/ResearchNP/JLEIC/USERS/trottar/OUTPUTS/pi_p_5on100.root"
# rootName="/home/trottar/ResearchNP/JLEIC/USERS/trottar/OUTPUTS/pi_n_5on100.root"
# rootName="/home/trottar/ResearchNP/JLEIC/USERS/trottar/OUTPUTS/k_lambda_5on100.root"
rootName="/home/trottar/ResearchNP/JLEIC/USERS/trottar/OUTPUTS/k_lambda_18on275.root"
tree = up.open(rootName)["Evnts"]
branch = r2p.pyBranch(tree)

pdf = matplotlib.backends.backend_pdf.PdfPages("fpiPlot_xpi.pdf")

# Define phyisics data
s_e = branch.findBranch("invts","s_e")
s_q = branch.findBranch("invts","s_q")
Q2 = branch.findBranch("invts","Q2")
xBj = branch.findBranch("invts","xBj")
t = -branch.findBranch("invts","tPrime")
y_D = branch.findBranch("invts","y_D")
nu = branch.findBranch("invts","nu")
TwoPdotk = branch.findBranch("invts","TwoPdotk")
TDIS_xbj = tree.array("TDIS_xbj")
sigma_dis = tree.array("sigma_dis")*(1e-5)
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
# fpi = tree.array("fpi")
fpi = 0.361*tree.array("f2N")
f2N = tree.array("f2N")
xpi = tree.array("xpi")
ypi = tree.array("ypi")
tpi = tree.array("tpi")
escat = tree.array("EScatRest")
pprz_inc = tree.array("pprz_inc")
ppiz_Lab = tree.array("ppiz_Lab")

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

Q2array =[10.,20.,30.,40.,50.,70.,80.,80.,100.]

cutDict = {}

i=0
for x in range(0,len(Q2array)) :
    Q2tmp = '{"Q2cut%i" : ((%0.1f < Q2) & (Q2 < %0.1f))}' % (i,Q2array[i]-2.5,Q2array[i]+2.5)
    ytmp = '{"ycut" : ((0.01 < y) & (y < 0.95))}'
    ttmp = '{"tcut" : ((-1.00 < t) & (t < 0.00))}'
    print '{"Q2cut%i" : ((%0.1f < Q2) & (Q2 < %0.1f))}' % (i,Q2array[i]-2.5,Q2array[i]+2.5)
    cutDict.update(eval(Q2tmp))
    cutDict.update(eval(ttmp))
    cutDict.update(eval(ytmp))
    i+=1
    
c = r2p.pyPlot(cutDict)

ycut1 = ["ycut"]
tcut1 = ["tcut"]
cut10 = ["Q2cut0","tcut","ycut"]
cut20 = ["Q2cut1","tcut","ycut"]
cut30 = ["Q2cut2","tcut","ycut"]
cut40 = ["Q2cut3","tcut","ycut"]
cut50 = ["Q2cut4","tcut","ycut"]
cut70 = ["Q2cut5","tcut","ycut"]
cut80 = ["Q2cut6","tcut","ycut"]
cut90 = ["Q2cut7","tcut","ycut"]
cut100 = ["Q2cut8","tcut","ycut"]

def sigmavxpi_Plot():
    
    f = plt.figure(figsize=(11.69,8.27))
    plt.style.use('classic')
    
    ax = f.add_subplot(451)
    xpiscat1 = ax.errorbar(c.applyCuts(xpi,cut10),c.applyCuts(tot_sigma,cut10),xerr=np.sqrt(c.applyCuts(xpi,cut10))/200,yerr=np.sqrt(c.applyCuts(tot_sigma,cut10))/200,fmt='o',label='$Q^2$=10 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=10 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    plt.ylabel('d$\sigma_{TDIS}$ ($nb/GeV^{2}$)')

    ax = f.add_subplot(452)
    xpiscat2 = ax.errorbar(c.applyCuts(xpi,cut20),c.applyCuts(tot_sigma,cut20),xerr=np.sqrt(c.applyCuts(xpi,cut20))/200,yerr=np.sqrt(c.applyCuts(tot_sigma,cut20))/200,fmt='o',label='$Q^2$=20 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=20 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    ax = f.add_subplot(453)
    xpiscat3 = ax.errorbar(c.applyCuts(xpi,cut30),c.applyCuts(tot_sigma,cut30),xerr=np.sqrt(c.applyCuts(xpi,cut30))/200,yerr=np.sqrt(c.applyCuts(tot_sigma,cut30))/200,fmt='o',label='$Q^2$=30 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    plt.title('d$\sigma_{TDIS}$ vs $x_\pi$', fontsize =20)
    
    ax = f.add_subplot(454)
    xpiscat4 = ax.errorbar(c.applyCuts(xpi,cut40),c.applyCuts(tot_sigma,cut40),xerr=np.sqrt(c.applyCuts(xpi,cut40))/200,yerr=np.sqrt(c.applyCuts(tot_sigma,cut40))/200,fmt='o',label='$Q^2$=40 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=40 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    ax = f.add_subplot(455)
    xpiscat5 = ax.errorbar(c.applyCuts(xpi,cut50),c.applyCuts(tot_sigma,cut50),xerr=np.sqrt(c.applyCuts(xpi,cut50))/200,yerr=np.sqrt(c.applyCuts(tot_sigma,cut50))/200,fmt='o',label='$Q^2$=50 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=50 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(456)
    xpiscat6 = ax.errorbar(c.applyCuts(xpi,cut70),c.applyCuts(tot_sigma,cut70),xerr=np.sqrt(c.applyCuts(xpi,cut70))/200,yerr=np.sqrt(c.applyCuts(tot_sigma,cut70))/200,fmt='o',label='$Q^2$=70 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=70 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))

    ax = f.add_subplot(457)
    xpiscat7 = ax.errorbar(c.applyCuts(xpi,cut80),c.applyCuts(tot_sigma,cut80),xerr=np.sqrt(c.applyCuts(xpi,cut80))/200,yerr=np.sqrt(c.applyCuts(tot_sigma,cut80))/200,fmt='o',label='$Q^2$=80 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=80 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(458)
    xpiscat8 = ax.errorbar(c.applyCuts(xpi,cut90),c.applyCuts(tot_sigma,cut90),xerr=np.sqrt(c.applyCuts(xpi,cut90))/200,yerr=np.sqrt(c.applyCuts(tot_sigma,cut90))/200,fmt='o',label='$Q^2$=90 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=90 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(459)
    xpiscat9 = ax.errorbar(c.applyCuts(xpi,cut100),c.applyCuts(tot_sigma,cut100),xerr=np.sqrt(c.applyCuts(xpi,cut100))/200,yerr=np.sqrt(c.applyCuts(tot_sigma,cut100))/200,fmt='o',label='$Q^2$=100 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=100 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    plt.xlabel('$x_\pi$')
    plt.tight_layout()

    plt.style.use('default')
    
def sigmavxBj_Plot():

    f = plt.figure(figsize=(11.69,8.27))
    plt.style.use('classic')   
    
    ax = f.add_subplot(451)
    # xBjscat1 = ax.errorbar(c.applyCuts(TDIS_xbj,cut10),c.applyCuts(tot_sigma,cut10),xerr=np.sqrt(c.applyCuts(TDIS_xbj,cut10))/200,yerr=np.sqrt(c.applyCuts(tot_sigma,cut10))/200,fmt='o',label='$Q^2$=10 $GeV^2$')
    xBjscat1 = ax.scatter(c.applyCuts(TDIS_xbj,cut10),c.applyCuts(tot_sigma,cut10),label='$Q^2$=10 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    plt.xlim(1e-3,1.0)
    plt.ylim(0.,7.)
    ax.text(0.65, 0.25, '$Q^2$=10 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))    
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    plt.title('d$\sigma_{TDIS}$ vs TDIS_xbj', fontsize =20)
    plt.ylabel('d$\sigma_{TDIS}$ ($nb/GeV^{2}$)')


    ax = f.add_subplot(452)
    # xBjscat2 = ax.errorbar(c.applyCuts(TDIS_xbj,cut20),c.applyCuts(tot_sigma,cut20),xerr=np.sqrt(c.applyCuts(TDIS_xbj,cut20))/200,yerr=np.sqrt(c.applyCuts(tot_sigma,cut20))/200,fmt='o',label='$Q^2$=20 $GeV^2$')
    xBjscat2 = ax.scatter(c.applyCuts(TDIS_xbj,cut20),c.applyCuts(tot_sigma,cut20),label='$Q^2$=20 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    plt.xlim(1e-3,1.0)
    plt.ylim(0.,7.)
    ax.text(0.65, 0.25, '$Q^2$=20 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))    
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))

    ax = f.add_subplot(453)
    # xBjscat3 = ax.errorbar(c.applyCuts(TDIS_xbj,cut30),c.applyCuts(tot_sigma,cut30),xerr=np.sqrt(c.applyCuts(TDIS_xbj,cut30))/200,yerr=np.sqrt(c.applyCuts(tot_sigma,cut30))/200,fmt='o',label='$Q^2$=30 $GeV^2$')
    xBjscat3 = ax.scatter(c.applyCuts(TDIS_xbj,cut30),c.applyCuts(tot_sigma,cut30),label='$Q^2$=30 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    plt.xlim(1e-3,1.0)
    plt.ylim(0.,7.)
    ax.text(0.65, 0.25, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))    
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    
    ax = f.add_subplot(454)
    # xBjscat4 = ax.errorbar(c.applyCuts(TDIS_xbj,cut40),c.applyCuts(tot_sigma,cut40),xerr=np.sqrt(c.applyCuts(TDIS_xbj,cut40))/200,yerr=np.sqrt(c.applyCuts(tot_sigma,cut40))/200,fmt='o',label='$Q^2$=40 $GeV^2$')
    xBjscat4 = ax.scatter(c.applyCuts(TDIS_xbj,cut40),c.applyCuts(tot_sigma,cut40),label='$Q^2$=40 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    plt.xlim(1e-3,1.0)
    plt.ylim(0.,7.)
    ax.text(0.65, 0.25, '$Q^2$=40 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))    
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    
    ax = f.add_subplot(455)
    # xBjscat5 = ax.errorbar(c.applyCuts(TDIS_xbj,cut50),c.applyCuts(tot_sigma,cut50),xerr=np.sqrt(c.applyCuts(TDIS_xbj,cut50))/200,yerr=np.sqrt(c.applyCuts(tot_sigma,cut50))/200,fmt='o',label='$Q^2$=50 $GeV^2$')
    xBjscat5 = ax.scatter(c.applyCuts(TDIS_xbj,cut50),c.applyCuts(tot_sigma,cut50),label='$Q^2$=50 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    plt.xlim(1e-3,1.0)
    plt.ylim(0.,.07)
    ax.text(0.65, 0.25, '$Q^2$=50 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))    
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))        
        
    ax = f.add_subplot(456)
    # xBjscat6 = ax.errorbar(c.applyCuts(TDIS_xbj,cut70),c.applyCuts(tot_sigma,cut70),xerr=np.sqrt(c.applyCuts(TDIS_xbj,cut70))/200,yerr=np.sqrt(c.applyCuts(tot_sigma,cut70))/200,fmt='o',label='$Q^2$=70 $GeV^2$')
    xBjscat6 = ax.scatter(c.applyCuts(TDIS_xbj,cut70),c.applyCuts(tot_sigma,cut70),label='$Q^2$=70 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,0.6)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=70 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    ax = f.add_subplot(457)

    # xBjscat7 = ax.errorbar(c.applyCuts(TDIS_xbj,cut80),c.applyCuts(tot_sigma,cut80),xerr=np.sqrt(c.applyCuts(TDIS_xbj,cut80))/200,yerr=np.sqrt(c.applyCuts(tot_sigma,cut80))/200,fmt='o',label='$Q^2$=80 $GeV^2$')
    xBjscat7 = ax.scatter(c.applyCuts(TDIS_xbj,cut80),c.applyCuts(tot_sigma,cut80),label='$Q^2$=80 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,0.6)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=80 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
        
    ax = f.add_subplot(458)
    # xBjscat8 = ax.errorbar(c.applyCuts(TDIS_xbj,cut90),c.applyCuts(tot_sigma,cut90),xerr=np.sqrt(c.applyCuts(TDIS_xbj,cut90))/200,yerr=np.sqrt(c.applyCuts(tot_sigma,cut90))/200,fmt='o',label='$Q^2$=90 $GeV^2$')
    xBjscat8 = ax.scatter(c.applyCuts(TDIS_xbj,cut90),c.applyCuts(tot_sigma,cut90),label='$Q^2$=90 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,0.6)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=90 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    ax = f.add_subplot(459)
    # xBjscat9 = ax.errorbar(c.applyCuts(TDIS_xbj,cut100),c.applyCuts(tot_sigma,cut100),xerr=np.sqrt(c.applyCuts(TDIS_xbj,cut100))/200,yerr=np.sqrt(c.applyCuts(tot_sigma,cut100))/200,fmt='o',label='$Q^2$=100 $GeV^2$')
    xBjscat9 = ax.scatter(c.applyCuts(TDIS_xbj,cut100),c.applyCuts(tot_sigma,cut100),label='$Q^2$=100 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,0.6)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=100 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    plt.xlabel('TDIS_xbj')
    plt.tight_layout()
    plt.setp(ax.get_xticklabels()[0], visible=False)
    plt.setp(ax.get_xticklabels()[-1], visible=False)
        
    # # Reset figure style, otherwise next plot will look weird
    plt.style.use('default')
    
    plt.close(f)

    return [xBjscat1,xBjscat2,xBjscat3,xBjscat4,xBjscat5,xBjscat6,xBjscat7,xBjscat8,xBjscat9]


def Lumi():
    [xBjscat1,xBjscat2,xBjscat3,xBjscat4,xBjscat5,xBjscat6,xBjscat7,xBjscat8,xBjscat9]=sigmavxBj_Plot()

    # Luminosity
    lumi = []
    uncern = []

    # xEvts,Q2Evts
    for j in range(0,8):

        tmp = 'xBjscat%i.get_offsets()' % (j+1)
        # print tmp
        binValue = eval(tmp)
        evts = len(binValue)
        
        xval = []
        sigval = []
    
        for i in range(0,evts):
            xval.append(binValue[i][0])
            sigval.append(binValue[i][1])
        
        binx = 0.1 # bin x size
        binQ2 = 10.0 # bin Q2 size
        print "---------------------------------"
        print "Events in bin: ",evts
        # print "\nBin value array: ", binValue
        print "\nbinx: ", binx, "\nbinQ2: ", binQ2
        
        avgSigma = np.average(sigval)

        print "\nAverage Sigma: ", avgSigma
        
        lumi = np.append(lumi,(evts/(avgSigma*(binQ2)*(binx))))
        
        # uncern = np.append(uncern,np.sqrt(avgLumi*avgSigma*(binQ2)*(binx))/evts)
        # uncern = np.append(uncern,np.sqrt(avgLumi*avgSigma*(binQ2)*(binx))/evts)
        
        print "Plot: ", j+1
        print "Luminosity: ", lumi[j]
        # print "Uncertainty: ", uncern[j]
        print "---------------------------------\n"

    avgLumi = np.average(lumi)
        
    print "\nLuminosity: ", lumi
    print "\nAverage Luminosity: ", avgLumi
    
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    scat1 = ax.hist(lumi,bins=c.setbin(lumi,20),label='Low',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.title('Luminosity', fontsize =20)
    plt.xlabel('Luminosity')

    plt.close(f)    

    return lumi

def namedtuple_to_str(t, field_widths=15):
    if isinstance(field_widths, int):
        field_widths = [field_widths] * len(t)
    field_pairs = ['{}={}'.format(field, getattr(t, field)) for field in t._fields]
    s = ' '.join('{{:{}}}'.format(w).format(f) for w,f in zip(field_widths, field_pairs))
    result = '{}( {} )'.format(type(t).__name__, s)
    return result

def piSF_GRV(x):

    # Quark charge
    qd = -1/3
    qd_bar = 1/3
    qu = 2/3
    qu_bar = -2/3
    qs = -1/3
    qs_bar = 1/3

    # x is included in contributions
    xu_contr = 1.377*(x**0.549)*(1+0.81*(np.sqrt(x))-4.36*x+19.4*(x**0.75))*((1-x)**3.027)
    xd_contr = 0.328*(x**0.366)*(1+1.114*(np.sqrt(x))-5.71*x+16.9*(x**0.75))*((1-x)**3.774)

    # 
    fit_f2pi = (qu**2)*xu_contr + (qd**2)*xd_contr
    
    return fit_f2pi

def piSF_BLFQ(x):
    
    # PDF pion valence quarks, a=b since m_u=m_d. a and b are parameters for fitting the PDF at fixed L_max (the basis resolution in the longitudinal direction)
    u_contr = []
    d_contr = []
    
    # Heavy flavor contributions
    s_contr = []
    
    # a and b are the extroplated values at L_max->inf (J. Lan et. al., arXiv preprint (2019) arXiv:1907.01509)
    a=b=0.5961

    # c and d are the extroplated values at L_max->inf (J. Lan et. al., arXiv preprint (2019) arXiv:1907.01509) but unlike pion they are different for strange quark contributions
    c=0.6337
    d=0.8546

    # Quark charge
    qd = -1/3
    qd_bar = 1/3
    qu = 2/3
    qu_bar = -2/3
    qs = -1/3
    qs_bar = 1/3

    for i,evt in enumerate(x):
        u_contr.append(((x[i]**a)*(1-x[i])**b)/sci.special.beta(a+1,b+1))
    for i,evt in enumerate(x):
        d_contr.append(((x[i]**a)*(1-x[i])**b)/sci.special.beta(a+1,b+1))        
    for i,evt in enumerate(x):
        s_contr.append(((x[i]**c)*(1-x[i])**d)/sci.special.beta(c+1,d+1))
    
    # f2pi is related to its PDF by the below equation, which is summed over quark flavor.
    fit_f2pi = []

    for i,evt in enumerate(x):
        fit_f2pi.append((qu**2)*(x[i])*u_contr[i]+(qd**2)*(x[i])*d_contr[i]+(qs**2)*(x[i])*s_contr[i])

    return fit_f2pi

def curve_fit(x, a, b):
    return a*np.log(b*x)

lumi = Lumi()

def pionPlots(Q2_inp):

    tot_lumi = lumi*(1e-6) # nano to femto
    
    lum10 = []
    lum100 = []

    if(Q2_inp==10):
        apply_cut=cut10
        Q2_val=10
    elif(Q2_inp==30):
        apply_cut=cut30
        Q2_val=30
    elif(Q2_inp==50):
        apply_cut=cut50
        Q2_val=50
    elif(Q2_inp==70):
        apply_cut=cut70
        Q2_val=70

    f = plt.figure(figsize=(11.69,8.27))
    plt.style.use('classic')

    ax = f.add_subplot(331)
    numBin1 = ax.hist(c.applyCuts(xpi,apply_cut),bins=c.setbin(xpi,200,0.,1.),histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_{\pi}$=(0.20,0.30)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.20,0.30)
    ax.text(0.65, 0.95, '$x_{\pi}$=(0.20,0.30)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts1 = []
    for val in numBin1[0]:
        if 0.20 < numBin1[1][i] < 0.30 :
            maxEvts1.append(numBin1[0][i])
        i+=1
    numEvts1 = np.sum(maxEvts1)
    lum10.append((10/(tot_lumi[1-1]))*numEvts1)
    lum100.append((100/(tot_lumi[1-1]))*numEvts1)
    
    plt.title('$x_{\pi}$ plots [$Q^2$ = 10 $GeV^2$]', fontsize =20)
    
    ax = f.add_subplot(332)
    numBin2 = ax.hist(c.applyCuts(xpi,apply_cut),bins=c.setbin(xpi,200,0.,1.),histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_{\pi}$=(0.30,0.40)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.30,0.40)
    ax.text(0.65, 0.95, '$x_{\pi}$=(0.30,0.40)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts2 = []
    for val in numBin2[0]:
        if 0.30 < numBin2[1][i] < 0.40 :
            maxEvts2.append(numBin2[0][i])
        i+=1
    numEvts2 = np.sum(maxEvts2)
    lum10.append((10/(tot_lumi[2-1]))*numEvts2)
    lum100.append((100/(tot_lumi[2-1]))*numEvts2)
    
    ax = f.add_subplot(333)
    numBin3 = ax.hist(c.applyCuts(xpi,apply_cut),bins=c.setbin(xpi,200,0.,1.),histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_{\pi}$=(0.40,0.50)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.40,0.50)
    ax.text(0.65, 0.95, '$x_{\pi}$=(0.40,0.50)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts3 = []
    for val in numBin3[0]:
        if 0.40 < numBin3[1][i] < 0.50 :
            maxEvts3.append(numBin3[0][i])
        i+=1
    numEvts3 = np.sum(maxEvts3)
    lum10.append((10/(tot_lumi[3-1]))*numEvts3)
    lum100.append((100/(tot_lumi[3-1]))*numEvts3)
    
    ax = f.add_subplot(334)
    numBin4 = ax.hist(c.applyCuts(xpi,apply_cut),bins=c.setbin(xpi,200,0.,1.),histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_{\pi}$=(0.50,0.60)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.50,0.60)
    ax.text(0.65, 0.95, '$x_{\pi}$=(0.50,0.60)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts4 = []
    for val in numBin4[0]:
        if 0.50 < numBin4[1][i] < 0.60 :
            maxEvts4.append(numBin4[0][i])
        i+=1
    numEvts4 = np.sum(maxEvts4)
    lum10.append((10/(tot_lumi[4-1]))*numEvts4)
    lum100.append((100/(tot_lumi[4-1]))*numEvts4)
    
    ax = f.add_subplot(335)
    numBin5 = ax.hist(c.applyCuts(xpi,apply_cut),bins=c.setbin(xpi,200,0.,1.),histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_{\pi}$=(0.60,0.70)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.60,0.70)
    ax.text(0.65, 0.95, '$x_{\pi}$=(0.60,0.70)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts5 = []
    for val in numBin5[0]:
        if 0.60 < numBin5[1][i] < 0.70 :
            maxEvts5.append(numBin5[0][i])
        i+=1
    numEvts5 = np.sum(maxEvts5)
    lum10.append((10/(tot_lumi[5-1]))*numEvts5)
    lum100.append((100/(tot_lumi[5-1]))*numEvts5)
    
    ax = f.add_subplot(336)
    numBin6 = ax.hist(c.applyCuts(xpi,apply_cut),bins=c.setbin(xpi,200,0.,1.),histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_{\pi}$=(0.70,0.80)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.70,0.80)
    ax.text(0.65, 0.95, '$x_{\pi}$=(0.70,0.80)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts6 = []
    for val in numBin6[0]:
        if 0.70 < numBin6[1][i] < 0.80 :
            maxEvts6.append(numBin6[0][i])
        i+=1
    numEvts6 = np.sum(maxEvts6)
    lum10.append((10/(tot_lumi[6-1]))*numEvts6)
    lum100.append((100/(tot_lumi[6-1]))*numEvts6)
    
    ax = f.add_subplot(337)
    numBin7 = ax.hist(c.applyCuts(xpi,apply_cut),bins=c.setbin(xpi,200,0.,1.),histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_{\pi}$=(0.80,0.90)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.80,0.90)
    ax.text(0.65, 0.95, '$x_{\pi}$=(0.80,0.90)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts7 = []
    for val in numBin7[0]:
        if 0.80 < numBin7[1][i] < 0.90 :
            maxEvts7.append(numBin7[0][i])
        i+=1
    numEvts7 = np.sum(maxEvts7)
    lum10.append((10/(tot_lumi[7-1]))*numEvts7)
    lum100.append((100/(tot_lumi[7-1]))*numEvts7)
    
    ax = f.add_subplot(338)
    numBin8 = ax.hist(c.applyCuts(xpi,apply_cut),bins=c.setbin(xpi,200,0.,1.),histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_{\pi}$=(0.90,1.00)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.90,1.00)
    ax.text(0.65, 0.95, '$x_{\pi}$=(0.90,1.00)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts8 = []
    for val in numBin8[0]:
        if 0.90 < numBin8[1][i] < 1.00 :
            maxEvts8.append(numBin8[0][i])
        i+=1
    numEvts8 = np.sum(maxEvts8)
    lum10.append((10/(tot_lumi[8-1]))*numEvts8)
    lum100.append((100/(tot_lumi[8-1]))*numEvts8)
    
    print "\n\nLumi10: ", lum10
    print "\n\nLumi100: ", lum100, "\n\n"

    fout = open('/home/trottar/ResearchNP/JLEIC/USERS/trottar/OUTPUTS/LuminosityTable.txt','a') 
    
    lum_tuple = namedtuple('lum_tuple',['binQ2','binxpi','lumi','nbin','nbin_10','nbin_100','uncern_10','uncern_100'])
    
    lum_data = []

    tot_data = []

    count = [1,2,3,4,5,6,7,8]
    
    xpiVal = [0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95]
    
    Q2Val = [Q2_inp] * len(count)

    uncern  = np.sqrt(lum10)/lum10
    uncern100  = np.sqrt(lum100)/lum100
    
    numEvts = []
    
    for i in range(0,len(uncern)) :
        tmp = 'numEvts%i' % count[i]
        numEvts.append(eval(tmp))

    for i in range(0,len(uncern)) :
        lum_data.append(lum_tuple(Q2Val[i],xpiVal[i],tot_lumi[i],numEvts[i],lum10[i],lum100[i],uncern[i], uncern100[i]))
        tot_data.append(lum_data[i])

    for t in tot_data:
        fout.write("%s\n" % namedtuple_to_str(t))

    fout.close()
    
    plt.style.use('default')
    plt.close(f)
    
    # Combined plots

    xpi_tot = []
    fpi_tot = []
    fpiuncern = []
    
    xpi1 = []
    fpi1 = []
    i=0
    for evt in c.applyCuts(xpi,apply_cut):
        if 0.240 < evt < 0.260 :
            xpi1.append(evt)
            fpi1.append(c.applyCuts(fpi,apply_cut)[i])
        i+=1
    if len(fpi1) > 1:
        xpi_tot.append(min(xpi1))
        fpi_tot.append(min(fpi1))
        # fpiuncern.append(min(fpi1)*uncern[0])
        fpiuncern.append(uncern[0])
    elif len(fpi1) == 1:
        xpi_tot.append((xpi1[0]))
        fpi_tot.append((fpi1[0]))
        # fpiuncern.append(min(fpi1)*uncern[0])
        fpiuncern.append(uncern[0])
    else:
        print "Invalid at ", xpi1
        
    xpi2 = []
    fpi2 = []
    i=0
    for evt in c.applyCuts(xpi,apply_cut):
        if 0.340 < evt < 0.360 :
            xpi2.append(evt)
            fpi2.append(c.applyCuts(fpi,apply_cut)[i])
        i+=1
    if len(fpi2) > 1:
        xpi_tot.append(min(xpi2))
        fpi_tot.append(min(fpi2))
        # fpiuncern.append(min(fpi2)*uncern[1])
        fpiuncern.append(uncern[1])
    elif len(fpi2) == 1:
        xpi_tot.append((xpi2[0]))
        fpi_tot.append((fpi2[0]))
        # fpiuncern.append(min(fpi2)*uncern[1])
        fpiuncern.append(uncern[1])
    else:
        print "Invalid at ", xpi2
        
    xpi3 = []
    fpi3 = []
    i=0
    for evt in c.applyCuts(xpi,apply_cut):
        if 0.440 < evt < 0.460 :
            xpi3.append(evt)
            fpi3.append(c.applyCuts(fpi,apply_cut)[i])
        i+=1
    if len(fpi3) > 1:
        xpi_tot.append(min(xpi3))
        fpi_tot.append(min(fpi3))
        # fpiuncern.append(min(fpi3)*uncern[2])
        fpiuncern.append(uncern[2])
    elif len(fpi3) == 1:
        xpi_tot.append((xpi3[0]))
        fpi_tot.append((fpi3[0]))
        # fpiuncern.append(uncern[2])
        fpiuncern.append(min(fpi3)*uncern[2])
    else:
        print "Invalid at ", xpi3
        
    xpi4 = []
    fpi4 = []
    i=0
    for evt in c.applyCuts(xpi,apply_cut):
        if 0.540 < evt < 0.560 :
            xpi4.append(evt)
            fpi4.append(c.applyCuts(fpi,apply_cut)[i])
        i+=1
    if len(fpi4) > 1:
        xpi_tot.append(min(xpi4))
        fpi_tot.append(min(fpi4))
        # fpiuncern.append(min(fpi4)*uncern[3])
        fpiuncern.append(uncern[3])
    elif len(fpi4) == 1:
        xpi_tot.append((xpi4[0]))
        fpi_tot.append((fpi4[0]))
        # fpiuncern.append(min(fpi4)*uncern[3])
        fpiuncern.append(uncern[3])
    else:
        print "Invalid at ", xpi4
        
    xpi5 = []
    fpi5 = []
    i=0
    for evt in c.applyCuts(xpi,apply_cut):
        if 0.640 < evt < 0.660 :
            xpi5.append(evt)
            fpi5.append(c.applyCuts(fpi,apply_cut)[i])
        i+=1
    if len(fpi5) > 1:
        xpi_tot.append(min(xpi5))
        fpi_tot.append(min(fpi5))
        # fpiuncern.append(min(fpi5)*uncern[4])
        fpiuncern.append(uncern[4])
    elif len(fpi5) == 1:
        xpi_tot.append((xpi5[0]))
        fpi_tot.append((fpi5[0]))
        # fpiuncern.append(uncern[4])
        fpiuncern.append(min(fpi5)*uncern[4])
    else:
        print "Invalid at ", xpi5
        
    xpi6 = []
    fpi6 = []
    i=0
    for evt in c.applyCuts(xpi,apply_cut):
        if 0.740 < evt < 0.760 :
            xpi6.append(evt)
            fpi6.append(c.applyCuts(fpi,apply_cut)[i])
        i+=1
    if len(fpi6) > 1:
        xpi_tot.append(min(xpi6))
        fpi_tot.append(min(fpi6))
        # fpiuncern.append(min(fpi6)*uncern[5])
        fpiuncern.append(uncern[5])
    elif len(fpi6) == 1:
        xpi_tot.append((xpi6[0]))
        fpi_tot.append((fpi6[0]))
        # fpiuncern.append(min(fpi6)*uncern[5])
        fpiuncern.append(uncern[5])
    else:
        print "Invalid at ", xpi6
        
    xpi7 = []
    fpi7 = []
    i=0
    for evt in c.applyCuts(xpi,apply_cut):
        if 0.840 < evt < 0.860 :
            xpi7.append(evt)
            fpi7.append(c.applyCuts(fpi,apply_cut)[i])
        i+=1
    if len(fpi7) > 1:
        xpi_tot.append(min(xpi7))
        fpi_tot.append(min(fpi7))
        # fpiuncern.append(min(fpi7)*uncern[6])
        fpiuncern.append(uncern[6])
    elif len(fpi7) == 1:
        xpi_tot.append((xpi7[0]))
        fpi_tot.append((fpi7[0]))
        # fpiuncern.append(uncern[6])
        fpiuncern.append(uncern[6])
    else:
        print "Invalid at ", xpi7
        
    xpi8 = []
    fpi8 = []
    i=0
    for evt in c.applyCuts(xpi,apply_cut):
        if 0.940 < evt < 0.960 :
            xpi8.append(evt)
            fpi8.append(c.applyCuts(fpi,apply_cut)[i])
        i+=1
    if len(fpi8) > 1:
        xpi_tot.append(min(xpi8))
        fpi_tot.append(min(fpi8))
        # fpiuncern.append(min(fpi8)*uncern[7])
        fpiuncern.append(uncern[7])
    elif len(fpi8) == 1:
        xpi_tot.append((xpi8[0]))
        fpi_tot.append((fpi8[0]))
        # fpiuncern.append(min(fpi8)*uncern[7])
        fpiuncern.append(uncern[7])
    else:
        print "Invalid at ", xpi8

    print "\n\nfpiuncern",fpiuncern,"\n\n"
    print "\n\nfpi",fpi_tot,"\n\n"
        
    x = np.sort(TDIS_xbj)
    # x = np.sort(xpi)
    
    # fit_f2pi = piSF_BLFQ(x)
    fit_f2pi = piSF_GRV(x)

    if(Q2_inp==10):
        BLFQ_Q2 = [[0.0009,0.002,0.004,0.0085,0.018],[0.52,0.45,0.40,0.32,0.25]]
        BLFQ_Q2_val=11
        GRV_Q2 = [[0.01,0.1,0.30],[0.25,0.20,0.18]]
        GRV_Q2_val=7
    elif(Q2_inp==30):
        BLFQ_Q2 = [[0.002,0.004,0.008,0.02,0.045],[0.6,0.48,0.35,0.30,0.25]]
        BLFQ_Q2_val=24
        GRV_Q2 = [[0.01,0.1,0.30],[0.40,0.20,0.16]]
        GRV_Q2_val=30
    elif(Q2_inp==50):
        BLFQ_Q2 = [0.008,0.019,0.038,0.075],[0.40,0.38,0.26,0.22]
        BLFQ_Q2_val=55
        GRV_Q2 = [0.01,0.1,0.30],[0.40,0.20,0.15]
        GRV_Q2_val=60
    elif(Q2_inp==70):
        BLFQ_Q2 = [0.023,0.045,0.09],[0.38,0.29,0.32]
        BLFQ_Q2_val=82
        GRV_Q2 = [0.01,0.1,0.30],[0.40,0.20,0.15]
        GRV_Q2_val=60

    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    uncernPlots10 = ax.errorbar(xpi_tot,fpi_tot,yerr=fpiuncern,fmt='o',label='$Q^2$ = %s $GeV^2$' % Q2_val,markersize='10')
    BLFQ_Plots10 = ax.errorbar(BLFQ_Q2[0],BLFQ_Q2[1],fmt='o',label='DESY-HERA-H1 \n $Q^2$ = %s $GeV^2$' % BLFQ_Q2_val,markersize='10')
    GRV_Plots10 = ax.plot(GRV_Q2[0],GRV_Q2[1],label='GRV fit \n $Q^2$ = %s $GeV^2$' % GRV_Q2_val,markersize='10')
    # fitPlot = ax.plot(x,fit_f2pi)
    # ax.fill_between(xpi_tot, fit_f2pi-error, fit_f2pi+error)
    ax.tick_params(labelsize=15)
    plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(0.01,1.0)
    plt.ylim(1e-6,0.5)
    plt.xlabel('$x_{\pi}$', fontsize =20)
    plt.ylabel('$F^{2}_{\pi}$', fontsize =20)
    plt.title('$F^{2}_{\pi}$ vs $x_{\pi}$', fontsize =20)
    leg = plt.legend(loc='center left', bbox_to_anchor=(0.8, 0.5), fontsize=10)
    leg.get_frame().set_alpha(1.)
    
    ####### Linear x, log y

    if(Q2_inp==10):
        GRV_Q2 = [[0.01,0.1,0.30],[0.25,0.20,0.18]]
        GRV_Q2_val=7
    elif(Q2_inp==30):
        GRV_Q2 = [[0.01,0.1,0.30],[0.40,0.20,0.16]]
        GRV_Q2_val=30
    elif(Q2_inp==50):
        GRV_Q2 = [[0.01,0.1,0.30],[0.40,0.20,0.15]]
        GRV_Q2_val=60
    elif(Q2_inp==70):
        GRV_Q2 = [[0.01,0.1,0.30],[0.40,0.20,0.15]]
        GRV_Q2_val=60
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    uncernPlots10 = ax.errorbar(xpi_tot,fpi_tot,yerr=fpiuncern,fmt='o',label='$Q^2$ = %s $GeV^2$' % Q2_val,markersize='10')
    GRV_Plots10 = ax.plot(GRV_Q2[0],GRV_Q2[1],label='GRV fit \n $Q^2$ = %s $GeV^2$' % GRV_Q2_val,markersize='10')
    ax.tick_params(labelsize=15)
    # fitPlot = ax.plot(x, fit_f2pi)
    # ax.fill_between(xpi_tot, fit_f2pi-error, fit_f2pi+error)
    # plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.1,1.0)
    # plt.ylim(1e-6,0.025)
    plt.ylim(1e-6,1.0)
    plt.xlabel('$x_{\pi}$', fontsize =20)
    plt.ylabel('$F^{2}_{\pi}$', fontsize =20)
    plt.title('$F^{2}_{\pi}$ vs $x_{\pi}$', fontsize =20)
    leg = plt.legend(loc='center left', bbox_to_anchor=(0.8, 0.5))
    leg.get_frame().set_alpha(1.)
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    f2PxbjPlot = ax.scatter(TDIS_xbj,f2N)
    plt.xlim(0.,1.)
    # plt.ylim(1e-6,0.5)
    ax.tick_params(labelsize=15)
    plt.xlabel('TDIS_xbj', fontsize =20)
    plt.ylabel('$F^{2}_{P}$', fontsize =20)
    plt.title('$F^{2}_{P}$ vs TDIS_xbj', fontsize =20)
    plt.close(f)
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    fpixbjPlot = ax.scatter(TDIS_xbj,fpi)
    plt.xlim(0.,1.)
    plt.ylim(1e-6,0.5)
    ax.tick_params(labelsize=15)
    plt.xlabel('TDIS_xbj', fontsize =20)
    plt.ylabel('$F^{2}_{\pi}$', fontsize =20)
    plt.title('$F^{2}_{\pi}$ vs TDIS_xbj', fontsize =20)
    plt.close(f)
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    fpixbjPlotLog = ax.scatter(TDIS_xbj,fpi)
    plt.xscale('log')
    plt.yscale('log')
    ax.tick_params(labelsize=15)
    plt.xlim(0.,1.)
    plt.ylim(1e-6,0.5)
    plt.xlabel('TDIS_xbj', fontsize =20)
    plt.ylabel('$F^{2}_{\pi}$', fontsize =20)
    plt.title('$F^{2}_{\pi}$ vs TDIS_xbj', fontsize =20)
    plt.close(f)
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    fpixpiPlot = ax.scatter(xpi,fpi)
    plt.xlim(0.,1)
    # plt.ylim(1e-6,0.5)
    ax.tick_params(labelsize=15)
    plt.xlabel('$x_{Bj}$', fontsize =20)
    plt.ylabel('$F^{2}_{\pi}$', fontsize =20)
    plt.title('$F^{2}_{\pi}$ vs $x_{Bj}$', fontsize =20)
    plt.close(f)
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    fpixpiPlotLog = ax.scatter(xpi,fpi)
    plt.xscale('log')
    # plt.yscale('log')
    ax.tick_params(labelsize=15)
    plt.xlim(0.,1.)
    # plt.ylim(1e-6,0.5)
    plt.xlabel('$x_{Bj}$', fontsize =20)
    plt.ylabel('$F^{2}_{\pi}$', fontsize =20)
    plt.title('$F^{2}_{\pi}$ vs $x_{Bj}$', fontsize =20)
    plt.close(f)    
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    xL_Plot = ax.hist(xL,bins=c.setbin(xL,200,0.,1.),label='all events',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.xlabel('$x_L$')
    plt.ylabel('$F^{2}_{\pi}$')
    plt.title('$F^{2}_{\pi}$ vs $x_L$', fontsize =20)
    plt.close(f)    
    
    phaseSpace = c.densityPlot(xpi, Q2, '$Q^2$ vs $x_{Bj}$','$x_{Bj}$','$Q^{2}$', 200, 200,  c, 0., 1.0, 0., 100.)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.,1.)
    plt.close(f)
    
    xpiPhase = xpi
    
    i=0
    tot_evts = []
    for evt in xpiPhase:
        if 0.1 < evt :
            tot_evts.append(evt)
    
    # Need to have numBin for each setting and the uncertainty will be the value uncern
    f = plt.figure(figsize=(11.69,8.27))
    plt.style.use('classic')    
    ax = f.add_subplot(331)    
    fpiscat1 = ax.errorbar(c.applyCuts(xpi,apply_cut),c.applyCuts(fpi,apply_cut),yerr=uncern[0],fmt='.',label='$x_{\pi}$=(0.20,0.30)')
    plt.subplots_adjust(hspace=0.3,wspace=0.45)
    ax.text(0.65, 0.95, '$x_{\pi}$=(0.20,0.30)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # plt.yscale('log')
    plt.xlim(0.20,0.30)
    # plt.ylim(1e-3,0.025)
    plt.xlabel('xpi')
    plt.ylabel('$F^{2}_{\pi}$')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    plt.title('$F^{2}_{\pi}$ vs xpi  [$Q^2$ = %s $GeV^2$]' % Q2_val, fontsize =20)
    
    ax = f.add_subplot(332)    
    fpiscat2 = ax.errorbar(c.applyCuts(xpi,apply_cut),c.applyCuts(fpi,apply_cut),yerr=uncern[1],fmt='.',label='$x_{\pi}$=(0.30,0.40)')
    plt.subplots_adjust(hspace=0.3,wspace=0.45)
    ax.text(0.65, 0.95, '$x_{\pi}$=(0.30,0.40)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # plt.yscale('log')
    plt.xlim(0.30,0.40)
    # plt.ylim(1e-3,0.025)
    plt.xlabel('xpi')
    plt.ylabel('$F^{2}_{\pi}$')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    
    ax = f.add_subplot(333)
    fpiscat3 = ax.errorbar(c.applyCuts(xpi,apply_cut),c.applyCuts(fpi,apply_cut),yerr=uncern[2],fmt='.',label='$x_{\pi}$=(0.40,0.50)')
    plt.subplots_adjust(hspace=0.3,wspace=0.45)
    ax.text(0.65, 0.95, '$x_{\pi}$=(0.40,0.50)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # plt.yscale('log')
    plt.xlim(0.40,0.50)
    # plt.ylim(1e-3,0.025)
    plt.xlabel('xpi')
    plt.ylabel('$F^{2}_{\pi}$')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    
    ax = f.add_subplot(334)
    fpiscat4 = ax.errorbar(c.applyCuts(xpi,apply_cut),c.applyCuts(fpi,apply_cut),yerr=uncern[3],fmt='.',label='$x_{\pi}$=(0.50,0.60)')
    plt.subplots_adjust(hspace=0.3,wspace=0.45)
    ax.text(0.65, 0.95, '$x_{\pi}$=(0.50,0.60)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # plt.yscale('log')
    plt.xlim(0.50,0.60)
    # plt.ylim(1e-3,0.025)
    plt.xlabel('xpi')
    plt.ylabel('$F^{2}_{\pi}$')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    
    ax = f.add_subplot(335)
    fpiscat5 = ax.errorbar(c.applyCuts(xpi,apply_cut),c.applyCuts(fpi,apply_cut),yerr=uncern[4],fmt='.',label='$x_{\pi}$=(0.60,0.70)')
    plt.subplots_adjust(hspace=0.3,wspace=0.45)
    ax.text(0.65, 0.95, '$x_{\pi}$=(0.60,0.70)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # plt.yscale('log')
    plt.xlim(0.60,0.70)
    # plt.ylim(1e-3,0.025)
    plt.xlabel('xpi')
    plt.ylabel('$F^{2}_{\pi}$')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    
    ax = f.add_subplot(336)
    fpiscat6 = ax.errorbar(c.applyCuts(xpi,apply_cut),c.applyCuts(fpi,apply_cut),yerr=uncern[5],fmt='.',label='$x_{\pi}$=(0.70,0.80)')
    plt.subplots_adjust(hspace=0.3,wspace=0.45)
    ax.text(0.65, 0.95, '$x_{\pi}$=(0.70,0.80)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # plt.yscale('log')
    plt.xlim(0.70,0.80)
    # plt.ylim(1e-3,0.025)
    plt.xlabel('xpi')
    plt.ylabel('$F^{2}_{\pi}$')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    
    ax = f.add_subplot(337)
    fpiscat7 = ax.errorbar(c.applyCuts(xpi,apply_cut),c.applyCuts(fpi,apply_cut),yerr=uncern[6],fmt='.',label='$x_{\pi}$=(0.80,0.90)')
    plt.subplots_adjust(hspace=0.3,wspace=0.45)
    ax.text(0.65, 0.95, '$x_{\pi}$=(0.80,0.90)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # plt.yscale('log')
    plt.xlim(0.80,0.90)
    # plt.ylim(1e-3,0.025)
    plt.xlabel('xpi')
    plt.ylabel('$F^{2}_{\pi}$')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    
    ax = f.add_subplot(338)
    fpiscat8 = ax.errorbar(c.applyCuts(xpi,apply_cut),c.applyCuts(fpi,apply_cut),yerr=uncern[7],fmt='.',label='$x_{\pi}$=(0.90,1.00)')
    plt.subplots_adjust(hspace=0.3,wspace=0.45)
    ax.text(0.65, 0.95, '$x_{\pi}$=(0.90,1.00)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # plt.yscale('log')
    plt.xlim(0.90,1.00)
    # plt.ylim(1e-3,0.025)
    plt.xlabel('xpi')
    plt.ylabel('$F^{2}_{\pi}$')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    
    # uncern plots

    uncern10 = []
    for i in range(0,8):
        uncern10.append(uncern[i])
        
    # uncern 100 plots
    
    uncern10 = []
    for i in range(0,8):
        uncern10.append(uncern100[i])
        
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    
    xpi1 = []
    fpi1 = []
    i=0
    for i in range(0,len(c.applyCuts(xpi,apply_cut))):
        if 0.20 < c.applyCuts(xpi,apply_cut)[i] < 0.30 :
            xpi1.append(c.applyCuts(xpi,apply_cut)[i])
            fpi1.append(c.applyCuts(fpi,apply_cut)[i])
        i+=1
    fpiPlot1 = ax.errorbar(xpi1,fpi1,yerr=uncern[0],fmt='.',label='$x_{\pi}$=(0.20,0.30)')
    plt.subplots_adjust(hspace=0.3,wspace=0.45)
    plt.yscale('log')
    # plt.xlim(0.20,0.30)
    # plt.ylim(1e-3,0.025)
    # plt.xlim(1e-3,1.)
    plt.xlabel('xpi')
    plt.ylabel('$F^{2}_{\pi}$')
    plt.title('$F^{2}_{\pi}$ vs xpi  [$Q^2$ = %s $GeV^2$]' % Q2_val, fontsize =20)
    
    xpi2 = []
    fpi2 = []
    i=0
    for i in range(0,len(c.applyCuts(xpi,apply_cut))):
        if 0.30 < c.applyCuts(xpi,apply_cut)[i] < 0.40 :
            xpi2.append(c.applyCuts(xpi,apply_cut)[i])
            fpi2.append(c.applyCuts(fpi,apply_cut)[i])
        i+=1
    fpiPlot2 = ax.errorbar(xpi2,fpi2,yerr=uncern[1],fmt='.',label='$x_{\pi}$=(0.30,0.40)')
    plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # plt.yscale('log')
    # plt.xlim(0.30,0.40)
    
    xpi3 = []
    fpi3 = []
    i=0
    for i in range(0,len(c.applyCuts(xpi,apply_cut))):
        if 0.40 < c.applyCuts(xpi,apply_cut)[i] < 0.50 :
            xpi3.append(c.applyCuts(xpi,apply_cut)[i])
            fpi3.append(c.applyCuts(fpi,apply_cut)[i])
        i+=1
    fpiPlot3 = ax.errorbar(xpi3,fpi3,yerr=uncern[2],fmt='.',label='$x_{\pi}$=(0.40,0.50)')
    plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # plt.yscale('log')
    # plt.xlim(0.40,0.50)
    
    xpi4 = []
    fpi4 = []
    i=0
    for i in range(0,len(c.applyCuts(xpi,apply_cut))):
        if 0.50 < c.applyCuts(xpi,apply_cut)[i] < 0.60 :
            xpi4.append(c.applyCuts(xpi,apply_cut)[i])
            fpi4.append(c.applyCuts(fpi,apply_cut)[i])
        i+=1
    fpiPlot4 = ax.errorbar(xpi4,fpi4,yerr=uncern[3],fmt='.',label='$x_{\pi}$=(0.50,0.60)')
    plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # plt.yscale('log')
    # plt.xlim(0.50,0.60)
    # plt.yscale('log')
    # plt.xlim(0.60,0.70)
    
    xpi5 = []
    fpi5 = []
    i=0
    for i in range(0,len(c.applyCuts(xpi,apply_cut))):
        if 0.60 < c.applyCuts(xpi,apply_cut)[i] < 0.70 :
            xpi5.append(c.applyCuts(xpi,apply_cut)[i])
            fpi5.append(c.applyCuts(fpi,apply_cut)[i])
        i+=1
    fpiPlot5 = ax.errorbar(xpi5,fpi5,yerr=uncern[4],fmt='.',label='$x_{\pi}$=(0.60,0.70)')
    plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # plt.yscale('log')
    # plt.xlim(0.60,0.70)
    
    xpi6 = []
    fpi6 = []
    i=0
    for i in range(0,len(c.applyCuts(xpi,apply_cut))):
        if 0.70 < c.applyCuts(xpi,apply_cut)[i] < 0.80 :
            xpi6.append(c.applyCuts(xpi,apply_cut)[i])
            fpi6.append(c.applyCuts(fpi,apply_cut)[i])
        i+=1
    fpiPlot6 = ax.errorbar(xpi6,fpi6,yerr=uncern[5],fmt='.',label='$x_{\pi}$=(0.70,0.80)')
    plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # plt.yscale('log')
    # plt.xlim(0.70,0.80)
    
    xpi7 = []
    fpi7 = []
    i=0
    for i in range(0,len(c.applyCuts(xpi,apply_cut))):
        if 0.80 < c.applyCuts(xpi,apply_cut)[i] < 0.90 :
            xpi7.append(c.applyCuts(xpi,apply_cut)[i])
            fpi7.append(c.applyCuts(fpi,apply_cut)[i])
        i+=1
    fpiPlot7 = ax.errorbar(xpi7,fpi7,yerr=uncern[6],fmt='.',label='$x_{\pi}$=(0.80,0.90)')
    plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # plt.yscale('log')
    # plt.xlim(0.80,0.90)
    
    xpi8 = []
    fpi8 = []
    i=0
    for i in range(0,len(c.applyCuts(xpi,apply_cut))):
        if 0.90 < c.applyCuts(xpi,apply_cut)[i] < 1.00 :
            xpi8.append(c.applyCuts(xpi,apply_cut)[i])
            fpi8.append(c.applyCuts(fpi,apply_cut)[i])
        i+=1
    fpiPlot8 = ax.errorbar(xpi8,fpi8,yerr=uncern[7],fmt='.',label='$x_{\pi}$=(0.90,1.00)')
    plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # plt.yscale('log')
    # plt.xlim(0.90,1.00)
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)
    
    # for f in xrange(1, plt.figure().number):
    #     pdf.savefig(f)
    # pdf.close()
    
def sigmavX(Q2_inp):

    if(Q2_inp==10):
        apply_cut=cut10
        Q2_val=10
        HERA_Q2_val=10
        xbj_HERA = np.array([1.61e-4,2.53e-4,4.00e-4,6.32e-4,1.02e-3,1.61e-3,2.53e-3,5.00e-3,2.10e-2])
        sigma_HERA = np.array([1.230,1.237,1.099,1.016,0.934,0.851,0.761,0.624,0.503])
        params,params_cov = optimize.curve_fit(curve_fit, xbj_HERA, sigma_HERA)
    elif(Q2_inp==30):
        apply_cut=cut30
        Q2_val=30
        HERA_Q2_val=27
        xbj_HERA = np.array([5.00e-4,8.00e-4,1.30e-3,2.10e-3,3.20e-3,5.00e-3,8.00e-3,1.30e-2])
        sigma_HERA = np.array([1.497,1.294,1.181,1.095,0.901,0.804,0.731,0.634])
        params,params_cov = optimize.curve_fit(curve_fit, xbj_HERA, sigma_HERA)
    elif(Q2_inp==50):
        apply_cut=cut50
        Q2_val=50
        HERA_Q2_val=45
        xbj_HERA = np.array([8.00e-4,1.30e-3,2.10e-3,3.20e-3,5.00e-3,8.00e-3,1.30e-2,2.10e-2])
        sigma_HERA = np.array([1.482,1.297,1.131,0.998,0.909,0.777,0.661,0.587])
        params,params_cov = optimize.curve_fit(curve_fit, xbj_HERA, sigma_HERA)   
    elif(Q2_inp==70):
        apply_cut=cut70
        Q2_val=70
        HERA_Q2_val=70
        xbj_HERA = np.array([1.30e-3,2.10e-3,3.20e-3,5.00e-3,8.00e-3,1.30e-2,2.10e-2,3.20e-2,5.00e-2])
        sigma_HERA = np.array([1.438,1.238,1.087,0.984,0.815,0.708,0.620,0.557,0.500])
        params,params_cov = optimize.curve_fit(curve_fit, xbj_HERA, sigma_HERA)

    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    plt.style.use('classic')    
        
    xBjscat = ax.scatter(c.applyCuts(TDIS_xbj,apply_cut),c.applyCuts(tot_sigma,apply_cut),label='$Q^2$=%s $GeV^2$' % Q2_val)
    HERAscat = ax.scatter(xbj_HERA,sigma_HERA,label='HERA $Q^2$=%s $GeV^2$' % HERA_Q2_val,color='r')
    fit = ax.plot(xbj_HERA,curve_fit(xbj_HERA, params[0], params[1]),color='r')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-5,1.0)
    plt.ylim(0.,1.5)
    # ax.text(0.45, 0.15, '$Q^2$=10 $GeV^2$', transform=ax.transAxes, fontsize=20, verticalalignment='top', horizontalalignment='right')
    plt.legend(loc='center left', bbox_to_anchor=(0.05, 0.30), fontsize=20)
    ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))

    plt.title('reduced d$\sigma$ vs $x_{Bj}$', fontsize =20)
    plt.ylabel('reduced d$\sigma$ ($10^-5$*$nb/GeV^{2}$)',fontsize=20)
    plt.xlabel('$x_{Bj}$',fontsize=20)

    plt.style.use('default')
        
    # for f in xrange(1, plt.figure().number):
    #     pdf.savefig(f)
    # pdf.close()
    
def main() :
    
    # sigmavxpi_Plot()
    for i in [10,30,50,70]:
        pionPlots(i)
        sigmavX(i)
    plt.show()
    
if __name__=='__main__': main()
