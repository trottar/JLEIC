#! /usr/bin/python

#
# Description:This will read in the array data file that contains all the leave histogram information
# ================================================================
# Time-stamp: 2019-04-08 08:07:33 trottar
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

from __future__ import division
import logging

# Gets rid of matplot logging DEBUG messages
plt_logger = logging.getLogger('matplotlib')
plt_logger.setLevel(logging.WARNING)

# Suppresses unwanted numpy warning
import warnings
import numpy as np
warnings.simplefilter(action='ignore', category=FutureWarning)

import seaborn as sns
import scipy as sci
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
from matplotlib import interactive
from matplotlib import colors
from matplotlib.ticker import MaxNLocator
from collections import namedtuple
from sys import path
import time,math,sys
# np.set_printoptions(threshold=sys.maxsize)

# My class function
sys.path.insert(0,'../../../Analysis/ROOTfiles/python')
from test_rootPy import pyPlot, pyCut

rootName = "TDISpion"
rootName_low = "TDISpion"
# rootName_low = "TDISpion_10-100"
# rootName_low = "TDISpion_10-100_OLD"
# rootName_low = "TDISpion_xpi"
rootName_high = "TDISpion_100-800"
# rootName_low = "TDISpion_med"
# rootName_high = "TDISpion_high"

tree1 = "T"

T1_arrkey =  "leafName"
T1_arrhist = "histData"

# pdf = matplotlib.backends.backend_pdf.PdfPages("%s.pdf" % rootName)
pdf = matplotlib.backends.backend_pdf.PdfPages("fpiPlot.pdf")
# pdf = matplotlib.backends.backend_pdf.PdfPages("fpiPlot_xbj.pdf")
# pdf = matplotlib.backends.backend_pdf.PdfPages("fpiPlot_xpi.pdf")
# pdf = matplotlib.backends.backend_pdf.PdfPages("Q2vxpi-xbj.pdf")
# pdf = matplotlib.backends.backend_pdf.PdfPages("sigmaPlots.pdf")

# Arguments for class function
p = pyPlot(rootName,tree1,T1_arrkey,T1_arrhist)
p1 = pyPlot(rootName_low,tree1,T1_arrkey,T1_arrhist)
p2 = pyPlot(rootName_high,tree1,T1_arrkey,T1_arrhist)

t_low = p1.lookup("tPrime")[0]
t_high = p2.lookup("tPrime")[0]

Q2_low = p1.lookup("Q2")[0]
Q2_high = p2.lookup("Q2")[0]

TDIS_xbj_low = p1.lookup("TDIS_xbj")[0]
TDIS_xbj_high = p2.lookup("TDIS_xbj")[0]

xpi_low = p1.lookup("xpi")[0]
# Look for xpi_low below xL if commented out
xpi_high = p2.lookup("xpi")[0]

# fpi_low = p1.lookup("fpi")[0]
# fpi_low = 0.361*p1.lookup("fpi")[0]
fpi_low = 0.361*p1.lookup("f2N")[0]
fpi_high = p2.lookup("fpi")[0]

f2N_low = p1.lookup("f2N")[0]
f2N_high = p2.lookup("f2N")[0]

# sigma_dis_low = p1.lookup("sigma_dis")[0]*(1e-5)
# sigma_dis_high = p2.lookup("sigma_dis")[0]*(1e-5)

sigma_dis_low = p1.lookup("sigma_dis")[0]
sigma_dis_high = p2.lookup("sigma_dis")[0]

pprz_inc = p1.lookup("pprz_inc")[0]
ppiz_Lab = p1.lookup("ppiz_Lab")[0]

pNz_Lab = pprz_inc-ppiz_Lab # Long neutron momentum
xL = pNz_Lab/pprz_inc # Frac of proton momentum
# sigma_tdis_low = p1.lookup("sigma_tdis")[0]
# sigma_tdis_high = p2.lookup("sigma_tdis")[0]

# xpi_low = TDIS_xbj_low/(1-xL)

# sigma_tdis_low = sigma_dis_low*(0.003/f2N_low)
# sigma_tdis_high = sigma_dis_high*(0.002/f2N_high)

sigma_tdis_low = sigma_dis_low*(fpi_low/f2N_low)
sigma_tdis_high = sigma_dis_high*(fpi_high/f2N_high)

y_low = p1.lookup("y")[0]
y_high = p2.lookup("y")[0]

yplus_low = 1+((1-y_low)*(1-y_low))
yplus_high = 1+((1-y_high)*(1-y_high))

# tot_sigma_low = (sigma_tdis_low)*((TDIS_xbj_low*(Q2_low*Q2_low)*(137)*(137))/(2*math.pi*yplus_low))
# tot_sigma_high = (sigma_tdis_high)*((TDIS_xbj_high*(Q2_high*Q2_high)*(137)*(137))/(2*math.pi*yplus_high))

# tot_sigma_low = (sigma_tdis_low)*(1e2)
# tot_sigma_high = (sigma_tdis_high)*(1e2)
tot_sigma_low = (sigma_tdis_low)
tot_sigma_high = (sigma_tdis_high)

# Define phyisics data
s_e = p.lookup("s_e")[0]
s_q = p.lookup("s_q")[0]
xBj = p.lookup("xBj")[0]
TDIS_xbj = p.lookup("TDIS_xbj")[0]
sigma_dis = p.lookup("sigma_dis")[0]*(1e-5)
TDIS_y = p.lookup("TDIS_y")[0]
ppix_Lab = p.lookup("ppix_Lab")[0] # Long pion momentum
ppiy_Lab = p.lookup("ppiy_Lab")[0]
ppiz_Lab = p.lookup("ppiz_Lab")[0]
EpiE_Lab = p.lookup("EpiE_Lab")[0]
pprx_inc = p.lookup("pprx_inc")[0] # Long proton momentum
ppry_inc = p.lookup("ppry_inc")[0]
pprz_inc = p.lookup("pprz_inc")[0]
EprE_inc = p.lookup("EprE_inc")[0]
y = p.lookup("y")[0]
Q2 = p.lookup("Q2")[0]
fpi = p.lookup("fpi")[0]
f2N = p.lookup("f2N")[0]
xpi = p.lookup("xpi")[0]
ypi = p.lookup("ypi")[0]
tpi = p.lookup("tpi")[0]
t = -p.lookup("tPrime")[0]
y_D = p.lookup("y_D")[0]
escat = p.lookup("EScatRest")[0]
nu = p.lookup("nu")[0]
TwoPdotk = p.lookup("TwoPdotk")[0]

# Q2_new = s_q/(xBj*y)
Q2_new = xBj/TwoPdotk
piQ2 = s_q/(xpi*ypi)
pNz_Lab = pprz_inc-ppiz_Lab # Long neutron momentum
ppiz_frac = ppiz_Lab/pprz_inc # Frac of proton momentum
pNz_frac = pNz_Lab/pprz_inc # Frac of proton momentum
yplus = 1+((1-y)*(1-y))
alt_sigma_dis = (1e1)*sigma_dis
tot_sigma = (sigma_dis)*((TDIS_xbj*(Q2*Q2)*(137)*(137))/(2*math.pi*yplus))
red_sig  = np.array([1.397,1.258,1.137,1.055,0.944,0.836,0.706,0.519])
red_y = np.array([0.657,0.416,0.263,0.163,0.103,0.066,0.033,0.008])
red_x = np.array([2.53e-4,4.0e-4,6.32e-4,1.02e-3,1.61e-3,2.52e-3,5.0e-3,2.10e-2])
red_Q2 = np.array([15,15,15,15,15,15,15,15])
red_yplus = 1+(1-red_y)*(1-red_y)
my_sigma = (red_sig)*((2*math.pi*red_yplus)/(red_x*(red_Q2*red_Q2)*(137)*(137)))

# Q2array_low =[10.,20.,30.,40.,50.,100.]
Q2array_low =[10.,30.,50.,70.]
# Q2array_low =[10.,16.,25.,40.,63,.,100.]

lowDict = {}

i=0
for x in range(0,len(Q2array_low)) :
    if i < (len(Q2array_low)-1):
        Q2tmp = '{"Q2cut%i" : ((%0.1f < Q2_low) & (Q2_low < %0.1f))}' % (i,Q2array_low[i]-2.5,Q2array_low[i]+2.5)
        # ttmp = '{"tcut" : (t_low < 0.05)}'
        ytmp = '{"ycut" : (0.05 < y_low)}'
        # ttmp = '{"tcut" : ((0.2 < t_low) & (t_low < 0.3))}'
        ttmp = '{"tcut" : (t_low < 1.0)}'
        # Q2fun = '{"Q2cut" : ((30 < Q2_low) & (Q2_low < 100))}'
        # xfun = '{"xbjcut" : ((0.03 < TDIS_xbj_low) & (TDIS_xbj_low < 0.1))}'
        print '{"Q2cut%i" : ((%0.1f < Q2_low) & (Q2_low < %0.1f))}' % (i,Q2array_low[i]-2.5,Q2array_low[i]+2.5)
    else:
        Q2tmp = '{"Q2cut%i" : ((%0.1f < Q2_low) & (Q2_low < %0.1f))}' % (i,Q2array_low[i]-2.5,Q2array_low[i]+2.5)
        # ttmp = '{"tcut" : (t_low < 0.05)}'
        ytmp = '{"ycut" : (0.05 < y_low)}'
        # ttmp = '{"tcut" : ((0.2 < t_low) & (t_low < 0.3))}'
        ttmp = '{"tcut" : (t_low < 1.0)}'
        # Q2fun = '{"Q2cut" : ((30 < Q2_low) & (Q2_low < 100))}'
        # xfun = '{"xbjcut" : ((0.03 < TDIS_xbj_low) & (TDIS_xbj_low < 0.1))}'
        print '{"Q2cut%i" : ((%0.1f < Q2_low) & (Q2_low < %0.1f))}' % (i,Q2array_low[i]-2.5,Q2array_low[i]+2.5)
    lowDict.update(eval(Q2tmp))
    # lowDict.update(eval(Q2fun))
    # lowDict.update(eval(xfun))
    lowDict.update(eval(ttmp))
    lowDict.update(eval(ytmp))
    i+=1

Q2array_high =[150.,200.,250.,300.,350.,400.,450.,500.,550.,600.,650.,700.,750.,800.]
    
highDict = {}
    
j=0
for x in range(0,len(Q2array_high)) :
    if j < (len(Q2array_high)-1):
        Q2tmp = '{"Q2cut%i" : ((%0.1f < Q2_high) & (Q2_high < %0.1f))}' % (i,Q2array_high[j]-2.5,Q2array_high[j]+2.5)
        ttmp = '{"tcut" : (t_high < 0.4)}'
        print '{"Q2cut%i" : ((%0.1f < Q2_high) & (Q2_high < %0.1f))}' % (i,Q2array_high[j]-2.5,Q2array_high[j]+2.5)
    else:
        Q2tmp = '{"Q2cut%i" : ((%0.1f < Q2_high) & (Q2_high < %0.1f))}' % (i,Q2array_high[j]-2.5,Q2array_high[j]+2.5)
        ttmp = '{"tcut" : (t_high < 0.4)}'
        print '{"Q2cut%i" : ((%0.1f < Q2_high) & (Q2_high < %0.1f))}' % (i,Q2array_high[j]-2.5,Q2array_high[j]+2.5)
    highDict.update(eval(Q2tmp))
    highDict.update(eval(ttmp))
    j+=1
    i+=1

sigmaPlotDict = {

    "Q2bin10" : ((10 < Q2_low) & (Q2_low < 16)),
    "Q2bin16" : ((16 < Q2_low) & (Q2_low < 25)),
    "Q2bin25" : ((25 < Q2_low) & (Q2_low < 40)),
    "Q2bin40" : ((40 < Q2_low) & (Q2_low < 63)),
    "Q2bin63" : ((63 < Q2_low) & (Q2_low < 100)),
    
}

c = pyCut(sigmaPlotDict)
clow = pyCut(lowDict)
chigh = pyCut(highDict)



def q2_Cut():

    cutQ2bin10 = ["Q2bin10"]
    cutQ2bin16 = ["Q2bin16"]
    cutQ2bin25 = ["Q2bin25"]
    cutQ2bin40 = ["Q2bin40"]
    cutQ2bin63 = ["Q2bin63"]
    ycut1 = ["ycut"]
    tcut1 = ["tcut"]
    # cut10 = ["Q2cut0","tcut"]
    # cut20 = ["Q2cut1","tcut"]
    # cut30 = ["Q2cut2","tcut"]
    # cut40 = ["Q2cut3","tcut"]
    # cut50 = ["Q2cut4","tcut"]
    # cut100 = ["Q2cut5","tcut"]
    
    cut10 = ["Q2cut0","tcut"]
    # cut10 = ["Q2cut","xbjcut"]
    cut30 = ["Q2cut1","tcut"]
    cut50 = ["Q2cut2","tcut"]
    cut70 = ["Q2cut3","tcut"]
    cut150 = ["Q2cut4","tcut"]
    cut200 = ["Q2cut5","tcut"]
    cut250 = ["Q2cut6","tcut"]
    cut300 = ["Q2cut7","tcut"]
    cut350 = ["Q2cut8","tcut"]
    cut400 = ["Q2cut9","tcut"]
    cut450 = ["Q2cut10","tcut"]
    cut500 = ["Q2cut11","tcut"]
    cut550 = ["Q2cut12","tcut"]
    cut600 = ["Q2cut13","tcut"]
    cut650 = ["Q2cut14","tcut"]
    cut700 = ["Q2cut15","tcut"]
    cut750 = ["Q2cut16","tcut"]
    cut800 = ["Q2cut17","tcut"]
        
    return[cutQ2bin10,cutQ2bin16,cutQ2bin25,cutQ2bin40,cutQ2bin63,ycut1,tcut1,cut10, cut30, cut50, cut70, cut150, cut200, cut250, cut300, cut350, cut400, cut450, cut500, cut550, cut600, cut650, cut700, cut750, cut800]
    
def densityPlot(x,y,title,xlabel,ylabel,binx,biny,xmin=None,xmax=None,ymin=None,ymax=None,cuts=None,figure=None,ax=None,layered=True):

    if cuts:
        xcut  = c.applyCuts(x,cuts)[0]
        ycut = c.applyCuts(y,cuts)[0]
    else:
        xcut = x
        ycut = y
    # if sub or figure:
    #     ax = figure.add_subplot(sub)
    if ax or figure:
        # ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
        print ""
    else:
        fig, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    if (xmin or xmax or ymin or ymax):
        # norm=colors.LogNorm() makes colorbar normed and logarithmic
        hist = ax.hist2d(xcut, ycut,bins=(p.setbin(x,binx,xmin,xmax)[0],p.setbin(y,biny,ymin,ymax)[0]))
    else:
        # norm=colors.LogNorm() makes colorbar normed and logarithmic
        hist = ax.hist2d(xcut, ycut,bins=(p.setbin(x,binx)[0],p.setbin(y,biny)[0]))
    if layered is True :
        # ax.set_title(title)
        # g = sns.jointplot(x=x, y=y, kind="hex", color="b")
        # g.ax_marg_x.hist(x, color="b", alpha=.6,bins=p.setbin(x,binx,xmin,xmax)[0])
        # g.ax_marg_y.hist(y, color="r", alpha=.6,orientation="horizontal",bins=p.setbin(y,biny,ymin,ymax)[0])
        # g.set_axis_labels(xlabel, ylabel)
        plt.colorbar(hist[3], ax=ax, spacing='proportional', label='Number of Events')
    # else :
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
    inputVal = [x,y]

    if (xmin or xmax or ymin or ymax):
        binVal = [p.setbin(x,binx,xmin,xmax)[0],p.setbin(y,biny,ymin,ymax)[0]]
    else:
        binVal = [p.setbin(x,binx)[0],p.setbin(y,biny)[0]]
    return binVal
    
def sigmavxpi_Plot():
    
    [cutQ2bin10,cutQ2bin16,cutQ2bin25,cutQ2bin40,cutQ2bin63,ycut1,tcut1,cut10, cut30, cut50, cut70, cut150, cut200, cut250, cut300, cut350, cut400, cut450, cut500, cut550, cut600, cut650, cut700, cut750, cut800] = q2_Cut()
    
    f = plt.figure(figsize=(11.69,8.27))
    plt.style.use('classic')
    
    ax = f.add_subplot(451)
    xpiscat1 = ax.errorbar(clow.applyCuts(xpi_low,cut10),clow.applyCuts(tot_sigma_low,cut10),xerr=np.sqrt(clow.applyCuts(xpi_low,cut10))/200,yerr=np.sqrt(clow.applyCuts(tot_sigma_low,cut10))/200,fmt='o',label='$Q^2$=10 $GeV^2$')
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
    xpiscat2 = ax.errorbar(clow.applyCuts(xpi_low,cut20),clow.applyCuts(tot_sigma_low,cut20),xerr=np.sqrt(clow.applyCuts(xpi_low,cut20))/200,yerr=np.sqrt(clow.applyCuts(tot_sigma_low,cut20))/200,fmt='o',label='$Q^2$=20 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=20 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    ax = f.add_subplot(453)
    xpiscat3 = ax.errorbar(clow.applyCuts(xpi_low,cut30),clow.applyCuts(tot_sigma_low,cut30),xerr=np.sqrt(clow.applyCuts(xpi_low,cut30))/200,yerr=np.sqrt(clow.applyCuts(tot_sigma_low,cut30))/200,fmt='o',label='$Q^2$=30 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    plt.title('d$\sigma_{TDIS}$ vs $x_\pi$', fontsize =20)
    
    ax = f.add_subplot(454)
    xpiscat4 = ax.errorbar(clow.applyCuts(xpi_low,cut40),clow.applyCuts(tot_sigma_low,cut40),xerr=np.sqrt(clow.applyCuts(xpi_low,cut40))/200,yerr=np.sqrt(clow.applyCuts(tot_sigma_low,cut40))/200,fmt='o',label='$Q^2$=40 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=40 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    ax = f.add_subplot(455)
    xpiscat5 = ax.errorbar(clow.applyCuts(xpi_low,cut50),clow.applyCuts(tot_sigma_low,cut50),xerr=np.sqrt(clow.applyCuts(xpi_low,cut50))/200,yerr=np.sqrt(clow.applyCuts(tot_sigma_low,cut50))/200,fmt='o',label='$Q^2$=50 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=50 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(456)
    xpiscat6 = ax.errorbar(clow.applyCuts(xpi_low,cut100),clow.applyCuts(tot_sigma_low,cut100),xerr=np.sqrt(clow.applyCuts(xpi_low,cut100))/200,yerr=np.sqrt(clow.applyCuts(tot_sigma_low,cut100))/200,fmt='o',label='$Q^2$=100 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=100 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))

    ax = f.add_subplot(457)
    xpiscat7 = ax.errorbar(chigh.applyCuts(xpi_high,cut150),chigh.applyCuts(tot_sigma_high,cut150),xerr=np.sqrt(chigh.applyCuts(xpi_high,cut150))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut150))/200,fmt='o',label='$Q^2$=150 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=150 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(458)
    xpiscat8 = ax.errorbar(chigh.applyCuts(xpi_high,cut200),chigh.applyCuts(tot_sigma_high,cut200),xerr=np.sqrt(chigh.applyCuts(xpi_high,cut200))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut200))/200,fmt='o',label='$Q^2$=200 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=200 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(459)
    xpiscat9 = ax.errorbar(chigh.applyCuts(xpi_high,cut250),chigh.applyCuts(tot_sigma_high,cut250),xerr=np.sqrt(chigh.applyCuts(xpi_high,cut250))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut250))/200,fmt='o',label='$Q^2$=250 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=250 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(4,5,10)
    xpiscat10 = ax.errorbar(chigh.applyCuts(xpi_high,cut300),chigh.applyCuts(tot_sigma_high,cut300),xerr=np.sqrt(chigh.applyCuts(xpi_high,cut300))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut300))/200,fmt='o',label='$Q^2$=300 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=300 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(4,5,11)
    xpiscat11 = ax.errorbar(chigh.applyCuts(xpi_high,cut350),chigh.applyCuts(tot_sigma_high,cut350),xerr=np.sqrt(chigh.applyCuts(xpi_high,cut350))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut350))/200,fmt='o',label='$Q^2$=350 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,0.8)
    ax.text(0.65, 0.25, '$Q^2$=350 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))

    ax = f.add_subplot(4,5,12)
    xpiscat12 = ax.errorbar(chigh.applyCuts(xpi_high,cut400),chigh.applyCuts(tot_sigma_high,cut400),xerr=np.sqrt(chigh.applyCuts(xpi_high,cut400))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut400))/200,fmt='o',label='$Q^2$=400 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,0.8)
    ax.text(0.65, 0.25, '$Q^2$=400 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(4,5,13)
    xpiscat13 = ax.errorbar(chigh.applyCuts(xpi_high,cut450),chigh.applyCuts(tot_sigma_high,cut450),xerr=np.sqrt(chigh.applyCuts(xpi_high,cut450))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut450))/200,fmt='o',label='$Q^2$=450 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,0.8)
    ax.text(0.65, 0.25, '$Q^2$=450 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(4,5,14)
    xpiscat14 = ax.errorbar(chigh.applyCuts(xpi_high,cut500),chigh.applyCuts(tot_sigma_high,cut500),xerr=np.sqrt(chigh.applyCuts(xpi_high,cut500))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut500))/200,fmt='o',label='$Q^2$=500 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,0.8)
    ax.text(0.65, 0.25, '$Q^2$=500 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(4,5,15)
    xpiscat15 = ax.errorbar(chigh.applyCuts(xpi_high,cut550),chigh.applyCuts(tot_sigma_high,cut550),xerr=np.sqrt(chigh.applyCuts(xpi_high,cut550))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut550))/2000,fmt='o',label='$Q^2$=550 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,0.8)
    ax.text(0.65, 0.25, '$Q^2$=550 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(4,5,16)
    xpiscat16 = ax.errorbar(chigh.applyCuts(xpi_high,cut600),chigh.applyCuts(tot_sigma_high,cut600),xerr=np.sqrt(chigh.applyCuts(xpi_high,cut600))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut600))/200,fmt='o',label='$Q^2$=600 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,0.6)
    ax.text(0.65, 0.25, '$Q^2$=600 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))

    ax = f.add_subplot(4,5,17)
    xpiscat17 = ax.errorbar(chigh.applyCuts(xpi_high,cut650),chigh.applyCuts(tot_sigma_high,cut650),xerr=np.sqrt(chigh.applyCuts(xpi_high,cut650))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut650))/200,fmt='o',label='$Q^2$=650 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,0.6)
    ax.text(0.65, 0.25, '$Q^2$=650 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    
    ax = f.add_subplot(4,5,18)
    xpiscat18 = ax.errorbar(chigh.applyCuts(xpi_high,cut700),chigh.applyCuts(tot_sigma_high,cut700),xerr=np.sqrt(chigh.applyCuts(xpi_high,cut700))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut700))/200,fmt='o',label='$Q^2$=700 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,0.6)
    ax.text(0.65, 0.25, '$Q^2$=700 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    
    ax = f.add_subplot(4,5,19)
    xpiscat19 = ax.errorbar(chigh.applyCuts(xpi_high,cut750),chigh.applyCuts(tot_sigma_high,cut750),xerr=np.sqrt(chigh.applyCuts(xpi_high,cut750))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut750))/200,fmt='o',label='$Q^2$=750 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,0.6)
    ax.text(0.65, 0.25, '$Q^2$=750 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    
    ax = f.add_subplot(4,5,20)
    xpiscat20 = ax.errorbar(chigh.applyCuts(xpi_high,cut800),chigh.applyCuts(tot_sigma_high,cut800),xerr=np.sqrt(chigh.applyCuts(xpi_high,cut800))/800,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut800))/800,fmt='o',label='$Q^2$=800 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,0.6)
    ax.text(0.65, 0.25, '$Q^2$=800 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    
    plt.xlabel('$x_\pi$')
    # plt.tight_layout()
    
def phaseSpace():

    plt.style.use('default')
    
    # Q2 vs xBj
    # hxbj_Q2_low = densityPlot(TDIS_xbj_low, Q2_low, '$Q^{2}$[10-100 GeV] vs TDIS_xbj[0.01-1.00]','TDIS_xbj','$Q^{2}$', 10, 20, 0.01, 1.00, 10., 100.)
    hxbj_Q2_low = densityPlot(TDIS_xbj_low, Q2_low, '$Q^{2}$[10-100 GeV] vs TDIS_xbj[0.01-1.00]','TDIS_xbj','$Q^{2}$', 10, 20, 0.035, 0.10, 35., 100.)
    plt.xscale('log')
    # plt.yscale('log')
    hxbj_Q2_high = densityPlot(TDIS_xbj_high, Q2_high, '$Q^{2}$[100-800 GeV] vs TDIS_xbj[0.1-1.0]','TDIS_xbj','$Q^{2}$', 10, 200, 0.1, 1.0, 100., 1000.)
    plt.xscale('log')
    # plt.yscale('log')

    # Q2 vs xBj chained
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    hxbj_Q2_lowFull = densityPlot(TDIS_xbj_low, Q2_low, '','TDIS_xbj','$Q^{2}$', 10, 20, 0.01, 1.00, 10., 100., figure=f,ax=ax)
    # hxbj_Q2_high = densityPlot(TDIS_xbj_high, Q2_high, '$Q^{2}$ vs TDIS_xbj (all sets)','TDIS_xbj','$Q^{2}$', 10, 200, 0.1, 1.0, 100., 1000., figure=f,ax=ax,layered=False)
    hxbj_Q2_highFull = densityPlot(TDIS_xbj_high, Q2_high, '$Q^{2}$ vs TDIS_xbj (all sets)','TDIS_xbj','$Q^{2}$', 10, 200, 0.1, 1.0, 100., 1000., figure=f,ax=ax)
    plt.xscale('log')
    plt.xlim(0.,1.)
    # plt.yscale('log')
    plt.ylim(0.,800.)
    
    # Q2 vs xpi
    hxpi_Q2_low = densityPlot(xpi_low, Q2_low, '$Q^{2}$[10-100 GeV] vs $x_{pi}$[0.01-1.00]','$x_{pi}$','$Q^{2}$', 10, 20, 0.01, 1., 10., 100.)
    plt.xscale('log')
    # plt.yscale('log')
    hxpi_Q2_high = densityPlot(xpi_high, Q2_high, '$Q^{2}$[100-800 GeV] vs $x_{pi}$[0.1-1.0]','$x_{pi}$','$Q^{2}$', 10, 200, 0.1, 1., 100., 800.)
    plt.xscale('log')
    # plt.yscale('log')

    # Q2 vs xpi chained
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    hxpi_Q2_low = densityPlot(xpi_low, Q2_low, '','$x_{pi}$','$Q^{2}$', 10, 20, 0.01, 1., 10., 100., figure=f,ax=ax)
    # hxpi_Q2_high = densityPlot(xpi_high, Q2_high, '$Q^{2}$ vs $x_{pi}$ (all sets)','$x_{pi}$','$Q^{2}$', 10, 200, 0.1, 1., 100., 800., figure=f,ax=ax,layered=False)
    hxpi_Q2_high = densityPlot(xpi_high, Q2_high, '$Q^{2}$ vs $x_{pi}$ (all sets)','$x_{pi}$','$Q^{2}$', 10, 200, 0.1, 1., 100., 800., figure=f,ax=ax)
    plt.xscale('log')
    plt.xlim(0.,1.)
    # plt.yscale('log')
    plt.ylim(0.,800.)

    
    # for f in xrange(1, plt.figure().number):
        # pdf.savefig(f)
    # pdf.close()

    # print len(hxbj_Q2_low[2])

    return hxbj_Q2_low,hxbj_Q2_high

def lumi():

    hxbj_Q2_low = phaseSpace()[0]
    hxbj_Q2_high = phaseSpace()[1]

    binx =  np.append(hxbj_Q2_low[0],hxbj_Q2_high[0])
    binx = np.sort(binx)
    binQ2 =  np.append(hxbj_Q2_low[1],hxbj_Q2_high[1])
    binQ2 = np.sort(binQ2)
    sigTDIS = np.append(tot_sigma_low,tot_sigma_high)
    tot_xbj = np.append(TDIS_xbj_low,TDIS_xbj_high)
    tot_Q2 = np.append(Q2_low,Q2_high)
    
    x_evts = len(tot_xbj)
    q2_evts = len(tot_Q2)
    print len(sigTDIS)

    # xEvts = np.array([])
    xEvts = []
    xval = []
    Q2Evts = []
    Q2val = []

    print x_evts,"Bin x: ", len(binx)
    print q2_evts,"Bin Q2: ", len(binQ2)

    print "\n", binx
    print "\n", binQ2
    
    print "Binning x..."
    for i in range(0,len(binx)-1) :
        p.progressBar(i,len(binx)-1,50)
        j=0
        for evt in tot_xbj:
            if 0 <= evt <= binx[0] or binx[i] <= evt <= binx[i+1]:
            # if binx[i] <= evt <= binx[i+1] or binx[len(binx)-1] < evt:
                xval.append([i,evt])
                j+=1
        for k in range(0,len(xval)):
            if xval[k][0] == i:
                xEvts.append([xval[k][1],j])
                # print xval[k][0], i
                # xEvts.append(evt)
            # else:
                # print evt, " not in range(",binx[i],",",binx[i+1],")"
    print "\nBinning Q2..."
    for i in range(0,len(binQ2)-1) :
        p.progressBar(i,len(binQ2)-1,50)
        j=0
        for evt in tot_Q2:
            if 0 <= evt <= binQ2[0] or binQ2[i] <= evt <= binQ2[i+1]:
            # if binQ2[i] <= evt <= binQ2[i+1] or binQ2[len(binQ2)-1] < evt:
                Q2val.append([i,evt])
                j+=1
        for k in range(0,len(Q2val)):
            if Q2val[k][0] == i:
                Q2Evts.append([Q2val[k][1],j])
                # print Q2val[k][0], i
                # Q2Evts.append(evt)
            # else:
                # print evt, " not in range(",binQ2[i],",",binQ2[i+1],")"

    # print "\n\n",xEvts
    # print "\n",Q2Evts

    print "\n\n",len(xEvts)
    print "\n",len(Q2Evts), "\n"

    x_binVal = []
    for i in range(0,len(xEvts)-1):
        if xEvts[i][1] != xEvts[i+1][1]:
            x_binVal.append(xEvts[i][1])
    print x_binVal
    print np.sum(x_binVal)

    Q2_binVal = []
    for i in range(0,len(Q2Evts)-1):
        if Q2Evts[i][1] != Q2Evts[i+1][1]:
            Q2_binVal.append(Q2Evts[i][1])
    print Q2_binVal
    print np.sum(Q2_binVal)

    print len(xval),len(xEvts)
    print len(Q2val),len(Q2Evts)
    
    return xval, Q2val

    
def sigmavxBj_Plot():

    [cutQ2bin10,cutQ2bin16,cutQ2bin25,cutQ2bin40,cutQ2bin63,ycut1,tcut1,cut10, cut30, cut50, cut70, cut150, cut200, cut250, cut300, cut350, cut400, cut450, cut500, cut550, cut600, cut650, cut700, cut750, cut800] = q2_Cut()
    
    f = plt.figure(figsize=(11.69,8.27))
    plt.style.use('classic')   
    
    ax = f.add_subplot(221)
    # xBjscat1 = ax.errorbar(clow.applyCuts(TDIS_xbj_low,cut10),clow.applyCuts(tot_sigma_low,cut10),xerr=np.sqrt(clow.applyCuts(TDIS_xbj_low,cut10))/200,yerr=np.sqrt(clow.applyCuts(tot_sigma_low,cut10))/200,fmt='o',label='$Q^2$=10 $GeV^2$')
    xBjscat1 = ax.scatter(clow.applyCuts(TDIS_xbj_low,cut10),clow.applyCuts(tot_sigma_low,cut10),label='$Q^2$=10 $GeV^2$')
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

    ax = f.add_subplot(222)
    # xBjscat2 = ax.errorbar(clow.applyCuts(TDIS_xbj_low,cut30),clow.applyCuts(tot_sigma_low,cut30),xerr=np.sqrt(clow.applyCuts(TDIS_xbj_low,cut30))/200,yerr=np.sqrt(clow.applyCuts(tot_sigma_low,cut30))/200,fmt='o',label='$Q^2$=30 $GeV^2$')
    xBjscat2 = ax.scatter(clow.applyCuts(TDIS_xbj_low,cut30),clow.applyCuts(tot_sigma_low,cut30),label='$Q^2$=30 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    plt.xlim(1e-3,1.0)
    plt.ylim(0.,7.)
    ax.text(0.65, 0.25, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))    
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    
    ax = f.add_subplot(223)
    # xBjscat3 = ax.errorbar(clow.applyCuts(TDIS_xbj_low,cut50),clow.applyCuts(tot_sigma_low,cut50),xerr=np.sqrt(clow.applyCuts(TDIS_xbj_low,cut50))/200,yerr=np.sqrt(clow.applyCuts(tot_sigma_low,cut50))/200,fmt='o',label='$Q^2$=50 $GeV^2$')
    xBjscat3 = ax.scatter(clow.applyCuts(TDIS_xbj_low,cut50),clow.applyCuts(tot_sigma_low,cut50),label='$Q^2$=50 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    plt.xlim(1e-3,1.0)
    plt.ylim(0.,.07)
    ax.text(0.65, 0.25, '$Q^2$=50 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))    
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    
    ax = f.add_subplot(224)
    # xBjscat4 = ax.errorbar(clow.applyCuts(TDIS_xbj_low,cut70),clow.applyCuts(tot_sigma_low,cut70),xerr=np.sqrt(clow.applyCuts(TDIS_xbj_low,cut70))/200,yerr=np.sqrt(clow.applyCuts(tot_sigma_low,cut70))/200,fmt='o',label='$Q^2$=70 $GeV^2$')
    xBjscat4 = ax.scatter(clow.applyCuts(TDIS_xbj_low,cut70),clow.applyCuts(tot_sigma_low,cut70),label='$Q^2$=70 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    plt.xlim(1e-3,1.0)
    plt.ylim(0.,.07)
    ax.text(0.65, 0.25, '$Q^2$=70 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))    
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    
    
    # ax = f.add_subplot(335)
    # # xBjscat5 = ax.errorbar(clow.applyCuts(TDIS_xbj_low,cut50),clow.applyCuts(tot_sigma_low,cut50),xerr=np.sqrt(clow.applyCuts(TDIS_xbj_low,cut50))/200,yerr=np.sqrt(clow.applyCuts(tot_sigma_low,cut50))/200,fmt='o',label='$Q^2$=50 $GeV^2$')
    # xBjscat5 = ax.scatter(clow.applyCuts(TDIS_xbj_low,cut50),clow.applyCuts(tot_sigma_low,cut50),label='$Q^2$=50 $GeV^2$')
    # plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    # plt.xlim(1e-3,0.6)
    # plt.ylim(0.,10.2)
    # ax.text(0.65, 0.25, '$Q^2$=50 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax = f.add_subplot(456)
    # # xBjscat6 = ax.errorbar(clow.applyCuts(TDIS_xbj_low,cut100),clow.applyCuts(tot_sigma_low,cut100),xerr=np.sqrt(clow.applyCuts(TDIS_xbj_low,cut100))/200,yerr=np.sqrt(clow.applyCuts(tot_sigma_low,cut100))/200,fmt='o',label='$Q^2$=100 $GeV^2$')
    # xBjscat6 = ax.scatter(clow.applyCuts(TDIS_xbj_low,cut100),clow.applyCuts(tot_sigma_low,cut100),label='$Q^2$=100 $GeV^2$')
    # plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    # plt.xlim(1e-3,0.6)
    # plt.ylim(0.,1.0)
    # ax.text(0.65, 0.25, '$Q^2$=100 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    # ax = f.add_subplot(457)
    # # xBjscat7 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut150),chigh.applyCuts(tot_sigma_high,cut150),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut150))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut150))/200,fmt='o',label='$Q^2$=150 $GeV^2$')
    # xBjscat7 = ax.scatter(chigh.applyCuts(TDIS_xbj_high,cut150),chigh.applyCuts(tot_sigma_high,cut150),label='$Q^2$=150 $GeV^2$')
    # plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    # plt.xlim(1e-3,0.6)
    # plt.ylim(0.,1.0)
    # ax.text(0.65, 0.25, '$Q^2$=150 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax = f.add_subplot(458)
    # # xBjscat8 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut200),chigh.applyCuts(tot_sigma_high,cut200),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut200))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut200))/200,fmt='o',label='$Q^2$=200 $GeV^2$')
    # xBjscat8 = ax.scatter(chigh.applyCuts(TDIS_xbj_high,cut200),chigh.applyCuts(tot_sigma_high,cut200),label='$Q^2$=200 $GeV^2$')
    # plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    # plt.xlim(1e-3,0.6)
    # plt.ylim(0.,1.0)
    # ax.text(0.65, 0.25, '$Q^2$=200 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax = f.add_subplot(459)
    # # xBjscat9 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut250),chigh.applyCuts(tot_sigma_high,cut250),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut250))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut250))/200,fmt='o',label='$Q^2$=250 $GeV^2$')
    # xBjscat9 = ax.scatter(chigh.applyCuts(TDIS_xbj_high,cut250),chigh.applyCuts(tot_sigma_high,cut250),label='$Q^2$=250 $GeV^2$')
    # plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    # plt.xlim(1e-3,0.6)
    # plt.ylim(0.,1.0)
    # ax.text(0.65, 0.25, '$Q^2$=250 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax = f.add_subplot(4,5,10)
    # # xBjscat10 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut300),chigh.applyCuts(tot_sigma_high,cut300),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut300))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut300))/200,fmt='o',label='$Q^2$=300 $GeV^2$')
    # xBjscat10 = ax.scatter(chigh.applyCuts(TDIS_xbj_high,cut300),chigh.applyCuts(tot_sigma_high,cut300),label='$Q^2$=300 $GeV^2$')
    # plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    # plt.xlim(1e-3,0.6)
    # plt.ylim(0.,1.0)
    # ax.text(0.65, 0.25, '$Q^2$=300 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax = f.add_subplot(4,5,11)
    # # xBjscat11 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut350),chigh.applyCuts(tot_sigma_high,cut350),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut350))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut350))/200,fmt='o',label='$Q^2$=350 $GeV^2$')
    # xBjscat11 = ax.scatter(chigh.applyCuts(TDIS_xbj_high,cut350),chigh.applyCuts(tot_sigma_high,cut350),label='$Q^2$=350 $GeV^2$')
    # plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    # plt.xlim(1e-3,0.6)
    # plt.ylim(0.,0.8)
    # ax.text(0.65, 0.25, '$Q^2$=350 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    # ax = f.add_subplot(4,5,12)
    # # xBjscat12 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut400),chigh.applyCuts(tot_sigma_high,cut400),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut400))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut400))/200,fmt='o',label='$Q^2$=400 $GeV^2$')
    # xBjscat12 = ax.scatter(chigh.applyCuts(TDIS_xbj_high,cut400),chigh.applyCuts(tot_sigma_high,cut400),label='$Q^2$=400 $GeV^2$')
    # plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    # plt.xlim(1e-3,0.6)
    # plt.ylim(0.,0.8)
    # ax.text(0.65, 0.25, '$Q^2$=400 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax = f.add_subplot(4,5,13)
    # # xBjscat13 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut450),chigh.applyCuts(tot_sigma_high,cut450),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut450))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut450))/200,fmt='o',label='$Q^2$=450 $GeV^2$')
    # xBjscat13 = ax.scatter(chigh.applyCuts(TDIS_xbj_high,cut450),chigh.applyCuts(tot_sigma_high,cut450),label='$Q^2$=450 $GeV^2$')
    # plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    # plt.xlim(1e-3,0.6)
    # plt.ylim(0.,0.8)
    # ax.text(0.65, 0.25, '$Q^2$=450 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax = f.add_subplot(4,5,14)
    # # xBjscat14 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut500),chigh.applyCuts(tot_sigma_high,cut500),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut500))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut500))/200,fmt='o',label='$Q^2$=500 $GeV^2$')
    # xBjscat14 = ax.scatter(chigh.applyCuts(TDIS_xbj_high,cut500),chigh.applyCuts(tot_sigma_high,cut500),label='$Q^2$=500 $GeV^2$')
    # plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    # plt.xlim(1e-3,0.6)
    # plt.ylim(0.,0.8)
    # ax.text(0.65, 0.25, '$Q^2$=500 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax = f.add_subplot(4,5,15)
    # # xBjscat15 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut550),chigh.applyCuts(tot_sigma_high,cut550),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut550))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut550))/2000,fmt='o',label='$Q^2$=550 $GeV^2$')
    # xBjscat15 = ax.scatter(chigh.applyCuts(TDIS_xbj_high,cut550),chigh.applyCuts(tot_sigma_high,cut550),label='$Q^2$=550 $GeV^2$')
    # plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    # plt.xlim(1e-3,0.6)
    # plt.ylim(0.,0.8)
    # ax.text(0.65, 0.25, '$Q^2$=550 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax = f.add_subplot(4,5,16)
    # # xBjscat16 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut600),chigh.applyCuts(tot_sigma_high,cut600),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut600))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut600))/200,fmt='o',label='$Q^2$=600 $GeV^2$')
    # xBjscat16 = ax.scatter(chigh.applyCuts(TDIS_xbj_high,cut600),chigh.applyCuts(tot_sigma_high,cut600),label='$Q^2$=600 $GeV^2$')
    # plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    # plt.xlim(1e-3,0.6)
    # plt.ylim(0.,0.6)
    # ax.text(0.65, 0.25, '$Q^2$=600 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    # ax = f.add_subplot(4,5,17)
    # # xBjscat17 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut650),chigh.applyCuts(tot_sigma_high,cut650),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut650))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut650))/200,fmt='o',label='$Q^2$=650 $GeV^2$')
    # xBjscat17 = ax.scatter(chigh.applyCuts(TDIS_xbj_high,cut650),chigh.applyCuts(tot_sigma_high,cut650),label='$Q^2$=650 $GeV^2$')
    # plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    # plt.xlim(1e-3,0.6)
    # plt.ylim(0.,0.6)
    # ax.text(0.65, 0.25, '$Q^2$=650 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    # ax = f.add_subplot(4,5,18)
    # # xBjscat18 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut700),chigh.applyCuts(tot_sigma_high,cut700),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut700))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut700))/200,fmt='o',label='$Q^2$=700 $GeV^2$')
    # xBjscat18 = ax.scatter(chigh.applyCuts(TDIS_xbj_high,cut700),chigh.applyCuts(tot_sigma_high,cut700),label='$Q^2$=700 $GeV^2$')
    # plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    # plt.xlim(1e-3,0.6)
    # plt.ylim(0.,0.6)
    # ax.text(0.65, 0.25, '$Q^2$=700 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    # ax = f.add_subplot(4,5,19)
    # # xBjscat19 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut750),chigh.applyCuts(tot_sigma_high,cut750),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut750))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut750))/200,fmt='o',label='$Q^2$=750 $GeV^2$')
    # xBjscat19 = ax.scatter(chigh.applyCuts(TDIS_xbj_high,cut750),chigh.applyCuts(tot_sigma_high,cut750),label='$Q^2$=750 $GeV^2$')
    # plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    # plt.xlim(1e-3,0.6)
    # plt.ylim(0.,0.6)
    # ax.text(0.65, 0.25, '$Q^2$=750 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    # ax = f.add_subplot(4,5,20)
    # # xBjscat20 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut800),chigh.applyCuts(tot_sigma_high,cut800),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut800))/800,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut800))/800,fmt='o',label='$Q^2$=800 $GeV^2$')
    # xBjscat20 = ax.scatter(chigh.applyCuts(TDIS_xbj_high,cut800),chigh.applyCuts(tot_sigma_high,cut800),label='$Q^2$=800 $GeV^2$')
    # plt.subplots_adjust(hspace=0.0,wspace=0.0)
    # plt.xscale('log')
    # plt.xlim(1e-3,0.6)
    # plt.ylim(0.,0.6)
    # ax.text(0.65, 0.25, '$Q^2$=800 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    
    plt.xlabel('TDIS_xbj')
    # plt.tight_layout()
    # plt.setp(ax.get_xticklabels()[0], visible=False)
    # plt.setp(ax.get_xticklabels()[-1], visible=False)
        
    # Reset figure style, otherwise next plot will look weird
    plt.style.use('default')

    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));    
    Q2scat1 = ax.scatter(clow.applyCuts(Q2_low,tcut1),clow.applyCuts(tot_sigma_low,tcut1))
    # plt.xscale('log')
    # plt.xlim(1e-3,1.0)
    # plt.ylim(0.,15.)
    plt.title('$d\sigma_{TDIS}$ vs $Q^2$', fontsize =20)
    plt.xlabel('$Q^2$')
    plt.ylabel('$d\sigma_{TDIS}$')

    

    # xEvts = lumi()
    # Q2Evts = lumi()
    
    # Luminosity
    lumi = []
    uncern = []
    allbinx = []

    # xEvts,Q2Evts
    for j in range(0,4):

        tmp = 'xBjscat%i.get_offsets()' % (j+1)
        # print tmp
        binValue = eval(tmp)
        evts = len(binValue)
        
        xval = []
        sigval = []
    
        for i in range(0,evts):
            xval.append(binValue[i][0])
            sigval.append(binValue[i][1])
        
        binx = 0.1
        # binx = max(xval)-min(xval)
        binQ2 = 1.0
        # print "---------------------------------"
        # print "Events in bin: ",evts, "\nBin value array: ", binValue
        # print "\nbinx: ", binx, "\nbinQ2: ", binQ2
        
        avgSigma = np.average(sigval)

        # print "\nAverage Sigma: ", avgSigma

        # Nevts = 
        
        lumi = np.append(lumi,(evts/(avgSigma*(binQ2)*(binx))))
        # lumi = np.append(lumi,(Nevts/(avgSigma*(binQ2)*(binx))))

        # print lumi, avgSigma, evts

        avgLumi = np.average(lumi)
        
        # uncern = np.append(uncern,np.sqrt(avgLumi*avgSigma*(binQ2)*(binx))/evts)
        # uncern = np.append(uncern,np.sqrt(avgLumi*avgSigma*(binQ2)*(binx))/evts)
        
        # print "Plot: ", j+1
        # print "Luminosity: ", lumi[j]
        # print "Uncertainty: ", uncern[j]
        # print "---------------------------------\n"

        allbinx.append(binx)
        
    print "\nbinx: ", allbinx
    print "\nLuminosity: ", lumi

    
    # f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    # scat1 = ax.hist(lumi,bins=p1.setbin(lumi,20)[0],label='Low',histtype='step', alpha=0.5, stacked=True, fill=True)
    # # scat1 = ax.hist(t_low,bins=p1.setbin(t_low,200,0.,1.0)[0],label='Low',histtype='step', alpha=0.5, stacked=True, fill=True)
    # # plt.title('$f_\pi$ vs TDIS_xbj', fontsize =20)
    # # plt.xlabel('TDIS_xbj')
    # # plt.ylabel('$f_\pi$')
    # leg = plt.legend(bbox_to_anchor=(0.95,0.3), loc="center right")
    # leg.get_frame().set_alpha(1.)

    return uncern,lumi

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


def pionPlots():

    [cutQ2bin10,cutQ2bin16,cutQ2bin25,cutQ2bin40,cutQ2bin63,ycut1,tcut1,cut10, cut30, cut50, cut70, cut150, cut200, cut250, cut300, cut350, cut400, cut450, cut500, cut550, cut600, cut650, cut700, cut750, cut800] = q2_Cut()
    
    lumi = sigmavxBj_Plot()[1]

    tot_lumi = []

    for i in range(0,4):
        for j in range (0,8):
            tot_lumi.append(0.031180)

    lum10 = []
    lum100 = []

    f = plt.figure(figsize=(11.69,8.27))
    plt.style.use('classic')

    ax = f.add_subplot(331)
    numBin10_1 = ax.hist(clow.applyCuts(xpi_low,cut10),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.20,0.30)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.20,0.30)
    ax.text(0.65, 0.95, '$x_\pi$=(0.20,0.30)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts10_1 = []
    for val in numBin10_1[0]:
        if 0.20 < numBin10_1[1][i] < 0.30 :
            maxEvts10_1.append(numBin10_1[0][i])
        i+=1
    numEvts10_1 = np.sum(maxEvts10_1)
    print (1e-6)*lumi[0], numEvts10_1
    lum10.append((10/(0.031180))*numEvts10_1)
    lum100.append((100/(0.031180))*numEvts10_1)
    
    plt.title('$x_{\pi}$ plots [$Q^2$ = 10 $GeV^2$]', fontsize =20)
    
    ax = f.add_subplot(332)
    numBin10_2 = ax.hist(clow.applyCuts(xpi_low,cut10),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.30,0.40)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.30,0.40)
    ax.text(0.65, 0.95, '$x_\pi$=(0.30,0.40)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts10_2 = []
    for val in numBin10_2[0]:
        if 0.30 < numBin10_2[1][i] < 0.40 :
            maxEvts10_2.append(numBin10_2[0][i])
        i+=1
    numEvts10_2 = np.sum(maxEvts10_2)
    print (1e-6)*lumi[0], numEvts10_2
    lum10.append((10/(0.031180))*numEvts10_2)
    lum100.append((100/(0.031180))*numEvts10_2)
    
    ax = f.add_subplot(333)
    numBin10_3 = ax.hist(clow.applyCuts(xpi_low,cut10),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.40,0.50)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.40,0.50)
    ax.text(0.65, 0.95, '$x_\pi$=(0.40,0.50)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts10_3 = []
    for val in numBin10_3[0]:
        if 0.40 < numBin10_3[1][i] < 0.50 :
            maxEvts10_3.append(numBin10_3[0][i])
        i+=1
    numEvts10_3 = np.sum(maxEvts10_3)
    print (1e-6)*lumi[0], numEvts10_3
    lum10.append((10/(0.031180))*numEvts10_3)
    lum100.append((100/(0.031180))*numEvts10_3)
    
    ax = f.add_subplot(334)
    numBin10_4 = ax.hist(clow.applyCuts(xpi_low,cut10),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.50,0.60)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.50,0.60)
    ax.text(0.65, 0.95, '$x_\pi$=(0.50,0.60)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts10_4 = []
    for val in numBin10_4[0]:
        if 0.50 < numBin10_4[1][i] < 0.60 :
            maxEvts10_4.append(numBin10_4[0][i])
        i+=1
    numEvts10_4 = np.sum(maxEvts10_4)
    print (1e-6)*lumi[0], numEvts10_4
    lum10.append((10/(0.031180))*numEvts10_4)
    lum100.append((100/(0.031180))*numEvts10_4)
    
    ax = f.add_subplot(335)
    numBin10_5 = ax.hist(clow.applyCuts(xpi_low,cut10),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.60,0.70)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.60,0.70)
    ax.text(0.65, 0.95, '$x_\pi$=(0.60,0.70)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts10_5 = []
    for val in numBin10_5[0]:
        if 0.60 < numBin10_5[1][i] < 0.70 :
            maxEvts10_5.append(numBin10_5[0][i])
        i+=1
    numEvts10_5 = np.sum(maxEvts10_5)
    print (1e-6)*lumi[0], numEvts10_5
    lum10.append((10/(0.031180))*numEvts10_5)
    lum100.append((100/(0.031180))*numEvts10_5)
    
    ax = f.add_subplot(336)
    numBin10_6 = ax.hist(clow.applyCuts(xpi_low,cut10),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.70,0.80)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.70,0.80)
    ax.text(0.65, 0.95, '$x_\pi$=(0.70,0.80)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts10_6 = []
    for val in numBin10_6[0]:
        if 0.70 < numBin10_6[1][i] < 0.80 :
            maxEvts10_6.append(numBin10_6[0][i])
        i+=1
    numEvts10_6 = np.sum(maxEvts10_6)
    print (1e-6)*lumi[0], numEvts10_6
    lum10.append((10/(0.031180))*numEvts10_6)
    lum100.append((100/(0.031180))*numEvts10_6)
    
    ax = f.add_subplot(337)
    numBin10_7 = ax.hist(clow.applyCuts(xpi_low,cut10),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.80,0.90)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.80,0.90)
    ax.text(0.65, 0.95, '$x_\pi$=(0.80,0.90)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts10_7 = []
    for val in numBin10_7[0]:
        if 0.80 < numBin10_7[1][i] < 0.90 :
            maxEvts10_7.append(numBin10_7[0][i])
        i+=1
    numEvts10_7 = np.sum(maxEvts10_7)
    print (1e-6)*lumi[0], numEvts10_7
    lum10.append((10/(0.031180))*numEvts10_7)
    lum100.append((100/(0.031180))*numEvts10_7)
    
    ax = f.add_subplot(338)
    numBin10_8 = ax.hist(clow.applyCuts(xpi_low,cut10),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.90,1.00)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.90,1.00)
    ax.text(0.65, 0.95, '$x_\pi$=(0.90,1.00)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts10_8 = []
    for val in numBin10_8[0]:
        if 0.90 < numBin10_8[1][i] < 1.00 :
            maxEvts10_8.append(numBin10_8[0][i])
        i+=1
    numEvts10_8 = np.sum(maxEvts10_8)
    print (1e-6)*lumi[0], numEvts10_8
    lum10.append((10/(0.031180))*numEvts10_8)
    lum100.append((100/(0.031180))*numEvts10_8)
    
    plt.style.use('default')
    plt.close(f)
    f = plt.figure(figsize=(11.69,8.27))
    plt.style.use('classic')

    ax = f.add_subplot(331)
    numBin30_1 = ax.hist(clow.applyCuts(xpi_low,cut30),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.20,0.30)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.20,0.30)
    ax.text(0.65, 0.95, '$x_\pi$=(0.20,0.30)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts30_1 = []
    for val in numBin30_1[0]:
        if 0.20 < numBin30_1[1][i] < 0.30 :
            maxEvts30_1.append(numBin30_1[0][i])
        i+=1
    numEvts30_1 = np.sum(maxEvts30_1)
    print (1e-6)*lumi[1], numEvts30_1
    lum10.append((10/(0.031180))*numEvts30_1)
    lum100.append((100/(0.031180))*numEvts30_1)
    
    plt.title('$x_{\pi}$ plots [$Q^2$ = 30 $GeV^2$]', fontsize =20)
    
    ax = f.add_subplot(332)
    numBin30_2 = ax.hist(clow.applyCuts(xpi_low,cut30),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.30,0.40)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.30,0.40)
    ax.text(0.65, 0.95, '$x_\pi$=(0.30,0.40)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts30_2 = []
    for val in numBin30_2[0]:
        if 0.30 < numBin30_2[1][i] < 0.40 :
            maxEvts30_2.append(numBin30_2[0][i])
        i+=1
    numEvts30_2 = np.sum(maxEvts30_2)
    print (1e-6)*lumi[1], numEvts30_2
    lum10.append((10/(0.031180))*numEvts30_2)
    lum100.append((100/(0.031180))*numEvts30_2)
    
    ax = f.add_subplot(333)
    numBin30_3 = ax.hist(clow.applyCuts(xpi_low,cut30),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.40,0.50)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.40,0.50)
    ax.text(0.65, 0.95, '$x_\pi$=(0.40,0.50)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts30_3 = []
    for val in numBin30_3[0]:
        if 0.40 < numBin30_3[1][i] < 0.50 :
            maxEvts30_3.append(numBin30_3[0][i])
        i+=1
    numEvts30_3 = np.sum(maxEvts30_3)
    print (1e-6)*lumi[1], numEvts30_3
    lum10.append((10/(0.031180))*numEvts30_3)
    lum100.append((100/(0.031180))*numEvts30_3)
    
    ax = f.add_subplot(334)
    numBin30_4 = ax.hist(clow.applyCuts(xpi_low,cut30),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.50,0.60)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.50,0.60)
    ax.text(0.65, 0.95, '$x_\pi$=(0.50,0.60)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts30_4 = []
    for val in numBin30_4[0]:
        if 0.50 < numBin30_4[1][i] < 0.60 :
            maxEvts30_4.append(numBin30_4[0][i])
        i+=1
    numEvts30_4 = np.sum(maxEvts30_4)
    print (1e-6)*lumi[1], numEvts30_4
    lum10.append((10/(0.031180))*numEvts30_4)
    lum100.append((100/(0.031180))*numEvts30_4)
    
    ax = f.add_subplot(335)
    numBin30_5 = ax.hist(clow.applyCuts(xpi_low,cut30),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.60,0.70)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.60,0.70)
    ax.text(0.65, 0.95, '$x_\pi$=(0.60,0.70)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts30_5 = []
    for val in numBin30_5[0]:
        if 0.60 < numBin30_5[1][i] < 0.70 :
            maxEvts30_5.append(numBin30_5[0][i])
        i+=1
    numEvts30_5 = np.sum(maxEvts30_5)
    print (1e-6)*lumi[1], numEvts30_5
    lum10.append((10/(0.031180))*numEvts30_5)
    lum100.append((100/(0.031180))*numEvts30_5)
    
    ax = f.add_subplot(336)
    numBin30_6 = ax.hist(clow.applyCuts(xpi_low,cut30),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.70,0.80)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.70,0.80)
    ax.text(0.65, 0.95, '$x_\pi$=(0.70,0.80)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts30_6 = []
    for val in numBin30_6[0]:
        if 0.70 < numBin30_6[1][i] < 0.80 :
            maxEvts30_6.append(numBin30_6[0][i])
        i+=1
    numEvts30_6 = np.sum(maxEvts30_6)
    print (1e-6)*lumi[1], numEvts30_6
    lum10.append((10/(0.031180))*numEvts30_6)
    lum100.append((100/(0.031180))*numEvts30_6)
    
    ax = f.add_subplot(337)
    numBin30_7 = ax.hist(clow.applyCuts(xpi_low,cut30),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.80,0.90)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.80,0.90)
    ax.text(0.65, 0.95, '$x_\pi$=(0.80,0.90)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts30_7 = []
    for val in numBin30_7[0]:
        if 0.80 < numBin30_7[1][i] < 0.90 :
            maxEvts30_7.append(numBin30_7[0][i])
        i+=1
    numEvts30_7 = np.sum(maxEvts30_7)
    print (1e-6)*lumi[1], numEvts30_7
    lum10.append((10/(0.031180))*numEvts30_7)
    lum100.append((100/(0.031180))*numEvts30_7)
    
    ax = f.add_subplot(338)
    numBin30_8 = ax.hist(clow.applyCuts(xpi_low,cut30),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.90,1.00)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.90,1.00)
    ax.text(0.65, 0.95, '$x_\pi$=(0.90,1.00)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts30_8 = []
    for val in numBin30_8[0]:
        if 0.90 < numBin30_8[1][i] < 1.00 :
            maxEvts30_8.append(numBin30_8[0][i])
        i+=1
    numEvts30_8 = np.sum(maxEvts30_8)
    print (1e-6)*lumi[1], numEvts30_8
    lum10.append((10/(0.031180))*numEvts30_8)
    lum100.append((100/(0.031180))*numEvts30_8)
    
    plt.style.use('default')
    plt.close(f)
    f = plt.figure(figsize=(11.69,8.27))
    plt.style.use('classic')
    
    ax = f.add_subplot(331)
    numBin50_1 = ax.hist(clow.applyCuts(xpi_low,cut50),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.20,0.30)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.20,0.30)
    ax.text(0.65, 0.95, '$x_\pi$=(0.20,0.30)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts50_1 = []
    for val in numBin50_1[0]:
        if 0.20 < numBin50_1[1][i] < 0.30 :
            maxEvts50_1.append(numBin50_1[0][i])
        i+=1
    numEvts50_1 = np.sum(maxEvts50_1)
    print (1e-6)*lumi[2], numEvts50_1
    lum10.append((10/(0.031180))*numEvts50_1)
    lum100.append((100/(0.031180))*numEvts50_1)
    
    plt.title('$x_{\pi}$ plots [$Q^2$ = 50 $GeV^2$]', fontsize =20)
    
    ax = f.add_subplot(332)
    numBin50_2 = ax.hist(clow.applyCuts(xpi_low,cut50),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.30,0.40)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.30,0.40)
    ax.text(0.65, 0.95, '$x_\pi$=(0.30,0.40)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts50_2 = []
    for val in numBin50_2[0]:
        if 0.30 < numBin50_2[1][i] < 0.40 :
            maxEvts50_2.append(numBin50_2[0][i])
        i+=1
    numEvts50_2 = np.sum(maxEvts50_2)
    print (1e-6)*lumi[2], numEvts50_2
    lum10.append((10/(0.031180))*numEvts50_2)
    lum100.append((100/(0.031180))*numEvts50_2)
    
    ax = f.add_subplot(333)
    numBin50_3 = ax.hist(clow.applyCuts(xpi_low,cut50),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.40,0.50)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.40,0.50)
    ax.text(0.65, 0.95, '$x_\pi$=(0.40,0.50)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts50_3 = []
    for val in numBin50_3[0]:
        if 0.40 < numBin50_3[1][i] < 0.50 :
            maxEvts50_3.append(numBin50_3[0][i])
        i+=1
    numEvts50_3 = np.sum(maxEvts50_3)
    print (1e-6)*lumi[2], numEvts50_3
    lum10.append((10/(0.031180))*numEvts50_3)
    lum100.append((100/(0.031180))*numEvts50_3)
    
    ax = f.add_subplot(334)
    numBin50_4 = ax.hist(clow.applyCuts(xpi_low,cut50),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.50,0.60)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.50,0.60)
    ax.text(0.65, 0.95, '$x_\pi$=(0.50,0.60)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts50_4 = []
    for val in numBin50_4[0]:
        if 0.50 < numBin50_4[1][i] < 0.60 :
            maxEvts50_4.append(numBin50_4[0][i])
        i+=1
    numEvts50_4 = np.sum(maxEvts50_4)
    print (1e-6)*lumi[2], numEvts50_4
    lum10.append((10/(0.031180))*numEvts50_4)
    lum100.append((100/(0.031180))*numEvts50_4)
    
    ax = f.add_subplot(335)
    numBin50_5 = ax.hist(clow.applyCuts(xpi_low,cut50),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.60,0.70)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.60,0.70)
    ax.text(0.65, 0.95, '$x_\pi$=(0.60,0.70)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts50_5 = []
    for val in numBin50_5[0]:
        if 0.60 < numBin50_5[1][i] < 0.70 :
            maxEvts50_5.append(numBin50_5[0][i])
        i+=1
    numEvts50_5 = np.sum(maxEvts50_5)
    print (1e-6)*lumi[2], numEvts50_5
    lum10.append((10/(0.031180))*numEvts50_5)
    lum100.append((100/(0.031180))*numEvts50_5)
    
    ax = f.add_subplot(336)
    numBin50_6 = ax.hist(clow.applyCuts(xpi_low,cut50),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.70,0.80)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.70,0.80)
    ax.text(0.65, 0.95, '$x_\pi$=(0.70,0.80)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts50_6 = []
    for val in numBin50_6[0]:
        if 0.70 < numBin50_6[1][i] < 0.80 :
            maxEvts50_6.append(numBin50_6[0][i])
        i+=1
    numEvts50_6 = np.sum(maxEvts50_6)
    print (1e-6)*lumi[2], numEvts50_6
    lum10.append((10/(0.031180))*numEvts50_6)
    lum100.append((100/(0.031180))*numEvts50_6)
    
    ax = f.add_subplot(337)
    numBin50_7 = ax.hist(clow.applyCuts(xpi_low,cut50),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.80,0.90)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.80,0.90)
    ax.text(0.65, 0.95, '$x_\pi$=(0.80,0.90)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts50_7 = []
    for val in numBin50_7[0]:
        if 0.80 < numBin50_7[1][i] < 0.90 :
            maxEvts50_7.append(numBin50_7[0][i])
        i+=1
    numEvts50_7 = np.sum(maxEvts50_7)
    print (1e-6)*lumi[2], numEvts50_7
    lum10.append((10/(0.031180))*numEvts50_7)
    lum100.append((100/(0.031180))*numEvts50_7)
    
    ax = f.add_subplot(338)
    numBin50_8 = ax.hist(clow.applyCuts(xpi_low,cut50),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.90,1.00)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.90,1.00)
    ax.text(0.65, 0.95, '$x_\pi$=(0.90,1.00)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts50_8 = []
    for val in numBin50_8[0]:
        if 0.90 < numBin50_8[1][i] < 1.00 :
            maxEvts50_8.append(numBin50_8[0][i])
        i+=1
    numEvts50_8 = np.sum(maxEvts50_8)
    print (1e-6)*lumi[2], numEvts50_8
    lum10.append((10/(0.031180))*numEvts50_8)
    lum100.append((100/(0.031180))*numEvts50_8)
    
    plt.style.use('default')
    plt.close(f)
    f = plt.figure(figsize=(11.69,8.27))
    plt.style.use('classic')
    
    ax = f.add_subplot(331)
    numBin70_1 = ax.hist(clow.applyCuts(xpi_low,cut70),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.20,0.30)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.20,0.30)
    ax.text(0.65, 0.95, '$x_\pi$=(0.20,0.30)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts70_1 = []
    for val in numBin70_1[0]:
        if 0.20 < numBin70_1[1][i] < 0.30 :
            maxEvts70_1.append(numBin70_1[0][i])
        i+=1
    numEvts70_1 = np.sum(maxEvts70_1)
    print (1e-6)*lumi[3], numEvts70_1
    lum10.append((10/(0.031180))*numEvts70_1)
    lum100.append((100/(0.031180))*numEvts70_1)
    
    plt.title('$x_{\pi}$ plots [$Q^2$ = 70 $GeV^2$]', fontsize =20)
    
    ax = f.add_subplot(332)
    numBin70_2 = ax.hist(clow.applyCuts(xpi_low,cut70),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.30,0.40)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.30,0.40)
    ax.text(0.65, 0.95, '$x_\pi$=(0.30,0.40)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts70_2 = []
    for val in numBin70_2[0]:
        if 0.30 < numBin70_2[1][i] < 0.40 :
            maxEvts70_2.append(numBin70_2[0][i])
        i+=1
    numEvts70_2 = np.sum(maxEvts70_2)
    print (1e-6)*lumi[3], numEvts70_2
    lum10.append((10/(0.031180))*numEvts70_2)
    lum100.append((100/(0.031180))*numEvts70_2)
    
    ax = f.add_subplot(333)
    numBin70_3 = ax.hist(clow.applyCuts(xpi_low,cut70),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.40,0.50)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.40,0.50)
    ax.text(0.65, 0.95, '$x_\pi$=(0.40,0.50)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts70_3 = []
    for val in numBin70_3[0]:
        if 0.40 < numBin70_3[1][i] < 0.50 :
            maxEvts70_3.append(numBin70_3[0][i])
        i+=1
    numEvts70_3 = np.sum(maxEvts70_3)
    print (1e-6)*lumi[3], numEvts70_3
    lum10.append((10/(0.031180))*numEvts70_3)
    lum100.append((100/(0.031180))*numEvts70_3)
    
    ax = f.add_subplot(334)
    numBin70_4 = ax.hist(clow.applyCuts(xpi_low,cut70),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.50,0.60)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.50,0.60)
    ax.text(0.65, 0.95, '$x_\pi$=(0.50,0.60)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts70_4 = []
    for val in numBin70_4[0]:
        if 0.50 < numBin70_4[1][i] < 0.60 :
            maxEvts70_4.append(numBin70_4[0][i])
        i+=1
    numEvts70_4 = np.sum(maxEvts70_4)
    print (1e-6)*lumi[3], numEvts70_4
    lum10.append((10/(0.031180))*numEvts70_4)
    lum100.append((100/(0.031180))*numEvts70_4)
    
    ax = f.add_subplot(335)
    numBin70_5 = ax.hist(clow.applyCuts(xpi_low,cut70),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.60,0.70)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.60,0.70)
    ax.text(0.65, 0.95, '$x_\pi$=(0.60,0.70)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts70_5 = []
    for val in numBin70_5[0]:
        if 0.60 < numBin70_5[1][i] < 0.70 :
            maxEvts70_5.append(numBin70_5[0][i])
        i+=1
    numEvts70_5 = np.sum(maxEvts70_5)
    print (1e-6)*lumi[3], numEvts70_5
    lum10.append((10/(0.031180))*numEvts70_5)
    lum100.append((100/(0.031180))*numEvts70_5)
    
    ax = f.add_subplot(336)
    numBin70_6 = ax.hist(clow.applyCuts(xpi_low,cut70),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.70,0.80)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.70,0.80)
    ax.text(0.65, 0.95, '$x_\pi$=(0.70,0.80)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts70_6 = []
    for val in numBin70_6[0]:
        if 0.70 < numBin70_6[1][i] < 0.80 :
            maxEvts70_6.append(numBin70_6[0][i])
        i+=1
    numEvts70_6 = np.sum(maxEvts70_6)
    print (1e-6)*lumi[3], numEvts70_6
    lum10.append((10/(0.031180))*numEvts70_6)
    lum100.append((100/(0.031180))*numEvts70_6)
    
    ax = f.add_subplot(337)
    numBin70_7 = ax.hist(clow.applyCuts(xpi_low,cut70),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.80,0.90)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.80,0.90)
    ax.text(0.65, 0.95, '$x_\pi$=(0.80,0.90)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts70_7 = []
    for val in numBin70_7[0]:
        if 0.80 < numBin70_7[1][i] < 0.90 :
            maxEvts70_7.append(numBin70_7[0][i])
        i+=1
    numEvts70_7 = np.sum(maxEvts70_7)
    print (1e-6)*lumi[3], numEvts70_7
    lum10.append((10/(0.031180))*numEvts70_7)
    lum100.append((100/(0.031180))*numEvts70_7)
    
    ax = f.add_subplot(338)
    numBin70_8 = ax.hist(clow.applyCuts(xpi_low,cut70),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$x_\pi$=(0.90,1.00)')
    plt.subplots_adjust(hspace=0.3,wspace=0.3)
    # plt.xscale('log')
    plt.xlim(0.90,1.00)
    ax.text(0.65, 0.95, '$x_\pi$=(0.90,1.00)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    plt.xlabel('$x_{\pi}$')
    plt.ylabel('$N^{bin}_{\pi}$')
    i=0
    j=0
    maxEvts70_8 = []
    for val in numBin70_8[0]:
        if 0.90 < numBin70_8[1][i] < 1.00 :
            maxEvts70_8.append(numBin70_8[0][i])
        i+=1
    numEvts70_8 = np.sum(maxEvts70_8)
    print (1e-6)*lumi[3], numEvts70_8
    lum10.append((10/(0.031180))*numEvts70_8)
    lum100.append((100/(0.031180))*numEvts70_8)
    
    # ax = f.add_subplot(335)
    # numBin5 = ax.hist(clow.applyCuts(xpi_low,cut50),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$Q^2$=50 $GeV^2$')
    # plt.subplots_adjust(hspace=0.0,wspace=0.3)
    # # plt.xscale('log')
    # plt.xlim(1e-6,1.)
    # ax.text(0.65, 0.95, '$Q^2$=50 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # plt.xlabel('$x_{\pi}$')
    # plt.ylabel('$N^{bin}_{\pi}$')
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    # i=0
    # j=0
    # # for val in numBin5[0]:
    # #     if max(numBin5[0]) == numBin5[0][i]:
    # #         print numBin5[0][i], j+1, numBin5[1][i], i
    # #         j+=1
    # #         tmp = i
    # #     i+=1
    # # maxEvts5 = numBin5[0][tmp]
    # # print lumi[4], maxEvts5
    # # lum10.append((10/(1e-6*lumi)[4])*maxEvts5)
    # # lum100.append((100/(1e-6*lumi)[4])*maxEvts5)
    # maxEvts5 = []
    # for val in numBin5[0]:
    #     if 1e-2 < numBin5[1][i] < 1e-1 :
    #         maxEvts5.append(numBin5[0][i])
    #         i+=1
    # numEvts5 = np.sum(maxEvts5)
    # print (1e-6)*lumi[4], numEvts5
    # lum10.append((10/(1e-6*lumi)[4])*numEvts5)
    # lum100.append((100/(1e-6*lumi)[4])*numEvts5)
    # Out of range
    # ax = f.add_subplot(336)    
    # numBin6 = ax.hist(clow.applyCuts(xpi_low,cut100),bins=p1.setbin(xpi_low,200,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$Q^2$=100 $GeV^2$')
    # plt.subplots_adjust(hspace=0.0,wspace=0.3)
    # # plt.xscale('log')
    # plt.xlim(1e-6,1.)
    # ax.text(0.65, 0.95, '$Q^2$=100 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # plt.xlabel('$x_{\pi}$')
    # plt.ylabel('$N^{bin}_{\pi}$')
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    # i=0
    # j=0
    # # for val in numBin6[0]:
    # #     if max(numBin6[0]) == numBin6[0][i]:
    # #         print numBin6[0][i], j+1, numBin6[1][i], i
    # #         j+=1
    # #         tmp = i
    # #     i+=1
    # # maxEvts6 = numBin6[0][tmp]
    # # print lumi[5], maxEvts6
    # # lum10.append((10/(1e-6*lumi)[5])*maxEvts6)
    # # lum100.append((100/(1e-6*lumi)[5])*maxEvts6)
    # maxEvts6 = []
    # for val in numBin6[0]:
    #     if 1e-2 < numBin6[1][i] < 1e-1 :
    #         maxEvts6.append(numBin6[0][i])
    #         i+=1
    # numEvts6 = np.sum(maxEvts6)
    # print (1e-6)*lumi[5], numEvts6
    # lum10.append((10/(1e-6*lumi)[5])*numEvts6)
    # lum100.append((100/(1e-6*lumi)[5])*numEvts6)

    fout = open('/home/trottar/ResearchNP/JLEIC/eic_SF/LuminosityTable.txt','w') 
    
    lum_tuple = namedtuple('lum_tuple',['binQ2','binxpi','lumi','nbin','nbin_10','nbin_100','uncern_10','uncern_100'])

    lum_data = []

    tot_data = []

    Q2Val = [10,10,10,10,10,10,10,10,30,30,30,30,30,30,30,30,50,50,50,50,50,50,50,50,70,70,70,70,70,70,70,70]

    count = [1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8]
    
    xpiVal = [0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95]

    uncern  = np.sqrt(lum10)/lum10
    uncern100  = np.sqrt(lum100)/lum100
    
    numEvts = []
    
    for i in range(0,len(uncern)) :
        tmp = 'numEvts%i_%i' % (Q2Val[i],count[i])
        numEvts.append(eval(tmp))

    for i in range(0,len(uncern)) :
        lum_data.append(lum_tuple(Q2Val[i],xpiVal[i],tot_lumi[i],numEvts[i],lum10[i],lum100[i],uncern[i], uncern100[i]))
        tot_data.append(lum_data[i])
        # print tot_data[i]

    for t in tot_data:
        # print namedtuple_to_str(t)
        fout.write("%s\n" % namedtuple_to_str(t))

    fout.close()
    
    plt.style.use('default')
    plt.close(f)
    
    # Combined plots

    xpi_tot10 = []
    fpi_tot10 = []
    fpiuncern10 = []
    
    xpi10_1 = []
    fpi10_1 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut10)[0]:
        # print clow.applyCuts(xpi_low,cut10)[0][i]
        if 0.240 < evt < 0.260 :
            xpi10_1.append(evt)
            fpi10_1.append(clow.applyCuts(fpi_low,cut10)[0][i])
            # print "good", i, xpi10_1
            # print xpi10_1
        i+=1
    if len(fpi10_1) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi10_1)
        xpi_tot10.append(min(xpi10_1))
        fpi_tot10.append(min(fpi10_1))
        fpiuncern10.append(uncern[0]*min(fpi10_1))
    if len(fpi10_1) == 1:
        # print "accpted: ", (fpi10_1)
        xpi_tot10.append((xpi10_1[0]))
        fpi_tot10.append((fpi10_1[0]))
        fpiuncern10.append(uncern[0]*min(fpi10_1))
    xpi10_2 = []
    fpi10_2 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut10)[0]:
        # print clow.applyCuts(xpi_low,cut10)[0][i]
        if 0.340 < evt < 0.360 :
            xpi10_2.append(evt)
            fpi10_2.append(clow.applyCuts(fpi_low,cut10)[0][i])
            # print "good", i, xpi10_2
            # print xpi10_2
        i+=1
    if len(fpi10_2) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi10_2)
        xpi_tot10.append(min(xpi10_2))
        fpi_tot10.append(min(fpi10_2))
        fpiuncern10.append(uncern[1]*min(fpi10_2))
    if len(fpi10_2) == 1:
        # print "accpted: ", (fpi10_2)
        xpi_tot10.append((xpi10_2[0]))
        fpi_tot10.append((fpi10_2[0]))
        fpiuncern10.append(uncern[1]*min(fpi10_2))
    xpi10_3 = []
    fpi10_3 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut10)[0]:
        # print clow.applyCuts(xpi_low,cut10)[0][i]
        if 0.440 < evt < 0.460 :
            xpi10_3.append(evt)
            fpi10_3.append(clow.applyCuts(fpi_low,cut10)[0][i])
            # print "good", i, xpi10_3
            # print xpi10_3
        i+=1
    if len(fpi10_3) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi10_3)
        xpi_tot10.append(min(xpi10_3))
        fpi_tot10.append(min(fpi10_3))
        fpiuncern10.append(uncern[2]*min(fpi10_3))
    if len(fpi10_3) == 1:
        # print "accpted: ", (fpi10_3)
        xpi_tot10.append((xpi10_3[0]))
        fpi_tot10.append((fpi10_3[0]))
        fpiuncern10.append(uncern[2]*min(fpi10_3))
    xpi10_4 = []
    fpi10_4 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut10)[0]:
        # print clow.applyCuts(xpi_low,cut10)[0][i]
        if 0.540 < evt < 0.560 :
            xpi10_4.append(evt)
            fpi10_4.append(clow.applyCuts(fpi_low,cut10)[0][i])
            # print "good", i, xpi10_4
            # print xpi10_4
        i+=1
    if len(fpi10_4) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi10_4)
        xpi_tot10.append(min(xpi10_4))
        fpi_tot10.append(min(fpi10_4))
        fpiuncern10.append(uncern[3]*min(fpi10_4))
    if len(fpi10_4) == 1:
        # print "accpted: ", (fpi10_4)
        xpi_tot10.append((xpi10_4[0]))
        fpi_tot10.append((fpi10_4[0]))
        fpiuncern10.append(uncern[3]*min(fpi10_4))
    xpi10_5 = []
    fpi10_5 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut10)[0]:
        # print clow.applyCuts(xpi_low,cut10)[0][i]
        if 0.640 < evt < 0.660 :
            xpi10_5.append(evt)
            fpi10_5.append(clow.applyCuts(fpi_low,cut10)[0][i])
            # print "good", i, xpi10_5
            # print xpi10_5
        i+=1
    if len(fpi10_5) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi10_5)
        xpi_tot10.append(min(xpi10_5))
        fpi_tot10.append(min(fpi10_5))
        fpiuncern10.append(uncern[4]*min(fpi10_5))
    if len(fpi10_5) == 1:
        # print "accpted: ", (fpi10_5)
        xpi_tot10.append((xpi10_5[0]))
        fpi_tot10.append((fpi10_5[0]))
        fpiuncern10.append(uncern[4]*min(fpi10_5))
    xpi10_6 = []
    fpi10_6 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut10)[0]:
        # print clow.applyCuts(xpi_low,cut10)[0][i]
        if 0.740 < evt < 0.760 :
            xpi10_6.append(evt)
            fpi10_6.append(clow.applyCuts(fpi_low,cut10)[0][i])
            # print "good", i, xpi10_6
            # print xpi10_6
        i+=1
    if len(fpi10_6) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi10_6)
        xpi_tot10.append(min(xpi10_6))
        fpi_tot10.append(min(fpi10_6))
        fpiuncern10.append(uncern[5]*min(fpi10_6))
    if len(fpi10_6) == 1:
        # print "accpted: ", (fpi10_6)
        xpi_tot10.append((xpi10_6[0]))
        fpi_tot10.append((fpi10_6[0]))
        fpiuncern10.append(uncern[5]*min(fpi10_6))
    xpi10_7 = []
    fpi10_7 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut10)[0]:
        # print clow.applyCuts(xpi_low,cut10)[0][i]
        if 0.840 < evt < 0.860 :
            xpi10_7.append(evt)
            fpi10_7.append(clow.applyCuts(fpi_low,cut10)[0][i])
            # print "good", i, xpi10_7
            # print xpi10_7
        i+=1
    if len(fpi10_7) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi10_7)
        xpi_tot10.append(min(xpi10_7))
        fpi_tot10.append(min(fpi10_7))
        fpiuncern10.append(uncern[6]*min(fpi10_7))
    if len(fpi10_7) == 1:
        # print "accpted: ", (fpi10_7)
        xpi_tot10.append((xpi10_7[0]))
        fpi_tot10.append((fpi10_7[0]))
        fpiuncern10.append(uncern[6]*min(fpi10_7))
    xpi10_8 = []
    fpi10_8 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut10)[0]:
        # print clow.applyCuts(xpi_low,cut10)[0][i]
        if 0.940 < evt < 0.960 :
            xpi10_8.append(evt)
            fpi10_8.append(clow.applyCuts(fpi_low,cut10)[0][i])
            # print "good", i, xpi10_8
            # print xpi10_8
        i+=1
    if len(fpi10_8) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi10_8)
        xpi_tot10.append(min(xpi10_8))
        fpi_tot10.append(min(fpi10_8))
        fpiuncern10.append(uncern[7]*min(fpi10_8))
    if len(fpi10_8) == 1:
        # print "accpted: ", (fpi10_8)
        xpi_tot10.append((xpi10_8[0]))
        fpi_tot10.append((fpi10_8[0]))
        fpiuncern10.append(uncern[7]*min(fpi10_8))

    # print xpi_tot10
    # print fpi_tot10

    xpi_tot30 = []
    fpi_tot30 = []
    fpiuncern30 = []
    
    xpi30_1 = []
    fpi30_1 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut30)[0]:
        # print clow.applyCuts(xpi_low,cut30)[0][i]
        if 0.240 < evt < 0.260 :
            xpi30_1.append(evt)
            fpi30_1.append(clow.applyCuts(fpi_low,cut30)[0][i])
            # print "good", i, xpi30_1
            # print xpi30_1
        i+=1
    if len(fpi30_1) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi30_1)
        xpi_tot30.append(min(xpi30_1))
        fpi_tot30.append(min(fpi30_1))
        fpiuncern30.append(uncern[8]*min(fpi30_1))
    if len(fpi30_1) == 1:
        # print "accpted: ", (fpi30_1)
        xpi_tot30.append((xpi30_1[0]))
        fpi_tot30.append((fpi30_1[0]))
        fpiuncern30.append(uncern[8]*min(fpi30_1))
    xpi30_2 = []
    fpi30_2 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut30)[0]:
        # print clow.applyCuts(xpi_low,cut30)[0][i]
        if 0.440 < evt < 0.460 :
            xpi30_2.append(evt)
            fpi30_2.append(clow.applyCuts(fpi_low,cut30)[0][i])
            # print "good", i, xpi30_2
            # print xpi30_2
        i+=1
    if len(fpi30_2) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi30_2)
        xpi_tot30.append(min(xpi30_2))
        fpi_tot30.append(min(fpi30_2))
        fpiuncern30.append(uncern[9]*min(fpi30_2))
    if len(fpi30_2) == 1:
        # print "accpted: ", (fpi30_2)
        xpi_tot30.append((xpi30_2[0]))
        fpi_tot30.append((fpi30_2[0]))
        fpiuncern30.append(uncern[9]*min(fpi30_2))
    xpi30_3 = []
    fpi30_3 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut30)[0]:
        # print clow.applyCuts(xpi_low,cut30)[0][i]
        if 0.340 < evt < 0.360 :
            xpi30_3.append(evt)
            fpi30_3.append(clow.applyCuts(fpi_low,cut30)[0][i])
            # print "good", i, xpi30_3
            # print xpi30_3
        i+=1
    if len(fpi30_3) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi30_3)
        xpi_tot30.append(min(xpi30_3))
        fpi_tot30.append(min(fpi30_3))
        fpiuncern30.append(uncern[10]*min(fpi30_3))
    if len(fpi30_3) == 1:
        # print "accpted: ", (fpi30_3)
        xpi_tot30.append((xpi30_3[0]))
        fpi_tot30.append((fpi30_3[0]))
        fpiuncern30.append(uncern[10]*min(fpi30_3))
    xpi30_4 = []
    fpi30_4 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut30)[0]:
        # print clow.applyCuts(xpi_low,cut30)[0][i]
        if 0.540 < evt < 0.560 :
            xpi30_4.append(evt)
            fpi30_4.append(clow.applyCuts(fpi_low,cut30)[0][i])
            # print "good", i, xpi30_4
            # print xpi30_4
        i+=1
    if len(fpi30_4) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi30_4)
        xpi_tot30.append(min(xpi30_4))
        fpi_tot30.append(min(fpi30_4))
        fpiuncern30.append(uncern[11]*min(fpi30_4))
    if len(fpi30_4) == 1:
        # print "accpted: ", (fpi30_4)
        xpi_tot30.append((xpi30_4[0]))
        fpi_tot30.append((fpi30_4[0]))
        fpiuncern30.append(uncern[11]*min(fpi30_4))
    xpi30_5 = []
    fpi30_5 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut30)[0]:
        # print clow.applyCuts(xpi_low,cut30)[0][i]
        if 0.640 < evt < 0.660 :
            xpi30_5.append(evt)
            fpi30_5.append(clow.applyCuts(fpi_low,cut30)[0][i])
            # print "good", i, xpi30_5
            # print xpi30_5
        i+=1
    if len(fpi30_5) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi30_5)
        xpi_tot30.append(min(xpi30_5))
        fpi_tot30.append(min(fpi30_5))
        fpiuncern30.append(uncern[12]*min(fpi30_5))
    if len(fpi30_5) == 1:
        # print "accpted: ", (fpi30_5)
        xpi_tot30.append((xpi30_5[0]))
        fpi_tot30.append((fpi30_5[0]))
        fpiuncern30.append(uncern[12]*min(fpi30_5))
    xpi30_6 = []
    fpi30_6 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut30)[0]:
        # print clow.applyCuts(xpi_low,cut30)[0][i]
        if 0.740 < evt < 0.760 :
            xpi30_6.append(evt)
            fpi30_6.append(clow.applyCuts(fpi_low,cut30)[0][i])
            # print "good", i, xpi30_6
            # print xpi30_6
        i+=1
    if len(fpi30_6) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi30_6)
        xpi_tot30.append(min(xpi30_6))
        fpi_tot30.append(min(fpi30_6))
        fpiuncern30.append(uncern[13]*min(fpi30_6))
    if len(fpi30_6) == 1:
        # print "accpted: ", (fpi30_6)
        xpi_tot30.append((xpi30_6[0]))
        fpi_tot30.append((fpi30_6[0]))
        fpiuncern30.append(uncern[13]*min(fpi30_6))
    xpi30_7 = []
    fpi30_7 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut30)[0]:
        # print clow.applyCuts(xpi_low,cut30)[0][i]
        if 0.840 < evt < 0.860 :
            xpi30_7.append(evt)
            fpi30_7.append(clow.applyCuts(fpi_low,cut30)[0][i])
            # print "good", i, xpi30_7
            # print xpi30_7
        i+=1
    if len(fpi30_7) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi30_7)
        xpi_tot30.append(min(xpi30_7))
        fpi_tot30.append(min(fpi30_7))
        fpiuncern30.append(uncern[14]*min(fpi30_7))
    if len(fpi30_7) == 1:
        # print "accpted: ", (fpi30_7)
        xpi_tot30.append((xpi30_7[0]))
        fpi_tot30.append((fpi30_7[0]))
        fpiuncern30.append(uncern[14]*min(fpi30_7))
    xpi30_8 = []
    fpi30_8 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut30)[0]:
        # print clow.applyCuts(xpi_low,cut30)[0][i]
        if 0.940 < evt < 0.960 :
            xpi30_8.append(evt)
            fpi30_8.append(clow.applyCuts(fpi_low,cut30)[0][i])
            # print "good", i, xpi30_8
            # print xpi30_8
        i+=1
    if len(fpi30_8) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi30_8)
        xpi_tot30.append(min(xpi30_8))
        fpi_tot30.append(min(fpi30_8))
        fpiuncern30.append(uncern[15]*min(fpi30_8))
    if len(fpi30_8) == 1:
        # print "accpted: ", (fpi30_8)
        xpi_tot30.append((xpi30_8[0]))
        fpi_tot30.append((fpi30_8[0]))
        fpiuncern30.append(uncern[15]*min(fpi30_8))

    # print xpi_tot30
    # print fpi_tot30

    xpi_tot50 = []
    fpi_tot50 = []
    fpiuncern50 = []
    
    xpi50_1 = []
    fpi50_1 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut50)[0]:
        # print clow.applyCuts(xpi_low,cut50)[0][i]
        if 0.240 < evt < 0.260 :
            xpi50_1.append(evt)
            fpi50_1.append(clow.applyCuts(fpi_low,cut50)[0][i])
            # print "good", i, xpi50_1
            # print xpi50_1
        i+=1
    if len(fpi50_1) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi50_1)
        xpi_tot50.append(min(xpi50_1))
        fpi_tot50.append(min(fpi50_1))
        fpiuncern50.append(uncern[16]*min(fpi50_1))
    if len(fpi50_1) == 1:
        # print "accpted: ", (fpi50_1)
        xpi_tot50.append((xpi50_1[0]))
        fpi_tot50.append((fpi50_1[0]))
        fpiuncern50.append(uncern[16]*min(fpi50_1))
    xpi50_2 = []
    fpi50_2 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut50)[0]:
        # print clow.applyCuts(xpi_low,cut50)[0][i]
        if 0.340 < evt < 0.360 :
            xpi50_2.append(evt)
            fpi50_2.append(clow.applyCuts(fpi_low,cut50)[0][i])
            # print "good", i, xpi50_2
            # print xpi50_2
        i+=1
    if len(fpi50_2) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi50_2)
        xpi_tot50.append(min(xpi50_2))
        fpi_tot50.append(min(fpi50_2))
        fpiuncern50.append(uncern[17]*min(fpi50_2))
    if len(fpi50_2) == 1:
        # print "accpted: ", (fpi50_2)
        xpi_tot50.append((xpi50_2[0]))
        fpi_tot50.append((fpi50_2[0]))
        fpiuncern50.append(uncern[17]*min(fpi50_2))
    xpi50_3 = []
    fpi50_3 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut50)[0]:
        # print clow.applyCuts(xpi_low,cut50)[0][i]
        if 0.440 < evt < 0.460 :
            xpi50_3.append(evt)
            fpi50_3.append(clow.applyCuts(fpi_low,cut50)[0][i])
            # print "good", i, xpi50_3
            # print xpi50_3
        i+=1
    if len(fpi50_3) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi50_3)
        xpi_tot50.append(min(xpi50_3))
        fpi_tot50.append(min(fpi50_3))
        fpiuncern50.append(uncern[18]*min(fpi50_3))
    if len(fpi50_3) == 1:
        # print "accpted: ", (fpi50_3)
        xpi_tot50.append((xpi50_3[0]))
        fpi_tot50.append((fpi50_3[0]))
        fpiuncern50.append(uncern[18]*min(fpi50_3))
    xpi50_4 = []
    fpi50_4 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut50)[0]:
        # print clow.applyCuts(xpi_low,cut50)[0][i]
        if 0.540 < evt < 0.560 :
            xpi50_4.append(evt)
            fpi50_4.append(clow.applyCuts(fpi_low,cut50)[0][i])
            # print "good", i, xpi50_4
            # print xpi50_4
        i+=1
    if len(fpi50_4) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi50_4)
        xpi_tot50.append(min(xpi50_4))
        fpi_tot50.append(min(fpi50_4))
        fpiuncern50.append(uncern[19]*min(fpi50_4))
    if len(fpi50_4) == 1:
        # print "accpted: ", (fpi50_4)
        xpi_tot50.append((xpi50_4[0]))
        fpi_tot50.append((fpi50_4[0]))
        fpiuncern50.append(uncern[19]*min(fpi50_4))
    xpi50_5 = []
    fpi50_5 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut50)[0]:
        # print clow.applyCuts(xpi_low,cut50)[0][i]
        if 0.640 < evt < 0.660 :
            xpi50_5.append(evt)
            fpi50_5.append(clow.applyCuts(fpi_low,cut50)[0][i])
            # print "good", i, xpi50_5
            # print xpi50_5
        i+=1
    if len(fpi50_5) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi50_5)
        xpi_tot50.append(min(xpi50_5))
        fpi_tot50.append(min(fpi50_5))
        fpiuncern50.append(uncern[20]*min(fpi50_5))
    if len(fpi50_5) == 1:
        # print "accpted: ", (fpi50_5)
        xpi_tot50.append((xpi50_5[0]))
        fpi_tot50.append((fpi50_5[0]))
        fpiuncern50.append(uncern[20]*min(fpi50_5))
    xpi50_6 = []
    fpi50_6 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut50)[0]:
        # print clow.applyCuts(xpi_low,cut50)[0][i]
        if 0.740 < evt < 0.760 :
            xpi50_6.append(evt)
            fpi50_6.append(clow.applyCuts(fpi_low,cut50)[0][i])
            # print "good", i, xpi50_6
            # print xpi50_6
        i+=1
    if len(fpi50_6) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi50_6)
        xpi_tot50.append(min(xpi50_6))
        fpi_tot50.append(min(fpi50_6))
        fpiuncern50.append(uncern[21]*min(fpi50_6))
    if len(fpi50_6) == 1:
        # print "accpted: ", (fpi50_6)
        xpi_tot50.append((xpi50_6[0]))
        fpi_tot50.append((fpi50_6[0]))
        fpiuncern50.append(uncern[21]*min(fpi50_6))
    xpi50_7 = []
    fpi50_7 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut50)[0]:
        # print clow.applyCuts(xpi_low,cut50)[0][i]
        if 0.840 < evt < 0.860 :
            xpi50_7.append(evt)
            fpi50_7.append(clow.applyCuts(fpi_low,cut50)[0][i])
            # print "good", i, xpi50_7
            # print xpi50_7
        i+=1
    if len(fpi50_7) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi50_7)
        xpi_tot50.append(min(xpi50_7))
        fpi_tot50.append(min(fpi50_7))
        fpiuncern50.append(uncern[22]*min(fpi50_7))
    if len(fpi50_7) == 1:
        # print "accpted: ", (fpi50_7)
        xpi_tot50.append((xpi50_7[0]))
        fpi_tot50.append((fpi50_7[0]))
        fpiuncern50.append(uncern[22]*min(fpi50_7))
    xpi50_8 = []
    fpi50_8 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut50)[0]:
        # print clow.applyCuts(xpi_low,cut50)[0][i]
        if 0.940 < evt < 0.960 :
            xpi50_8.append(evt)
            fpi50_8.append(clow.applyCuts(fpi_low,cut50)[0][i])
            # print "good", i, xpi50_8
            # print xpi50_8
        i+=1
    if len(fpi50_8) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi50_8)
        xpi_tot50.append(min(xpi50_8))
        fpi_tot50.append(min(fpi50_8))
        fpiuncern50.append(uncern[23]*min(fpi50_8))
    if len(fpi50_8) == 1:
        # print "accpted: ", (fpi50_8)
        xpi_tot50.append((xpi50_8[0]))
        fpi_tot50.append((fpi50_8[0]))
        fpiuncern50.append(uncern[23]*min(fpi50_8))

    # print xpi_tot50
    # print fpi_tot50

    xpi_tot70 = []
    fpi_tot70 = []
    fpiuncern70 = []
    
    xpi70_1 = []
    fpi70_1 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut70)[0]:
        # print clow.applyCuts(xpi_low,cut70)[0][i]
        if 0.240 < evt < 0.260 :
            xpi70_1.append(evt)
            fpi70_1.append(clow.applyCuts(fpi_low,cut70)[0][i])
            # print "good", i, xpi70_1
            # print xpi70_1
        i+=1
    if len(fpi70_1) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi70_1)
        xpi_tot70.append(min(xpi70_1))
        fpi_tot70.append(min(fpi70_1))
        fpiuncern70.append(uncern[24]*min(fpi70_1))
    if len(fpi70_1) == 1:
        # print "accpted: ", (fpi70_1)
        xpi_tot70.append((xpi70_1[0]))
        fpi_tot70.append((fpi70_1[0]))
        fpiuncern70.append(uncern[24]*min(fpi70_1))
    xpi70_2 = []
    fpi70_2 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut70)[0]:
        # print clow.applyCuts(xpi_low,cut70)[0][i]
        if 0.340 < evt < 0.360 :
            xpi70_2.append(evt)
            fpi70_2.append(clow.applyCuts(fpi_low,cut70)[0][i])
            # print "good", i, xpi70_2
            # print xpi70_2
        i+=1
    if len(fpi70_2) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi70_2)
        xpi_tot70.append(min(xpi70_2))
        fpi_tot70.append(min(fpi70_2))
        fpiuncern70.append(uncern[25]*min(fpi70_2))
    if len(fpi70_2) == 1:
        # print "accpted: ", (fpi70_2)
        xpi_tot70.append((xpi70_2[0]))
        fpi_tot70.append((fpi70_2[0]))
        fpiuncern70.append(uncern[25]*min(fpi70_2))
    xpi70_3 = []
    fpi70_3 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut70)[0]:
        # print clow.applyCuts(xpi_low,cut70)[0][i]
        if 0.440 < evt < 0.460 :
            xpi70_3.append(evt)
            fpi70_3.append(clow.applyCuts(fpi_low,cut70)[0][i])
            # print "good", i, xpi70_3
            # print xpi70_3
        i+=1
    if len(fpi70_3) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi70_3)
        xpi_tot70.append(min(xpi70_3))
        fpi_tot70.append(min(fpi70_3))
        fpiuncern70.append(uncern[26]*min(fpi70_3))
    if len(fpi70_3) == 1:
        # print "accpted: ", (fpi70_3)
        xpi_tot70.append((xpi70_3[0]))
        fpi_tot70.append((fpi70_3[0]))
        fpiuncern70.append(uncern[26]*min(fpi70_3))
    xpi70_4 = []
    fpi70_4 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut70)[0]:
        # print clow.applyCuts(xpi_low,cut70)[0][i]
        if 0.540 < evt < 0.560 :
            xpi70_4.append(evt)
            fpi70_4.append(clow.applyCuts(fpi_low,cut70)[0][i])
            # print "good", i, xpi70_4
            # print xpi70_4
        i+=1
    if len(fpi70_4) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi70_4)
        xpi_tot70.append(min(xpi70_4))
        fpi_tot70.append(min(fpi70_4))
        fpiuncern70.append(uncern[27]*min(fpi70_4))
    if len(fpi70_4) == 1:
        # print "accpted: ", (fpi70_4)
        xpi_tot70.append((xpi70_4[0]))
        fpi_tot70.append((fpi70_4[0]))
        fpiuncern70.append(uncern[27]*min(fpi70_4))
    xpi70_5 = []
    fpi70_5 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut70)[0]:
        # print clow.applyCuts(xpi_low,cut70)[0][i]
        if 0.640 < evt < 0.660 :
            xpi70_5.append(evt)
            fpi70_5.append(clow.applyCuts(fpi_low,cut70)[0][i])
            # print "good", i, xpi70_5
            # print xpi70_5
        i+=1
    if len(fpi70_5) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi70_5)
        xpi_tot70.append(min(xpi70_5))
        fpi_tot70.append(min(fpi70_5))
        fpiuncern70.append(uncern[28]*min(fpi70_5))
    if len(fpi70_5) == 1:
        # print "accpted: ", (fpi70_5)
        xpi_tot70.append((xpi70_5[0]))
        fpi_tot70.append((fpi70_5[0]))
        fpiuncern70.append(uncern[28]*min(fpi70_5))
    xpi70_6 = []
    fpi70_6 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut70)[0]:
        # print clow.applyCuts(xpi_low,cut70)[0][i]
        if 0.740 < evt < 0.760 :
            xpi70_6.append(evt)
            fpi70_6.append(clow.applyCuts(fpi_low,cut70)[0][i])
            # print "good", i, xpi70_6
            # print xpi70_6
        i+=1
    if len(fpi70_6) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi70_6)
        xpi_tot70.append(min(xpi70_6))
        fpi_tot70.append(min(fpi70_6))
        fpiuncern70.append(uncern[29]*min(fpi70_6))
    if len(fpi70_6) == 1:
        # print "accpted: ", (fpi70_6)
        xpi_tot70.append((xpi70_6[0]))
        fpi_tot70.append((fpi70_6[0]))
        fpiuncern70.append(uncern[29]*min(fpi70_6))
    xpi70_7 = []
    fpi70_7 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut70)[0]:
        # print clow.applyCuts(xpi_low,cut70)[0][i]
        if 0.840 < evt < 0.860 :
            xpi70_7.append(evt)
            fpi70_7.append(clow.applyCuts(fpi_low,cut70)[0][i])
            # print "good", i, xpi70_7
            # print xpi70_7
        i+=1
    if len(fpi70_7) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi70_7)
        xpi_tot70.append(min(xpi70_7))
        fpi_tot70.append(min(fpi70_7))
        fpiuncern70.append(uncern[30]*min(fpi70_7))
    if len(fpi70_7) == 1:
        # print "accpted: ", (fpi70_7)
        xpi_tot70.append((xpi70_7[0]))
        fpi_tot70.append((fpi70_7[0]))
        fpiuncern70.append(uncern[30]*min(fpi70_7))
    xpi70_8 = []
    fpi70_8 = []
    i=0
    for evt in clow.applyCuts(xpi_low,cut70)[0]:
        # print clow.applyCuts(xpi_low,cut70)[0][i]
        if 0.940 < evt < 0.960 :
            xpi70_8.append(evt)
            fpi70_8.append(clow.applyCuts(fpi_low,cut70)[0][i])
            # print "good", i, xpi70_8
            # print xpi70_8
        i+=1
    if len(fpi70_8) > 1:
        # print "multiple"
        # print "accpted: ", min(fpi70_8)
        xpi_tot70.append(min(xpi70_8))
        fpi_tot70.append(min(fpi70_8))
        fpiuncern70.append(uncern[31]*min(fpi70_8))
    if len(fpi70_8) == 1:
        # print "accpted: ", (fpi70_8)
        xpi_tot70.append((xpi70_8[0]))
        fpi_tot70.append((fpi70_8[0]))
        fpiuncern70.append(uncern[31]*min(fpi70_8))

    # print xpi_tot70
    # print fpi_tot70
    
    print "xpi\n",(xpi_tot10),"fpi\n",(fpi_tot10), "uncern\n",(fpiuncern10)
    print len(xpi_tot10),len(fpi_tot10), len(fpiuncern10)
    print len(xpi_tot30),len(fpi_tot30), len(fpiuncern30)
    print len(xpi_tot50),len(fpi_tot50), len(fpiuncern50)
    print len(xpi_tot70),len(fpi_tot70), len(fpiuncern70)

    # counts 0 to 1 over intervals of 1e-3
    # x = np.linspace(0.,1.,64397,endpoint=True)
    # x = np.sort(TDIS_xbj_low)
    x = np.sort(xpi_low)

    Q2 = Q2_low

    # fit_f2pi_10 = piSF_BLFQ(x)
    fit_f2pi_10 = piSF_GRV(x)
    
    print "->", len(fit_f2pi_10)
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    uncernPlots10_10 = ax.errorbar(xpi_tot10,fpi_tot10,yerr=fpiuncern10,fmt='.',label='$Q^2$ = 10 $GeV^2$')
    BLFQ_Plots10_10 = ax.errorbar([0.0009,0.002,0.004,0.0085,0.018],[0.52,0.45,0.40,0.32,0.25],fmt='.',label='DESY-HERA-H1 \n $Q^2$ = 11 $GeV^2$')
    GRV_Plots10_10 = ax.plot([0.01,0.1,0.30],[0.25,0.20,0.18],label='GRV fit \n $Q^2$ = 7 $GeV^2$')
    # fitPlot_10 = ax.plot(x,fit_f2pi_10)
    # ax.fill_between(xpi_tot10, fit_f2pi_10-error, fit_f2pi_10+error)
    plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(0.01,1.0)
    plt.ylim(1e-6,0.4)
    plt.xlabel('$x_\pi$')
    plt.ylabel('$F^{2}_{\pi}$')
    plt.title('$F^{2}_{\pi}$ vs $x_\pi$', fontsize =20)
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    uncernPlots10_30 = ax.errorbar(xpi_tot30,fpi_tot30,yerr=fpiuncern30,fmt='.',label='$Q^2$ = 30 $GeV^2$')
    BLFQ_Plots10_30 = ax.errorbar([0.002,0.004,0.008,0.02,0.045],[0.6,0.48,0.35,0.30,0.25],fmt='.',label='DESY-HERA-H1 \n $Q^2$ = 24 $GeV^2$')
    GRV_Plots10_30 = ax.plot([0.01,0.1,0.30],[0.40,0.20,0.16],label='GRV fit \n $Q^2$ = 30 $GeV^2$')
    plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(0.01,1.0)
    plt.ylim(1e-6,0.4)
    plt.xlabel('$x_\pi$')
    plt.ylabel('$F^{2}_{\pi}$')
    plt.title('$F^{2}_{\pi}$ vs $x_\pi$', fontsize =20)
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    uncernPlots10_50 = ax.errorbar(xpi_tot50,fpi_tot50,yerr=fpiuncern50,fmt='.',label='$Q^2$ = 50 $GeV^2$')
    BLFQ_Plots10_50 = ax.errorbar([0.008,0.019,0.038,0.075],[0.40,0.38,0.26,0.22],fmt='.',label='DESY-HERA-H1 \n $Q^2$ = 55 $GeV^2$')
    GRV_Plots10_50 = ax.plot([0.01,0.1,0.30],[0.40,0.20,0.15],label='GRV fit \n $Q^2$ = 60 $GeV^2$')
    plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(0.01,1.0)
    plt.ylim(1e-6,0.4)
    plt.xlabel('$x_\pi$')
    plt.ylabel('$F^{2}_{\pi}$')
    plt.title('$F^{2}_{\pi}$ vs $x_\pi$', fontsize =20)
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    uncernPlots10_70 = ax.errorbar(xpi_tot70,fpi_tot70,yerr=fpiuncern70,fmt='.',label='$Q^2$ = 70 $GeV^2$')
    BLFQ_Plots10_70 = ax.errorbar([0.023,0.045,0.09],[0.38,0.29,0.32],fmt='.',label='DESY-HERA-H1 \n $Q^2$ = 82 $GeV^2$')
    GRV_Plots10_70 = ax.plot([0.01,0.1,0.30],[0.40,0.20,0.15],label='GRV fit \n $Q^2$ = 60 $GeV^2$')
    plt.xscale('log')
    # plt.yscale('log')
    plt.ylim(1e-6,0.4)
    plt.xlim(0.01,1.0)
    plt.xlabel('$x_\pi$')
    plt.ylabel('$F^{2}_{\pi}$')
    plt.title('$F^{2}_{\pi}$ vs $x_\pi$', fontsize =20)
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)

    ####### Linear x, log y
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    uncernPlots10_10 = ax.errorbar(xpi_tot10,fpi_tot10,yerr=fpiuncern10,fmt='.',label='$Q^2$ = 10 $GeV^2$')
    GRV_Plots10_10 = ax.plot([0.01,0.1,0.30],[0.25,0.20,0.18],label='GRV fit \n $Q^2$ = 7 $GeV^2$')
    # fitPlot_10 = ax.plot(x, fit_f2pi_10)
    # ax.fill_between(xpi_tot10, fit_f2pi_10-error, fit_f2pi_10+error)
    # plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.1,1.0)
    plt.ylim(1e-6,0.025)
    plt.xlabel('$x_\pi$')
    plt.ylabel('$F^{2}_{\pi}$')
    plt.title('$F^{2}_{\pi}$ vs $x_\pi$', fontsize =20)
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    uncernPlots10_30 = ax.errorbar(xpi_tot30,fpi_tot30,yerr=fpiuncern30,fmt='.',label='$Q^2$ = 30 $GeV^2$')
    GRV_Plots10_30 = ax.plot([0.01,0.1,0.30],[0.40,0.20,0.16],label='GRV fit \n $Q^2$ = 30 $GeV^2$')
    # plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.1,1.0)
    plt.ylim(1e-6,0.025)
    plt.xlabel('$x_\pi$')
    plt.ylabel('$F^{2}_{\pi}$')
    plt.title('$F^{2}_{\pi}$ vs $x_\pi$', fontsize =20)
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    uncernPlots10_50 = ax.errorbar(xpi_tot50,fpi_tot50,yerr=fpiuncern50,fmt='.',label='$Q^2$ = 50 $GeV^2$')
    GRV_Plots10_50 = ax.plot([0.01,0.1,0.30],[0.40,0.20,0.15],label='GRV fit \n $Q^2$ = 60 $GeV^2$')
    # plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.1,1.0)
    plt.ylim(1e-6,0.025)
    plt.xlabel('$x_\pi$')
    plt.ylabel('$F^{2}_{\pi}$')
    plt.title('$F^{2}_{\pi}$ vs $x_\pi$', fontsize =20)
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    uncernPlots10_70 = ax.errorbar(xpi_tot70,fpi_tot70,yerr=fpiuncern70,fmt='.',label='$Q^2$ = 70 $GeV^2$')
    GRV_Plots10_70 = ax.plot([0.01,0.1,0.30],[0.40,0.20,0.15],label='GRV fit \n $Q^2$ = 60 $GeV^2$')
    # plt.xscale('log')
    plt.yscale('log')
    plt.ylim(1e-6,0.025)
    plt.xlim(0.1,1.0)
    plt.xlabel('$x_\pi$')
    plt.ylabel('$F^{2}_{\pi}$')
    plt.title('$F^{2}_{\pi}$ vs $x_\pi$', fontsize =20)
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    f2PxbjPlot = ax.scatter(TDIS_xbj_low,f2N_low)
    plt.xlim(0.,1.)
    # plt.ylim(1e-6,0.5)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('$F^{2}_{P}$')
    plt.title('$F^{2}_{P}$ vs TDIS_xbj', fontsize =20)
    plt.close(f)
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    fpixbjPlot = ax.scatter(TDIS_xbj_low,fpi_low)
    plt.xlim(0.,1.)
    plt.ylim(1e-6,0.5)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('$F^{2}_{\pi}$')
    plt.title('$F^{2}_{\pi}$ vs TDIS_xbj', fontsize =20)
    plt.close(f)
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    fpixbjPlotLog = ax.scatter(TDIS_xbj_low,fpi_low)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.,1.)
    plt.ylim(1e-6,0.5)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('$F^{2}_{\pi}$')
    plt.title('$F^{2}_{\pi}$ vs TDIS_xbj', fontsize =20)
    plt.close(f)
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    fpixpiPlot = ax.scatter(xpi_low,fpi_low)
    plt.xlim(0.,1.)
    # plt.ylim(1e-6,0.5)
    plt.xlabel('$x_\pi$')
    plt.ylabel('$F^{2}_{\pi}$')
    plt.title('$F^{2}_{\pi}$ vs $x_\pi$', fontsize =20)
    plt.close(f)
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    fpixpiPlotLog = ax.scatter(xpi_low,fpi_low)
    plt.xscale('log')
    # plt.yscale('log')
    plt.xlim(0.,1.)
    # plt.ylim(1e-6,0.5)
    plt.xlabel('$x_\pi$')
    plt.ylabel('$F^{2}_{\pi}$')
    plt.title('$F^{2}_{\pi}$ vs $x_\pi$', fontsize =20)
    plt.close(f)
    
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    xL_Plot = ax.hist(xL,bins=p1.setbin(xL,200,0.,1.)[0],label='all events',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.xlabel('$x_L$')
    plt.ylabel('$F^{2}_{\pi}$')
    plt.title('$F^{2}_{\pi}$ vs $x_L$', fontsize =20)
    plt.close(f)
    
    
    phaseSpace10_10 = densityPlot(xpi_low, Q2_low, '$Q^2$ vs $x_\pi$','$x_\pi$','$Q^{2}$', 200, 200, 0., 1.0, 0., 100.)
    # phaseSpace10_10 = ax.scatter(clow.applyCuts(TDIS_xbj_low,tcut1),clow.applyCuts(Q2_low,tcut1))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.,1.)
    plt.close(f)
    
    xpiPhase = xpi_low
    
    print len(xpiPhase)
    
    i=0
    tot_evts = []
    for evt in xpiPhase:
        if 0.1 < evt :
            tot_evts.append(evt)
    print len(tot_evts)
    
    # # Need to have numBin for each setting and the uncertainty will be the value uncern
    # f = plt.figure(figsize=(11.69,8.27))
    # plt.style.use('classic')
    
    # ax = f.add_subplot(331)    
    # fpiscat10_1 = ax.errorbar(clow.applyCuts(xpi_low,cut10),clow.applyCuts(fpi_low,cut10),yerr=uncern[0],fmt='.',label='$x_\pi$=(0.20,0.30)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.20,0.30)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.20,0.30)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))

    # plt.title('$F^{2}_{\pi}$ vs xpi  [$Q^2$ = 10 $GeV^2$]', fontsize =20)
    
    # ax = f.add_subplot(332)    
    # fpiscat10_2 = ax.errorbar(clow.applyCuts(xpi_low,cut10),clow.applyCuts(fpi_low,cut10),yerr=uncern[1],fmt='.',label='$x_\pi$=(0.30,0.40)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.30,0.40)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.30,0.40)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(333)
    # fpiscat10_3 = ax.errorbar(clow.applyCuts(xpi_low,cut10),clow.applyCuts(fpi_low,cut10),yerr=uncern[2],fmt='.',label='$x_\pi$=(0.40,0.50)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.40,0.50)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.40,0.50)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(334)
    # fpiscat10_4 = ax.errorbar(clow.applyCuts(xpi_low,cut10),clow.applyCuts(fpi_low,cut10),yerr=uncern[3],fmt='.',label='$x_\pi$=(0.50,0.60)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.50,0.60)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.50,0.60)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(335)
    # fpiscat10_5 = ax.errorbar(clow.applyCuts(xpi_low,cut10),clow.applyCuts(fpi_low,cut10),yerr=uncern[4],fmt='.',label='$x_\pi$=(0.60,0.70)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.60,0.70)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.60,0.70)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(336)
    # fpiscat10_6 = ax.errorbar(clow.applyCuts(xpi_low,cut10),clow.applyCuts(fpi_low,cut10),yerr=uncern[5],fmt='.',label='$x_\pi$=(0.70,0.80)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.70,0.80)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.70,0.80)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(337)
    # fpiscat10_7 = ax.errorbar(clow.applyCuts(xpi_low,cut10),clow.applyCuts(fpi_low,cut10),yerr=uncern[6],fmt='.',label='$x_\pi$=(0.80,0.90)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.80,0.90)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.80,0.90)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(338)
    # fpiscat10_8 = ax.errorbar(clow.applyCuts(xpi_low,cut10),clow.applyCuts(fpi_low,cut10),yerr=uncern[7],fmt='.',label='$x_\pi$=(0.90,1.00)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.90,1.00)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.90,1.00)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # plt.style.use('default')
    # f = plt.figure(figsize=(11.69,8.27))
    # plt.style.use('classic')
    
    # ax = f.add_subplot(331)    
    # fpiscat30_1 = ax.errorbar(clow.applyCuts(xpi_low,cut30),clow.applyCuts(fpi_low,cut30),yerr=uncern[8],fmt='.',label='$x_\pi$=(0.20,0.30)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.20,0.30)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.20,0.30)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))

    # plt.title('$F^{2}_{\pi}$ vs xpi  [$Q^2$ = 30 $GeV^2$]', fontsize =20)
    
    # ax = f.add_subplot(332)    
    # fpiscat30_2 = ax.errorbar(clow.applyCuts(xpi_low,cut30),clow.applyCuts(fpi_low,cut30),yerr=uncern[9],fmt='.',label='$x_\pi$=(0.30,0.40)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.30,0.40)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.30,0.40)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(333)
    # fpiscat30_3 = ax.errorbar(clow.applyCuts(xpi_low,cut30),clow.applyCuts(fpi_low,cut30),yerr=uncern[10],fmt='.',label='$x_\pi$=(0.40,0.50)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.40,0.50)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.40,0.50)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(334)
    # fpiscat30_4 = ax.errorbar(clow.applyCuts(xpi_low,cut30),clow.applyCuts(fpi_low,cut30),yerr=uncern[11],fmt='.',label='$x_\pi$=(0.50,0.60)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.50,0.60)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.50,0.60)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(335)
    # fpiscat30_5 = ax.errorbar(clow.applyCuts(xpi_low,cut30),clow.applyCuts(fpi_low,cut30),yerr=uncern[12],fmt='.',label='$x_\pi$=(0.60,0.70)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.60,0.70)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.60,0.70)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(336)
    # fpiscat30_6 = ax.errorbar(clow.applyCuts(xpi_low,cut30),clow.applyCuts(fpi_low,cut30),yerr=uncern[13],fmt='.',label='$x_\pi$=(0.70,0.80)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.70,0.80)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.70,0.80)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(337)
    # fpiscat30_7 = ax.errorbar(clow.applyCuts(xpi_low,cut30),clow.applyCuts(fpi_low,cut30),yerr=uncern[14],fmt='.',label='$x_\pi$=(0.80,0.90)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.80,0.90)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.80,0.90)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(338)
    # fpiscat30_8 = ax.errorbar(clow.applyCuts(xpi_low,cut30),clow.applyCuts(fpi_low,cut30),yerr=uncern[15],fmt='.',label='$x_\pi$=(0.90,1.00)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.90,1.00)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.90,1.00)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    
    # plt.style.use('default')
    # f = plt.figure(figsize=(11.69,8.27))
    # plt.style.use('classic')

    # ax = f.add_subplot(331)    
    # fpiscat50_1 = ax.errorbar(clow.applyCuts(xpi_low,cut50),clow.applyCuts(fpi_low,cut50),yerr=uncern[16],fmt='.',label='$x_\pi$=(0.20,0.30)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.20,0.30)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.20,0.30)
    # # plt.ylim(1e-3,0.025)
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))

    # plt.title('$F^{2}_{\pi}$ vs xpi  [$Q^2$ = 50 $GeV^2$]', fontsize =20)
    
    # ax = f.add_subplot(332)    
    # fpiscat50_2 = ax.errorbar(clow.applyCuts(xpi_low,cut50),clow.applyCuts(fpi_low,cut50),yerr=uncern[17],fmt='.',label='$x_\pi$=(0.30,0.40)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.30,0.40)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.30,0.40)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(333)
    # fpiscat50_3 = ax.errorbar(clow.applyCuts(xpi_low,cut50),clow.applyCuts(fpi_low,cut50),yerr=uncern[18],fmt='.',label='$x_\pi$=(0.40,0.50)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.40,0.50)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.40,0.50)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(334)
    # fpiscat50_4 = ax.errorbar(clow.applyCuts(xpi_low,cut50),clow.applyCuts(fpi_low,cut50),yerr=uncern[19],fmt='.',label='$x_\pi$=(0.50,0.60)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.50,0.60)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.50,0.60)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(335)
    # fpiscat50_5 = ax.errorbar(clow.applyCuts(xpi_low,cut50),clow.applyCuts(fpi_low,cut50),yerr=uncern[20],fmt='.',label='$x_\pi$=(0.60,0.70)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.60,0.70)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.60,0.70)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(336)
    # fpiscat50_6 = ax.errorbar(clow.applyCuts(xpi_low,cut50),clow.applyCuts(fpi_low,cut50),yerr=uncern[21],fmt='.',label='$x_\pi$=(0.70,0.80)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.70,0.80)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.70,0.80)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(337)
    # fpiscat50_7 = ax.errorbar(clow.applyCuts(xpi_low,cut50),clow.applyCuts(fpi_low,cut50),yerr=uncern[22],fmt='.',label='$x_\pi$=(0.80,0.90)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.80,0.90)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.80,0.90)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(338)
    # fpiscat50_8 = ax.errorbar(clow.applyCuts(xpi_low,cut50),clow.applyCuts(fpi_low,cut50),yerr=uncern[23],fmt='.',label='$x_\pi$=(0.90,1.00)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.90,1.00)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.90,1.00)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    
    
    # plt.style.use('default')
    # f = plt.figure(figsize=(11.69,8.27))
    # plt.style.use('classic')
    
    # ax = f.add_subplot(331)    
    # fpiscat70_1 = ax.errorbar(clow.applyCuts(xpi_low,cut70),clow.applyCuts(fpi_low,cut70),yerr=uncern[24],fmt='.',label='$x_\pi$=(0.20,0.30)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.20,0.30)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.20,0.30)
    # # plt.ylim(1e-3,0.025)
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))

    # plt.title('$F^{2}_{\pi}$ vs xpi  [$Q^2$ = 70 $GeV^2$]', fontsize =20)
    
    # ax = f.add_subplot(332)    
    # fpiscat70_2 = ax.errorbar(clow.applyCuts(xpi_low,cut70),clow.applyCuts(fpi_low,cut70),yerr=uncern[25],fmt='.',label='$x_\pi$=(0.30,0.40)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.30,0.40)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.30,0.40)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(333)
    # fpiscat70_3 = ax.errorbar(clow.applyCuts(xpi_low,cut70),clow.applyCuts(fpi_low,cut70),yerr=uncern[26],fmt='.',label='$x_\pi$=(0.40,0.50)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.40,0.50)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.40,0.50)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(334)
    # fpiscat70_4 = ax.errorbar(clow.applyCuts(xpi_low,cut70),clow.applyCuts(fpi_low,cut70),yerr=uncern[27],fmt='.',label='$x_\pi$=(0.50,0.60)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.50,0.60)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.50,0.60)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(335)
    # fpiscat70_5 = ax.errorbar(clow.applyCuts(xpi_low,cut70),clow.applyCuts(fpi_low,cut70),yerr=uncern[28],fmt='.',label='$x_\pi$=(0.60,0.70)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.60,0.70)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.60,0.70)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(336)
    # fpiscat70_6 = ax.errorbar(clow.applyCuts(xpi_low,cut70),clow.applyCuts(fpi_low,cut70),yerr=uncern[29],fmt='.',label='$x_\pi$=(0.70,0.80)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.70,0.80)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.70,0.80)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(337)
    # fpiscat70_7 = ax.errorbar(clow.applyCuts(xpi_low,cut70),clow.applyCuts(fpi_low,cut70),yerr=uncern[30],fmt='.',label='$x_\pi$=(0.80,0.90)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.80,0.90)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.80,0.90)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))    

    # ax = f.add_subplot(338)
    # fpiscat70_8 = ax.errorbar(clow.applyCuts(xpi_low,cut70),clow.applyCuts(fpi_low,cut70),yerr=uncern[31],fmt='.',label='$x_\pi$=(0.90,1.00)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # ax.text(0.65, 0.95, '$x_\pi$=(0.90,1.00)', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # # plt.yscale('log')
    # plt.xlim(0.90,1.00)
    # # plt.ylim(1e-3,0.025)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # # ax.xaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_formatter(plt.NullFormatter())
    # # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))

    # uncern plots

    uncern10 = []
    uncern30 = []
    uncern50 = []
    uncern70 = []
    for i in range(0,8):
        uncern10.append(uncern[i])
    for i in range(8,16):
        uncern30.append(uncern[i])
    for i in range(16,24):
        uncern50.append(uncern[i])
    for i in range(24,32):
        uncern70.append(uncern[i])
    
    # uncern 100 plots
    
    uncern10 = []
    uncern30 = []
    uncern50 = []
    uncern70 = []
    for i in range(0,8):
        uncern10.append(uncern100[i])
    for i in range(8,16):
        uncern30.append(uncern100[i])
    for i in range(16,24):
        uncern50.append(uncern100[i])
    for i in range(24,32):
        uncern70.append(uncern100[i])
    
    # uncern100 plots
    
    # plt.style.use('default')
    # f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    
    # uncern10_1 = ax.scatter(0.25 ,uncern100[0],label='$x_\pi$=(0.20,0.30)')
    # uncern10_2 = ax.scatter(0.35 ,uncern100[1],label='$x_\pi$=(0.30,0.40)')
    # uncern10_3 = ax.scatter(0.45 ,uncern100[2],label='$x_\pi$=(0.40,0.50)')
    # uncern10_4 = ax.scatter(0.55 ,uncern100[3],label='$x_\pi$=(0.50,0.60)')
    # uncern10_5 = ax.scatter(0.65 ,uncern100[4],label='$x_\pi$=(0.60,0.70)')
    # uncern10_6 = ax.scatter(0.75 ,uncern100[5],label='$x_\pi$=(0.70,0.80)')
    # uncern10_7 = ax.scatter(0.85 ,uncern100[6],label='$x_\pi$=(0.80,0.90)')
    # uncern10_8 = ax.scatter(0.95 ,uncern100[7],label='$x_\pi$=(0.90,1.00)')
    # plt.yscale('log')
    # plt.xlabel('xpi')
    # plt.ylabel('$\delta^{100}_{F^{2}_{\pi}}$')
    # plt.title('$\delta^{100}_{F^{2}_{\pi}}$ vs xpi  [$Q^2$ = 10 $GeV^2$]', fontsize =20)
    # leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # leg.get_frame().set_alpha(1.)

    # plt.style.use('default')
    # f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    
    # uncern30_1 = ax.scatter(0.25 ,uncern100[8],label='$x_\pi$=(0.20,0.30)')
    # uncern30_2 = ax.scatter(0.35 ,uncern100[9],label='$x_\pi$=(0.30,0.40)')
    # uncern30_3 = ax.scatter(0.45 ,uncern100[10],label='$x_\pi$=(0.40,0.50)')
    # uncern30_4 = ax.scatter(0.55 ,uncern100[11],label='$x_\pi$=(0.50,0.60)')
    # uncern30_5 = ax.scatter(0.65 ,uncern100[12],label='$x_\pi$=(0.60,0.70)')
    # uncern30_6 = ax.scatter(0.75 ,uncern100[13],label='$x_\pi$=(0.70,0.80)')
    # uncern30_7 = ax.scatter(0.85 ,uncern100[14],label='$x_\pi$=(0.80,0.90)')
    # uncern30_8 = ax.scatter(0.95 ,uncern100[15],label='$x_\pi$=(0.90,1.00)')
    # plt.yscale('log')
    # plt.xlabel('xpi')
    # plt.ylabel('$\delta^{100}_{F^{2}_{\pi}}$')
    # plt.title('$\delta^{100}_{F^{2}_{\pi}}$ vs xpi  [$Q^2$ = 30 $GeV^2$]', fontsize =20)
    # leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # leg.get_frame().set_alpha(1.)

    # plt.style.use('default')
    # f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    
    # uncern50_1 = ax.scatter(0.25 ,uncern100[16],label='$x_\pi$=(0.20,0.30)')
    # uncern50_2 = ax.scatter(0.35 ,uncern100[17],label='$x_\pi$=(0.30,0.40)')
    # uncern50_3 = ax.scatter(0.45 ,uncern100[18],label='$x_\pi$=(0.40,0.50)')
    # uncern50_4 = ax.scatter(0.55 ,uncern100[19],label='$x_\pi$=(0.50,0.60)')
    # uncern50_5 = ax.scatter(0.65 ,uncern100[20],label='$x_\pi$=(0.60,0.70)')
    # uncern50_6 = ax.scatter(0.75 ,uncern100[21],label='$x_\pi$=(0.70,0.80)')
    # uncern50_7 = ax.scatter(0.85 ,uncern100[22],label='$x_\pi$=(0.80,0.90)')
    # uncern50_8 = ax.scatter(0.95 ,uncern100[23],label='$x_\pi$=(0.90,1.00)')
    # plt.yscale('log')
    # plt.xlabel('xpi')
    # plt.ylabel('$\delta^{100}_{F^{2}_{\pi}}$')
    # plt.title('$\delta^{100}_{F^{2}_{\pi}}$ vs xpi  [$Q^2$ = 50 $GeV^2$]', fontsize =20)
    # leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # leg.get_frame().set_alpha(1.)

    # plt.style.use('default')
    # f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    
    # uncern70_1 = ax.scatter(0.25 ,uncern100[24],label='$x_\pi$=(0.20,0.30)')
    # uncern70_2 = ax.scatter(0.35 ,uncern100[25],label='$x_\pi$=(0.30,0.40)')
    # uncern70_3 = ax.scatter(0.45 ,uncern100[26],label='$x_\pi$=(0.40,0.50)')
    # uncern70_4 = ax.scatter(0.55 ,uncern100[27],label='$x_\pi$=(0.50,0.60)')
    # uncern70_5 = ax.scatter(0.65 ,uncern100[28],label='$x_\pi$=(0.60,0.70)')
    # uncern70_6 = ax.scatter(0.75 ,uncern100[29],label='$x_\pi$=(0.70,0.80)')
    # uncern70_7 = ax.scatter(0.85 ,uncern100[30],label='$x_\pi$=(0.80,0.90)')
    # uncern70_8 = ax.scatter(0.95 ,uncern100[31],label='$x_\pi$=(0.90,1.00)')
    # plt.yscale('log')
    # plt.xlabel('xpi')
    # plt.ylabel('$\delta^{100}_{F^{2}_{\pi}}$')
    # plt.title('$\delta^{100}_{F^{2}_{\pi}}$ vs xpi  [$Q^2$ = 70 $GeV^2$]', fontsize =20)
    # leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # leg.get_frame().set_alpha(1.)
                        
    # plt.style.use('default')
    # f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));    
    # xpi10_1 = []
    # fpi10_1 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut10)[0])):
    #     # print clow.applyCuts(xpi_low,cut10)[0][i]
    #     if 0.20 < clow.applyCuts(xpi_low,cut10)[0][i] < 0.30 :
    #         xpi10_1.append(clow.applyCuts(xpi_low,cut10)[0][i])
    #         fpi10_1.append(clow.applyCuts(fpi_low,cut10)[0][i])
    #         # print xpi10_1
    #     i+=1
    # fpiPlot10_1 = ax.errorbar(xpi10_1,fpi10_1,yerr=uncern[0],fmt='.',label='$x_\pi$=(0.20,0.30)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # plt.yscale('log')
    # # plt.xlim(0.20,0.30)
    # # plt.ylim(1e-3,0.025)
    # # plt.xlim(1e-3,1.)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # plt.title('$F^{2}_{\pi}$ vs xpi  [$Q^2$ = 10 $GeV^2$]', fontsize =20)
    # xpi10_2 = []
    # fpi10_2 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut10)[0])):
    #     # print clow.applyCuts(xpi_low,cut10)[0][i]
    #     if 0.30 < clow.applyCuts(xpi_low,cut10)[0][i] < 0.40 :
    #         xpi10_2.append(clow.applyCuts(xpi_low,cut10)[0][i])
    #         fpi10_2.append(clow.applyCuts(fpi_low,cut10)[0][i])
    #         # print xpi10_2
    #     i+=1
    # fpiPlot10_2 = ax.errorbar(xpi10_2,fpi10_2,yerr=uncern[1],fmt='.',label='$x_\pi$=(0.30,0.40)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.30,0.40)
    # xpi10_3 = []
    # fpi10_3 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut10)[0])):
    #     # print clow.applyCuts(xpi_low,cut10)[0][i]
    #     if 0.40 < clow.applyCuts(xpi_low,cut10)[0][i] < 0.50 :
    #         xpi10_3.append(clow.applyCuts(xpi_low,cut10)[0][i])
    #         fpi10_3.append(clow.applyCuts(fpi_low,cut10)[0][i])
    #         # print xpi10_3
    #     i+=1
    # fpiPlot10_3 = ax.errorbar(xpi10_3,fpi10_3,yerr=uncern[2],fmt='.',label='$x_\pi$=(0.40,0.50)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.40,0.50)
    # xpi10_4 = []
    # fpi10_4 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut10)[0])):
    #     # print clow.applyCuts(xpi_low,cut10)[0][i]
    #     if 0.50 < clow.applyCuts(xpi_low,cut10)[0][i] < 0.60 :
    #         xpi10_4.append(clow.applyCuts(xpi_low,cut10)[0][i])
    #         fpi10_4.append(clow.applyCuts(fpi_low,cut10)[0][i])
    #         # print xpi10_4
    #     i+=1
    # fpiPlot10_4 = ax.errorbar(xpi10_4,fpi10_4,yerr=uncern[3],fmt='.',label='$x_\pi$=(0.50,0.60)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.50,0.60)
    # # plt.yscale('log')
    # # plt.xlim(0.60,0.70)
    # xpi10_5 = []
    # fpi10_5 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut10)[0])):
    #     # print clow.applyCuts(xpi_low,cut10)[0][i]
    #     if 0.60 < clow.applyCuts(xpi_low,cut10)[0][i] < 0.70 :
    #         xpi10_5.append(clow.applyCuts(xpi_low,cut10)[0][i])
    #         fpi10_5.append(clow.applyCuts(fpi_low,cut10)[0][i])
    #         # print xpi10_5
    #     i+=1
    # fpiPlot10_5 = ax.errorbar(xpi10_5,fpi10_5,yerr=uncern[4],fmt='.',label='$x_\pi$=(0.60,0.70)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.60,0.70)        
    # xpi10_6 = []
    # fpi10_6 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut10)[0])):
    #     # print clow.applyCuts(xpi_low,cut10)[0][i]
    #     if 0.70 < clow.applyCuts(xpi_low,cut10)[0][i] < 0.80 :
    #         xpi10_6.append(clow.applyCuts(xpi_low,cut10)[0][i])
    #         fpi10_6.append(clow.applyCuts(fpi_low,cut10)[0][i])
    #         # print xpi10_6
    #     i+=1
    # fpiPlot10_6 = ax.errorbar(xpi10_6,fpi10_6,yerr=uncern[5],fmt='.',label='$x_\pi$=(0.70,0.80)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.70,0.80)    
    # xpi10_7 = []
    # fpi10_7 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut10)[0])):
    #     # print clow.applyCuts(xpi_low,cut10)[0][i]
    #     if 0.80 < clow.applyCuts(xpi_low,cut10)[0][i] < 0.90 :
    #         xpi10_7.append(clow.applyCuts(xpi_low,cut10)[0][i])
    #         fpi10_7.append(clow.applyCuts(fpi_low,cut10)[0][i])
    #         # print xpi10_7
    #     i+=1
    # fpiPlot10_7 = ax.errorbar(xpi10_7,fpi10_7,yerr=uncern[6],fmt='.',label='$x_\pi$=(0.80,0.90)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.80,0.90)    
    # xpi10_8 = []
    # fpi10_8 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut10)[0])):
    #     # print clow.applyCuts(xpi_low,cut10)[0][i]
    #     if 0.90 < clow.applyCuts(xpi_low,cut10)[0][i] < 1.00 :
    #         xpi10_8.append(clow.applyCuts(xpi_low,cut10)[0][i])
    #         fpi10_8.append(clow.applyCuts(fpi_low,cut10)[0][i])
    #         # print xpi10_8
    #     i+=1
    # fpiPlot10_8 = ax.errorbar(xpi10_8,fpi10_8,yerr=uncern[7],fmt='.',label='$x_\pi$=(0.90,1.00)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.90,1.00)
    # leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # leg.get_frame().set_alpha(1.)
    # plt.style.use('default')
    # f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    # xpi30_1 = []
    # fpi30_1 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut30)[0])):
    #     # print clow.applyCuts(xpi_low,cut30)[0][i]
    #     if 0.20 < clow.applyCuts(xpi_low,cut30)[0][i] < 0.30 :
    #         xpi30_1.append(clow.applyCuts(xpi_low,cut30)[0][i])
    #         fpi30_1.append(clow.applyCuts(fpi_low,cut30)[0][i])
    #         # print xpi30_1
    #     i+=1
    # fpiPlot30_1 = ax.errorbar(xpi30_1,fpi30_1,yerr=uncern[8],fmt='.',label='$x_\pi$=(0.20,0.30)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # plt.yscale('log')
    # # plt.xlim(0.20,0.30)
    # # plt.ylim(1e-3,0.025)
    # # plt.xlim(1e-3,1.)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # plt.title('$F^{2}_{\pi}$ vs xpi  [$Q^2$ = 30 $GeV^2$]', fontsize =20)
    # xpi30_2 = []
    # fpi30_2 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut30)[0])):
    #     # print clow.applyCuts(xpi_low,cut30)[0][i]
    #     if 0.30 < clow.applyCuts(xpi_low,cut30)[0][i] < 0.40 :
    #         xpi30_2.append(clow.applyCuts(xpi_low,cut30)[0][i])
    #         fpi30_2.append(clow.applyCuts(fpi_low,cut30)[0][i])
    #         # print xpi30_2
    #     i+=1
    # fpiPlot30_2 = ax.errorbar(xpi30_2,fpi30_2,yerr=uncern[9],fmt='.',label='$x_\pi$=(0.30,0.40)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.30,0.40)
    # xpi30_3 = []
    # fpi30_3 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut30)[0])):
    #     # print clow.applyCuts(xpi_low,cut30)[0][i]
    #     if 0.40 < clow.applyCuts(xpi_low,cut30)[0][i] < 0.50 :
    #         xpi30_3.append(clow.applyCuts(xpi_low,cut30)[0][i])
    #         fpi30_3.append(clow.applyCuts(fpi_low,cut30)[0][i])
    #         # print xpi30_3
    #     i+=1
    # fpiPlot30_3 = ax.errorbar(xpi30_3,fpi30_3,yerr=uncern[10],fmt='.',label='$x_\pi$=(0.40,0.50)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.40,0.50)
    # xpi30_4 = []
    # fpi30_4 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut30)[0])):
    #     # print clow.applyCuts(xpi_low,cut30)[0][i]
    #     if 0.50 < clow.applyCuts(xpi_low,cut30)[0][i] < 0.60 :
    #         xpi30_4.append(clow.applyCuts(xpi_low,cut30)[0][i])
    #         fpi30_4.append(clow.applyCuts(fpi_low,cut30)[0][i])
    #         # print xpi30_4
    #     i+=1
    # fpiPlot30_4 = ax.errorbar(xpi30_4,fpi30_4,yerr=uncern[11],fmt='.',label='$x_\pi$=(0.50,0.60)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.50,0.60)
    # # plt.yscale('log')
    # # plt.xlim(0.60,0.70)
    # xpi30_5 = []
    # fpi30_5 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut30)[0])):
    #     # print clow.applyCuts(xpi_low,cut30)[0][i]
    #     if 0.60 < clow.applyCuts(xpi_low,cut30)[0][i] < 0.70 :
    #         xpi30_5.append(clow.applyCuts(xpi_low,cut30)[0][i])
    #         fpi30_5.append(clow.applyCuts(fpi_low,cut30)[0][i])
    #         # print xpi30_5
    #     i+=1
    # fpiPlot30_5 = ax.errorbar(xpi30_5,fpi30_5,yerr=uncern[12],fmt='.',label='$x_\pi$=(0.60,0.70)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.60,0.70)    
    # xpi30_6 = []
    # fpi30_6 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut30)[0])):
    #     # print clow.applyCuts(xpi_low,cut30)[0][i]
    #     if 0.70 < clow.applyCuts(xpi_low,cut30)[0][i] < 0.80 :
    #         xpi30_6.append(clow.applyCuts(xpi_low,cut30)[0][i])
    #         fpi30_6.append(clow.applyCuts(fpi_low,cut30)[0][i])
    #         # print xpi30_6
    #     i+=1
    # fpiPlot30_6 = ax.errorbar(xpi30_6,fpi30_6,yerr=uncern[13],fmt='.',label='$x_\pi$=(0.70,0.80)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.70,0.80)    
    # xpi30_7 = []
    # fpi30_7 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut30)[0])):
    #     # print clow.applyCuts(xpi_low,cut30)[0][i]
    #     if 0.80 < clow.applyCuts(xpi_low,cut30)[0][i] < 0.90 :
    #         xpi30_7.append(clow.applyCuts(xpi_low,cut30)[0][i])
    #         fpi30_7.append(clow.applyCuts(fpi_low,cut30)[0][i])
    #         # print xpi30_7
    #     i+=1
    # fpiPlot30_7 = ax.errorbar(xpi30_7,fpi30_7,yerr=uncern[14],fmt='.',label='$x_\pi$=(0.80,0.90)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.80,0.90)
    # xpi30_8 = []
    # fpi30_8 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut30)[0])):
    #     # print clow.applyCuts(xpi_low,cut30)[0][i]
    #     if 0.90 < clow.applyCuts(xpi_low,cut30)[0][i] < 1.00 :
    #         xpi30_8.append(clow.applyCuts(xpi_low,cut30)[0][i])
    #         fpi30_8.append(clow.applyCuts(fpi_low,cut30)[0][i])
    #         # print xpi30_8
    #     i+=1
    # fpiPlot30_8 = ax.errorbar(xpi30_8,fpi30_8,yerr=uncern[15],fmt='.',label='$x_\pi$=(0.90,1.00)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.90,1.00)
    # leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # leg.get_frame().set_alpha(1.)
    # plt.style.use('default')
    # f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    # xpi50_1 = []
    # fpi50_1 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut50)[0])):
    #     # print clow.applyCuts(xpi_low,cut50)[0][i]
    #     if 0.20 < clow.applyCuts(xpi_low,cut50)[0][i] < 0.30 :
    #         xpi50_1.append(clow.applyCuts(xpi_low,cut50)[0][i])
    #         fpi50_1.append(clow.applyCuts(fpi_low,cut50)[0][i])
    #         # print xpi50_1
    #     i+=1
    # fpiPlot50_1 = ax.errorbar(xpi50_1,fpi50_1,yerr=uncern[16],fmt='.',label='$x_\pi$=(0.20,0.30)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # plt.yscale('log')
    # # plt.xlim(0.20,0.30)
    # # plt.ylim(1e-3,0.025)
    # # plt.xlim(1e-3,1.)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # plt.title('$F^{2}_{\pi}$ vs xpi  [$Q^2$ = 50 $GeV^2$]', fontsize =20)
    # xpi50_2 = []
    # fpi50_2 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut50)[0])):
    #     # print clow.applyCuts(xpi_low,cut50)[0][i]
    #     if 0.30 < clow.applyCuts(xpi_low,cut50)[0][i] < 0.40 :
    #         xpi50_2.append(clow.applyCuts(xpi_low,cut50)[0][i])
    #         fpi50_2.append(clow.applyCuts(fpi_low,cut50)[0][i])
    #         # print xpi50_2
    #     i+=1
    # fpiPlot50_2 = ax.errorbar(xpi50_2,fpi50_2,yerr=uncern[17],fmt='.',label='$x_\pi$=(0.30,0.40)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.30,0.40)
    # xpi50_3 = []
    # fpi50_3 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut50)[0])):
    #     # print clow.applyCuts(xpi_low,cut50)[0][i]
    #     if 0.40 < clow.applyCuts(xpi_low,cut50)[0][i] < 0.50 :
    #         xpi50_3.append(clow.applyCuts(xpi_low,cut50)[0][i])
    #         fpi50_3.append(clow.applyCuts(fpi_low,cut50)[0][i])
    #         # print xpi50_3
    #     i+=1
    # fpiPlot50_3 = ax.errorbar(xpi50_3,fpi50_3,yerr=uncern[18],fmt='.',label='$x_\pi$=(0.40,0.50)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.40,0.50)
    # xpi50_4 = []
    # fpi50_4 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut50)[0])):
    #     # print clow.applyCuts(xpi_low,cut50)[0][i]
    #     if 0.50 < clow.applyCuts(xpi_low,cut50)[0][i] < 0.60 :
    #         xpi50_4.append(clow.applyCuts(xpi_low,cut50)[0][i])
    #         fpi50_4.append(clow.applyCuts(fpi_low,cut50)[0][i])
    #         # print xpi50_4
    #     i+=1
    # fpiPlot50_4 = ax.errorbar(xpi50_4,fpi50_4,yerr=uncern[19],fmt='.',label='$x_\pi$=(0.50,0.60)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.50,0.60)
    # # plt.yscale('log')
    # # plt.xlim(0.60,0.70)
    # xpi50_5 = []
    # fpi50_5 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut50)[0])):
    #     # print clow.applyCuts(xpi_low,cut50)[0][i]
    #     if 0.60 < clow.applyCuts(xpi_low,cut50)[0][i] < 0.70 :
    #         xpi50_5.append(clow.applyCuts(xpi_low,cut50)[0][i])
    #         fpi50_5.append(clow.applyCuts(fpi_low,cut50)[0][i])
    #         # print xpi50_5
    #     i+=1
    # fpiPlot50_5 = ax.errorbar(xpi50_5,fpi50_5,yerr=uncern[20],fmt='.',label='$x_\pi$=(0.60,0.70)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.60,0.70)        
    # xpi50_6 = []
    # fpi50_6 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut50)[0])):
    #     # print clow.applyCuts(xpi_low,cut50)[0][i]
    #     if 0.70 < clow.applyCuts(xpi_low,cut50)[0][i] < 0.80 :
    #         xpi50_6.append(clow.applyCuts(xpi_low,cut50)[0][i])
    #         fpi50_6.append(clow.applyCuts(fpi_low,cut50)[0][i])
    #         # print xpi50_6
    #     i+=1
    # fpiPlot50_6 = ax.errorbar(xpi50_6,fpi50_6,yerr=uncern[21],fmt='.',label='$x_\pi$=(0.70,0.80)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.70,0.80)    
    # xpi50_7 = []
    # fpi50_7 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut50)[0])):
    #     # print clow.applyCuts(xpi_low,cut50)[0][i]
    #     if 0.80 < clow.applyCuts(xpi_low,cut50)[0][i] < 0.90 :
    #         xpi50_7.append(clow.applyCuts(xpi_low,cut50)[0][i])
    #         fpi50_7.append(clow.applyCuts(fpi_low,cut50)[0][i])
    #         # print xpi50_7
    #     i+=1
    # fpiPlot50_7 = ax.errorbar(xpi50_7,fpi50_7,yerr=uncern[22],fmt='.',label='$x_\pi$=(0.80,0.90)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.80,0.90)
    # xpi50_8 = []
    # fpi50_8 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut50)[0])):
    #     # print clow.applyCuts(xpi_low,cut50)[0][i]
    #     if 0.90 < clow.applyCuts(xpi_low,cut50)[0][i] < 1.00 :
    #         xpi50_8.append(clow.applyCuts(xpi_low,cut50)[0][i])
    #         fpi50_8.append(clow.applyCuts(fpi_low,cut50)[0][i])
    #         # print xpi50_8
    #     i+=1
    # fpiPlot50_8 = ax.errorbar(xpi50_8,fpi50_8,yerr=uncern[23],fmt='.',label='$x_\pi$=(0.90,1.00)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.90,1.00)
    # leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # leg.get_frame().set_alpha(1.)
    # plt.style.use('default')
    # f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    # xpi70_1 = []
    # fpi70_1 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut70)[0])):
    #     # print clow.applyCuts(xpi_low,cut70)[0][i]
    #     if 0.20 < clow.applyCuts(xpi_low,cut70)[0][i] < 0.30 :
    #         xpi70_1.append(clow.applyCuts(xpi_low,cut70)[0][i])
    #         fpi70_1.append(clow.applyCuts(fpi_low,cut70)[0][i])
    #         # print xpi70_1
    #     i+=1
    # fpiPlot70_1 = ax.errorbar(xpi70_1,fpi70_1,yerr=uncern[24],fmt='.',label='$x_\pi$=(0.20,0.30)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # plt.yscale('log')
    # # plt.xlim(0.20,0.30)
    # # plt.ylim(1e-3,0.025)
    # # plt.xlim(1e-3,1.)
    # plt.xlabel('xpi')
    # plt.ylabel('$F^{2}_{\pi}$')
    # plt.title('$F^{2}_{\pi}$ vs xpi  [$Q^2$ = 70 $GeV^2$]', fontsize =20)
    # xpi70_2 = []
    # fpi70_2 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut70)[0])):
    #     # print clow.applyCuts(xpi_low,cut70)[0][i]
    #     if 0.30 < clow.applyCuts(xpi_low,cut70)[0][i] < 0.40 :
    #         xpi70_2.append(clow.applyCuts(xpi_low,cut70)[0][i])
    #         fpi70_2.append(clow.applyCuts(fpi_low,cut70)[0][i])
    #         # print xpi70_2
    #     i+=1
    # fpiPlot70_2 = ax.errorbar(xpi70_2,fpi70_2,yerr=uncern[25],fmt='.',label='$x_\pi$=(0.30,0.40)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.30,0.40)
    # xpi70_3 = []
    # fpi70_3 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut70)[0])):
    #     # print clow.applyCuts(xpi_low,cut70)[0][i]
    #     if 0.40 < clow.applyCuts(xpi_low,cut70)[0][i] < 0.50 :
    #         xpi70_3.append(clow.applyCuts(xpi_low,cut70)[0][i])
    #         fpi70_3.append(clow.applyCuts(fpi_low,cut70)[0][i])
    #         # print xpi70_3
    #     i+=1
    # fpiPlot70_3 = ax.errorbar(xpi70_3,fpi70_3,yerr=uncern[26],fmt='.',label='$x_\pi$=(0.40,0.50)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.40,0.50)
    # xpi70_4 = []
    # fpi70_4 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut70)[0])):
    #     # print clow.applyCuts(xpi_low,cut70)[0][i]
    #     if 0.50 < clow.applyCuts(xpi_low,cut70)[0][i] < 0.60 :
    #         xpi70_4.append(clow.applyCuts(xpi_low,cut70)[0][i])
    #         fpi70_4.append(clow.applyCuts(fpi_low,cut70)[0][i])
    #         # print xpi70_4
    #     i+=1
    # fpiPlot70_4 = ax.errorbar(xpi70_4,fpi70_4,yerr=uncern[27],fmt='.',label='$x_\pi$=(0.50,0.60)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.50,0.60)
    # # plt.yscale('log')
    # # plt.xlim(0.60,0.70)
    # xpi70_5 = []
    # fpi70_5 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut70)[0])):
    #     # print clow.applyCuts(xpi_low,cut70)[0][i]
    #     if 0.60 < clow.applyCuts(xpi_low,cut70)[0][i] < 0.70 :
    #         xpi70_5.append(clow.applyCuts(xpi_low,cut70)[0][i])
    #         fpi70_5.append(clow.applyCuts(fpi_low,cut70)[0][i])
    #         # print xpi70_5
    #     i+=1
    # fpiPlot70_5 = ax.errorbar(xpi70_5,fpi70_5,yerr=uncern[28],fmt='.',label='$x_\pi$=(0.60,0.70)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.60,0.70)        
    # xpi70_6 = []
    # fpi70_6 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut70)[0])):
    #     # print clow.applyCuts(xpi_low,cut70)[0][i]
    #     if 0.70 < clow.applyCuts(xpi_low,cut70)[0][i] < 0.80 :
    #         xpi70_6.append(clow.applyCuts(xpi_low,cut70)[0][i])
    #         fpi70_6.append(clow.applyCuts(fpi_low,cut70)[0][i])
    #         # print xpi70_6
    #     i+=1
    # fpiPlot70_6 = ax.errorbar(xpi70_6,fpi70_6,yerr=uncern[29],fmt='.',label='$x_\pi$=(0.70,0.80)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.70,0.80)
    # xpi70_7 = []
    # fpi70_7 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut70)[0])):
    #     # print clow.applyCuts(xpi_low,cut70)[0][i]
    #     if 0.80 < clow.applyCuts(xpi_low,cut70)[0][i] < 0.90 :
    #         xpi70_7.append(clow.applyCuts(xpi_low,cut70)[0][i])
    #         fpi70_7.append(clow.applyCuts(fpi_low,cut70)[0][i])
    #         # print xpi70_7
    #     i+=1
    # fpiPlot70_7 = ax.errorbar(xpi70_7,fpi70_7,yerr=uncern[30],fmt='.',label='$x_\pi$=(0.80,0.90)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.80,0.90)    
    # xpi70_8 = []
    # fpi70_8 = []
    # i=0
    # for i in range(0,len(clow.applyCuts(xpi_low,cut70)[0])):
    #     # print clow.applyCuts(xpi_low,cut70)[0][i]
    #     if 0.90 < clow.applyCuts(xpi_low,cut70)[0][i] < 1.00 :
    #         xpi70_8.append(clow.applyCuts(xpi_low,cut70)[0][i])
    #         fpi70_8.append(clow.applyCuts(fpi_low,cut70)[0][i])
    #         # print xpi70_8
    #     i+=1
    # fpiPlot70_8 = ax.errorbar(xpi70_8,fpi70_8,yerr=uncern[31],fmt='.',label='$x_\pi$=(0.90,1.00)')
    # plt.subplots_adjust(hspace=0.3,wspace=0.45)
    # # plt.yscale('log')
    # # plt.xlim(0.90,1.00)
    # leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    # leg.get_frame().set_alpha(1.)
    
    # Reset figure style, otherwise next plot will look weird
    plt.style.use('default')
    
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    tscat = ax.hist(t_low,bins=p1.setbin(t_low,200,0.,1.0)[0],label='Low',histtype='step', alpha=0.5, stacked=True, fill=True)
    tscat1 = ax.hist(clow.applyCuts(t_low,tcut1),bins=p1.setbin(t_low,200,0.,1.0)[0],label='Low cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.title('t Distribution', fontsize =20)
    plt.xlabel('t')
    plt.ylabel('Number of Events')
    leg = plt.legend(bbox_to_anchor=(0.95,0.3), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.close(f)
    
    for f in xrange(1, plt.figure().number):
        pdf.savefig(f)
    pdf.close()

def sigmavX():
    
    # Cross-section vs xBj
    hTDIS_xbj_tot_sigma_low = densityPlot(TDIS_xbj_low, tot_sigma_low, 'reduced $\sigma_{dis}$[1-8 GeV] vs TDIS_xbj[0.001-0.008]','TDIS_xbj','reduced $\sigma_{dis}$', 20, 20, 0.001, 0.008)
    plt.xscale('log')
    plt.yscale('log')
    hTDIS_xbj_tot_sigma_mid = densityPlot(TDIS_xbj_mid, tot_sigma_mid, 'reduced $\sigma_{dis}$[10-80 GeV] vs TDIS_xbj[0.01-0.08]','TDIS_xbj','reduced $\sigma_{dis}$', 20, 20, 0.01, 0.08)
    plt.xscale('log')
    plt.yscale('log')
    hTDIS_xbj_tot_sigma_high = densityPlot(TDIS_xbj_high, tot_sigma_high, 'reduced $\sigma_{dis}$[100-800 GeV] vs TDIS_xbj[0.1-0.8]','TDIS_xbj','reduced $\sigma_{dis}$', 20, 20, 0.1, 0.8)
    plt.xscale('log')
    plt.yscale('log')

    # Cross-section vs xBj chained
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    hTDIS_xbj_tot_sigma_low = densityPlot(TDIS_xbj_low, tot_sigma_low, '','TDIS_xbj','reduced $\sigma_{dis}$', 20, 20, 0.001, 0.008, figure=f,ax=ax)
    hTDIS_xbj_tot_sigma_mid = densityPlot(TDIS_xbj_mid, tot_sigma_mid, '','TDIS_xbj','reduced $\sigma_{dis}$', 20, 20, 0.01, 0.08, figure=f,ax=ax,layered=False)
    hTDIS_xbj_tot_sigma_high = densityPlot(TDIS_xbj_high, tot_sigma_high, 'reduced $\sigma_{dis}$ vs TDIS_xbj (all sets)','TDIS_xbj','reduced $\sigma_{dis}$', 20, 20, 0.1, 0.8, figure=f,ax=ax,layered=False)
    plt.xscale('log')
    plt.xlim(0.,1.)
    plt.yscale('log')
    plt.ylim(0.,0.15)
    
    # Cross-section vs xpi
    hxpi_tot_sigma_low = densityPlot(xpi_low, tot_sigma_low, 'reduced $\sigma_{dis}$[1-8 GeV] vs $x_{pi}$[0.001-.008]','$x_{pi}$','reduced $\sigma_{dis}$', 20, 20, 0.001, 0.008)
    plt.xscale('log')
    plt.yscale('log')
    hxpi_tot_sigma_mid = densityPlot(xpi_mid, tot_sigma_mid, 'reduced $\sigma_{dis}$[10-80 GeV] vs $x_{pi}$[0.01-0.08]','$x_{pi}$','reduced $\sigma_{dis}$', 20, 20, 0.01, 0.08)
    plt.xscale('log')
    plt.yscale('log')
    hxpi_tot_sigma_high = densityPlot(xpi_high, tot_sigma_high, 'reduced $\sigma_{dis}$[100-800 GeV] vs $x_{pi}$[0.1-0.8]','$x_{pi}$','reduced $\sigma_{dis}$', 20, 20, 0.1, 0.8)
    plt.xscale('log')
    plt.yscale('log')

    # Cross-section vs xpi chained
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    hxpi_tot_sigma_low = densityPlot(xpi_low, tot_sigma_low, '','$x_{pi}$','reduced $\sigma_{dis}$', 20, 20, 0.001, 0.008, figure=f,ax=ax)
    hxpi_tot_sigma_mid = densityPlot(xpi_mid, tot_sigma_mid, '','$x_{pi}$','reduced $\sigma_{dis}$', 20, 20, 0.01, 0.08, figure=f,ax=ax,layered=False)
    hxpi_tot_sigma_high = densityPlot(xpi_high, tot_sigma_high, 'reduced $\sigma_{dis}$ vs $x_{pi}$ (all sets)','$x_{pi}$','reduced $\sigma_{dis}$', 20, 20, 0.1, 0.8, figure=f,ax=ax,layered=False)
    plt.xscale('log')
    plt.xlim(0.,1.)
    plt.yscale('log')
    plt.ylim(0.,0.15)

    for f in xrange(1, plt.figure().number):
        pdf.savefig(f)
    pdf.close()

def yCutPlots():

    [cutQ2bin10,cutQ2bin16,cutQ2bin25,cutQ2bin40,cutQ2bin63,ycut1,tcut1,cut10, cut30, cut50, cut70, cut150, cut200, cut250, cut300, cut350, cut400, cut450, cut500, cut550, cut600, cut650, cut700, cut750, cut800] = q2_Cut()

    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    Q2hist = ax.hist(Q2_low,bins=p1.setbin(Q2_low,200,0.,100.0)[0],histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.xscale('log')
    plt.title('$Q^2$ Distribution', fontsize =20)
    plt.xlabel('$Q^2$')
    plt.ylabel('Number of Events')
    # leg = plt.legend(bbox_to_anchor=(0.95,0.3), loc="center right")
    # leg.get_frame().set_alpha(1.)
    
    phaseSpace10_10 = densityPlot(xpi_low, Q2_low, '$Q^2$ vs $x_\pi$','$x_\pi$','$Q^{2}$', 200, 200, 0., 1.0, 0., 100.)
    # phaseSpace10_10 = ax.scatter(clow.applyCuts(TDIS_xbj_low,tcut1),clow.applyCuts(Q2_low,tcut1))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.,1.)

    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    q2Plot_1 = ax.scatter(xpi_low,Q2_low)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.,1.)
    # plt.ylim(1e-6,0.5)
    plt.xlabel('$x_\pi$')
    plt.ylabel('$Q^2$')
    plt.title('$Q^2$ vs $x_\pi$', fontsize =20)

    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    q2Plot_2 = ax.scatter(TDIS_xbj,Q2_low)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.,1.)
    # plt.ylim(1e-6,0.5)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('$Q^2$')
    plt.title('$Q^2$ vs TDIS_xbj', fontsize =20)

    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    q2Plot_4 = ax.scatter(clow.applyCuts(TDIS_xbj,ycut1),clow.applyCuts(Q2_low,ycut1))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.,1.)
    # plt.ylim(1e-6,0.5)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('$Q^2$')
    plt.title('$Q^2$ vs TDIS_xbj [y > 0.05]', fontsize =20)
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    q2Plot_3 = ax.scatter(clow.applyCuts(xpi_low,ycut1),clow.applyCuts(Q2_low,ycut1))
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.,1.)
    # plt.ylim(1e-6,0.5)
    plt.xlabel('$x_\pi$')
    plt.ylabel('$Q^2$')
    plt.title('$Q^2$ vs $x_\pi$ [y > 0.05]', fontsize =20)

    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    q2Plot_5 = ax.scatter(clow.applyCuts(xpi_low,ycut1),clow.applyCuts(Q2_low,ycut1))
    plt.yscale('log')
    plt.xlim(0.,1.)
    # plt.ylim(1e-6,0.5)
    plt.xlabel('$x_\pi$')
    plt.ylabel('$Q^2$')
    plt.title('$Q^2$ vs $x_\pi$ [y > 0.05]', fontsize =20)

    xpi_ycut = clow.applyCuts(xpi_low,ycut1)[0]
    Q2_ycut = clow.applyCuts(Q2_low,ycut1)[0]
    Q2xpi_ycut = [xpi_ycut,clow.applyCuts(Q2_low,ycut1)[0]]
    
    binxpi = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    binQ2 = [10,20,30,40,50,60,70,80,90,100]
    
    numEvt_xpi = []
    numEvt_Q2 = []
    numEvt = []

    # Step 7 in email, bin only in xpi(0.1)
    for i in range(0,len(binxpi)-1):
        tmp=[]
        for evt in xpi_ycut:
            if binxpi[i] <= evt <= binxpi[i+1]:
                tmp.append(evt)
        numEvt_xpi.append([binxpi[i],len(tmp)])
        # numEvt_xpi.append(len(tmp))
    # print numEvt_xpi

    for i in range(0,len(binQ2)-1):
        tmp=[]
        for evt in Q2_ycut:
            if binQ2[i] <= evt <= binQ2[i+1]:
                tmp.append(evt)
                # print binQ2[i], "<", evt, "<", binQ2[i+1]
        numEvt_Q2.append([binQ2[i],len(tmp)])
        # numEvt_Q2.append(len(tmp))
    # print numEvt_Q2

    fout = open('/home/trottar/ResearchNP/JLEIC/eic_SF/xpiBin.txt','w') 
    
    xpi_tuple = namedtuple('xpi_tuple',['binxpi','nbin'])
    
    xpi_data = []
    tot_xpiData = []
    
    for i in range(0,len(numEvt_xpi)) :
        xpi_data.append(xpi_tuple(numEvt_xpi[i][0],numEvt_xpi[i][1]))
        tot_xpiData.append(xpi_data[i])
        # print tot_data[i]

    for t in tot_xpiData:
        # print namedtuple_to_str(t)
        fout.write("%s\n" % namedtuple_to_str(t))

    fout.close()
    
    # Step 8 in email, bin in xpi and Q2
    for i in range(0,len(binQ2)-1):
        for j in range(0,len(binxpi)-1):
            tmp=[]
            for k in range(0,len(Q2xpi_ycut[0])):
                if binxpi[j] <= Q2xpi_ycut[0][k] <= binxpi[j+1] and binQ2[i] <= Q2xpi_ycut[1][k] <= binQ2[i+1]:
                    tmp.append(Q2xpi_ycut[0][k])
            numEvt.append([binQ2[i],binxpi[j],len(tmp)])
            # numEvt.append(len(tmp))
    print numEvt

    fout = open('/home/trottar/ResearchNP/JLEIC/eic_SF/Q2-xpi_Bin.txt','w') 
    
    Q2xpi_tuple = namedtuple('Q2xpi_tuple',['binQ2','binxpi','nbin'])

    Q2xpi_data = []
    tot_Q2xpiData = []
    
    for i in range(0,len(numEvt)) :
        Q2xpi_data.append(Q2xpi_tuple(numEvt[i][0],numEvt[i][1],numEvt[i][2]))
        tot_Q2xpiData.append(Q2xpi_data[i])
        # print tot_data[i]

    for t in tot_Q2xpiData:
        # print namedtuple_to_str(t)
        fout.write("%s\n" % namedtuple_to_str(t))

    fout.close()

    for f in xrange(1, plt.figure().number):
        pdf.savefig(f)
    pdf.close()

def sigmaPlot():

    [cutQ2bin10,cutQ2bin16,cutQ2bin25,cutQ2bin40,cutQ2bin63,ycut1,tcut1,cut10, cut30, cut50, cut70, cut150, cut200, cut250, cut300, cut350, cut400, cut450, cut500, cut550, cut600, cut650, cut700, cut750, cut800] = q2_Cut()

    binQ2 = [10,16,25,40,63,100]

    numEvt_Q2 = []
    bin_sigma = []
    
    for i in range(0,len(binQ2)-1):
        tmp1=[]
        tmp2=[]
        for evt in Q2_low:
            if binQ2[i] <= evt <= binQ2[i+1]:
                tmp1.append(evt)
                tmp2.append(tot_sigma_low[i])
                # print binQ2[i], "<", evt, "<", binQ2[i+1]
        # numEvt_Q2.append([binQ2[i],len(tmp1)])
        numEvt_Q2.append(len(tmp1))
        bin_sigma.append(sum(tmp2)/numEvt_Q2[i])
    print numEvt_Q2
    print bin_sigma

    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    sigmaQ2 = ax.scatter([10,16,25,40,63],bin_sigma,label='$Q^2$=10 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.title('$d\sigma^{Q2}_{TDIS}$ vs $Q^2$', fontsize =20)
    plt.xlabel('$Q^2$')
    plt.ylabel('$d\sigma^{Q2}_{TDIS}$')
    # ax.text(0.65, 0.25, '$Q^2$=10 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')

    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    numEvtQ2 = ax.scatter([10,16,25,40,63],numEvt_Q2,label='$Q^2$=10 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.title('NumEvts vs $Q^2$', fontsize =20)
    plt.xlabel('$Q^2$')
    plt.ylabel('NumEvts')
    # ax.text(0.65, 0.25, '$Q^2$=10 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')

    binxbj = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

    numEvt_xbj = []
    bin_sigma = []
    
    for i in range(0,len(binxbj)-1):
        tmp1=[]
        tmp2=[]
        for evt in c.applyCuts(TDIS_xbj_low,cutQ2bin10)[0]:
            if binxbj[i] <= evt <= binxbj[i+1]:
                tmp1.append(evt)
                tmp2.append(c.applyCuts(tot_sigma_low,cutQ2bin10)[0][i])
                # print binxbj[i], "<", evt, "<", binxbj[i+1]
        # numEvt_xbj.append([binxbj[i],len(tmp1)])
        numEvt_xbj.append(len(tmp1))
        bin_sigma.append(sum(tmp2)/numEvt_xbj[i])
    print numEvt_xbj
    print bin_sigma
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    sigmaQ2_1 = ax.scatter([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],bin_sigma,label='$Q^2$=10-16 $GeV^2$')
    plt.title('$d\sigma^{Q2,x}_{TDIS}$ vs TDIS_xbj', fontsize =20)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('$d\sigma^{Q2,x}_{TDIS}$')
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)

    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    numEvt_1 = ax.scatter([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],numEvt_xbj,label='$Q^2$=10-16 $GeV^2$')
    plt.title('NumEvts vs TDIS_xbj', fontsize =20)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('NumEvts')
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)

    numEvt_xbj = []
    bin_sigma = []
    
    for i in range(0,len(binxbj)-1):
        tmp1=[]
        tmp2=[]
        for evt in c.applyCuts(TDIS_xbj_low,cutQ2bin16)[0]:
            if binxbj[i] <= evt <= binxbj[i+1]:
                tmp1.append(evt)
                tmp2.append(c.applyCuts(tot_sigma_low,cutQ2bin16)[0][i])
                # print binxbj[i], "<", evt, "<", binxbj[i+1]
        # numEvt_xbj.append([binxbj[i],len(tmp1)])
        numEvt_xbj.append(len(tmp1))
        bin_sigma.append(sum(tmp2)/numEvt_xbj[i])
    print numEvt_xbj
    print bin_sigma
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    sigmaQ2_2 = ax.scatter([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],bin_sigma,label='$Q^2$=16-25 $GeV^2$')
    plt.title('$d\sigma^{Q2,x}_{TDIS}$ vs TDIS_xbj', fontsize =20)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('$d\sigma^{Q2,x}_{TDIS}$')
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)

    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    numEvt_2 = ax.scatter([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],numEvt_xbj,label='$Q^2$=16-25 $GeV^2$')
    plt.title('NumEvts vs TDIS_xbj', fontsize =20)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('NumEvts')
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)

    numEvt_xbj = []
    bin_sigma = []
    
    for i in range(0,len(binxbj)-1):
        tmp1=[]
        tmp2=[]
        for evt in c.applyCuts(TDIS_xbj_low,cutQ2bin25)[0]:
            if binxbj[i] <= evt <= binxbj[i+1]:
                tmp1.append(evt)
                tmp2.append(c.applyCuts(tot_sigma_low,cutQ2bin25)[0][i])
                # print binxbj[i], "<", evt, "<", binxbj[i+1]
        # numEvt_xbj.append([binxbj[i],len(tmp1)])
        numEvt_xbj.append(len(tmp1))
        bin_sigma.append(sum(tmp2)/numEvt_xbj[i])
    print numEvt_xbj
    print bin_sigma
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    sigmaQ2_3 = ax.scatter([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],bin_sigma,label='$Q^2$=25-40 $GeV^2$')
    plt.title('$d\sigma^{Q2,x}_{TDIS}$ vs TDIS_xbj', fontsize =20)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('$d\sigma^{Q2,x}_{TDIS}$')
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)

    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    numEvt_3 = ax.scatter([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],numEvt_xbj,label='$Q^2$=25-40 $GeV^2$')
    plt.title('NumEvts vs TDIS_xbj', fontsize =20)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('NumEvts')
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)

    numEvt_xbj = []
    bin_sigma = []
    
    for i in range(0,len(binxbj)-1):
        tmp1=[]
        tmp2=[]
        for evt in c.applyCuts(TDIS_xbj_low,cutQ2bin40)[0]:
            if binxbj[i] <= evt <= binxbj[i+1]:
                tmp1.append(evt)
                tmp2.append(c.applyCuts(tot_sigma_low,cutQ2bin40)[0][i])
                # print binxbj[i], "<", evt, "<", binxbj[i+1]
        # numEvt_xbj.append([binxbj[i],len(tmp1)])
        numEvt_xbj.append(len(tmp1))
        bin_sigma.append(sum(tmp2)/numEvt_xbj[i])
    print numEvt_xbj
    print bin_sigma
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    sigmaQ2_4 = ax.scatter([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],bin_sigma,label='$Q^2$=40-63 $GeV^2$')
    plt.title('$d\sigma^{Q2,x}_{TDIS}$ vs TDIS_xbj', fontsize =20)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('$d\sigma^{Q2,x}_{TDIS}$')
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)

    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    numEvt_4 = ax.scatter([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],numEvt_xbj,label='$Q^2$=40-63 $GeV^2$')
    plt.title('NumEvts vs TDIS_xbj', fontsize =20)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('NumEvts')
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)

    numEvt_xbj = []
    bin_sigma = []
    
    for i in range(0,len(binxbj)-1):
        tmp1=[]
        tmp2=[]
        for evt in c.applyCuts(TDIS_xbj_low,cutQ2bin63)[0]:
            if binxbj[i] <= evt <= binxbj[i+1]:
                tmp1.append(evt)
                tmp2.append(c.applyCuts(tot_sigma_low,cutQ2bin63)[0][i])
                # print binxbj[i], "<", evt, "<", binxbj[i+1]
        # numEvt_xbj.append([binxbj[i],len(tmp1)])
        numEvt_xbj.append(len(tmp1))
        bin_sigma.append(sum(tmp2)/numEvt_xbj[i])
    print numEvt_xbj
    print bin_sigma
    
    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    sigmaQ2_5 = ax.scatter([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],bin_sigma,label='$Q^2$=63-100 $GeV^2$')
    plt.title('$d\sigma^{Q2,x}_{TDIS}$ vs TDIS_xbj', fontsize =20)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('$d\sigma^{Q2,x}_{TDIS}$')
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)

    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    numEvt_5 = ax.scatter([0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],numEvt_xbj,label='$Q^2$=63-100 $GeV^2$')
    plt.title('NumEvts vs TDIS_xbj', fontsize =20)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('NumEvts')
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)

    binxbj = [0.2,0.3]

    numEvt_xbj = []
    bin_sigma = []
    
    for i in range(0,len(binxbj)-1):
        tmp1=[]
        tmp2=[]
        for evt in c.applyCuts(TDIS_xbj_low,cutQ2bin10)[0]:
            if binxbj[i] <= evt <= binxbj[i+1]:
                tmp1.append(evt)
                tmp2.append(c.applyCuts(tot_sigma_low,cutQ2bin10)[0][i])
                # print binxbj[i], "<", evt, "<", binxbj[i+1]
        # numEvt_xbj.append([binxbj[i],len(tmp1)])
        numEvt_xbj.append(len(tmp1))
        bin_sigma.append(tmp2[0])
    print len(numEvt_xbj), numEvt_xbj
    print len(bin_sigma), bin_sigma

    plt.style.use('default')
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    sigmaQ2_1 = ax.scatter(bin_sigma,numEvt_xbj,label='$Q^2$=10-16 $GeV^2$')
    plt.ylim(0,500)
    plt.title('NumEvts vs $d\sigma_{TDIS}$', fontsize =20)
    plt.xlabel('$d\sigma_{TDIS}$')
    plt.ylabel('NumEvts')
    leg = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    leg.get_frame().set_alpha(1.)

    
    for f in xrange(1, plt.figure().number):
        pdf.savefig(f)
    pdf.close()
    
def main() :

    # sigmavxBj_Plot()
    # sigmavxpi_Plot()
    # phaseSpace()
    # lumi()
    pionPlots()
    # sigmavX()
    # yCutPlots()
    # sigmaPlot()
    plt.show()
    
if __name__=='__main__': main()
