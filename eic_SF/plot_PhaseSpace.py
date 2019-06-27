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

import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
from matplotlib import interactive
from matplotlib import colors
from matplotlib.ticker import MaxNLocator
from sys import path
import time,math,sys
# np.set_printoptions(threshold=sys.maxsize)

# My class function
sys.path.insert(0,'../../../Analysis/ROOTfiles/python')
from test_rootPy import pyPlot, pyCut

rootName = "TDISpion"
rootName_low = "TDISpion_10-100"
rootName_high = "TDISpion_100-800"
# rootName_low = "TDISpion_med"
# rootName_high = "TDISpion_high"

tree1 = "T"

T1_arrkey =  "leafName"
T1_arrhist = "histData"

# pdf = matplotlib.backends.backend_pdf.PdfPages("%s.pdf" % rootName)
pdf = matplotlib.backends.backend_pdf.PdfPages("sigmaPlot.pdf")

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
xpi_high = p2.lookup("xpi")[0]

fpi_low = p1.lookup("fpi")[0]
fpi_high = p2.lookup("fpi")[0]

f2N_low = p1.lookup("f2N")[0]
f2N_high = p2.lookup("f2N")[0]

# sigma_dis_low = p1.lookup("sigma_dis")[0]*(1e-5)
# sigma_dis_high = p2.lookup("sigma_dis")[0]*(1e-5)

sigma_dis_low = p1.lookup("sigma_dis")[0]
sigma_dis_high = p2.lookup("sigma_dis")[0]

# sigma_tdis_low = p1.lookup("sigma_tdis")[0]
# sigma_tdis_high = p2.lookup("sigma_tdis")[0]

sigma_tdis_low = sigma_dis_low*abs(0.003/f2N_low)
sigma_tdis_high = sigma_dis_high*abs(0.002/f2N_high)

# sigma_tdis_low = sigma_dis_low*(fpi_low/f2N_low)
# sigma_tdis_high = sigma_dis_high*(fpi_high/f2N_high)

y_low = p1.lookup("y")[0]
y_high = p2.lookup("y")[0]

yplus_low = 1+((1-y_low)*(1-y_low))
yplus_high = 1+((1-y_high)*(1-y_high))

# tot_sigma_low = (sigma_tdis_low)*((TDIS_xbj_low*(Q2_low*Q2_low)*(137)*(137))/(2*math.pi*yplus_low))
# tot_sigma_high = (sigma_tdis_high)*((TDIS_xbj_high*(Q2_high*Q2_high)*(137)*(137))/(2*math.pi*yplus_high))

tot_sigma_low = (sigma_dis_low)
tot_sigma_high = (sigma_dis_high)

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

Q2array_low =[10.,20.,30.,40.,50.,100.]

lowDict = {}

i=0
for x in range(0,len(Q2array_low)) :
    if i < (len(Q2array_low)-1):
        Q2tmp = '{"Q2cut%i" : ((%0.1f < Q2_low) & (Q2_low < %0.1f))}' % (i,Q2array_low[i]-1.0,Q2array_low[i]+1.0)
        print '{"Q2cut%i" : ((%0.1f < Q2_low) & (Q2_low < %0.1f))}' % (i,Q2array_low[i]-1.0,Q2array_low[i]+1.0)
    else:
        Q2tmp = '{"Q2cut%i" : ((%0.1f < Q2_low) & (Q2_low < %0.1f))}' % (i,Q2array_low[i]-1.0,Q2array_low[i]+1.0)
        print '{"Q2cut%i" : ((%0.1f < Q2_low) & (Q2_low < %0.1f))}' % (i,Q2array_low[i]-1.0,Q2array_low[i]+1.0)
    lowDict.update(eval(Q2tmp))
    i+=1

Q2array_high =[150.,200.,250.,300.,350.,400.,450.,500.,550.,600.,650.,700.,750.,800.]
    
highDict = {}
    
j=0
for x in range(0,len(Q2array_high)) :
    if j < (len(Q2array_high)-1):
        Q2tmp = '{"Q2cut%i" : ((%0.1f < Q2_high) & (Q2_high < %0.1f))}' % (i,Q2array_high[j]-1.0,Q2array_high[j]+1.0)
        print '{"Q2cut%i" : ((%0.1f < Q2_high) & (Q2_high < %0.1f))}' % (i,Q2array_high[j]-1.0,Q2array_high[j]+1.0)
    else:
        Q2tmp = '{"Q2cut%i" : ((%0.1f < Q2_high) & (Q2_high < %0.1f))}' % (i,Q2array_high[j]-1.0,Q2array_high[j]+1.0)
        print '{"Q2cut%i" : ((%0.1f < Q2_high) & (Q2_high < %0.1f))}' % (i,Q2array_high[j]-1.0,Q2array_high[j]+1.0)
    highDict.update(eval(Q2tmp))
    j+=1
    i+=1

clow = pyCut(lowDict)
chigh = pyCut(highDict)

def q2_Cut():

    cut10 = ["Q2cut0"]
    cut20 = ["Q2cut1"]
    cut30 = ["Q2cut2"]
    cut40 = ["Q2cut3"]
    cut50 = ["Q2cut4"]
    cut100 = ["Q2cut5"]
    cut150 = ["Q2cut6"]
    cut200 = ["Q2cut7"]
    cut250 = ["Q2cut8"]
    cut300 = ["Q2cut9"]
    cut350 = ["Q2cut10"]
    cut400 = ["Q2cut11"]
    cut450 = ["Q2cut12"]
    cut500 = ["Q2cut13"]
    cut550 = ["Q2cut14"]
    cut600 = ["Q2cut15"]
    cut650 = ["Q2cut16"]
    cut700 = ["Q2cut17"]
    cut750 = ["Q2cut18"]
    cut800 = ["Q2cut19"]
        
    return[cut10, cut20, cut30, cut40, cut50, cut100, cut150, cut200, cut250, cut300, cut350, cut400, cut450, cut500, cut550, cut600, cut650, cut700, cut750, cut800]
    
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
        plt.colorbar(hist[3], ax=ax, spacing='proportional', label='Number of Events')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

def sigmavxBj_Plot():
    
    [cut10,cut20,cut30,cut40,cut50,cut100,cut150,cut200,cut250,cut300,cut350,cut400,cut450,cut500,cut550,cut600,cut650,cut700,cut750,cut800] = q2_Cut()
    
    f = plt.figure(figsize=(11.69,8.27))
    plt.style.use('classic')
    
    ax = f.add_subplot(451)
    # xBjscat1 = ax.errorbar(clow.applyCuts(TDIS_xbj_low,cut10),clow.applyCuts(tot_sigma_low,cut10),xerr=np.sqrt(clow.applyCuts(TDIS_xbj_low,cut10))/200,yerr=np.sqrt(clow.applyCuts(tot_sigma_low,cut10))/200,fmt='o',label='$Q^2$=10 $GeV^2$')
    xBjscat1 = ax.scatter(clow.applyCuts(TDIS_xbj_low,cut10),clow.applyCuts(tot_sigma_low,cut10),label='$Q^2$=10 $GeV^2$')
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
    xBjscat2 = ax.errorbar(clow.applyCuts(TDIS_xbj_low,cut20),clow.applyCuts(tot_sigma_low,cut20),xerr=np.sqrt(clow.applyCuts(TDIS_xbj_low,cut20))/200,yerr=np.sqrt(clow.applyCuts(tot_sigma_low,cut20))/200,fmt='o',label='$Q^2$=20 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=20 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    ax = f.add_subplot(453)
    xBjscat3 = ax.errorbar(clow.applyCuts(TDIS_xbj_low,cut30),clow.applyCuts(tot_sigma_low,cut30),xerr=np.sqrt(clow.applyCuts(TDIS_xbj_low,cut30))/200,yerr=np.sqrt(clow.applyCuts(tot_sigma_low,cut30))/200,fmt='o',label='$Q^2$=30 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=30 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    plt.title('d$\sigma_{TDIS}$ vs TDIS_xbj', fontsize =20)
    
    ax = f.add_subplot(454)
    xBjscat4 = ax.errorbar(clow.applyCuts(TDIS_xbj_low,cut40),clow.applyCuts(tot_sigma_low,cut40),xerr=np.sqrt(clow.applyCuts(TDIS_xbj_low,cut40))/200,yerr=np.sqrt(clow.applyCuts(tot_sigma_low,cut40))/200,fmt='o',label='$Q^2$=40 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=40 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    
    ax = f.add_subplot(455)
    xBjscat5 = ax.errorbar(clow.applyCuts(TDIS_xbj_low,cut50),clow.applyCuts(tot_sigma_low,cut50),xerr=np.sqrt(clow.applyCuts(TDIS_xbj_low,cut50))/200,yerr=np.sqrt(clow.applyCuts(tot_sigma_low,cut50))/200,fmt='o',label='$Q^2$=50 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.2)
    ax.text(0.65, 0.25, '$Q^2$=50 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(456)
    xBjscat6 = ax.errorbar(clow.applyCuts(TDIS_xbj_low,cut100),clow.applyCuts(tot_sigma_low,cut100),xerr=np.sqrt(clow.applyCuts(TDIS_xbj_low,cut100))/200,yerr=np.sqrt(clow.applyCuts(tot_sigma_low,cut100))/200,fmt='o',label='$Q^2$=100 $GeV^2$')
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
    xBjscat7 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut150),chigh.applyCuts(tot_sigma_high,cut150),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut150))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut150))/200,fmt='o',label='$Q^2$=150 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=150 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(458)
    xBjscat8 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut200),chigh.applyCuts(tot_sigma_high,cut200),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut200))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut200))/200,fmt='o',label='$Q^2$=200 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=200 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(459)
    xBjscat9 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut250),chigh.applyCuts(tot_sigma_high,cut250),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut250))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut250))/200,fmt='o',label='$Q^2$=250 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=250 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(4,5,10)
    xBjscat10 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut300),chigh.applyCuts(tot_sigma_high,cut300),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut300))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut300))/200,fmt='o',label='$Q^2$=300 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,1.0)
    ax.text(0.65, 0.25, '$Q^2$=300 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(4,5,11)
    xBjscat11 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut350),chigh.applyCuts(tot_sigma_high,cut350),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut350))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut350))/200,fmt='o',label='$Q^2$=350 $GeV^2$')
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
    xBjscat12 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut400),chigh.applyCuts(tot_sigma_high,cut400),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut400))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut400))/200,fmt='o',label='$Q^2$=400 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,0.8)
    ax.text(0.65, 0.25, '$Q^2$=400 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(4,5,13)
    xBjscat13 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut450),chigh.applyCuts(tot_sigma_high,cut450),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut450))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut450))/200,fmt='o',label='$Q^2$=450 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,0.8)
    ax.text(0.65, 0.25, '$Q^2$=450 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(4,5,14)
    xBjscat14 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut500),chigh.applyCuts(tot_sigma_high,cut500),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut500))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut500))/200,fmt='o',label='$Q^2$=500 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,0.8)
    ax.text(0.65, 0.25, '$Q^2$=500 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(4,5,15)
    xBjscat15 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut550),chigh.applyCuts(tot_sigma_high,cut550),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut550))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut550))/2000,fmt='o',label='$Q^2$=550 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,0.8)
    ax.text(0.65, 0.25, '$Q^2$=550 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())

    ax = f.add_subplot(4,5,16)
    xBjscat16 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut600),chigh.applyCuts(tot_sigma_high,cut600),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut600))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut600))/200,fmt='o',label='$Q^2$=600 $GeV^2$')
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
    xBjscat17 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut650),chigh.applyCuts(tot_sigma_high,cut650),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut650))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut650))/200,fmt='o',label='$Q^2$=650 $GeV^2$')
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
    xBjscat18 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut700),chigh.applyCuts(tot_sigma_high,cut700),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut700))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut700))/200,fmt='o',label='$Q^2$=700 $GeV^2$')
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
    xBjscat19 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut750),chigh.applyCuts(tot_sigma_high,cut750),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut750))/200,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut750))/200,fmt='o',label='$Q^2$=750 $GeV^2$')
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
    # xBjscat20 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,cut800),chigh.applyCuts(tot_sigma_high,cut800),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,cut800))/800,yerr=np.sqrt(chigh.applyCuts(tot_sigma_high,cut800))/800,fmt='o',label='$Q^2$=800 $GeV^2$')
    xBjscat20 = ax.scatter(chigh.applyCuts(TDIS_xbj_high,cut800),chigh.applyCuts(tot_sigma_high,cut800),label='$Q^2$=800 $GeV^2$')
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    plt.xscale('log')
    plt.xlim(1e-3,1.)
    plt.ylim(0.,0.6)
    ax.text(0.65, 0.25, '$Q^2$=800 $GeV^2$', transform=ax.transAxes, fontsize=10, verticalalignment='top', horizontalalignment='right')
    # ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_locator(MaxNLocator(prune='both'))
    # ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    
    plt.xlabel('TDIS_xbj')
    # plt.tight_layout()
    # plt.setp(ax.get_xticklabels()[0], visible=False)
    # plt.setp(ax.get_xticklabels()[-1], visible=False)
    
    # leg = plt.legend(bbox_to_anchor=(0.95,0.7), loc 
    # leg = plt.legend(bbox_to_anchor=(0.95,0.3), loc="center right")
    # leg.get_frame().set_alpha(1.)
    # plt.xscale('log')
    # plt.ylim(0.,15e-4)
    # plt.yscale('log')
    # print "tot_sigma.add_subplot(6,6,1,totsigscat1.get_offsets()
    # print "tot_sigma Q^2=35.0:\n",totsigscat3.get_offsets()
    # print "tot_sigma Q^2=60.0:\n",totsigscat4.get_offsets()
    # print "tot_sigma Q^2=120.0:\n",totsigscat5.get_offsets()
    # print "tot_sigma Q^2=250.0:\n",totsigscat6.get_offsets()
    
    binValue = xBjscat1.get_offsets()
    evts = len(binValue)

    xval = np.array([])
    sigval = np.array([])
    
    for i in range(0,evts):
        xval = np.append(xval,binValue[i][0])
        sigval = np.append(sigval,binValue[i][1])
        
    binx = max(xval)-min(xval)
    binQ2 = 1.0
    print "---------------------------------"
    print "Events in bin: ",evts, "\nBin value array: ", binValue
    print "\nbinx: ", binx, "\nbinQ2: ", binQ2

    # Reset figure style, otherwise next plot will look weird
    plt.style.use('default')

    avgSigma = np.average(sigval)

    print "\nAverage Sigma: ", avgSigma

    # lumi1 = evts/(clow.applyCuts(tot_sigma_low,cut100)[0]*(binQ2)*(binx))
    lumi1 = evts/(avgSigma*(binQ2)*(binx))

    print "\nLuminosity: ", lumi1*(1e-6), "nb"
    print "---------------------------------\n"
    
    # Luminosity
    # f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    # scat1 = ax.hist(lumi1,bins=p1.setbin(lumi1,200)[0],label='Low',histtype='step', alpha=0.5, stacked=True, fill=True)
    # # scat1 = ax.hist(t_low,bins=p1.setbin(t_low,200,0.,1.0)[0],label='Low',histtype='step', alpha=0.5, stacked=True, fill=True)
    # # plt.title('$f_\pi$ vs TDIS_xbj', fontsize =20)
    # # plt.xlabel('TDIS_xbj')
    # # plt.ylabel('$f_\pi$')
    # leg = plt.legend(bbox_to_anchor=(0.95,0.3), loc="center right")
    # leg.get_frame().set_alpha(1.)
    
    
def sigmavxpi_Plot():
    
    [cut10,cut20,cut30,cut40,cut50,cut100,cut150,cut200,cut250,cut300,cut350,cut400,cut450,cut500,cut550,cut600,cut650,cut700,cut750,cut800] = q2_Cut()
    
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
    
    # leg = plt.legend(bbox_to_anchor=(0.95,0.7), loc 
    # leg = plt.legend(bbox_to_anchor=(0.95,0.3), loc="center right")
    # leg.get_frame().set_alpha(1.)
    # plt.xscale('log')
    # plt.ylim(0.,15e-4)
    # plt.yscale('log')
    # print "tot_sigma.add_subplot(6,6,1,totsigscat1.get_offsets()
    # print "tot_sigma Q^2=15.0:\n",totsigscat2.get_offsets()
    # print "tot_sigma Q^2=35.0:\n",totsigscat3.get_offsets()
    # print "tot_sigma Q^2=60.0:\n",totsigscat4.get_offsets()
    # print "tot_sigma Q^2=120.0:\n",totsigscat5.get_offsets()
    # print "tot_sigma Q^2=250.0:\n",totsigscat6.get_offsets()
    
def phaseSpace():

    plt.style.use('default')
    
    # Q2 vs xBj
    hxbj_Q2_low = densityPlot(TDIS_xbj_low, Q2_low, '$Q^{2}$[10-100 GeV] vs TDIS_xbj[0.01-1.00]','TDIS_xbj','$Q^{2}$', 10, 20, 0.01, 1.00, 10., 100.)
    plt.xscale('log')
    # plt.yscale('log')
    hxbj_Q2_high = densityPlot(TDIS_xbj_high, Q2_high, '$Q^{2}$[100-800 GeV] vs TDIS_xbj[0.1-1.0]','TDIS_xbj','$Q^{2}$', 10, 200, 0.1, 1.0, 100., 1000.)
    plt.xscale('log')
    # plt.yscale('log')

    # Q2 vs xBj chained
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    hxbj_Q2_low = densityPlot(TDIS_xbj_low, Q2_low, '','TDIS_xbj','$Q^{2}$', 10, 20, 0.01, 1.00, 10., 100., figure=f,ax=ax)
    hxbj_Q2_high = densityPlot(TDIS_xbj_high, Q2_high, '$Q^{2}$ vs TDIS_xbj (all sets)','TDIS_xbj','$Q^{2}$', 10, 200, 0.1, 1.0, 100., 1000., figure=f,ax=ax,layered=False)
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
    hxpi_Q2_high = densityPlot(xpi_high, Q2_high, '$Q^{2}$ vs $x_{pi}$ (all sets)','$x_{pi}$','$Q^{2}$', 10, 200, 0.1, 1., 100., 800., figure=f,ax=ax,layered=False)
    plt.xscale('log')
    plt.xlim(0.,1.)
    # plt.yscale('log')
    plt.ylim(0.,800.)

    
    for f in xrange(1, plt.figure().number):
        pdf.savefig(f)
    pdf.close()

def pionPlots():

    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    fpiscat1 = ax.errorbar(TDIS_xbj_low,fpi_low,xerr=np.sqrt(TDIS_xbj_low)/200,yerr=np.sqrt(fpi_low)/200,fmt='o',label='Low')
    fpiscat2 = ax.errorbar(TDIS_xbj_high,fpi_high,xerr=np.sqrt(TDIS_xbj_high)/200,yerr=np.sqrt(fpi_high)/200,fmt='o',label='High')
    plt.title('$f_\pi$ vs TDIS_xbj', fontsize =20)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('$f_\pi$')
    leg = plt.legend(bbox_to_anchor=(0.95,0.3), loc="center right")
    leg.get_frame().set_alpha(1.)

    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    tscat1 = ax.hist(t_low,bins=p1.setbin(t_low,200,0.,1.0)[0],label='Low',histtype='step', alpha=0.5, stacked=True, fill=True)
    tscat2 = ax.hist(t_high,bins=p2.setbin(t_high,200,0.,1.0)[0],label='High',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.title('t Distribution', fontsize =20)
    plt.xlabel('t')
    plt.ylabel('Number of Events')
    leg = plt.legend(bbox_to_anchor=(0.95,0.3), loc="center right")
    leg.get_frame().set_alpha(1.)

    # plt.xscale('log')
    # plt.xlim(1e-3,1.)
    # plt.ylim(0.,0.6)
    
    # for f in xrange(1, plt.figure().number):
    #     pdf.savefig(f)
    # pdf.close()    

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
    
def main() :

    sigmavxBj_Plot()
    sigmavxpi_Plot()
    phaseSpace()
    # pionPlots()
    # sigmavX()
    plt.show()
    
if __name__=='__main__': main()
