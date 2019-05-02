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
from sys import path
import time,math,sys
# np.set_printoptions(threshold=sys.maxsize)

# My class function
sys.path.insert(0,'../../../Analysis/ROOTfiles/python')
from test_rootPy import pyPlot, pyCut

rootName = "TDISpion_80k"
# rootName = "TDISpion_low"
# rootName = "TDISpion_med"
# rootName = "TDISpion_high"
# rootName = "TDISpion"
# rootName = "chain_20k"


tree1 = "T"

T1_arrkey =  "leafName"
T1_arrhist = "histData"

pdf = matplotlib.backends.backend_pdf.PdfPages("%s.pdf" % rootName)

# Arguments for class function
p = pyPlot(rootName,tree1,T1_arrkey,T1_arrhist) 

# Define phyisics data
xBj = p.lookup("xBj")[0]
TDIS_xbj = p.lookup("TDIS_xbj")[0]
sigma_dis = p.lookup("sigma_dis")[0]*(1e-6)
TDIS_y = p.lookup("TDIS_y")[0]
y = p.lookup("y")[0]
Q2 = p.lookup("Q2")[0]
f2N = p.lookup("f2N")[0]
tpi = p.lookup("tpi")[0]
t = -p.lookup("tPrime")[0]
y_D = p.lookup("y_D")[0]
escat = p.lookup("EScatRest")[0]
nu = p.lookup("nu")[0]
TwoPdotk = p.lookup("TwoPdotk")[0]
yplus = 1+((1-y)*(1-y))
alt_sigma_dis = (1e1)*sigma_dis
tot_sigma = (sigma_dis)*((TDIS_xbj*(Q2*Q2)*(137)*(137))/(2*math.pi*yplus))
red_sig  = np.array([1.397,1.258,1.137,1.055,0.944,0.836,0.706,0.519])
red_y = np.array([0.657,0.416,0.263,0.163,0.103,0.066,0.033,0.008])
red_x = np.array([2.53e-4,4.0e-4,6.32e-4,1.02e-3,1.61e-3,2.52e-3,5.0e-3,2.10e-2])
red_Q2 = np.array([15,15,15,15,15,15,15,15])
red_yplus = 1+(1-red_y)*(1-red_y)
my_sigma = (red_sig)*((2*math.pi*red_yplus)/(red_x*(red_Q2*red_Q2)*(137)*(137)))

# arrPlot = arrPlot[(arrCut > low) & (arrCut < high)]
# To call item, cutDict.get(key,"Leaf name not found")
cutDict = {
    "xcut1" : ((1.03e-4 < TDIS_xbj) & (TDIS_xbj < 2.1e-2)),
    "xcut2" : ((2.53e-4 < TDIS_xbj) & (TDIS_xbj < 2.1e-2)),
    "xcut3" : ((8.0e-4 < TDIS_xbj) & (TDIS_xbj < 3.2e-2)),
    "xcut4" : ((1.3e-3 < TDIS_xbj) & (TDIS_xbj < 5.0e-2)),
    "xcut5" : ((2.1e-3 < TDIS_xbj) & (TDIS_xbj < 1.8e-1)),
    "xcut6" : ((5.0e-3 < TDIS_xbj) & (TDIS_xbj < 1.8e-1)),
    "Q2cut1" : ((6.4 < Q2) & (Q2 < 6.6)),
    "Q2cut2" : ((14.9 < Q2) & (Q2 < 15.1)),
    "Q2cut3" : ((34.9 < Q2) & (Q2 < 35.1)),
    "Q2cut4" : ((59.9 < Q2) & (Q2 < 60.1)),
    "Q2cut5" : ((119.9 < Q2) & (Q2 < 120.1)),
    "Q2cut6" : ((249.9 < Q2) & (Q2 < 250.1)),

}

# Please find below the bins to make for Q2, xpi. I also include y here. For this, please run the simulation for 5 GeV electrons and 100 GeV protons. This should give a value for s= 4x 5 x 100 = 2000 GeV^2. Note that s is related to Q2 as s = Q2 x (x_B)/y. We will assume an implicit cut on y < 0.7 for now, which means Q2 x (xB) < 1400. We can relax that later. 

# xarray   =[0.001,0.002,0.002,0.004,0.004,0.004,0.004,0.006,0.006,0.006,0.006,0.008,0.008,0.008,0.008,0.008]
# Q2array  =[1.0,1.0,2.0,1.0,2.0,3.0,4.0,1.5,3.0,4.5,6.0,1.0,2.0,4.0,6.0,8.0]
# yarray   =[0.8,0.4,0.8,0.2,0.4,0.6,0.8,0.2,0.4,0.6,0.8,0.1,0.2,0.4,0.6,0.8]

# xarray   =[0.01,0.01,0.01,0.01,0.01,0.02,0.02,0.02,0.02,0.02,0.04,0.04,0.04,0.04,0.04,0.06,0.06,0.06,0.06,0.06,0.08,0.08,0.08,0.08,0.08]
# Q2array  =[1.25,2.5,5.0,7.5,10.0,2.5,5.0,10.0,15.0,20.0,5.0,10.0,20.0,30.0,40.0,7.5,15.0,30.0,45.0,60.0,10.0,20.0,40.0,60.0,80.0]
# yarray   =[0.1,0.2,0.4,0.6,0.8,0.1,0.2,0.4,0.6,0.8,0.1,0.2,0.4,0.6,0.8,0.1,0.2,0.4,0.6,0.8,0.1,0.2,0.4,0.6,0.8]

xarray   =[0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,0.6,0.6,0.6,0.6,0.8,0.8,0.8,0.8,0.8]
Q2array  =[12.5,25.0,50.0,75.0,100.0,25.0,50.0,100.0,150.0,200.0,50.0,100.0,200.0,300.0,400.0,75.0,150.0,300.0,450.0,600.0,100.0,200.0,400.0,600.0,800.0]
yarray   =[0.1,0.2,0.4,0.6,0.8,0.1,0.2,0.4,0.6,0.8,0.1,0.2,0.4,0.6,0.8,0.1,0.2,0.4,0.6,0.8,0.1,0.2,0.4,0.6,0.8]

cutBinDict = {}

i=0
for x in range(0,len(xarray)) :
    if i < (len(xarray)-1):
        xtmp = '{"xcut%ia" : ((%0.4f  < TDIS_xbj) & (TDIS_xbj < %0.4f))}' % (i,xarray[i]-0.0001,xarray[i]+0.0001)
        ytmp = '{"ycut%ia" : ((%0.1f < y) & (y < %0.1f))}' % (i,yarray[i]-0.1,yarray[i]+0.1)
        Q2tmp = '{"Q2cut%ia" : ((%0.1f < Q2) & (Q2 < %0.1f))}' % (i,Q2array[i]-0.1,Q2array[i]+0.1)
    else:
        xtmp = '{"xcut%ia" : ((%0.4f  < TDIS_xbj) & (TDIS_xbj < %0.4f))}' % (i,xarray[i]-0.0001,xarray[i]+0.0001)
        ytmp = '{"ycut%ia" : ((%0.1f < y) & (y < %0.1f))}' % (i,yarray[i]-0.1,yarray[i]+0.1)
        Q2tmp = '{"Q2cut%ia" : ((%0.1f < Q2) & (Q2 < %0.1f))}' % (i,Q2array[i]-0.1,Q2array[i]+0.1)
    cutBinDict.update(eval(xtmp))
    cutBinDict.update(eval(ytmp))
    cutBinDict.update(eval(Q2tmp))
    i+=1

c = pyCut(cutDict)

cbin = pyCut(cutBinDict)

def sigmaDIS_Cut():
    
    cuts1 = ["xcut1", "Q2cut1"]
    cuts2 = ["xcut2", "Q2cut2"]
    cuts3 = ["xcut3", "Q2cut3"]
    cuts4 = ["xcut4", "Q2cut4"]
    cuts5 = ["xcut5", "Q2cut5"]
    cuts6 = ["xcut6", "Q2cut6"]


    return[cuts1,cuts2,cuts3,cuts4,cuts5,cuts6]

def sigmaBin_Cut():

    cuts1 = ["xcut1a", "ycut1a", "Q2cut1a"]
    cuts2 = ["xcut2a", "ycut2a", "Q2cut2a"]
    cuts3 = ["xcut3a", "ycut3a", "Q2cut3a"]
    cuts4 = ["xcut4a", "ycut4a", "Q2cut4a"]
    cuts5 = ["xcut5a", "ycut5a", "Q2cut5a"]
    cuts6 = ["xcut6a", "ycut6a", "Q2cut6a"]
    cuts7 = ["xcut7a", "ycut7a", "Q2cut7a"]
    cuts8 = ["xcut8a", "ycut8a", "Q2cut8a"]
    cuts9 = ["xcut9a", "ycut9a", "Q2cut9a"]
    cuts10 = ["xcut10a", "ycut10a", "Q2cut10a"]
    cuts11 = ["xcut11a", "ycut11a", "Q2cut11a"]
    cuts12 = ["xcut12a", "ycut12a", "Q2cut12a"]
    cuts13 = ["xcut13a", "ycut13a", "Q2cut13a"]
    cuts14 = ["xcut14a", "ycut14a", "Q2cut14a"]
    cuts15 = ["xcut15a", "ycut15a", "Q2cut15a"]
                                                    

    return[cuts1,cuts2,cuts3,cuts4,cuts5,cuts6,cuts7,cuts8,cuts9,cuts10,cuts11,cuts12,cuts13,cuts14,cuts15]

def densityPlot(x,y,title,xlabel,ylabel,binx,biny,xmin=None,xmax=None,ymin=None,ymax=None,cuts=None,figure=None,sub=None):

    if cuts:
        xcut  = c.applyCuts(x,cuts)[0]
        ycut = c.applyCuts(y,cuts)[0]
    else:
        xcut = x
        ycut = y
    if sub or figure:
        ax = figure.add_subplot(sub)
    else:
        fig, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    if (xmin or xmax or ymin or ymax):
        hist = ax.hist2d(xcut, ycut,bins=(p.setbin(x,binx,xmin,xmax)[0],p.setbin(y,biny,ymin,ymax)[0]), norm=colors.LogNorm())
    else:
        hist = ax.hist2d(xcut, ycut,bins=(p.setbin(x,binx)[0],p.setbin(y,biny)[0]), norm=colors.LogNorm())
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
def sigmaDIS_Plot():

    [cuts1,cuts2,cuts3,cuts4,cuts5,cuts6] = sigmaDIS_Cut()

    # print 'Q2 \n:', Q2
    # print 'TwoPdotk \n:', TwoPdotk
    # print 'TDIS_xbj \n:', TDIS_xbj, np.average(TDIS_xbj)
    # print 'xBj \n:', xBj, np.average(xBj)
    # print 'TDIS_xbj*TwoPdotk \n:', TDIS_xbj*TwoPdotk
    # print 'y \n:', y, np.average(y)
    # print 'y_D \n:', y_D, np.average(y_D)

    # print 'tot_sigma \n:', tot_sigma, np.average(tot_sigma)

    f = plt.figure(figsize=(11.69,8.27))
    f.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.7, wspace=0.35, hspace=0.45)
    ax = f.add_subplot(331)
    hsigma_dis = ax.hist(sigma_dis,bins=p.setbin(sigma_dis,200,0.,10.0)[0],histtype='step', alpha=0.5, stacked=True, fill=True )
    plt.title("d$\sigma_{dis}$")
    ax = f.add_subplot(332)
    totsigplot = ax.hist(tot_sigma,bins=p.setbin(tot_sigma,200,0.,100.0)[0],histtype='step', alpha=0.5, stacked=True, fill=True);
    plt.title('Reduced $\sigma_{dis}$',fontsize =16)
    ax = f.add_subplot(333)
    f2Nplot = ax.hist(f2N,bins=p.setbin(f2N,200,0.,100.0)[0],histtype='step', alpha=0.5, stacked=True, fill=True);
    plt.title('$f^{2}_{N}$',fontsize =16)
    ax = f.add_subplot(334)
    hxbj = ax.hist(TDIS_xbj,bins=p.setbin(TDIS_xbj,200,0.,2.0)[0],histtype='step', alpha=0.5, stacked=True, fill=True )
    plt.title("TDIS_xbj")
    ax = f.add_subplot(335)
    hy = ax.hist(y,bins=p.setbin(y,200,0.,2.0)[0],histtype='step', alpha=0.5, stacked=True, fill=True )
    plt.title("y")
    ax = f.add_subplot(336)
    h1mmiss = ax.hist(yplus,bins=p.setbin(yplus,200,0.002,2.0)[0],histtype='step', alpha=0.5, stacked=True, fill=True )
    plt.title("$Y_+$")

    f = plt.figure(figsize=(11.69,8.27))
    f.tight_layout()
    hQ2_sigma_dis = densityPlot(Q2, sigma_dis, 'd$\sigma_{dis}$ vs $Q^2$','$Q^2$','d$\sigma_{dis}$', 200, 200, 0., 350., 0., 1e-4, figure=f, sub=131)
    hxbj_sigma_dis = densityPlot(TDIS_xbj, sigma_dis, 'd$\sigma_{dis}$ vs TDIS_xbj','TDIS_xbj','d$\sigma_{dis}$', 200, 200, 0., 1., 0., 1e-4, figure=f, sub=132)
    hy_sigma_dis = densityPlot(y, sigma_dis, 'd$\sigma_{dis}$ vs y','y','d$\sigma_{dis}$', 200, 200, 0., 1., 0., 1e-4, figure=f, sub=133)
    
    hxbj_Q2 = densityPlot(TDIS_xbj, Q2, '$Q^{2}$ vs TDIS_xbj','TDIS_xbj','$Q^{2}$', 200, 200, 0., 1., 0., 350.)

    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    mysigscat = ax.scatter(red_x,my_sigma,label='Zeus')
    sigscat2 = ax.scatter(c.applyCuts(TDIS_xbj,cuts2)[0],c.applyCuts(alt_sigma_dis,cuts2)[0],label='$Q^2$=15.0 $GeV^2$')
    plt.title('d$\sigma_{DIS}$ vs TDIS_xbj', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.2,0.3), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('d$\sigma_{DIS}$ (fb)')
    plt.xlim(1e-4,1)
    plt.xscale('log')
    plt.yscale('log')
    
    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    sigscat1 = ax.scatter(c.applyCuts(TDIS_xbj,cuts1),c.applyCuts(sigma_dis,cuts1),label='$Q^2$=6.5 $GeV^2$')
    sigscat2 = ax.scatter(c.applyCuts(TDIS_xbj,cuts2),c.applyCuts(sigma_dis,cuts2),label='$Q^2$=15.0 $GeV^2$')
    sigscat3 = ax.scatter(c.applyCuts(TDIS_xbj,cuts3),c.applyCuts(sigma_dis,cuts3),label='$Q^2$=35.0 $GeV^2$')
    sigscat4 = ax.scatter(c.applyCuts(TDIS_xbj,cuts4),c.applyCuts(sigma_dis,cuts4),label='$Q^2$=60.0 $GeV^2$')
    sigscat5 = ax.scatter(c.applyCuts(TDIS_xbj,cuts5),c.applyCuts(sigma_dis,cuts5),label='$Q^2$=120.0 $GeV^2$')
    sigscat6 = ax.scatter(c.applyCuts(TDIS_xbj,cuts6),c.applyCuts(sigma_dis,cuts6),label='$Q^2$=250.0 $GeV^2$')
    plt.title('d$\sigma_{DIS}$ vs TDIS_xbj', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.2,0.3), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('d$\sigma_{DIS}$ (nb)')
    plt.xlim(1e-4,1)
    plt.xscale('log')
    # plt.yscale('log')
    print "sigma_dis Q^2=6.5:\n",sigscat1.get_offsets()
    print "sigma_dis Q^2=15.0:\n",sigscat2.get_offsets()
    print "sigma_dis Q^2=35.0:\n",sigscat3.get_offsets()
    print "sigma_dis Q^2=60.0:\n",sigscat4.get_offsets()
    print "sigma_dis Q^2=120.0:\n",sigscat5.get_offsets()
    print "sigma_dis Q^2=250.0:\n",sigscat6.get_offsets()

    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    f2Nscat1 = ax.scatter(c.applyCuts(TDIS_xbj,cuts1),c.applyCuts(f2N,cuts1),label='$Q^2$=6.5 $GeV^2$')
    f2Nscat2 = ax.scatter(c.applyCuts(TDIS_xbj,cuts2),c.applyCuts(f2N,cuts2),label='$Q^2$=15.0 $GeV^2$')
    f2Nscat3 = ax.scatter(c.applyCuts(TDIS_xbj,cuts3),c.applyCuts(f2N,cuts3),label='$Q^2$=35.0 $GeV^2$')
    f2Nscat4 = ax.scatter(c.applyCuts(TDIS_xbj,cuts4),c.applyCuts(f2N,cuts4),label='$Q^2$=60.0 $GeV^2$')
    f2Nscat5 = ax.scatter(c.applyCuts(TDIS_xbj,cuts5),c.applyCuts(f2N,cuts5),label='$Q^2$=120.0 $GeV^2$')
    f2Nscat6 = ax.scatter(c.applyCuts(TDIS_xbj,cuts6),c.applyCuts(f2N,cuts6),label='$Q^2$=250.0 $GeV^2$')
    plt.title('$F^{N}_{2}$ vs TDIS_xbj', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.2,0.3), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('$F^{N}_{2}$')
    plt.xlim(1e-4,1)
    plt.xscale('log')
    # plt.ylim(0.2,0.6)
    # plt.yscale('log')
    print "f2N Q^2=6.5:\n",f2Nscat1.get_offsets()
    print "f2N Q^2=15.0:\n",f2Nscat2.get_offsets()
    print "f2N Q^2=35.0:\n",f2Nscat3.get_offsets()
    print "f2N Q^2=60.0:\n",f2Nscat4.get_offsets()
    print "f2N Q^2=120.0:\n",f2Nscat5.get_offsets()
    print "f2N Q^2=250.0:\n",f2Nscat6.get_offsets()

    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    # totsigscat = ax.scatter(TDIS_xbj,tot_sigma,label='No cuts')
    totsigscat1 = ax.scatter(c.applyCuts(TDIS_xbj,cuts1),c.applyCuts(tot_sigma,cuts1),label='$Q^2$=6.5 $GeV^2$')
    totsigscat2 = ax.scatter(c.applyCuts(TDIS_xbj,cuts2),c.applyCuts(tot_sigma,cuts2),label='$Q^2$=15.0 $GeV^2$')
    totsigscat3 = ax.scatter(c.applyCuts(TDIS_xbj,cuts3),c.applyCuts(tot_sigma,cuts3),label='$Q^2$=35.0 $GeV^2$')
    totsigscat4 = ax.scatter(c.applyCuts(TDIS_xbj,cuts4),c.applyCuts(tot_sigma,cuts4),label='$Q^2$=60.0 $GeV^2$')
    totsigscat5 = ax.scatter(c.applyCuts(TDIS_xbj,cuts5),c.applyCuts(tot_sigma,cuts5),label='$Q^2$=120.0 $GeV^2$')
    totsigscat6 = ax.scatter(c.applyCuts(TDIS_xbj,cuts6),c.applyCuts(tot_sigma,cuts6),label='$Q^2$=250.0 $GeV^2$')
    plt.title('Reduced $\sigma_{DIS}$ vs TDIS_xbj', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.2,0.3), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('Reduced $\sigma_{DIS}$ ($nb^{-1}$)')
    plt.xlim(1e-4,1)
    plt.xscale('log')
    # plt.ylim(0.,15e-4)
    # plt.yscale('log')
    print "tot_sigma Q^2=6.5:\n",totsigscat1.get_offsets()
    print "tot_sigma Q^2=15.0:\n",totsigscat2.get_offsets()
    print "tot_sigma Q^2=35.0:\n",totsigscat3.get_offsets()
    print "tot_sigma Q^2=60.0:\n",totsigscat4.get_offsets()
    print "tot_sigma Q^2=120.0:\n",totsigscat5.get_offsets()
    print "tot_sigma Q^2=250.0:\n",totsigscat6.get_offsets()
    
    # f, ax = plt.subplots()
    # thist1 = ax.hist(c.applyCuts(t,cuts1),bins=p.setbin(t,100,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$Q^2$=6.5 $GeV^2$' )
    # thist2 = ax.hist(c.applyCuts(t,cuts2),bins=p.setbin(t,100,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$Q^2$=15.0 $GeV^2$' )
    # thist3 = ax.hist(c.applyCuts(t,cuts3),bins=p.setbin(t,100,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$Q^2$=35.0 $GeV^2$' )
    # plt.title("t cuts", fontsize =16)
    # leg = plt.legend(bbox_to_anchor=(0.6,0.5), loc="center right")
    # leg.get_frame().set_alpha(1.)
    
    for f in xrange(1, plt.figure().number):
        pdf.savefig(f)
    pdf.close()

def sigmaBin_Plot():

    [cuts1,cuts2,cuts3,cuts4,cuts5,cuts6,cuts7,cuts8,cuts9,cuts10,cuts11,cuts12,cuts13,cuts14,cuts15] = sigmaBin_Cut()

    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    sigbin1 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts1)[0],cbin.applyCuts(sigma_dis,cuts1)[0],label='cut1')
    sigbin2 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts2)[0],cbin.applyCuts(sigma_dis,cuts2)[0],label='cut2')
    sigbin3 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts3)[0],cbin.applyCuts(sigma_dis,cuts3)[0],label='cut3')
    sigbin4 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts4)[0],cbin.applyCuts(sigma_dis,cuts4)[0],label='cut4')
    sigbin5 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts5)[0],cbin.applyCuts(sigma_dis,cuts5)[0],label='cut5')
    sigbin6 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts6)[0],cbin.applyCuts(sigma_dis,cuts6)[0],label='cut6')
    sigbin7 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts7)[0],cbin.applyCuts(sigma_dis,cuts7)[0],label='cut7')
    sigbin8 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts8)[0],cbin.applyCuts(sigma_dis,cuts8)[0],label='cut8')
    sigbin9 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts9)[0],cbin.applyCuts(sigma_dis,cuts9)[0],label='cut9')
    sigbin10 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts10)[0],cbin.applyCuts(sigma_dis,cuts10)[0],label='cut10')
    sigbin11 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts11)[0],cbin.applyCuts(sigma_dis,cuts11)[0],label='cut11')
    sigbin12 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts12)[0],cbin.applyCuts(sigma_dis,cuts12)[0],label='cut12')
    sigbin13 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts13)[0],cbin.applyCuts(sigma_dis,cuts13)[0],label='cut13')
    sigbin14 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts14)[0],cbin.applyCuts(sigma_dis,cuts14)[0],label='cut14')
    sigbin15 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts15)[0],cbin.applyCuts(sigma_dis,cuts15)[0],label='cut15')
    plt.title('d$\sigma_{DIS}$ vs TDIS_xbj', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.2,0.3), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('d$\sigma_{DIS}$ (fb)')
    plt.xlim(1e-4,1)
    plt.xscale('log')
    plt.yscale('log')
    print "sigma_dis cut1:\n",sigbin1.get_offsets()
    print "sigma_dis cut2:\n",sigbin2.get_offsets()
    print "sigma_dis cut3:\n",sigbin3.get_offsets()
    print "sigma_dis cut4:\n",sigbin4.get_offsets()
    print "sigma_dis cut5:\n",sigbin5.get_offsets()
    print "sigma_dis cut6:\n",sigbin6.get_offsets()
    print "sigma_dis cut7:\n",sigbin7.get_offsets()
    print "sigma_dis cut8:\n",sigbin8.get_offsets()
    print "sigma_dis cut9:\n",sigbin9.get_offsets()
    print "sigma_dis cut10:\n",sigbin10.get_offsets()
    print "sigma_dis cut11:\n",sigbin11.get_offsets()
    print "sigma_dis cut12:\n",sigbin12.get_offsets()
    print "sigma_dis cut13:\n",sigbin13.get_offsets()
    print "sigma_dis cut14:\n",sigbin14.get_offsets()
    print "sigma_dis cut15:\n",sigbin15.get_offsets()
    
    for f in xrange(1, plt.figure().number):
        pdf.savefig(f)
    pdf.close()    
    
def main() :

    sigmaBin_Plot()
    # sigmaDIS_Plot()
    # recreateLeaves()
    plt.show()
    
if __name__=='__main__': main()
