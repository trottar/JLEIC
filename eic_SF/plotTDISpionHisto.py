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
from test_rootPy import pyPlot

rootName = "TDISpion_80k"
# rootName = "TDISpion"

tree1 = "T"

T1_arrkey =  "leafName"
T1_arrhist = "histData"

# Arguments for class function
p = pyPlot(rootName,tree1,T1_arrkey,T1_arrhist) 

pdf = matplotlib.backends.backend_pdf.PdfPages("%s.pdf" % rootName)

def cut(key):

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
        "Q2cut6" : ((249.9 < Q2) & (Q2 < 250.1))

    }
    
    return[cutDict.get(key,"Leaf name not found")]

def applyCuts(leaf,cuts=None):

    if cuts:
        tmp = leaf
        tmp = eval(p.applyCuts(leaf,cuts)[0])
    else:
        tmp = leaf
    
    return[tmp]
    

def selectCut():
    
    cuts1 = ["xcut1", "Q2cut1"]
    cuts2 = ["xcut2", "Q2cut2"]
    cuts3 = ["xcut3", "Q2cut3"]
    cuts4 = ["xcut4", "Q2cut4"]
    cuts5 = ["xcut5", "Q2cut5"]
    cuts6 = ["xcut6", "Q2cut6"]


    return[cuts1,cuts2,cuts3,cuts4,cuts5,cuts6]

def densityPlot(x,y,title,xlabel,ylabel,binx,biny,xmin=None,xmax=None,ymin=None,ymax=None,cuts=None,figure=None,sub=None):

    if cuts:
        xcut  = applyCuts(x,cuts)[0]
        ycut = applyCuts(y,cuts)[0]
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

# Define phyisics data
xBj = p.lookup("xBj")[0]
TDIS_xbj = p.lookup("TDIS_xbj")[0]
sigma_dis = p.lookup("sigma_dis")[0]
TDIS_y = p.lookup("TDIS_y")[0]
y = p.lookup("y")[0]
Q2 = p.lookup("Q2")[0]
f2N = p.lookup("f2N")[0]
tpi = p.lookup("tpi")[0]
t = -p.lookup("tPrime")[0]
y_D = p.lookup("y_D")[0]
# Yplus = p.lookup("Yplus")[0]
# sigma_dis = (Q2/y)*p.lookup("sigma_dis")[0]
# invts.y_D  = invts.Q2/(invts.xBj*invts.TwoPdotk);
TwoPdotk = p.lookup("TwoPdotk")[0]
# y = Q2/(TDIS_xbj*TwoPdotk)
yplus = 1+((1-y)*(1-y))
tot_sigma = (sigma_dis)*((TDIS_xbj*(Q2*Q2)*(137)*(137))/(2*math.pi*yplus))
# tot_sigma = (f2N)*((TDIS_xbj*(Q2*Q2)*(137)*(137))/(2*math.pi*yplus))

def customPlots():

    [cuts1,cuts2,cuts3,cuts4,cuts5,cuts6] = selectCut()

    # print 'Q2 \n:', Q2
    # print 'TwoPdotk \n:', TwoPdotk
    # print 'TDIS_xbj \n:', TDIS_xbj, np.average(TDIS_xbj)
    # print 'xBj \n:', xBj, np.average(xBj)
    # print 'TDIS_xbj*TwoPdotk \n:', TDIS_xbj*TwoPdotk
    # print 'y \n:', y, np.average(y)
    # print 'y_D \n:', y_D, np.average(y_D)

    # print 'tot_sigma \n:', tot_sigma, np.average(tot_sigma)

    f = plt.figure(figsize=(11.69,8.27))
    f.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.7, wspace=0.3, hspace=0.3)
    ax = f.add_subplot(331)
    hsigma_dis = ax.hist(sigma_dis,bins=p.setbin(sigma_dis,200,0.,10.0)[0],histtype='step', alpha=0.5, stacked=True, fill=True )
    plt.title("d$\sigma_{dis}$")
    ax = f.add_subplot(332)
    totsigplot = ax.hist(tot_sigma,bins=p.setbin(tot_sigma,200,0.,100.0)[0],histtype='step', alpha=0.5, stacked=True, fill=True);
    plt.title('Total $\sigma_{dis}$',fontsize =16)
    ax = f.add_subplot(333)
    totsigplot = ax.hist(f2N,bins=p.setbin(f2N,200,0.,100.0)[0],histtype='step', alpha=0.5, stacked=True, fill=True);
    plt.title('Total $\sigma_{dis}$',fontsize =16)
    ax = f.add_subplot(334)
    hxbj = ax.hist(TDIS_xbj,bins=p.setbin(TDIS_xbj,200,0.,2.0)[0],histtype='step', alpha=0.5, stacked=True, fill=True )
    plt.title("TDIS_xbj")
    ax = f.add_subplot(335)
    # hy_D = ax.hist(y_D,bins=p.setbin(y_D,200,0.,2.0)[0],histtype='step', stacked=True, fill=True )
    hy = ax.hist(y,bins=p.setbin(y,200,0.,2.0)[0],histtype='step', alpha=0.5, stacked=True, fill=True )
    plt.title("y")
    ax = f.add_subplot(336)
    h1mmiss = ax.hist(yplus,bins=p.setbin(yplus,200,0.002,2.0)[0],histtype='step', alpha=0.5, stacked=True, fill=True )
    plt.title("$Y_+$")

    f = plt.figure(figsize=(11.69,8.27))
    f.tight_layout()
    hQ2_sigma_dis = densityPlot(Q2, sigma_dis, 'd$\sigma_{dis}$ vs $Q^2$','$Q^2$','d$\sigma_{dis}$', 200, 200, 0., 350., 0., 10., figure=f, sub=131)
    hxbj_sigma_dis = densityPlot(TDIS_xbj, sigma_dis, 'd$\sigma_{dis}$ vs TDIS_xbj','TDIS_xbj','d$\sigma_{dis}$', 200, 200, 0., 1., 0., 10., figure=f, sub=132)
    hy_sigma_dis = densityPlot(y, sigma_dis, 'd$\sigma_{dis}$ vs y','y','d$\sigma_{dis}$', 200, 200, 0., 1., 0., 10., figure=f, sub=133)

    hxbj_Q2 = densityPlot(TDIS_xbj, Q2, '$Q^{2}$ vs TDIS_xbj','TDIS_xbj','$Q^{2}$', 200, 200, 0., 1., 0., 350.)
    
    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    sigscat1 = ax.scatter(applyCuts(TDIS_xbj,cuts1),applyCuts(sigma_dis,cuts1),label='$Q^2$=6.5 $GeV^2$')
    sigscat2 = ax.scatter(applyCuts(TDIS_xbj,cuts2),applyCuts(sigma_dis,cuts2),label='$Q^2$=15.0 $GeV^2$')
    sigscat3 = ax.scatter(applyCuts(TDIS_xbj,cuts3),applyCuts(sigma_dis,cuts3),label='$Q^2$=35.0 $GeV^2$')
    sigscat4 = ax.scatter(applyCuts(TDIS_xbj,cuts4),applyCuts(sigma_dis,cuts4),label='$Q^2$=60.0 $GeV^2$')
    sigscat5 = ax.scatter(applyCuts(TDIS_xbj,cuts5),applyCuts(sigma_dis,cuts5),label='$Q^2$=120.0 $GeV^2$')
    sigscat6 = ax.scatter(applyCuts(TDIS_xbj,cuts6),applyCuts(sigma_dis,cuts6),label='$Q^2$=250.0 $GeV^2$')
    plt.title('d$\sigma_{DIS}$ vs TDIS_xbj', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.4,0.5), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('d$\sigma_{DIS}$ (nb)')
    plt.xlim(1e-4,1)
    plt.xscale('log')
    # plt.ylim(0.,10)
    # plt.yscale('log')
    print "sigma_dis Q^2=6.5:\n",sigscat1.get_offsets()
    print "sigma_dis Q^2=15.0:\n",sigscat2.get_offsets()
    print "sigma_dis Q^2=35.0:\n",sigscat3.get_offsets()
    print "sigma_dis Q^2=60.0:\n",sigscat4.get_offsets()
    print "sigma_dis Q^2=120.0:\n",sigscat5.get_offsets()
    print "sigma_dis Q^2=250.0:\n",sigscat6.get_offsets()

    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    f2Nscat1 = ax.scatter(applyCuts(TDIS_xbj,cuts1),applyCuts(f2N,cuts1),label='$Q^2$=6.5 $GeV^2$')
    f2Nscat2 = ax.scatter(applyCuts(TDIS_xbj,cuts2),applyCuts(f2N,cuts2),label='$Q^2$=15.0 $GeV^2$')
    f2Nscat3 = ax.scatter(applyCuts(TDIS_xbj,cuts3),applyCuts(f2N,cuts3),label='$Q^2$=35.0 $GeV^2$')
    f2Nscat4 = ax.scatter(applyCuts(TDIS_xbj,cuts4),applyCuts(f2N,cuts4),label='$Q^2$=60.0 $GeV^2$')
    f2Nscat5 = ax.scatter(applyCuts(TDIS_xbj,cuts5),applyCuts(f2N,cuts5),label='$Q^2$=120.0 $GeV^2$')
    f2Nscat6 = ax.scatter(applyCuts(TDIS_xbj,cuts6),applyCuts(f2N,cuts6),label='$Q^2$=250.0 $GeV^2$')
    plt.title('$F^{N}_{2}$ vs TDIS_xbj', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.4,0.5), loc="center right")
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
    totsigscat1 = ax.scatter(applyCuts(TDIS_xbj,cuts1),applyCuts(tot_sigma,cuts1),label='$Q^2$=6.5 $GeV^2$')
    totsigscat2 = ax.scatter(applyCuts(TDIS_xbj,cuts2),applyCuts(tot_sigma,cuts2),label='$Q^2$=15.0 $GeV^2$')
    totsigscat3 = ax.scatter(applyCuts(TDIS_xbj,cuts3),applyCuts(tot_sigma,cuts3),label='$Q^2$=35.0 $GeV^2$')
    totsigscat4 = ax.scatter(applyCuts(TDIS_xbj,cuts4),applyCuts(tot_sigma,cuts4),label='$Q^2$=60.0 $GeV^2$')
    totsigscat5 = ax.scatter(applyCuts(TDIS_xbj,cuts5),applyCuts(tot_sigma,cuts5),label='$Q^2$=120.0 $GeV^2$')
    totsigscat6 = ax.scatter(applyCuts(TDIS_xbj,cuts6),applyCuts(tot_sigma,cuts6),label='$Q^2$=250.0 $GeV^2$')
    plt.title('Total $\sigma_{DIS}$ vs TDIS_xbj', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.4,0.5), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('Total $\sigma_{DIS}$ ($nb^{-1}$)')
    plt.xlim(1e-4,1)
    plt.xscale('log')
    # plt.ylim(0.,10)
    # plt.yscale('log')
    print "tot_sigma Q^2=6.5:\n",totsigscat1.get_offsets()
    print "tot_sigma Q^2=15.0:\n",totsigscat2.get_offsets()
    print "tot_sigma Q^2=35.0:\n",totsigscat3.get_offsets()
    print "tot_sigma Q^2=60.0:\n",totsigscat4.get_offsets()
    print "tot_sigma Q^2=120.0:\n",totsigscat5.get_offsets()
    print "tot_sigma Q^2=250.0:\n",totsigscat6.get_offsets()
    
    # f, ax = plt.subplots()
    # thist1 = ax.hist(applyCuts(t,cuts1),bins=p.setbin(t,100,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$Q^2$=6.5 $GeV^2$' )
    # thist2 = ax.hist(applyCuts(t,cuts2),bins=p.setbin(t,100,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$Q^2$=15.0 $GeV^2$' )
    # thist3 = ax.hist(applyCuts(t,cuts3),bins=p.setbin(t,100,0.,1.)[0],histtype='step', alpha=0.5, stacked=True, fill=True,label='$Q^2$=35.0 $GeV^2$' )
    # plt.title("t cuts", fontsize =16)
    # leg = plt.legend(bbox_to_anchor=(0.6,0.5), loc="center right")
    # leg.get_frame().set_alpha(1.)
    
    for f in xrange(1, plt.figure().number):
        pdf.savefig(f)
    pdf.close()

def main() :

    customPlots()
    # recreateLeaves()
    plt.show()
    
if __name__=='__main__': main()
