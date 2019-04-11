#! /usr/bin/python

#
# Description:This will read in the array data file that contains all the leave histogram information
# ================================================================
# Time-stamp: "2019-04-11 16:13:54 trottar"
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

import matplotlib.pyplot as plt
from matplotlib import interactive
from matplotlib import colors
import numpy as np
import math
import sys

# rootName = "TDISpion_80k"
rootName = "TDISpion"

tree1 = "Evnts"

T1_arrkey =  "leafName"
T1_arrhist = "histData"

# Retrieves the array data file and creates new leaves from this
def pullArray():
    
    data = np.load("%s.npz" % rootName)

    T1 = data[T1_arrkey]
    T1_hist = data[T1_arrhist]

    return[T1,T1_hist]

# Creates a dictionary that stores leaf names with the corresponding array
def dictionary():

    [T1,T1_hist] = pullArray()
    
    T1_leafdict = dict(zip(T1, T1_hist))

    return[T1_leafdict]

# A quick way to look up a leaf name for its array
def lookup(key):

    [T1_leafdict] = dictionary()
    
    T_leafFound = T1_leafdict.get(key,"Leaf name not found")

    return[T_leafFound]

# Recreates the histograms of the root file
def recreateLeaves():

    [T1_leafdict] = dictionary()

    binwidth = 1.0
    
    i=1
    print("Looing at TTree %s" % tree1)
    print("Enter n to see next plot and q to exit program\n")
    for key,arr in T1_leafdict.items():
        # print key, -
        if (np.all(arr == 0.)):
            print("Histogram %s: Empty array" % key)
        elif ( 2. > len(arr)) :
            print("Histogram %s: Only one element" % key)
        else:
            binwidth = (abs(arr).max())/100
            plt.figure()
            plt.hist(arr,bins=np.arange(min(arr), max(arr) + binwidth, binwidth),histtype='step', stacked=True, fill=False )
            plt.title(key, fontsize =16)
            foutname = 'fig_'+key+'.png'
            i+=1

    print("\nTTree %s completed" % tree1)

def cut(cut,plot,low,high):
    
    [T1_leafdict] = dictionary()

    # arrCut = T1_leafdict[cut]
    # arrPlot = T1_leafdict[plot]

    arrCut = cut
    arrPlot = plot
    
    arrPlot = arrPlot[(arrCut > low) & (arrCut < high)]

    return[arrPlot]

def cutRecursive(lastCut,newcut,plot,low,high):
    
    [T1_leafdict] = dictionary()

    # arrLast = lastCut
    # arrCut = T1_leafdict[newcut]
    # arrPlot = T1_leafdict[plot]

    arrLast = lastCut
    arrCut = newcut
    arrPlot = plot
    
    arrPlot = arrPlot[(arrCut > low) & (arrCut < high)]

    arrNew = np.intersect1d(arrLast,arrPlot)

    return[arrNew]

# Wider the binsize the fewer bins
def setbin(leaf,binsize):

    
    binwidth = (abs(leaf).max()-abs(leaf).min())/binsize
    
    bins = np.arange(min(leaf), max(leaf) + binwidth, binwidth)

    return[bins]
    
def densityPlot(x,y,title,xlabel,ylabel,binx,biny):

    fig, ax = plt.subplots(tight_layout=True)
    hist = ax.hist2d(x, y,bins=(setbin(x,binx)[0],setbin(y,biny)[0]), norm=colors.LogNorm())
    plt.title(title, fontsize =16)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
# Can call arrays to create your own plots
def customPlots():

    xBj = lookup("xBj")[0]
    TDIS_xbj = lookup("TDIS_xbj")[0]
    sigma_dis = lookup("sigma_dis")[0]
    Q2 = lookup("Q2")[0]
    f2N = lookup("f2N")[0]
    tpi = lookup("tpi")[0]
    y_D = lookup("y_D")[0]
    # invts.y_D  = invts.Q2/(invts.xBj*invts.TwoPdotk);
    TwoPdotk = lookup("TwoPdotk")[0]
    y = Q2/(TDIS_xbj*TwoPdotk)
    # Yplus = lookup("Yplus")[0]+(1e-15)
    yplus = 1+((1-y)*(1-y))

    print 'Q2 \n:', Q2
    print 'TwoPdotk \n:', TwoPdotk
    print 'TDIS_xbj \n:', TDIS_xbj, np.sum(TDIS_xbj)
    print 'xBj \n:', xBj, np.sum(xBj)
    print 'TDIS_xbj*TwoPdotk \n:', TDIS_xbj*TwoPdotk
    print 'y \n:', y
    print 'y_D \n:', y_D, np.sum(y_D)

    y1 = 0.008
    x1 = 2.10e-2
    Q1 = 15.0
    sig1 = 0.519

    print '\n\n\n', ((2*math.pi*(1/137)*(1/137)*(1+((1-y1)*(1-y1))))/(x1*(Q1*Q1*Q1*Q1)))*(sig1), '\n\n\n'

    # Cuts
    xlow = [1.04e-4,2.53e-4,8.0e-4]
    xhigh = [2.2e-2,2.2e-2,2.2e-2]

    Q2low = [6.4,14.9,34.9]
    Q2high = [6.6,15.1,35.1]
    
    f, ax = plt.subplots()
    sigplot = ax.hist(sigma_dis,bins=setbin(sigma_dis,100)[0],histtype='step', stacked=True, fill=True)
    plt.title('d$sigma_{dis}$',fontsize =16)

    f, ax = plt.subplots()
    yplot = ax.hist(y,bins=setbin(y,100)[0],histtype='step', stacked=True, fill=True)
    plt.title('y',fontsize =16)
    
    f, ax = plt.subplots()
    yplusplot = ax.hist(yplus,bins=setbin(yplus,100)[0],histtype='step', stacked=True, fill=True)
    plt.title('$Y_{+}$',fontsize =16)
    # plt.xlim(1e-4,1e-1)
    # plt.xscale('log')
    # plt.ylim(0.,10)
    # plt.yscale('log')

    cut1_TDIS_xbj = cut(Q2,TDIS_xbj,Q2low[0],Q2high[0])[0]
    cut2_TDIS_xbj = cut(Q2,TDIS_xbj,Q2low[1],Q2high[1])[0]
    cut3_TDIS_xbj = cut(Q2,TDIS_xbj,Q2low[2],Q2high[2])[0]
    
    # Density plot
    densityPlot(TDIS_xbj,Q2,'$Q^{2}$ vs TDIS_xbj','TDIS_xbj','$Q^{2}$',100,100)
    
    cut1_sigma_dis = cut(Q2,sigma_dis,Q2low[0],Q2high[0])[0]
    cut2_sigma_dis = cut(Q2,sigma_dis,Q2low[1],Q2high[1])[0]
    cut3_sigma_dis = cut(Q2,sigma_dis,Q2low[2],Q2high[2])[0]
    
    f, ax = plt.subplots()
    sigscat1 = ax.scatter(cut1_TDIS_xbj,cut1_sigma_dis,label='$Q^2$=6.5 $GeV^2$')
    sigscat2 = ax.scatter(cut2_TDIS_xbj,cut2_sigma_dis,label='$Q^2$=15.0 $GeV^2$')
    sigscat3 = ax.scatter(cut3_TDIS_xbj,cut3_sigma_dis,label='$Q^2$=35.0 $GeV^2$')
    plt.title('d$\sigma_{DIS}$ vs TDIS_xbj', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.4,0.5), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('d$\sigma_{DIS}$ (nb)')
    # plt.xlim(1e-4,1e-1)
    plt.xscale('log')
    # plt.ylim(0.,10)
    plt.yscale('log')
    print "sigma_dis 1:\n",sigscat1.get_offsets()

    cut1_f2N = cut(Q2,f2N,Q2low[0],Q2high[0])[0]
    cut2_f2N = cut(Q2,f2N,Q2low[1],Q2high[1])[0]
    cut3_f2N = cut(Q2,f2N,Q2low[2],Q2high[2])[0]

    f, ax = plt.subplots()
    f2Nscat1 = ax.scatter(cut1_TDIS_xbj,cut1_f2N,label='$Q^2$=6.5 $GeV^2$')
    f2Nscat2 = ax.scatter(cut2_TDIS_xbj,cut2_f2N,label='$Q^2$=15.0 $GeV^2$')
    f2Nscat3 = ax.scatter(cut3_TDIS_xbj,cut3_f2N,label='$Q^2$=35.0 $GeV^2$')
    plt.title('$F^{N}_{2}$ vs TDIS_xbj', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.4,0.5), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('$F^{N}_{2}$')
    plt.xlim(1e-3,1e-1)
    plt.xscale('log')
    plt.ylim(0.2,0.6)
    # plt.yscale('log')
    print "f2N $Q^2$=6.5:\n",f2Nscat1.get_offsets()
    print "f2N $Q^2$=15.0:\n",f2Nscat2.get_offsets()
    print "f2N $Q^2$=35.0:\n",f2Nscat3.get_offsets()

    tot_sigma = (sigma_dis)*((TDIS_xbj*(Q2*Q2*Q2*Q2)*(137)*(137))/(2*math.pi*yplus))
    # tot_sigma = (f2N)*((TDIS_xbj*(Q2*Q2*Q2*Q2)*(137)*(137))/(2*math.pi*yplus))
    
    f, ax = plt.subplots()
    totsigplot = ax.hist(tot_sigma,bins=setbin(tot_sigma,100)[0],histtype='step', stacked=True, fill=True)
    plt.title('$sigma_{dis}$',fontsize =16)
    
    cut1a_tot_sigma = cut(Q2,tot_sigma,Q2low[0],Q2high[0])[0]
    cut1_tot_sigma = cutRecursive(cut1a_tot_sigma,TDIS_xbj,tot_sigma,xlow[0],xhigh[0])[0]
    cut2a_tot_sigma = cut(Q2,tot_sigma,Q2low[1],Q2high[1])[0]
    cut2_tot_sigma = cutRecursive(cut2a_tot_sigma,TDIS_xbj,tot_sigma,xlow[1],xhigh[1])[0]
    cut3a_tot_sigma = cut(Q2,tot_sigma,Q2low[2],Q2high[2])[0]
    cut3_tot_sigma = cutRecursive(cut3a_tot_sigma,TDIS_xbj,tot_sigma,xlow[2],xhigh[2])[0]
    
    cut1b_TDIS_xbj = cutRecursive(cut1_TDIS_xbj,TDIS_xbj,TDIS_xbj,xlow[0],xhigh[0])[0]
    cut2b_TDIS_xbj = cutRecursive(cut2_TDIS_xbj,TDIS_xbj,TDIS_xbj,xlow[1],xhigh[1])[0]
    cut3b_TDIS_xbj = cutRecursive(cut3_TDIS_xbj,TDIS_xbj,TDIS_xbj,xlow[2],xhigh[2])[0]
    
    f, ax = plt.subplots()
    # totsigscat = ax.scatter(TDIS_xbj,tot_sigma,label='No cuts')
    totsigscat1 = ax.scatter(cut1_TDIS_xbj,cut1a_tot_sigma,label='$Q^2$=6.5 $GeV^2$')
    totsigscat2 = ax.scatter(cut2_TDIS_xbj,cut2a_tot_sigma,label='$Q^2$=15.0 $GeV^2$')
    totsigscat3 = ax.scatter(cut3_TDIS_xbj,cut3a_tot_sigma,label='$Q^2$=35.0 $GeV^2$')
    plt.title('$\sigma_{DIS}$ vs TDIS_xbj', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.4,0.5), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('$\sigma_{DIS}$ ($nb^{-1}$)')
    # plt.xlim(1e-4,1e-1)
    plt.xscale('log')
    # plt.ylim(0.,10)
    plt.yscale('log')
    print "tot_sigma $Q^2$=6.5:\n",totsigscat1.get_offsets()
    print "tot_sigma $Q^2$=15.0:\n",totsigscat2.get_offsets()
    print "tot_sigma $Q^2$=35.0:\n",totsigscat3.get_offsets()


    # Histogram
    # cut1_tpi = cut(Q2,tpi,Q2low[0],Q2high[0])[0]
    # cut2_tpi = cut(Q2,tpi,Q2low[1],Q2high[1])[0]
    # cut3_tpi = cut(Q2,tpi,Q2low[2],Q2high[2])[0]
    # cut4a_tpi = cut(Q2,tpi,1.0,2.0)[0]
    # cut4_tpi = cutRecursive(cut4a_tpi,TDIS_xbj,tpi,0.001,0.01)[0]

    # f, ax = plt.subplots()
    # tpihist4a = ax.hist(cut4a_tpi,bins=setbin(cut4a_tpi,100)[0],histtype='step', stacked=True, fill=True,label='1.0 $GeV^2$ < $Q^2$  < 2.0 $GeV^2$' )
    # tpihist1 = ax.hist(cut1_tpi,bins=setbin(cut1_tpi,100)[0],histtype='step', stacked=True, fill=True,label='$Q^2$=6.5 $GeV^2$' )
    # tpihist2 = ax.hist(cut2_tpi,bins=setbin(cut2_tpi,100)[0],histtype='step', stacked=True, fill=True,label='$Q^2$=15.0 $GeV^2$' )
    # tpihist3 = ax.hist(cut3_tpi,bins=setbin(cut3_tpi,100)[0],histtype='step', stacked=True, fill=True,label='$Q^2$=35.0 $GeV^2$' )
    # tpihist4 = ax.hist(cut4_tpi,bins=setbin(cut4_tpi,100)[0],histtype='step', stacked=True, fill=True,label='0.001 < TDIS_xbj  < 0.01\n 1.0 $GeV^2$ < $Q^2$  < 2.0 $GeV^2$' )
    # plt.title("$t_{\pi}$ cuts", fontsize =16)
    # leg = plt.legend(bbox_to_anchor=(0.6,0.5), loc="center right")
    # leg.get_frame().set_alpha(1.)
    
def main() :

    customPlots()
    # recreateLeaves()
    plt.show()
    
if __name__=='__main__': main()
