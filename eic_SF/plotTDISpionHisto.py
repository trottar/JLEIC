#! /usr/bin/python

#
# Description:This will read in the array data file that contains all the leave histogram information
# ================================================================
# Time-stamp: "2019-04-08 08:05:34 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

import logging

# Gets rid of matplot logging DEBUG messages
plt_logger = logging.getLogger('matplotlib')
plt_logger.setLevel(logging.WARNING)

from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TTreeReader, TTreeReaderValue
from ROOT import gROOT
from rootpy.interactive import wait
import matplotlib.pyplot as plt
from matplotlib import interactive
from matplotlib import colors
import numpy as np
import sys

rootName = "TDISpion_80k"
# rootName = "TDISpion"

tree1 = "Evnts"

T1_arrkey =  "leafName1"
T1_arrhist = "histData1"

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
    
    T1_leafdict.get(key,"Leaf name not found")

    return[T1_leafdict.get(key,"Leaf name not found")]

# Recreates the histograms of the root file
def recreateLeaves():

    [T1_leafdict] = dictionary()

    binwidth = 1.0
    
    i=1
    print("Looing at TTree %s" % tree1)
    print("Enter n to see next plot and q to exit program\n")
    for key,arr in T1_leafdict.items():
        # print key, "->", arr)
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

    [T1,T1_hist] = pullArray()
    
    [T1_leafdict] = dictionary()

    arrCut = T1_leafdict[cut]
    arrPlot = T1_leafdict[plot]
    
    arrPlot = arrPlot[(arrCut > low) & (arrCut < high)]

    return[arrPlot]

def cutRecursive(lastCut,newcut,plot,low,high):

    [T1,T1_hist] = pullArray()
    
    [T1_leafdict] = dictionary()

    arrLast = lastCut
    arrCut = T1_leafdict[newcut]
    arrPlot = T1_leafdict[plot]
    
    arrPlot = arrPlot[(arrCut > low) & (arrCut < high)]

    arrNew = np.intersect1d(arrLast,arrPlot)

    return[arrNew]    
    
def densityPlot(x,y,title,xlabel,ylabel):

    fig, ax = plt.subplots(tight_layout=True)
    
    hist = ax.hist2d(x, y,bins=40, norm=colors.LogNorm())
    plt.title(title, fontsize =16)
    plt.xlabel(title)
    plt.ylabel(title)
    
    # Can call arrays to create your own plots
def customPlots():

    [T1_leafdict] = dictionary()

    TDIS_xbj = lookup("TDIS_xbj")
    sigma_dis = lookup("sigma_dis")
    Q2 = lookup("Q2")

    # Density plot
    densityPlot(TDIS_xbj[0],Q2[0],'$Q^{2}$ vs TDIS_xbj','TDIS_xbj','$Q^{2}$')

    # Scatter plot
    cut1_sigma_dis = cut("Q2","sigma_dis",6.9,7.1)[0]
    cut1_TDIS_xbj = cut("Q2","TDIS_xbj",6.9,7.1)[0]

    cut2_sigma_dis = cut("Q2","sigma_dis",14.9,15.1)[0]
    cut2_TDIS_xbj = cut("Q2","TDIS_xbj",14.9,15.1)[0]

    cut3_sigma_dis = cut("Q2","sigma_dis",29.9,30.1)[0]
    cut3_TDIS_xbj = cut("Q2","TDIS_xbj",29.9,30.1)[0]
    
    plt.figure()
    plt.scatter(cut1_TDIS_xbj,cut1_sigma_dis,label='$Q^2$=7.0 $GeV^2$')
    plt.scatter(cut2_TDIS_xbj,cut2_sigma_dis,label='$Q^2$=15.0 $GeV^2$')
    plt.scatter(cut3_TDIS_xbj,cut3_sigma_dis,label='$Q^2$=30.0 $GeV^2$')
    plt.title('$\sigma_{DIS}$ vs TDIS_xbj', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.4,0.5), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('$\sigma_{DIS}$ (nb)')
    plt.xlim(1e-4,1e-1)
    plt.xscale('log')
    # plt.ylim(0.,10)
    plt.yscale('log')

    # Scatter plot
    cut1_f2N = cut("Q2","f2N",6.9,7.1)[0]
    cut2_f2N = cut("Q2","f2N",14.9,15.1)[0]
    cut3_f2N = cut("Q2","f2N",29.9,30.1)[0]
    
    plt.figure()
    plt.scatter(cut1_TDIS_xbj,cut1_f2N,label='$Q^2$=7.0 $GeV^2$')
    plt.scatter(cut2_TDIS_xbj,cut2_f2N,label='$Q^2$=15.0 $GeV^2$')
    plt.scatter(cut3_TDIS_xbj,cut3_f2N,label='$Q^2$=30.0 $GeV^2$')
    plt.title('$F^{N}_{2}$ vs TDIS_xbj', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.4,0.5), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('$F^{N}_{2}$')
    plt.xlim(1e-4,1e-1)
    plt.xscale('log')
    # plt.ylim(0.,10)
    # plt.yscale('log')

    # Histogram
    cut1_tpi = cut("Q2","tpi",6.9,7.1)[0]
    cut2_tpi = cut("Q2","tpi",14.9,15.1)[0]
    cut3_tpi = cut("Q2","tpi",29.9,30.1)[0]
    cut4a_tpi = cut("Q2","tpi",1.0,2.0)[0]
    cut4_tpi = cutRecursive(cut4a_tpi,"TDIS_xbj","tpi",0.001,0.01)[0]

    binwidth = (abs(cut1_tpi).max())/100
    plt.figure()
    plt.hist(cut4a_tpi,bins=np.arange(min(cut1_tpi), max(cut1_tpi) + binwidth, binwidth),histtype='step', stacked=True, fill=True,label='1.0 $GeV^2$ < $Q^2$  < 2.0 $GeV^2$' )
    plt.hist(cut1_tpi,bins=np.arange(min(cut1_tpi), max(cut1_tpi) + binwidth, binwidth),histtype='step', stacked=True, fill=True,label='$Q^2$=7.0 $GeV^2$' )
    plt.hist(cut2_tpi,bins=np.arange(min(cut1_tpi), max(cut1_tpi) + binwidth, binwidth),histtype='step', stacked=True, fill=True,label='$Q^2$=15.0 $GeV^2$' )
    plt.hist(cut3_tpi,bins=np.arange(min(cut1_tpi), max(cut1_tpi) + binwidth, binwidth),histtype='step', stacked=True, fill=True,label='$Q^2$=30.0 $GeV^2$' )
    plt.hist(cut4_tpi,bins=np.arange(min(cut1_tpi), max(cut1_tpi) + binwidth, binwidth),histtype='step', stacked=True, fill=True,label='0.001 < TDIS_xbj  < 0.01\n 1.0 $GeV^2$ < $Q^2$  < 2.0 $GeV^2$' )
    plt.title("$t_{\pi}$ cuts", fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.6,0.5), loc="center right")
    leg.get_frame().set_alpha(1.)
    
def main() :

    customPlots()
    # recreateLeaves()
    plt.show()
    
if __name__=='__main__': main()
