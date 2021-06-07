#! /usr/bin/python

#
# Description:This will read in the array data file that contains all the leave histogram information
# include...
# ================================================================
# Time-stamp: "2019-08-14 21:10:09 trottar"
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

from ROOT import TFile, TH1F
import scipy as sc
from scipy import stats, optimize, interpolate
import matplotlib.pyplot as plt
from matplotlib import interactive
from matplotlib import colors
import uproot as up
import time, math, sys

# garbage collector
import gc
gc.collect()

class pyDict(dict):
    
    def __init__(self,inputTree):
        self.inputTree = inputTree
            
class pyPlot(pyDict):
    
    def __init__(self, cutDict=None):
        self.cutDict = cutDict

    def setbin(self,plot,numbin,xmin=None,xmax=None):
        
        if (xmin or xmax):
            leaf = self.fixBin(plot,plot,xmin,xmax)
        else:
            leaf = plot
            
        binwidth = (abs(leaf).max()-abs(leaf).min())/numbin
        
        bins = np.arange(min(leaf), max(leaf) + binwidth, binwidth)

        return bins

    def fixBin(self,cut,plot,low,high):
            
        arrCut = cut
        arrPlot = plot
        arrPlot = arrPlot[(arrCut > low) & (arrCut < high)]

        return arrPlot

    def read_dict(self,f):

        cutDict = {}
        cut_new = ()
        for line in f:
            if "#" in line:
                continue
            else:
                line  = line.split("=")
                cuts = line[1]
                cutName = {line[0].rstrip() : cuts}
                cutDict.update(cutName)
        return cutDict

    # Create a working dictionary for cuts
    def w_dict(self,cuts):

        inputDict = self.cutDict
        subDict = inputDict[cuts]
        subDict = subDict.split(",")
        cut_arr = [evt for evt in subDict]
        return cut_arr
            
    def cut(self,key,cuts=None):

        if cuts:
            inputDict = self.cutDict
            subDict = inputDict[cuts]
            value = subDict.get(key,"Leaf name not found")
            return value
        else:
            return self.cutDict.get(key,"Leaf name not found")

    # Old version of apply cuts
    def applyCuts(self,leaf,cuts=None,either=False):
        
        if cuts:
            if either==True:
                tmp = leaf
                applycut = 'tmp['
                i=0
                while i < (len(cuts)-1):
                    applycut += 'self.cut("%s") | ' % cuts[i]
                    i+=1
                applycut += 'self.cut("%s")]' % cuts[len(cuts)-1]
                tmp = eval(applycut)
            else:
                tmp = leaf
                applycut = 'tmp['
                i=0
                while i < (len(cuts)-1):
                    applycut += 'self.cut("%s") & ' % cuts[i]
                    i+=1
                applycut += 'self.cut("%s")]' % cuts[len(cuts)-1]
                tmp = eval(applycut)
        else:
            print('No cuts applied to %s' % leaf)
            tmp = leaf
            
        return tmp

    # New version of applying cuts
    def add_cut(self,arr, cuts):
        
        arr_cut = arr  
        applycut = "arr_cut["
        inputDict = self.cutDict
        subDict = inputDict[cuts]
        for i,(key,val) in enumerate(subDict.items()):
            if i == len(subDict)-1:
                applycut += 'self.cut("%s","%s")]' % (key,cuts)
            else:
                applycut += 'self.cut("%s","%s") & ' % (key,cuts)
        arr_cut = eval(applycut)        
        return arr_cut


    def progressBar(self,value, endvalue, bar_length):

        percent = float(value) / endvalue
        arrow = '=' * int(round(percent * bar_length)-1) + '>'
        spaces = ' ' * (bar_length - len(arrow))
        
        sys.stdout.write(" \r[{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
        sys.stdout.flush()

    # Recreates the histograms of the root file
    def recreateLeaves(self):
                
        binwidth = 1.0
    
        i=1
        print("Looing at TTree %s" % self.tree1)
        print("Enter n to see next plot and q to exit program\n")
        # for key,arr in self.T1_leafdict.dictionary().items():
        for key,arr in self.T1_leafdict.items():
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

        print("\nTTree %s completed" % self.tree1)
