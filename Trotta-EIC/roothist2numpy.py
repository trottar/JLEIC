#! /usr/bin/python

#
# Description: A light code that will convert the leaves of a ROOT file into arrays which can be easily manipulated and plotted in python
# ================================================================
# Time-stamp: "2019-08-07 19:31:09 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

from root_numpy import root2array, tree2array
import ROOT
import numpy as np
import sys,time,os

# rootName = sys.argv[1]
rootName = "TDISpion"

inputROOT = "%s.root" % rootName

tree1 = "Evnts"

# The critical function which takes the leaf names and stores the histogram values as arrays
def pullRootFiles():
    
    start = time. time()

    # Opens root file and gets tree info
    rfile = ROOT.TFile(inputROOT)
    Tree1 = rfile.Get(tree1)
    l1 = Tree1.GetListOfLeaves()
    
    # Stores the leaf names
    T1 = [l1.At(i).GetName() for i in range(0,l1.GetEntries()) if not "." in l1.At(i).GetName()]

    print("\nConverting root file to numpy arrays...\n")

    # Convert a TTree in a ROOT file into a NumPy structured array
    T1_hist = root2array(inputROOT, tree1, T1)

    # Matches elements to leaf array
    T1_hist = zip(*T1_hist)

    end = time. time()

    print("\nTime to pull root file: %0.1f seconds" % (end-start))

    return[T1_hist,T1]
    
# Saves the arrays to a *.npz file for future python use
def sendArraytoFile():
    
    [T1_hist,T1] = pullRootFiles()
    
    np.savez_compressed(rootName, leafName=T1, histData=T1_hist)

def main() :

    sendArraytoFile()
    
if __name__=='__main__': main()
