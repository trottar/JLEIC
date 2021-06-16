#! /usr/bin/python

#
# Description: A light code that will convert the leaves of a ROOT file into arrays which can be easily manipulated and plotted in python
# ================================================================
# Time-stamp: "2019-08-14 15:26:57 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

import uproot as up
import numpy as np
import sys,time,os

# rootName = sys.argv[1]
rootName = "TDISpion"

inputROOT = "%s.root" % rootName

tree1 = "Evnts"

start = time. time()

# The critical function which takes the leaf names and stores the histogram values as arrays
def pullRootFiles():

    # Opens root file and gets tree info
    rfile = up.open(inputROOT)
    print rfile.keys()
    Tree1 = rfile[tree1]
    print len(Tree1.keys())

    print Tree1.allkeys()
    
    # Stores the leaf names
    
    print("\nConverting root file to numpy arrays...\n")
    
    # Convert a TTree in a ROOT file into a NumPy structured array
    T1_array = Tree1.arrays()

    Tree1.show()
    
    T1_array = T1_array.items()

    T1 = []
    T1_hist = []
    
    # Matches elements to leaf array
    for i,evts in enumerate(T1_array):
        if type(T1_array[i][1][0]) is np.dtype('>f8').type: # Check if float (i.e. histogram)
            T1.append(T1_array[i][0])
            T1_hist.append(T1_array[i][1])
    
    # Need to get invts branch separately
    T1_invt =  Tree1.array("invts")
    T1_invt  = zip(*T1_invt) # Match elements to proper leaves
    
    T1_invName = []
    tmp = Tree1["invts"].interpretation.fromdtype.descr # Converts branch names that are a dtype object to a list
    
    for name,typ in tmp:
        T1_invName.append(name)

    for i,evts in enumerate(T1_invt):
        T1.append(T1_invName[i])
        T1_hist.append(T1_invt[i])
    
    # for i,evts in enumerate(T1_hist):
    #     print T1[i], np.array(T1_hist[i])
    
    return[T1_hist,T1]
    
# Saves the arrays to a *.npz file for future python use
def sendArraytoFile():
    [T1_hist,T1] = pullRootFiles()    
    # np.savez_compressed(rootName, leafName=T1, histData=T1_hist)
    np.savez(rootName, leafName=T1, histData=T1_hist)

def main() :
    sendArraytoFile()
    # pullRootFiles()
    
if __name__=='__main__': main()

end = time. time()
print("\nTime to pull root file: %0.1f seconds" % (end-start))
