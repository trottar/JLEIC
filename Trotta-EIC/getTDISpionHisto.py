#! /usr/bin/python

#
# Descriptions:A light code that will convert the leaves of a ROOT file into arrays which can be easily manipulated and plotted in python
# ================================================================
# Time-stamp: "2019-07-27 04:26:55 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TTreeReader, TTreeReaderValue
from ROOT import gROOT
from rootpy.interactive import wait
import numpy as np
import concurrent.futures
import sys,time,os,multiprocessing

# rootName = sys.argv[1]
rootName = "TDISpion"

inputROOT = "%s.root" % rootName

tree1 = "Evnts"

# Makes a progess bar in terminal
def progressBar(value, endvalue , bar_length):

    percent = float(value) / endvalue
    arrow = '=' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))
    
    sys.stdout.write(" \r[{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()

# Test if a value is numeric
def is_numeric(obj):
    
    attrs = ['__add__', '__sub__', '__mul__', '__div__', '__pow__']
    
    return all(hasattr(obj, attr) for attr in attrs)

# Opens root file and then prints out the tree and leaf information
def getTree():

    f = TFile.Open(inputROOT,"read")
    if not f or f.IsZombie() :
        exit

    Tree1 = f.Get(tree1)

    # print("Tree %s" % tree1)
    # print("="*78)
    # Tree1.Print()
    # print("="*78)
    # print("\n")

    return[f,Tree1]

# Grabs the leaf names and stores them
def getLeaves():

    [f,Tree1] = getTree()

    l1 = Tree1.GetListOfLeaves()
    
    # Stores the leaf names
    T1 = [l1.At(i).GetName() for i in range(0,l1.GetEntries())]
    
    return[f,T1]

# The critical function which takes the leaf names and stores the histogram values as arrays
def pullRootFiles():

    [f,T1] = getLeaves()
    
    hist = np.empty(shape=(1,1))
    
    # Begin timing the act of grabbing root files and converting to arrays
    start = time. time()

    print("Converting all root files to numpy arrays...")
    
    T1_array = range(len(T1))
    
    # Multiprocessing for multiple processes at once
    pool = multiprocessing.Pool(processes=6)
    process_result = pool.map(loopRoot,T1_array)
    pool.close()
    pool.join()
    
    for j,evt in enumerate(T1):
        hist = np.append(hist,[j,np.asarray(process_result[j])])

    # Close root file, no longer needed. Conversion complete!
    f.Close()

    # T1_hist is of form [0 array(histogram data) ... N array(histogram data)] where the 0 to N  elements are the leave numbers
    # and histogram data elements are an array containing the data from the leaves
    
    T1_hist = []

    # odd numbers -- 0 to N elements
    # even numbers -- array elements

    # for n in range(1,len(hist)):
    for k,evt in enumerate(hist):
        if (float(k+1)/2).is_integer(): # even
            T1_hist.append(hist[k+1])
        
    # T1[(0,len(T1))] will give you the leaf names
    # T1_hist[(0,len(T1))] will give you the histogram data for each leaf
            
    end = time. time()
    
    print("\nTime to pull root file: %0.1f seconds" % (end-start))
    
    return[T1_hist,T1]

# def loopRoot(key,ttree,index):
def loopRoot(index):

    [f,T1] = getLeaves()
    
    # TTree1
    T1string = "Tree1 = f.%s" % tree1
    exec(T1string)

    # print "\nProcess %s working" % os.getpid() ## Debugging multiprocessing
    tmp = [getattr(e,T1[index]) for i,e in enumerate(Tree1) if is_numeric(getattr(e,T1[index]))]
    # print "\nProcess %s done" % os.getpid() ## Debugging multiprocessing

    return tmp

# Saves the arrays to a *.npz file for future python use
def sendArraytoFile():
    
    [T1_hist,T1] = pullRootFiles()
    
    # for i,evt in enumerate(T1): ## DEBUGGING
        # print '\n', i, '\n', T1[i], T1_hist[i][0]
        
    np.savez_compressed(rootName, leafName=T1, histData=T1_hist)

def main() :

    sendArraytoFile()
    
if __name__=='__main__': main()
