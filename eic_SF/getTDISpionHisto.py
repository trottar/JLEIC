#! /usr/bin/python

#
# Descriptions:A light code that will convert the leaves of a ROOT file into arrays which can be easily manipulated and plotted in python
# ================================================================
# Time-stamp: "2019-04-08 05:51:43 trottar"
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
import sys

rootName = "TDISpion_80k"
# rootName = "TDISpion"

inputROOT = "%s.root" % rootName

tree1 = "Evnts"

def progressBar(value, endvalue, bar_length):

        percent = float(value) / endvalue
        arrow = '=' * int(round(percent * bar_length)-1) + '>'
        spaces = ' ' * (bar_length - len(arrow))

        sys.stdout.write(" \r[{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
        sys.stdout.flush()

def getTree():

    f = TFile.Open(inputROOT,"read")
    if not f or f.IsZombie() :
        exit

    Tree1 = f.Get(tree1)

    print("Tree %s" % tree1)
    print("="*78)
    Tree1.Print()
    print("="*78)
    print("\n")

    return[f,Tree1]

def getLeaves():

    [f,Tree1] = getTree()
    
    T1 = np.array([])

    l1 = Tree1.GetListOfLeaves()
    for i in range(0,l1.GetEntries()) :
        T1_hist = l1.At(i)
        T1 = np.append(T1,str(T1_hist.GetName()))
        
    return[f,T1]

def pullRootFiles():

    [f,T1] = getLeaves()

    # TTree1
    T1string = "Tree1 = f.%s" % tree1
    exec(T1string)
    
    hist1 = np.empty(shape=(1,1))
    
    j=0
    for n in np.nditer(T1):
        # progressBar(j,len(T1),70)
        print("\n%s" % str(T1[j]))
        th = []
        i=0
        for e in Tree1:
            if  '.' in str(T1[j]):
                progressBar(i,Tree1.GetEntries(),50)
                th.append(0.)
                # print("%i::Hist:%s" % (i,str(hist1[i])))
            elif 'TLorentzVector' in str(getattr(e,T1[j])):
                progressBar(i,Tree1.GetEntries(),50)
                th.append(0.)
                # print("%i::Hist:%s" % (i,str(hist2[i])))
            else:
                progressBar(i,Tree1.GetEntries(),50)
                th.append(getattr(e,T1[j]))
                # print("%i::Hist:%s" % (i,str(hist1[i])))
            i+=1
        hist1 = np.append(hist1,[j,np.asarray(th)])
        j+=1
        
    # hist is of form [0 array(histogram data) ... N array(histogram data)] where the 0 to N  elements are the leave numbers and histogram data elements are an array containing the data from the leaves
    # T1 contains strings with the names for each leaf (e.g. element 0 = 'e_Inc.')
    # The goal now is to match the 0 to N elements with the strings in T1 to match the string value with there corresponding histogram data
    # hist[0]=0.0, hist[1]=0, hist[2]=[0 0 ... 0] ,hist[3]=1
    
    # ->Therefore...
    # odd numbers -- 0 to N elements
    # even numbers -- array elements

    # T1_hist = np.zeros(shape=(1,1))

    # T1_hist = np.array([])

    T1_hist = []

    k=1
    l=0
    m=0
    for n in range(1,len(hist1)):
        if (float((k-1))/2).is_integer(): # odd
            print("\n%s:" % T1[l]) ##for debugging
            l+=1
        elif (float(k)/2).is_integer(): # even
            T1_hist.append(hist1[k])
            print("%s" % str(T1_hist[m])) ##for debugging
            m+=1
        k+=1

    # This above is correct
    # T1[(0,len(T1))] will give you the leaf names
    # T1_hist[(0,len(T1))] will give you the histogram data for each leaf

    # Next step is to copy this to some type of data file for future use.
    # I would like all the arrays to be saved into the same data file but we will see how practical that will be.

    # From there we can figure out the best way to make cuts and then finally extend  this to the hcana replay files.
    # This works pretty quickly so hopefully it is not bogged down too much by the number of events

    return[T1,T1_hist]

def sendArraytoFile():

    [T1,T1_hist] = pullRootFiles()

    np.savez(rootName, leafName1=T1, histData1=T1_hist)
    

def main() :

    sendArraytoFile()
    
if __name__=='__main__': main()
