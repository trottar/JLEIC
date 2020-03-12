#! /usr/bin/python

#
# Description:Chain a series of ROOT2PY data file (*.npz) into one data file for future use 
# ================================================================
# Time-stamp: 2019-04-10 21:14:43 trottar
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

import numpy as np
import sys

numDatFiles = 2

chainName = "chain_8005-8010"

T1_arrkey =  "leafName"
T1_arrhist = "histData"

def progressBar(value, endvalue, bar_length):

        percent = float(value) / endvalue
        arrow = '=' * int(round(percent * bar_length)-1) + '>'
        spaces = ' ' * (bar_length - len(arrow))

        sys.stdout.write(" \r[{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
        sys.stdout.flush()


# Retrieves the array data file and creates new leaves from this
def pullArray():

    rootDict = {
         0 : "/home/trottar/ResearchNP/JLEIC/eic_SF/TDISpion_1-4",
         1 : "/home/trottar/ResearchNP/JLEIC/eic_SF/TDISpion_4-8",
    }

    dictT1 = {}

    dictT1_hist = {}

    j=0
    print("Gathering data files...")
    for i in range(0,numDatFiles):
        progressBar(j,numDatFiles, 50)
        rootName = rootDict[i]
        data = np.load("%s.npz" % rootName)
        
        dictT1["T1_%i" % i] = data[T1_arrkey]
        dictT1_hist["T1__hist_%i" % i] = data[T1_arrhist]
        j+=1

    return[dictT1,dictT1_hist]

# Creates a dictionary that stores leaf names with the corresponding array
def mergeKeys():

    [dictT1,dictT1_hist] = pullArray()

    T1_leafdicttmp = np.array([])
    
    for i in range(0,numDatFiles):
        T1_leafdicttmp = np.append(T1_leafdicttmp, dict(zip(dictT1.get("T1_%i" % i, "Leaf not found"),dictT1_hist.get("T1__hist_%i" % i, "Leaf not found"))))
        T1_leafdictmerge = [dict(zip(dictT1.get("T1_%i" % i, "Leaf not found"), dictT1_hist.get("T1__hist_%i" % i, "Leaf not found"))), dict(zip(dictT1.get("T1_%i" % i, "Leaf not found"), dictT1_hist.get("T1__hist_%i" % i, "Leaf not found")))]
        T1_leafdictnew = {}
        for k in dict(zip(dictT1.get("T1_%i" % i, "Leaf not found"), dictT1_hist.get("T1__hist_%i" % i, "Leaf not found"))).iterkeys():
                  T1_leafdictnew[k] = tuple(T1_leafdictnew[k] for T1_leafdictnew in T1_leafdictmerge)
    
    return[T1_leafdictnew,dictT1]

def dictionary():

    [T1_leafdictnew,dictT1] = mergeKeys()

    T1_leafdict = {}
    
    for key,arr in T1_leafdictnew.items():
        tmp = []
        for i in range(0,numDatFiles):
             tmp = np.append(tmp, T1_leafdictnew[key][i])
        T1_leafdict[key] = tmp
    
    return[T1_leafdict]

def sendArraytoFile():

    [T1_leafdict] = dictionary()
    
    T1 = np.array([])
    T1_hist = []

    for key,arr in T1_leafdict.items():
        T1 = np.append(T1,key)
        T1_hist.append(arr)

    print("\n\nLoading chained data to file, this may take a few minutes.")
    np.savez(chainName, leafName=T1, histData=T1_hist)
        
def main() :

    sendArraytoFile()
    
if __name__=='__main__': main()	    
