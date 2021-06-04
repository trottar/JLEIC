#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2021-06-03 20:58:22 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

import pandas as pd
import numpy as np
import uproot as up
import awkward as ak
import sys

import matplotlib.pyplot as plt

kinematics = sys.argv[1]
rootName="../../OUTPUTS/%s.root" % kinematics

tName = "Evnts"
print("Reading tree {} in file {}".format(tName, rootName))
tdata = up.open(rootName)[tName]

bnames = tdata.keys()
print ("Tree has the following branches:")
print ("  [{}]".format(', '.join(bnames)))

tdict = {}
print("Tree has the following data in it's branches:")
for i,val in enumerate(bnames):
    print("{0} = {1}".format(val,tdata[val].array()))
    tdict[val] = tdata[val].array()

print("Tree dictionary of all branches:")
print("{}".format(tdict))

def findKey(key='invts'):
    """
    findKey(key='invts')
    
    """
    if key=='invts':
        # [key for key in tdict][0] finds the first key in the dictionary
        return list(tdict[[k for k in tdict][0]])
    else:
        # [k for k in tdict if k==key] finds the key matching the input
        return list(tdict[[k for k in tdict if k==key][0]])

numEvts = len(findKey()) 
print("The number of events simulated is {}".format(numEvts))
#print(findKey.__doc__)

"""
Need to figure out how to properly access the 'invts' branch
"""
#table = ak.Table(tdict["invts"])
#print(table)

TDIS_xbj_raw = findKey('TDIS_xbj')
fpi_raw = findKey('fpi')
print("Average (no bin)",np.average(TDIS_xbj_raw))

bins =  np.linspace(0.,1.,1000) # 0.001 bin size
print("Bins will be", bins)

# Bins data weighted by TDIS_xbj
TDIS_xbj_bin = (np.histogram(TDIS_xbj_raw, bins, weights=TDIS_xbj_raw)[0] / np.histogram(TDIS_xbj_raw, bins)[0])
fpi_bin = (np.histogram(fpi_raw, bins, weights=TDIS_xbj_raw)[0] / np.histogram(fpi_raw, bins)[0])
print(TDIS_xbj_bin)
print(fpi_bin)

print("Average (no bin)",np.average(TDIS_xbj_bin))

plt.scatter(TDIS_xbj_bin,fpi_bin)
plt.show()

'''
Next step is to restructure MC root files then bin everything
'''