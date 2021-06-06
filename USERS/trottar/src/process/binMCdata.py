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
import awkward as ak
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import dict as d

numEvts = len(d.findKey()) 
print("The number of events simulated is {}".format(numEvts))
#print(d.findKey.__doc__)

"""
Need to figure out how to properly access the 'invts' branch
"""
#table = ak.Table(tdict["invts"])
#print(table)

TDIS_xbj_raw = d.findKey('TDIS_xbj')
fpi_raw = d.findKey('fpi')
Q2_raw = d.findKey('TDIS_Q2')
t_raw = d.findKey('TDIS_t')
xL_raw = d.findKey('xL')
print("Average (no bin)",np.average(TDIS_xbj_raw))

xbinwidth = 0.1
xbins =  np.linspace(0.,1.,1/xbinwidth) 
qbinwidth = 10
#qbins =  [60,120,240,480]
qbins =  np.linspace(0.,1000.,100)
print("xBj Bins will be", xbins)
print("Q2 Bins will be", qbins)

def binBool(rawdata,bindata,binwidth):
    booldata = []
    for i,r in enumerate(rawdata):
        for b in bindata:
            if (b-binwidth/20 < r < b+binwidth/20):
                booldata.append(True)
                break
        if i+1 != len(booldata):
            booldata.append(False)
    return booldata

xbinVal = np.array(binBool(TDIS_xbj_raw,xbins,xbinwidth))
qbinVal = np.array(binBool(Q2_raw,qbins,qbinwidth))

TDIS_xbj_xbin = np.array(TDIS_xbj_raw)
Q2_xbin = np.array(Q2_raw)
fpi_xbin = np.array(fpi_raw)
t_xbin = np.array(t_raw)
xL_xbin = np.array(xL_raw)
TDIS_xbj_xbin[~xbinVal] = 0.
Q2_xbin[~xbinVal] = 0.
fpi_xbin[~xbinVal] = 0.
t_xbin[~xbinVal] = 0.
xL_xbin[~xbinVal] = 0.

TDIS_xbj_qbin = np.array(TDIS_xbj_xbin)
Q2_qbin = np.array(Q2_xbin)
fpi_qbin = np.array(fpi_xbin)
t_qbin = np.array(t_xbin)
xL_qbin = np.array(xL_xbin)
TDIS_xbj_qbin[~qbinVal] = 0.
Q2_qbin[~qbinVal] = 0.
fpi_qbin[~qbinVal] = 0.
t_qbin[~qbinVal] = 0.
xL_qbin[~qbinVal] = 0.

'''
Next step is to restructure MC root files then bin everything
'''