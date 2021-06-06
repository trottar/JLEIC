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
qbins =  [60,120,240,480]
#qbins =  np.linspace(0.,1000.,100)
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

def binData(lst, binType):
    arr_bin  = np.array(lst)
    arr_bin[~binType] = 0.
    return arr_bin

TDIS_xbj_xbin = binData(TDIS_xbj_raw, xbinVal)
Q2_xbin = binData(Q2_raw, xbinVal)
fpi_xbin = binData(fpi_raw, xbinVal)
t_xbin = binData(t_raw, xbinVal)
xL_xbin = binData(xL_raw, xbinVal)

TDIS_xbj_qbin = binData(TDIS_xbj_xbin, qbinVal)
Q2_qbin = binData(Q2_xbin, qbinVal)
fpi_qbin = binData(fpi_xbin, qbinVal)
t_qbin = binData(t_xbin, qbinVal)
xL_qbin = binData(xL_xbin, qbinVal)

TDIS_xbj_qbin = np.trim_zeros(TDIS_xbj_qbin)
Q2_qbin = np.trim_zeros(Q2_qbin)
fpi_qbin = np.trim_zeros(fpi_qbin)
t_qbin = np.trim_zeros(t_qbin)
xL_qbin = np.trim_zeros(xL_qbin)

'''
Next step is to restructure MC root files then bin everything
'''