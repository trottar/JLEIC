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
import csv

import dict as d

numEvts = len(d.findKey()) 
print("The number of events simulated is {}".format(numEvts))
#print(d.findKey.__doc__)

TDIS_xbj_raw = d.findKey('TDIS_xbj')
fpi_raw = d.findKey('fpi')
Q2_raw = d.findKey('TDIS_Q2')
t_raw = d.findKey('TDIS_t')
xL_raw = d.findKey('xL')
y_raw = d.findKey('TDIS_y')
sigma_dis_raw = d.findKey("sigma_dis")
f2N_raw = d.findKey("f2N")
xpi_raw = d.findKey("xpi")
ypi_raw = d.findKey("ypi")
tpi_raw = d.findKey("tpi")

xbinwidth = 0.01
qbinwidth = 10
tbinwidth = 0.1
#xbins =  np.linspace(0.,1.,100) 
xbins  = np.arange(xbinwidth/2,1.,xbinwidth).tolist()
qbins =  np.arange(qbinwidth/2,1000.,qbinwidth).tolist()
tbins =  np.arange(tbinwidth/2,1.,tbinwidth).tolist()
print("xBj Bins will be", xbins)
print("Q2 Bins will be", qbins)
print("t Bins will be", tbins)

def binBool(rawdata,bindata,binwidth):
    booldata = []
    for i,r in enumerate(rawdata):
        for b in bindata:
            if (b-(binwidth/20) < r < b+(binwidth/20)):
                booldata.append(True)
                break
        if i+1 != len(booldata):
            booldata.append(False)
    return booldata

xbinVal = np.array(binBool(TDIS_xbj_raw,xbins,xbinwidth))
qbinVal = np.array(binBool(Q2_raw,qbins,qbinwidth))
tbinVal = np.array(binBool(t_raw,qbins,tbinwidth))

def binData(lst, binType):
    arr_bin  = np.array(lst)
    arr_bin[~binType] = 0.
    return arr_bin

TDIS_xbj_xbin = binData(TDIS_xbj_raw, xbinVal)
Q2_xbin = binData(Q2_raw, xbinVal)
fpi_xbin = binData(fpi_raw, xbinVal)
t_xbin = binData(t_raw, xbinVal)
xL_xbin = binData(xL_raw, xbinVal)
y_xbin = binData(y_raw, xbinVal)
sigma_dis_xbin = binData(sigma_dis_raw, xbinVal)*(1e-5) # used to compare to HERA cross-section
f2N_xbin = binData(f2N_raw, xbinVal)
xpi_xbin = binData(xpi_raw, xbinVal)
ypi_xbin = binData(ypi_raw, xbinVal)
tpi_xbin = binData(tpi_raw, xbinVal)

TDIS_xbj_qbin = binData(TDIS_xbj_xbin, qbinVal)
Q2_qbin = binData(Q2_xbin, qbinVal)
fpi_qbin = binData(fpi_xbin, qbinVal)
t_qbin = binData(t_xbin, qbinVal)
xL_qbin = binData(xL_xbin, qbinVal)
y_qbin = binData(y_xbin, qbinVal)
sigma_dis_qbin = binData(sigma_dis_xbin, qbinVal)
f2N_qbin = binData(f2N_xbin, qbinVal)
xpi_qbin = binData(xpi_xbin, qbinVal)
ypi_qbin = binData(ypi_xbin, qbinVal)
tpi_qbin = binData(tpi_xbin, qbinVal)

TDIS_xbj_qbin = np.trim_zeros(TDIS_xbj_qbin)
Q2_qbin = np.trim_zeros(Q2_qbin)
fpi_qbin = np.trim_zeros(fpi_qbin)
t_qbin = np.trim_zeros(t_qbin)
xL_qbin = np.trim_zeros(xL_qbin)
y_qbin = np.trim_zeros(y_qbin)
sigma_dis_qbin = np.trim_zeros(sigma_dis_qbin)
f2N_qbin = np.trim_zeros(f2N_qbin)
xpi_qbin = np.trim_zeros(xpi_qbin)
ypi_qbin = np.trim_zeros(ypi_qbin)
tpi_qbin = np.trim_zeros(tpi_qbin)

TDIS_xbj_tbin = binData(TDIS_xbj_xbin, tbinVal)
Q2_tbin = binData(Q2_xbin, tbinVal)
fpi_tbin = binData(fpi_xbin, tbinVal)
t_tbin = binData(t_xbin, tbinVal)
xL_tbin = binData(xL_xbin, tbinVal)
y_tbin = binData(y_xbin, tbinVal)
sigma_dis_tbin = binData(sigma_dis_xbin, tbinVal)
f2N_tbin = binData(f2N_xbin, tbinVal)
xpi_tbin = binData(xpi_xbin, tbinVal)
ypi_tbin = binData(ypi_xbin, tbinVal)
tpi_tbin = binData(tpi_xbin, tbinVal)

TDIS_xbj_tbin = np.trim_zeros(TDIS_xbj_tbin)
Q2_tbin = np.trim_zeros(Q2_tbin)
fpi_tbin = np.trim_zeros(fpi_tbin)
t_tbin = np.trim_zeros(t_tbin)
xL_tbin = np.trim_zeros(xL_tbin)
y_tbin = np.trim_zeros(y_tbin)
sigma_dis_tbin = np.trim_zeros(sigma_dis_tbin)
f2N_tbin = np.trim_zeros(f2N_tbin)
xpi_tbin = np.trim_zeros(xpi_tbin)
ypi_tbin = np.trim_zeros(ypi_tbin)
tpi_tbin = np.trim_zeros(tpi_tbin)

d.tdict['TDIS_xbj'] = TDIS_xbj_qbin
d.tdict['TDIS_Q2'] = Q2_qbin
d.tdict['fpi'] = fpi_qbin
d.tdict['TDIS_t'] = t_qbin
d.tdict['xL'] = xL_qbin
d.tdict['TDIS_y'] = y_qbin
d.tdict['sigma_dis'] = sigma_dis_qbin
d.tdict['f2N'] = f2N_qbin
d.tdict['xpi'] = xpi_qbin
d.tdict['ypi'] = ypi_qbin
d.tdict['tpi'] = tpi_qbin

with open('datafiles/test.csv', 'w') as f:
    w = csv.writer(f)
    w.writerow(d.tdict.keys())
    w.writerows(zip(*d.tdict.values()))