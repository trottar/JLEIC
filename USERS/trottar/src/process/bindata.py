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
print("\n\nThe number of events simulated is {}\n\n".format(numEvts))
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
qbinwidth = 10.
tbinwidth = 0.01
xLbinwidth = 0.01

xbins  = np.arange(xbinwidth/2,1.,xbinwidth).tolist()
qbins =  np.arange(qbinwidth/2,1000.,qbinwidth).tolist()
tbins =  np.arange(tbinwidth/2,1.,tbinwidth).tolist()
xLbins =  np.arange(xLbinwidth/2,1.,xLbinwidth).tolist()

print("xBj Bins will be", xbins)
print("Q2 Bins will be", qbins)
print("t Bins will be", tbins)
print("xL Bins will be", xLbins,"\n\n")

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
tbinVal = np.array(binBool(t_raw,tbins,tbinwidth))
xLbinVal = np.array(binBool(xL_raw,xLbins,xLbinwidth))

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

'''
TDIS_xbj_bin = np.trim_zeros(TDIS_xbj_qbin)
Q2_bin = np.trim_zeros(Q2_qbin)
fpi_bin = np.trim_zeros(fpi_qbin)
t_bin = np.trim_zeros(t_qbin)
xL_bin = np.trim_zeros(xL_qbin)
y_bin = np.trim_zeros(y_qbin)
sigma_dis_bin = np.trim_zeros(sigma_dis_qbin)
f2N_bin = np.trim_zeros(f2N_qbin)
xpi_bin = np.trim_zeros(xpi_qbin)
ypi_bin = np.trim_zeros(ypi_qbin)
tpi_bin = np.trim_zeros(tpi_qbin)
'''
'''
TDIS_xbj_tbin = binData(TDIS_xbj_qbin, tbinVal)
Q2_tbin = binData(Q2_qbin, tbinVal)
fpi_tbin = binData(fpi_qbin, tbinVal)
t_tbin = binData(t_qbin, tbinVal)
xL_tbin = binData(xL_qbin, tbinVal)
y_tbin = binData(y_qbin, tbinVal)
sigma_dis_tbin = binData(sigma_dis_qbin, tbinVal)
f2N_tbin = binData(f2N_qbin, tbinVal)
xpi_tbin = binData(xpi_qbin, tbinVal)
ypi_tbin = binData(ypi_qbin, tbinVal)
tpi_tbin = binData(tpi_qbin, tbinVal)

TDIS_xbj_tbin = np.trim_zeros(TDIS_xbj_tbin)
Q2_bin = np.trim_zeros(Q2_tbin)
fpi_bin = np.trim_zeros(fpi_tbin)
t_bin = np.trim_zeros(t_tbin)
xL_bin = np.trim_zeros(xL_tbin)
y_bin = np.trim_zeros(y_tbin)
sigma_dis_bin = np.trim_zeros(sigma_dis_tbin)
f2N_bin = np.trim_zeros(f2N_tbin)
xpi_bin = np.trim_zeros(xpi_tbin)
ypi_bin = np.trim_zeros(ypi_tbin)
tpi_bin = np.trim_zeros(tpi_tbin)
'''
#'''
TDIS_xbj_xLbin = binData(TDIS_xbj_qbin, xLbinVal)
Q2_xLbin = binData(Q2_qbin, xLbinVal)
fpi_xLbin = binData(fpi_qbin, xLbinVal)
t_xLbin = binData(t_qbin, xLbinVal)
xL_xLbin = binData(xL_qbin, xLbinVal)
y_xLbin = binData(y_qbin, xLbinVal)
sigma_dis_xLbin = binData(sigma_dis_qbin, xLbinVal)
f2N_xLbin = binData(f2N_qbin, xLbinVal)
xpi_xLbin = binData(xpi_qbin, xLbinVal)
ypi_xLbin = binData(ypi_qbin, xLbinVal)
tpi_xLbin = binData(tpi_qbin, xLbinVal)

TDIS_xbj_bin = np.trim_zeros(TDIS_xbj_xLbin)
Q2_bin = np.trim_zeros(Q2_xLbin)
fpi_bin = np.trim_zeros(fpi_xLbin)
t_bin = np.trim_zeros(t_xLbin)
xL_bin = np.trim_zeros(xL_xLbin)
y_bin = np.trim_zeros(y_xLbin)
sigma_dis_bin = np.trim_zeros(sigma_dis_xLbin)
f2N_bin = np.trim_zeros(f2N_xLbin)
xpi_bin = np.trim_zeros(xpi_xLbin)
ypi_bin = np.trim_zeros(ypi_xLbin)
tpi_bin = np.trim_zeros(tpi_xLbin)
#'''

if (len(fpi_qbin) == 0):
    print("Error: {} is null".format("fpi_qbin"))
else:
    d.tdict['TDIS_xbj'] = TDIS_xbj_bin
    d.tdict['TDIS_Q2'] = Q2_bin
    d.tdict['fpi'] = fpi_bin
    d.tdict['TDIS_t'] = t_bin
    d.tdict['xL'] = xL_bin
    d.tdict['TDIS_y'] = y_bin
    d.tdict['sigma_dis'] = sigma_dis_bin
    d.tdict['f2N'] = f2N_bin
    d.tdict['xpi'] = xpi_bin
    d.tdict['ypi'] = ypi_bin
    d.tdict['tpi'] = tpi_bin

with open('datafiles/{0}{1}{2}{3}_{4}.csv'.format('x%0.3f' % xbinwidth,'q%0.1f' % qbinwidth,'t%0.3f' % tbinwidth,'xL%0.3f' % xLbinwidth,d.rootName.strip("../../OUTPUTS/").strip(".root")), 'w') as f:
    w = csv.writer(f)
    w.writerow(d.tdict.keys())
    w.writerows(zip(*d.tdict.values()))