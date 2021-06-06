#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2021-06-05 23:49:10 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import MaxNLocator

import sys
from sys import path
sys.path.insert(0,'./cuts/')
import cuts as c

import bindata as data

def densityPlot(x,y,title,xlabel,ylabel,binx,biny,
                    xmin=None,xmax=None,ymin=None,ymax=None,cuts=None,fig=None,ax=None,layered=True):

    if ax or fig:
        print("")
    else:
        fig, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))

    # norm=colors.LogNorm() makes colorbar normed and logarithmic
    hist = ax.hist2d(x, y,bins=(binx,biny),norm=colors.LogNorm())
    if layered is True :
        plt.colorbar(hist[3], ax=ax, spacing='proportional', label='Number of Events')

    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)

    inputVal = [x,y]
        
    return fig

# Create cut dictionary
cutDict = {}

'''
xbj = data.TDIS_xbj_raw
Q2 = data.Q2_raw
y = data.y_raw
t = data.t_raw
fpi = data.fpi_raw
xpi = data.TDIS_xbj_raw/(1-data.xL_raw)

xbj = data.TDIS_xbj_xbin
Q2 = data.Q2_xbin
y = data.y_xbin
t = data.t_xbin
fpi = data.fpi_xbin
xpi = data.TDIS_xbj_xbin/(1-data.xL_xbin)

xbj = data.TDIS_xbj_qbin
Q2 = data.Q2_qbin
y = data.y_qbin
t = data.t_qbin
fpi = data.fpi_qbin
xpi = data.TDIS_xbj_qbin/(1-data.xL_qbin)
'''
xbj = data.TDIS_xbj_tbin
Q2 = data.Q2_tbin
y = data.y_tbin
t = data.t_tbin
fpi = data.fpi_tbin
xpi = data.TDIS_xbj_tbin/(1-data.xL_tbin)

Q2binarray = [7,15,30,60,120,240,480,1000]
for i,x in enumerate(Q2binarray) :
    Q2tmp = '{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,Q2binarray[i]-data.qbinwidth/20,Q2binarray[i]+data.qbinwidth/20)
    print('{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,Q2binarray[i]-data.qbinwidth/20,Q2binarray[i]+data.qbinwidth/20))
    cutDict.update(eval(Q2tmp))

# xarray = np.arange(0.05,1.0,0.1).tolist() # table, large x (usual)
xarray = np.arange(data.xbinwidth/2,1.0,data.xbinwidth).tolist() # table, small x 
for i,x in enumerate(xarray):
    xtmp = '{"xcut%i" : ((%0.4f <= xpi) & (xpi <= %0.4f))}' % (i,xarray[i]-data.xbinwidth/2,xarray[i]+data.xbinwidth/2)
    print('{"xcut%i" : ((%0.4f <= xpi) & (xpi <= %0.4f))}' % (i,xarray[i]-data.xbinwidth/2,xarray[i]+data.xbinwidth/2))
    cutDict.update(eval(xtmp))

ytmp = '{"ycut" : ((0.01 <= y) & (y <= 0.95))}'
ttmp = '{"tcut" : ((-1.00 <= t) & (t <= 0.00))}'
cutDict.update(eval(ttmp))
cutDict.update(eval(ytmp))
cut = c.pyPlot(cutDict)    

ycut1 = ["ycut"]
tcut1 = ["tcut"]
#cut_t_y = ["xcut0","Q2cut3","tcut","ycut"] # Q2= 60 GeV^2
cut_t_y = ["xcut0","tcut","ycut"] # Q2= 60 GeV^2

def fpivxpi_Plot():
    
    f = plt.figure(figsize=(11.69,8.27))
    plt.rcParams.update({'font.size': 15})
    plt.style.use('classic')
    
    ax = f.add_subplot(221)
    xpiscat4 = ax.scatter(cut.applyCuts(xpi,cut_t_y),cut.applyCuts(fpi,cut_t_y),label='$Q^2$=60 $GeV^2$')
    #xpiscat4a = ax.scatter(xpi,fpi,label='$Q^2$=60 $GeV^2$')
    plt.plot([0.001,0.01,0.1],[1.2,0.45,0.25], label="GRV fit",color="y")
    plt.xscale('log')
    plt.ylim(-0.1,0.3)
    plt.xlim(1e-3,1.)
    ax.text(0.25, 0.65, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.set_yticks([-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3])
    ax.set_xticks([1e-2,1e-1])
    
    plt.ylabel('$F^{\pi}_{2}$', fontsize=20)
    
    '''
    ax = f.add_subplot(222)
    xpiscat5 = ax.errorbar(cut.applyCuts(xpi,cut120),cut.applyCuts(fpi,cut120),yerr=np.sqrt(cut.applyCuts(lumi,cut120))/cut.applyCuts(lumi,cut120),fmt='.',label='$Q^2$=120 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.5,0.25,0.15], label="GRV fit",color="y")
    plt.xscale('log')
    plt.ylim(-0.1,0.3)
    plt.xlim(1e-3,1.)
    ax.text(0.25, 0.65, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.set_yticks([-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3])
    ax.set_xticks([1e-2,1e-1])
    
    ax = f.add_subplot(223)
    xpiscat6 = ax.errorbar(cut.applyCuts(xpi,cut240),cut.applyCuts(fpi,cut240),yerr=np.sqrt(cut.applyCuts(lumi,cut240))/cut.applyCuts(lumi,cut240),fmt='.',label='$Q^2$=240 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit",color="y")
    plt.xscale('log')
    plt.ylim(-0.1,0.3)
    plt.xlim(1e-3,1.)
    ax.text(0.25, 0.65, '$Q^2$=240 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.set_yticks([-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3])
    ax.set_xticks([1e-2,1e-1])

    ax = f.add_subplot(224)
    xpiscat7 = ax.errorbar(cut.applyCuts(xpi,cut480),cut.applyCuts(fpi,cut480),yerr=np.sqrt(cut.applyCuts(lumi,cut480))/cut.applyCuts(lumi,cut480),fmt='.',label='$Q^2$=480 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit",color="y")
    plt.xscale('log')
    plt.ylim(-0.1,0.3)
    plt.xlim(1e-3,1.)
    ax.text(0.25, 0.65, '$Q^2$=480 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.set_yticks([-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3])
    ax.set_xticks([1e-2,1e-1])
    '''

    plt.xlabel('$x_\pi$', fontsize=20)    
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.0,wspace=0.0)

    plt.style.use('default')

fpivxpi_Plot()

fig = plt.figure(figsize=(17,12),facecolor='silver')

ax = fig.add_subplot(331)
#plt.scatter(data.t_qbin,data.fpi_qbin)
densityPlot(data.t_qbin,data.fpi_qbin, '$fpi$ vs $t$','$t$','$fpi$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(332)
#plt.scatter(TDIS_xbj_xbin,fpi_xbin)
densityPlot(data.TDIS_xbj_xbin,data.fpi_xbin, '$fpi$ vs $x$','$x$','$fpi$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(333)
#plt.scatter(TDIS_xbj_qbin,fpi_qbin)
densityPlot(data.TDIS_xbj_qbin,data.fpi_qbin, '$fpi$ vs $x$','$x$','$fpi$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(334)
densityPlot(data.TDIS_xbj_qbin,data.xL_qbin, '$xL$ vs $x$','$x$','$xL$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(335)
densityPlot(xbj,t, '$t$ vs $x$','$x$','$t$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(336)
densityPlot(data.xL_qbin,data.t_qbin, '$t$ vs $xL$','$xL$','$t$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(337)
densityPlot(data.TDIS_xbj_xbin,data.Q2_xbin, '$Q^2$ vs $x$','$x$','$Q^{2}$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(338)
densityPlot(data.TDIS_xbj_qbin,data.Q2_qbin, '$Q^2$ vs $x$','$x$','$Q^{2}$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(339)
densityPlot(cut.applyCuts(xbj,cut_t_y),cut.applyCuts(Q2,cut_t_y), '$Q^2$ vs $x$','$x$','$Q^{2}$', 200, 200, ax=ax, fig=fig)
#plt.scatter(cut.applyCuts(xbj,cut_t_y),cut.applyCuts(Q2,cut_t_y))

plt.tight_layout()
plt.show()

'''
Leftover code

#densityPlot(TDIS_xbj_raw,xL_raw, '$xL$ vs $x$','$x$','$xL$', 200, 200)
#plt.show()

#densityPlot(TDIS_xbj_raw,t_raw, '$t$ vs $x$','$x$','$t$', 200, 200)
#plt.show()

#densityPlot(xL_raw,t_raw, '$t$ vs $xL$','$xL$','$t$', 200, 200)
#plt.show()

#plt.scatter(TDIS_xbj_raw,fpi_raw)
#plt.show()

#plt.scatter(TDIS_xbj_raw,Q2_raw)
#plt.show()

#densityPlot(TDIS_xbj_raw,Q2_raw, '$Q^2$ vs $x$','$x$','$Q^{2}$', 200, 200)
#plt.show()

# Bins data weighted by Q2
Q2_qbin = (np.histogram(Q2_raw, qbins, weights=Q2_raw)[0] / np.histogram(Q2_raw, qbins)[0])
TDIS_xbj_qbin = (np.histogram(TDIS_xbj_raw, qbins, weights=Q2_raw)[0] / np.histogram(TDIS_xbj_raw, qbins)[0])
fpi_qbin = (np.histogram(fpi_raw, qbins, weights=Q2_raw)[0] / np.histogram(fpi_raw, qbins)[0])
print("\n\n",TDIS_xbj_qbin)
print(fpi_qbin)
print(Q2_qbin,"\n\n")

#plt.scatter(TDIS_xbj_qbin,fpi_qbin)
#plt.show()

#plt.scatter(TDIS_xbj_qbin,Q2_qbin)
#plt.show()

# Bins data weighted by TDIS_xbj
Q2_xbin = (np.histogram(Q2_raw, xbins, weights=TDIS_xbj_raw)[0] / np.histogram(Q2_raw, xbins)[0])
TDIS_xbj_xbin = (np.histogram(TDIS_xbj_raw, xbins, weights=TDIS_xbj_raw)[0] / np.histogram(TDIS_xbj_raw, xbins)[0])
t_xbin = (np.histogram(t_raw, xbins,weights=TDIS_xbj_raw)[0] / np.histogram(t_raw, xbins)[0])
fpi_xbin = (np.histogram(fpi_raw, xbins, weights=TDIS_xbj_raw)[0] / np.histogram(fpi_raw, xbins)[0])
print(TDIS_xbj_xbin)
print(fpi_xbin)
print(t_xbin)
print(Q2_xbin)
'''