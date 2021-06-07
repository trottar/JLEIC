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

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import MaxNLocator
from scipy.interpolate import griddata

import sys
from sys import path
sys.path.insert(0,'./cuts/')
import cuts as c

df = pd.read_csv(r'./datafiles/test.csv')
print(df)

xbj = df['TDIS_xbj']
Q2 = df['TDIS_Q2']
fpi = df['fpi']
t = df['TDIS_t']
xL = df['xL']
y = df['TDIS_y']
sigma_dis = df['sigma_dis']
f2N = df['f2N']
xpi = df['xpi']
#xpi = xbj/(1.-xL)
ypi = df['ypi']
tpi = df['tpi']


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

qbinwidth = 10
xbinwidth = 0.01
tbinwidth = 0.1

Q2binarray = [7,15,30,60,120,240,480,1000]
#Q2binarray = np.arange(qbinwidth/2,1000.,qbinwidth).tolist()
for i,x in enumerate(Q2binarray) :
    Q2tmp = '{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,Q2binarray[i]-qbinwidth/2,Q2binarray[i]+qbinwidth/2)
    print('{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,Q2binarray[i]-qbinwidth/2,Q2binarray[i]+qbinwidth/2))
    cutDict.update(eval(Q2tmp))

xarray = np.arange(xbinwidth/2,1.0,xbinwidth).tolist()
for i,x in enumerate(xarray):
    xtmp = '{"xcut%i" : ((%0.4f <= xbj) & (xbj <= %0.4f))}' % (i,xarray[i]-xbinwidth/2,xarray[i]+xbinwidth/2)
    print('{"xcut%i" : ((%0.4f <= xbj) & (xbj <= %0.4f))}' % (i,xarray[i]-xbinwidth/2,xarray[i]+xbinwidth/2))
    cutDict.update(eval(xtmp))

'''
tarray = np.arange(tbinwidth/2,1.0,tbinwidth).tolist()
for i,t in enumerate(tarray):
    ttmp = '{"tcut%i" : ((%0.4f <= t) & (t <= %0.4f))}' % (i,tarray[i]-tbinwidth/2,tarray[i]+tbinwidth/2)
    print('{"tcut%i" : ((%0.4f <= t) & (t <= %0.4f))}' % (i,tarray[i]-tbinwidth/2,tarray[i]+tbinwidth/2))
    cutDict.update(eval(ttmp))
'''

ytmp = '{"ycut" : ((0.01 <= y) & (y <= 0.95))}'
cutDict.update(eval(ytmp))
cut = c.pyPlot(cutDict)    

ycut1 = ["ycut"]

cut_q = ["xcut2","Q2cut3","ycut"] # Q2= 60 GeV^2

cut_x2_q3 = ["xcut2","Q2cut3","ycut"] # Q2= 60 GeV^2
cut_x4_q3 = ["xcut4","Q2cut3","ycut"] # Q2= 60 GeV^2
cut_x6_q3= ["xcut6","Q2cut3","ycut"] # Q2= 60 GeV^2
cut_x8_q3 = ["xcut8","Q2cut3","ycut"] # Q2= 60 GeV^2

cut_x2_q4 = ["xcut2","Q2cut4","ycut"] # Q2= 120 GeV^2
cut_x4_q4 = ["xcut4","Q2cut4","ycut"] # Q2= 120 GeV^2
cut_x6_q4= ["xcut6","Q2cut4","ycut"] # Q2= 120 GeV^2
cut_x8_q4 = ["xcut8","Q2cut4","ycut"] # Q2= 120 GeV^2

cut_x2_q5 = ["xcut2","Q2cut5","ycut"] # Q2= 240 GeV^2
cut_x4_q5 = ["xcut4","Q2cut5","ycut"] # Q2= 240 GeV^2
cut_x6_q5= ["xcut6","Q2cut5","ycut"] # Q2= 240 GeV^2
cut_x8_q5 = ["xcut8","Q2cut5","ycut"] # Q2= 240 GeV^2

cut_x2_q6 = ["xcut2","Q2cut6","ycut"] # Q2= 480 GeV^2
cut_x4_q6 = ["xcut4","Q2cut6","ycut"] # Q2= 480 GeV^2
cut_x6_q6= ["xcut6","Q2cut6","ycut"] # Q2= 480 GeV^2
cut_x8_q6 = ["xcut8","Q2cut6","ycut"] # Q2= 480 GeV^2

def F2pi(xpi, Q2):
    points,values=np.load('./../../analysis/interpGrids/xpiQ2.npy'),np.load('./../../analysis/interpGrids/F2pi.npy')
    F2pi=lambda xpi,Q2: griddata(points,values,(np.log10(xpi),np.log10(Q2)))
    return F2pi(xpi,Q2)

print("~",F2pi(xpi,Q2))

# Calculate cross-section using Patrick's interpolate grid
def ds_dxdQ2dxLdt(x, xL,t):
    points60,values60=np.load('./../../analysis/xsec/pointsxsec60.npy'),np.load('./../../analysis/xsec/valuesxsec60.npy')
    points120,values120=np.load('./../../analysis/xsec/pointsxsec120.npy'),np.load('./../../analysis/xsec/valuesxsec120.npy')
    points240,values240=np.load('./../../analysis/xsec/pointsxsec240.npy'),np.load('./../../analysis/xsec/valuesxsec240.npy')
    points480,values480=np.load('./../../analysis/xsec/pointsxsec480.npy'),np.load('./../../analysis/xsec/valuesxsec480.npy')
    sigma60=lambda x,xL,t: griddata(points60,values60,(x,xL,t))
    sigma120=lambda x,xL,t: griddata(points120,values120,(x,xL,t))
    sigma240=lambda x,xL,t: griddata(points240,values240,(x,xL,t))
    sigma480=lambda x,xL,t: griddata(points480,values480,(x,xL,t))

    return [sigma60(x,xL,t),sigma120(x,xL,t),sigma240(x,xL,t),sigma480(x,xL,t)]

print("~",ds_dxdQ2dxLdt(xbj,xL,t))

def fpivxpi_Plot():
    
    f = plt.figure(figsize=(11.69,8.27))
    plt.rcParams.update({'font.size': 15})
    plt.style.use('classic')
    
    ax = f.add_subplot(221)
    #xpiscat3 = ax.scatter(cut.applyCuts(xpi,cut_q),cut.applyCuts(fpi,cut_q),label='$Q^2$=?? $GeV^2$',c='red')
    xpiscat3 = ax.scatter(cut.applyCuts(xpi,cut_x2_q3),F2pi(cut.applyCuts(xpi,cut_x2_q3),cut.applyCuts(Q2,cut_x2_q3)),label='$Q^2$=60 $GeV^2$, $x_{\pi}$ = 0.02-0.03')
    xpiscat3 = ax.scatter(cut.applyCuts(xpi,cut_x4_q3),F2pi(cut.applyCuts(xpi,cut_x4_q3),cut.applyCuts(Q2,cut_x4_q3)),label='$Q^2$=60 $GeV^2$, $x_{\pi}$ = 0.04-0.05')
    xpiscat3 = ax.scatter(cut.applyCuts(xpi,cut_x6_q3),F2pi(cut.applyCuts(xpi,cut_x6_q3),cut.applyCuts(Q2,cut_x6_q3)),label='$Q^2$=60 $GeV^2$, $x_{\pi}$ = 0.06-0.07')
    xpiscat3 = ax.scatter(cut.applyCuts(xpi,cut_x8_q3),F2pi(cut.applyCuts(xpi,cut_x8_q3),cut.applyCuts(Q2,cut_x8_q3)),label='$Q^2$=60 $GeV^2$, $x_{\pi}$ = 0.08-0.09')
    
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
    
    ax = f.add_subplot(222)
    xpiscat4 = ax.scatter(cut.applyCuts(xpi,cut_x2_q4),F2pi(cut.applyCuts(xpi,cut_x2_q4),cut.applyCuts(Q2,cut_x2_q4)),label='$Q^2$=60 $GeV^2$, $x_{\pi}$ = 0.02-0.03')
    xpiscat4 = ax.scatter(cut.applyCuts(xpi,cut_x4_q4),F2pi(cut.applyCuts(xpi,cut_x4_q4),cut.applyCuts(Q2,cut_x4_q4)),label='$Q^2$=60 $GeV^2$, $x_{\pi}$ = 0.04-0.05')
    xpiscat4 = ax.scatter(cut.applyCuts(xpi,cut_x6_q4),F2pi(cut.applyCuts(xpi,cut_x6_q4),cut.applyCuts(Q2,cut_x6_q4)),label='$Q^2$=60 $GeV^2$, $x_{\pi}$ = 0.06-0.07')
    xpiscat4 = ax.scatter(cut.applyCuts(xpi,cut_x8_q4),F2pi(cut.applyCuts(xpi,cut_x8_q4),cut.applyCuts(Q2,cut_x8_q4)),label='$Q^2$=60 $GeV^2$, $x_{\pi}$ = 0.08-0.09')
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
    xpiscat5 = ax.scatter(cut.applyCuts(xpi,cut_x2_q5),F2pi(cut.applyCuts(xpi,cut_x2_q5),cut.applyCuts(Q2,cut_x2_q5)),label='$Q^2$=60 $GeV^2$, $x_{\pi}$ = 0.02-0.03')
    xpiscat5 = ax.scatter(cut.applyCuts(xpi,cut_x4_q5),F2pi(cut.applyCuts(xpi,cut_x4_q5),cut.applyCuts(Q2,cut_x4_q5)),label='$Q^2$=60 $GeV^2$, $x_{\pi}$ = 0.04-0.05')
    xpiscat5 = ax.scatter(cut.applyCuts(xpi,cut_x6_q5),F2pi(cut.applyCuts(xpi,cut_x6_q5),cut.applyCuts(Q2,cut_x6_q5)),label='$Q^2$=60 $GeV^2$, $x_{\pi}$ = 0.06-0.07')
    xpiscat5 = ax.scatter(cut.applyCuts(xpi,cut_x8_q5),F2pi(cut.applyCuts(xpi,cut_x8_q5),cut.applyCuts(Q2,cut_x8_q5)),label='$Q^2$=60 $GeV^2$, $x_{\pi}$ = 0.08-0.09')
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit",color="y")
    plt.xscale('log')
    plt.ylim(-0.1,0.3)
    plt.xlim(1e-3,1.)
    ax.text(0.25, 0.65, '$Q^2$=240 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.set_yticks([-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3])
    ax.set_xticks([1e-2,1e-1])

    ax = f.add_subplot(224)
    xpiscat6 = ax.scatter(cut.applyCuts(xpi,cut_x2_q6),F2pi(cut.applyCuts(xpi,cut_x2_q6),cut.applyCuts(Q2,cut_x2_q6)),label='$Q^2$=60 $GeV^2$, $x_{\pi}$ = 0.02-0.03')
    xpiscat6 = ax.scatter(cut.applyCuts(xpi,cut_x4_q6),F2pi(cut.applyCuts(xpi,cut_x4_q6),cut.applyCuts(Q2,cut_x4_q6)),label='$Q^2$=60 $GeV^2$, $x_{\pi}$ = 0.04-0.05')
    xpiscat6 = ax.scatter(cut.applyCuts(xpi,cut_x6_q6),F2pi(cut.applyCuts(xpi,cut_x6_q6),cut.applyCuts(Q2,cut_x6_q6)),label='$Q^2$=60 $GeV^2$, $x_{\pi}$ = 0.06-0.07')
    xpiscat6 = ax.scatter(cut.applyCuts(xpi,cut_x8_q6),F2pi(cut.applyCuts(xpi,cut_x8_q6),cut.applyCuts(Q2,cut_x8_q6)),label='$Q^2$=60 $GeV^2$, $x_{\pi}$ = 0.08-0.09')
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit",color="y")
    plt.xscale('log')
    plt.ylim(-0.1,0.3)
    plt.xlim(1e-3,1.)
    ax.text(0.25, 0.65, '$Q^2$=480 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.set_yticks([-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3])
    ax.set_xticks([1e-2,1e-1])

    plt.xlabel('$x_\pi$', fontsize=20)    
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.0,wspace=0.0)

    plt.style.use('default')

fig = plt.figure(figsize=(17,12),facecolor='silver')

ax = fig.add_subplot(331)
#plt.scatter(t,fpi)
densityPlot(t,fpi, '','$t$','$fpi$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(332)
#plt.scatter(xbj,fpi)
densityPlot(xbj,fpi, '','$x$','$fpi$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(333)
densityPlot(xbj,xL, '','$x$','$xL$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(334)
densityPlot(xbj,t, '','$x$','$t$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(335)
densityPlot(xL,t, '','$xL$','$t$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(336)
densityPlot(xbj,Q2, '','$x$','$Q^{2}$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(337)
densityPlot(cut.applyCuts(xbj,cut_q),cut.applyCuts(Q2,cut_q), '','$x$','$Q^{2}$', 200, 200, ax=ax, fig=fig)
#plt.scatter(cut.applyCuts(xbj,cut_t_y),cut.applyCuts(Q2,cut_t_y))

ax = fig.add_subplot(338)
fpivxpi_Plot()

ax = fig.add_subplot(339)


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
Q2 = (np.histogram(Q2_raw, qbins, weights=Q2_raw)[0] / np.histogram(Q2_raw, qbins)[0])
TDIS_xbj = (np.histogram(TDIS_xbj_raw, qbins, weights=Q2_raw)[0] / np.histogram(TDIS_xbj_raw, qbins)[0])
fpi = (np.histogram(fpi_raw, qbins, weights=Q2_raw)[0] / np.histogram(fpi_raw, qbins)[0])
print("\n\n",TDIS_xbj)
print(fpi)
print(Q2,"\n\n")

#plt.scatter(TDIS_xbj,fpi)
#plt.show()

#plt.scatter(TDIS_xbj,Q2)
#plt.show()

# Bins data weighted by TDIS_xbj
Q2 = (np.histogram(Q2_raw, xbins, weights=TDIS_xbj_raw)[0] / np.histogram(Q2_raw, xbins)[0])
TDIS_xbj = (np.histogram(TDIS_xbj_raw, xbins, weights=TDIS_xbj_raw)[0] / np.histogram(TDIS_xbj_raw, xbins)[0])
t = (np.histogram(t_raw, xbins,weights=TDIS_xbj_raw)[0] / np.histogram(t_raw, xbins)[0])
fpi = (np.histogram(fpi_raw, xbins, weights=TDIS_xbj_raw)[0] / np.histogram(fpi_raw, xbins)[0])
print(TDIS_xbj)
print(fpi)
print(t)
print(Q2)
'''