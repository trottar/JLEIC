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

#df = pd.read_csv(r'./datafiles/test.csv')
df = pd.read_csv(r'./datafiles/x0.010q10.0t0.010xL0.010_pi_n_10on135_x0.001-1.000_q1.0-1000.0.csv') # xL bin, no t bin
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
lumi = df['tot_int_lumi']


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
        
    return fig

# Create cut dictionary
cutDict = {}

xbinwidth = 0.01
qbinwidth = 10
tbinwidth = 0.01
xLbinwidth = 0.01

#xlmin,xlmax = 0.8000,0.8100
xlmin,xlmax = 0.8400,0.8500

qbinarray = [7,15,30,60,120,240,480,1000]
#qbinarray = np.arange(qbinwidth/2,1000.,qbinwidth).tolist()
for i,q in enumerate(qbinarray) :
    qtmp = '{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,qbinarray[i]-qbinwidth/2,qbinarray[i]+qbinwidth/2)
    print('{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,qbinarray[i]-qbinwidth/2,qbinarray[i]+qbinwidth/2))
    cutDict.update(eval(qtmp))

xarray = np.arange(xbinwidth/2,1.0,xbinwidth).tolist()
for i,x in enumerate(xarray):
    xtmp = '{"xcut%i" : ((%0.4f <= xbj) & (xbj <= %0.4f))}' % (i,xarray[i]-xbinwidth/2,xarray[i]+xbinwidth/2)
    print('{"xcut%i" : ((%0.4f <= xbj) & (xbj <= %0.4f))}' % (i,xarray[i]-xbinwidth/2,xarray[i]+xbinwidth/2))
    cutDict.update(eval(xtmp))

tarray = np.arange(tbinwidth/2,1.0,tbinwidth).tolist()
for i,tval in enumerate(tarray):
    ttmp = '{"tcut%i" : ((%0.4f <= t) & (t <= %0.4f))}' % (i,tarray[i]-tbinwidth/2,tarray[i]+tbinwidth/2)
    print('{"tcut%i" : ((%0.4f <= t) & (t <= %0.4f))}' % (i,tarray[i]-tbinwidth/2,tarray[i]+tbinwidth/2))
    cutDict.update(eval(ttmp))

xLarray = np.arange(xLbinwidth/2,1.0,xLbinwidth).tolist()
for i,x in enumerate(xLarray):
    xLtmp = '{"xLcut%i" : ((%0.4f <= xL) & (xL <= %0.4f))}' % (i,xLarray[i]-xLbinwidth/2,xLarray[i]+xLbinwidth/2)
    print('{"xLcut%i" : ((%0.4f <= xL) & (xL <= %0.4f))}' % (i,xLarray[i]-xLbinwidth/2,xLarray[i]+xLbinwidth/2))
    cutDict.update(eval(xLtmp))

ytmp = '{"ycut" : ((0.01 <= y) & (y <= 0.95))}'
cutDict.update(eval(ytmp))
cut = c.pyPlot(cutDict)    

ycut1 = ["ycut"]

#cut_q = ["Q2cut3","xcut2","ycut"] # Q2= 60 GeV^2
#cut_q = ["Q2cut3","ycut"] # Q2= 60 GeV^2
cut_q = ["xLcut90","xcut1","ycut"] # Q2= 60 GeV^2

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

cut7 = ["Q2cut0","ycut"]
cut15 = ["Q2cut1","ycut"]
cut30 = ["Q2cut2","ycut"]
'''
cut60 = ["Q2cut3","xLcut80","ycut"]
cut120 = ["Q2cut4","xLcut80","ycut"]
cut240 = ["Q2cut5","xLcut80","ycut"]
cut480 = ["Q2cut6","xLcut80","ycut"]
cut1000 = ["Q2cut7","xLcut80","ycut"]
'''
cut60 = ["Q2cut3","xLcut85","ycut"]
cut120 = ["Q2cut4","xLcut85","ycut"]
cut240 = ["Q2cut5","xLcut85","ycut"]
cut480 = ["Q2cut6","xLcut85","ycut"]
cut1000 = ["Q2cut7","xLcut85","ycut"]
def F2pi(xpi, Q2):
    points,values=np.load('./../../analysis/interpGrids/xpiQ2.npy'),np.load('./../../analysis/interpGrids/F2pi.npy')
    F2pi=lambda xpi,Q2: griddata(points,values,(np.log10(xpi),np.log10(Q2)))
    return F2pi(xpi,Q2)

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

def fpivxpi_Plot():
    
    f = plt.figure(figsize=(11.69,8.27))
    plt.rcParams.update({'font.size': 15})
    plt.style.use('classic')
    
    ax = f.add_subplot(221)
    xpiscat4 = ax.errorbar(cut.applyCuts(xpi,cut60),cut.applyCuts(fpi,cut60),yerr=np.sqrt(cut.applyCuts(lumi,cut60))/cut.applyCuts(lumi,cut60),fmt='.',label='$Q^2$=60 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
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

    plt.xlabel('$x_\pi$', fontsize=20)    
    plt.title("{0} $\leq$ xL $\leq$ {1}".format(xlmin,xlmax))
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.0,wspace=0.0)

    plt.style.use('default')


fig = plt.figure(figsize=(17,12),facecolor='silver')

ax = fig.add_subplot(331)
#plt.scatter(t,fpi)
densityPlot(t,fpi, '','$t$','$fpi$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(332)
densityPlot(xbj,fpi, '','$x$','$fpi$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(333)
denplt = densityPlot(cut.applyCuts(xbj,cut60),cut.applyCuts(xL,cut60), '','$x$','$xL$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(334)
densityPlot(xbj,t, '','$x$','$t$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(335)
densityPlot(xL,t, '','$xL$','$t$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(336)
densityPlot(xbj,Q2, '','$x$','$Q^{2}$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(337)
densityPlot(cut.applyCuts(xbj,cut60),cut.applyCuts(Q2,cut60), '','$x$','$Q^{2}$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(338)
fpivxpi_Plot()

ax = fig.add_subplot(339)
plt.plot(np.average(cut.applyCuts(xbj,cut60)),np.average(cut.applyCuts(xL,cut60)))
plt.title("{0}{1}".format.(np.average(cut.applyCuts(xbj,cut60)),np.average(cut.applyCuts(xL,cut60))))
print("~~~~",np.average(cut.applyCuts(xbj,cut60)),np.average(cut.applyCuts(xL,cut60)))


plt.tight_layout()
plt.show()