#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2021-10-20 11:52:27 trottar"
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
sys.path.insert(1,'./src/process/cuts/') # Note: this is relative to the bash script NOT this python script!
import cuts as c

kinematics = sys.argv[1]

xbinwidth = float(sys.argv[2])
qbinwidth = float(sys.argv[3])
tbinwidth = float(sys.argv[4])
xLbinwidth = float(sys.argv[5])

#xLselection, tselection = "xLcut8","tcut1"
xLselection, tselection = "",""

df = pd.read_csv(r'./src/process/datafiles/x{0:0.3f}q{1:0.1f}t{2:0.3f}xL{3:0.3f}_{4}.csv'.format(xbinwidth,qbinwidth,tbinwidth,xLbinwidth,kinematics)) # xL bin, no t bin
print(df)

xbj = df['TDIS_xbj']
Q2 = df['TDIS_Q2']
#fpi = df['fpi']
t = df['TDIS_t']
xL = df['xL']
y = df['TDIS_y']
sigma_tdis = df['sigma_tdis']
f2N = df['f2N']
xpi = df['xpi']
xpi2 = xbj/(1.-xL)
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
cutnameDict= { "" : "All xL and t"}

qbinarray = [7.0,15.0,30.0,60.0,120.0,240.0,480.0,1000.0]+np.arange(qbinwidth/2,6,qbinwidth).tolist()
#qbinarray = np.arange(qbinwidth/2,1000.,qbinwidth).tolist()
for i,q in enumerate(qbinarray) :
    qtmp = '{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,qbinarray[i]-qbinwidth/2,qbinarray[i]+qbinwidth/2)
    cutDict.update(eval(qtmp))
    cutnameDict.update({"Q2cut{}".format(i) : "(({0:0.1f} <= Q2) & (Q2 <= {1:0.1f}))".format(qbinarray[i]-qbinwidth/2,qbinarray[i]+qbinwidth/2)})

xarray = np.arange(xbinwidth/2,1.0,xbinwidth).tolist()
for i,x in enumerate(xarray):
    xtmp = '{"xcut%i" : ((%0.4f <= xbj) & (xbj <= %0.4f))}' % (i,xarray[i]-xbinwidth/2,xarray[i]+xbinwidth/2)
    cutDict.update(eval(xtmp))

tarray = np.arange(tbinwidth/2,1.0,tbinwidth).tolist()
for i,tval in enumerate(tarray):
    ttmp = '{"tcut%i" : ((-%0.4f >= t) & (t >= -%0.4f))}' % (i,tarray[i]-tbinwidth/2,tarray[i]+tbinwidth/2)
    cutDict.update(eval(ttmp))
    cutnameDict.update({"tcut{}".format(i) : "((-{0:0.1f} >= t) & (t >= -{1:0.1f}))".format(tarray[i]-tbinwidth/2,tarray[i]+tbinwidth/2)})

xLarray = np.arange(xLbinwidth/2,1.0,xLbinwidth).tolist()
for i,x in enumerate(xLarray):
    xLtmp = '{"xLcut%i" : ((%0.4f <= xL) & (xL <= %0.4f))}' % (i,xLarray[i]-xLbinwidth/2,xLarray[i]+xLbinwidth/2)
    cutDict.update(eval(xLtmp)) 
    cutnameDict.update({"xLcut{}".format(i) : "(({0:0.1f} <= xL) & (xL <= {1:0.1f}))".format(xLarray[i]-xLbinwidth/2,xLarray[i]+xLbinwidth/2)})
    
ytmp = '{"ycut" : ((0.01 <= y) & (y <= 0.95))}'
cutDict.update(eval(ytmp))
cut = c.pyPlot(cutDict)    

cutname_table = pd.DataFrame(cutnameDict, columns=cutnameDict.keys(), index=[0])
cutname_table = cutname_table.reindex(sorted(cutname_table.columns), axis=1)
print("\n\nBin selection...")
print(cutname_table.transpose(),"\n\n")

xLcutName = cutname_table[xLselection].to_string(index=False).replace("<=","$\leq$").replace("xL) & (xL","xL").strip("(").strip(")")
tcutName = cutname_table[tselection].to_string(index=False).replace(">=","$\geq$").replace("t) & (t","t").strip("(").strip(")")

ycut1 = ["ycut"]

if xLselection != "" or tselection != "":
    cut7 = ["Q2cut0","ycut","{}".format(xLselection),"{}".format(tselection)]
    cut15 = ["Q2cut1","ycut","{}".format(xLselection),"{}".format(tselection)]
    cut30 = ["Q2cut2","ycut","{}".format(xLselection),"{}".format(tselection)]
    cut60 = ["Q2cut3","ycut","{}".format(xLselection),"{}".format(tselection)]
    cut120 = ["Q2cut4","ycut","{}".format(xLselection),"{}".format(tselection)]
    cut240 = ["Q2cut5","ycut","{}".format(xLselection),"{}".format(tselection)]
    cut480 = ["Q2cut6","ycut","{}".format(xLselection),"{}".format(tselection)]
    cut1000 = ["Q2cut7","ycut","{}".format(xLselection),"{}".format(tselection)]
else:
    cut1 = ["Q2cut8","ycut"]
    cut2 = ["Q2cut9","ycut"]
    cut3 = ["Q2cut10","ycut"]
    cut4 = ["Q2cut11","ycut"]
    cut5 = ["Q2cut12","ycut"]
    cut6 = ["Q2cut13","ycut"]
    cut7 = ["Q2cut0","ycut"]
    cut15 = ["Q2cut1","ycut"]
    cut30 = ["Q2cut2","ycut"]
    cut60 = ["Q2cut3","ycut"]
    cut120 = ["Q2cut4","ycut"]
    cut240 = ["Q2cut5","ycut"]
    cut480 = ["Q2cut6","ycut"]
    cut1000 = ["Q2cut7","ycut"]

def F2pi(xpi, Q2):
    points,values=np.load('./analysis/interpGrids/xpiQ2.npy'),np.load('./analysis/interpGrids/F2pi.npy')
    F2pi=lambda xpi,Q2: griddata(points,values,(np.log10(xpi),np.log10(Q2)))
    return F2pi(xpi,Q2)

# Calculate cross-section using Patrick's interpolate grid
def ds_dxdQ2dxLdt(x, xL,t):
    points60,values60=np.load('./analysis/xsec/pointsxsec60.npy'),np.load('./analysis/xsec/valuesxsec60.npy')
    points120,values120=np.load('./analysis/xsec/pointsxsec120.npy'),np.load('./analysis/xsec/valuesxsec120.npy')
    points240,values240=np.load('./analysis/xsec/pointsxsec240.npy'),np.load('./analysis/xsec/valuesxsec240.npy')
    points480,values480=np.load('./analysis/xsec/pointsxsec480.npy'),np.load('./analysis/xsec/valuesxsec480.npy')
    sigma60=lambda x,xL,t: griddata(points60,values60,(x,xL,t))
    sigma120=lambda x,xL,t: griddata(points120,values120,(x,xL,t))
    sigma240=lambda x,xL,t: griddata(points240,values240,(x,xL,t))
    sigma480=lambda x,xL,t: griddata(points480,values480,(x,xL,t))

    return [sigma60(x,xL,t),sigma120(x,xL,t),sigma240(x,xL,t),sigma480(x,xL,t)]

fpi = F2pi(xpi, Q2)
dsigma = ds_dxdQ2dxLdt(xpi, xL,t)

def dsigma_Plot():
    
    f = plt.figure(figsize=(11.69,8.27))
    plt.rcParams.update({'font.size': 15})
    plt.style.use('classic')

    ax = f.add_subplot(221)
    xpiscat4 = ax.errorbar(cut.applyCuts(xpi,cut60),cut.applyCuts(dsigma[0],cut60),yerr=np.sqrt(cut.applyCuts(lumi,cut60))/cut.applyCuts(lumi,cut60),fmt='.',label='$Q^2$=60 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    plt.xscale('log')
    #plt.ylim(0.,0.3)
    plt.xlim(1e-2,1.)
    ax.text(0.25, 0.65, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    #ax.set_yticks([0.0,0.1,0.2,0.3])
    ax.set_xticks([1.,1e-1])
    
    plt.ylabel(r'$\frac{d\sigma}{dxdQ^2dx_Ldt}$', fontsize=20)
    
    ax = f.add_subplot(222)
    xpiscat5 = ax.errorbar(cut.applyCuts(xpi,cut120),cut.applyCuts(dsigma[1],cut120),yerr=np.sqrt(cut.applyCuts(lumi,cut120))/cut.applyCuts(lumi,cut120),fmt='.',label='$Q^2$=120 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    plt.xscale('log')
    #plt.ylim(0.,0.3)
    plt.xlim(1e-2,1.)
    ax.text(0.25, 0.65, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    #ax.set_yticks([0.1,0.2,0.3])
    ax.set_xticks([1.,1e-1])

    plt.title("{0}\n{1}".format(xLcutName,tcutName))
    
    ax = f.add_subplot(223)
    xpiscat6 = ax.errorbar(cut.applyCuts(xpi,cut240),cut.applyCuts(dsigma[2],cut240),yerr=np.sqrt(cut.applyCuts(lumi,cut240))/cut.applyCuts(lumi,cut240),fmt='.',label='$Q^2$=240 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    plt.xscale('log')
    #plt.ylim(0.,0.3)
    plt.xlim(1e-2,1.)
    ax.text(0.25, 0.65, '$Q^2$=240 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    #ax.set_yticks([0.0,0.1,0.2,0.3])
    ax.set_xticks([1e-2,1.,1e-1])

    ax = f.add_subplot(224)
    xpiscat7 = ax.errorbar(cut.applyCuts(xpi,cut480),cut.applyCuts(dsigma[3],cut480),yerr=np.sqrt(cut.applyCuts(lumi,cut480))/cut.applyCuts(lumi,cut480),fmt='.',label='$Q^2$=480 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    plt.xscale('log')
    #plt.ylim(0.,0.3)
    plt.xlim(1e-2,1.)
    ax.text(0.25, 0.65, '$Q^2$=480 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    #ax.set_yticks([0.0,0.1,0.2,0.3])
    ax.set_xticks([1.,1e-1])

    plt.xlabel('$x_\pi$', fontsize=20)    
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.0,wspace=0.0)

    plt.style.use('default')


def fpivxpi_Plot():
    
    f = plt.figure(figsize=(11.69,8.27))
    plt.rcParams.update({'font.size': 15})
    plt.style.use('classic')
    
    ax = f.add_subplot(221)
    xpiscat4 = ax.errorbar(cut.applyCuts(xpi,cut60),cut.applyCuts(fpi,cut60),yerr=np.sqrt(cut.applyCuts(lumi,cut60))/cut.applyCuts(lumi,cut60),fmt='.',label='$Q^2$=60 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    plt.plot([0.001,0.01,0.1],[1.2,0.45,0.25], label="GRV fit",color="y")
    plt.xscale('log')
    plt.ylim(0.,0.3)
    plt.xlim(1e-2,1.)
    ax.text(0.25, 0.65, '$Q^2$=60 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.set_yticks([0.0,0.1,0.2,0.3])
    ax.set_xticks([1.,1e-1])
    
    plt.ylabel('$F^{\pi}_{2}$', fontsize=20)
    
    ax = f.add_subplot(222)
    xpiscat5 = ax.errorbar(cut.applyCuts(xpi,cut120),cut.applyCuts(fpi,cut120),yerr=np.sqrt(cut.applyCuts(lumi,cut120))/cut.applyCuts(lumi,cut120),fmt='.',label='$Q^2$=120 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.5,0.25,0.15], label="GRV fit",color="y")
    plt.xscale('log')
    plt.ylim(0.,0.3)
    plt.xlim(1e-2,1.)
    ax.text(0.25, 0.65, '$Q^2$=120 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.set_yticks([0.1,0.2,0.3])
    ax.set_xticks([1.,1e-1])
    
    plt.title("{0}\n{1}".format(xLcutName,tcutName))

    ax = f.add_subplot(223)
    xpiscat6 = ax.errorbar(cut.applyCuts(xpi,cut240),cut.applyCuts(fpi,cut240),yerr=np.sqrt(cut.applyCuts(lumi,cut240))/cut.applyCuts(lumi,cut240),fmt='.',label='$Q^2$=240 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit",color="y")
    plt.xscale('log')
    plt.ylim(0.,0.3)
    plt.xlim(1e-2,1.)
    ax.text(0.25, 0.65, '$Q^2$=240 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_locator(MaxNLocator(prune='both'))
    ax.set_yticks([0.0,0.1,0.2,0.3])
    ax.set_xticks([1e-2,1.,1e-1])

    ax = f.add_subplot(224)
    xpiscat7 = ax.errorbar(cut.applyCuts(xpi,cut480),cut.applyCuts(fpi,cut480),yerr=np.sqrt(cut.applyCuts(lumi,cut480))/cut.applyCuts(lumi,cut480),fmt='.',label='$Q^2$=480 $GeV^2$',ecolor='cyan',capsize=2, capthick=2)
    plt.plot([0.01,0.1,0.3],[0.55,0.25,0.15], label="GRV fit",color="y")
    plt.xscale('log')
    plt.ylim(0.,0.3)
    plt.xlim(1e-2,1.)
    ax.text(0.25, 0.65, '$Q^2$=480 $GeV^2$', transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
    ax.yaxis.set_major_formatter(plt.NullFormatter())
    ax.set_yticks([0.0,0.1,0.2,0.3])
    ax.set_xticks([1.,1e-1])

    plt.xlabel('$x_\pi$', fontsize=20)    
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.0,wspace=0.0)

    plt.style.use('default')

def phaseSpace_Plots():

    fig = plt.figure(figsize=(17,12),facecolor='silver')

    ax = fig.add_subplot(331)
    plt.scatter(cut.applyCuts(t,cut60),cut.applyCuts(fpi,cut60))
    #densityPlot(t,fpi, '','$t$','$fpi$', 200, 200, ax=ax, fig=fig)

    ax = fig.add_subplot(332)
    #densityPlot(xbj,fpi, '','$x$','$fpi$', 200, 200, ax=ax, fig=fig)

    ax = fig.add_subplot(333)
    denplt = densityPlot(xbj,xL, '','$x$','$xL$', 200, 200, ax=ax, fig=fig)

    ax = fig.add_subplot(334)
    densityPlot(xbj,t, '','$x$','$t$', 200, 200, ax=ax, fig=fig)

    ax = fig.add_subplot(335)
    densityPlot(xL,t, '','$xL$','$t$', 200, 200, ax=ax, fig=fig)

    ax = fig.add_subplot(336)
    densityPlot(xbj,Q2, '','$x$','$Q^{2}$', 200, 200, ax=ax, fig=fig)

    ax = fig.add_subplot(337)
    densityPlot(cut.applyCuts(xbj,cut60),cut.applyCuts(t,cut60), '','$x$','t', 200, 200, ax=ax, fig=fig)

    ax = fig.add_subplot(338)
    plt.scatter([np.average(cut.applyCuts(xbj,cut60))],[np.average(cut.applyCuts(xL,cut60))],label='$Q^2$=60 $GeV^2$')
    plt.scatter([np.average(cut.applyCuts(xbj,cut120))],[np.average(cut.applyCuts(xL,cut120))],label='$Q^2$=120 $GeV^2$')
    plt.scatter([np.average(cut.applyCuts(xbj,cut240))],[np.average(cut.applyCuts(xL,cut240))],label='$Q^2$=240 $GeV^2$')
    plt.scatter([np.average(cut.applyCuts(xbj,cut480))],[np.average(cut.applyCuts(xL,cut480))],label='$Q^2$=480 $GeV^2$')
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('xL')
    #plt.title("x = {0:0.3f}, xL = {1:0.3f}".format(np.average(cut.applyCuts(xbj,cut60)),np.average(cut.applyCuts(xL,cut60))))
    #print("~~~~",np.average(cut.applyCuts(xbj,cut60)),np.average(cut.applyCuts(xL,cut60)))

    ax = fig.add_subplot(339)

    plt.tight_layout()

def TheoryTable():

    def dict2df(inp_d):
        out_t = pd.DataFrame(inp_d, columns=inp_d.keys())
        out_t = out_t.reindex(sorted(out_t.columns), axis=1)
        return out_t

    if xLselection == "" or tselection == "":
        cut1Dict = {

            "Q2-1" : cut.applyCuts(Q2,cut1),
            "xbj-1" : cut.applyCuts(xbj,cut1),
            "xpi-1" : cut.applyCuts(xpi,cut1),
            "xL-1" : cut.applyCuts(xL,cut1),
            "y-1" : cut.applyCuts(y,cut1),
            "lumi-1" : cut.applyCuts(lumi,cut1),
        }        

        cut1Table = dict2df(cut1Dict)

        cut2Dict = {

            "Q2-2" : cut.applyCuts(Q2,cut2),
            "xbj-2" : cut.applyCuts(xbj,cut2),
            "xpi-2" : cut.applyCuts(xpi,cut2),
            "xL-2" : cut.applyCuts(xL,cut2),
            "y-2" : cut.applyCuts(y,cut2),
            "lumi-2" : cut.applyCuts(lumi,cut2),
        }        

        cut2Table = dict2df(cut2Dict)

        cut3Dict = {

            "Q2-3" : cut.applyCuts(Q2,cut3),
            "xbj-3" : cut.applyCuts(xbj,cut3),
            "xpi-3" : cut.applyCuts(xpi,cut3),
            "xL-3" : cut.applyCuts(xL,cut3),
            "y-3" : cut.applyCuts(y,cut3),
            "lumi-3" : cut.applyCuts(lumi,cut3),
        }        

        cut3Table = dict2df(cut3Dict)

        cut4Dict = {

            "Q2-4" : cut.applyCuts(Q2,cut4),
            "xbj-4" : cut.applyCuts(xbj,cut4),
            "xpi-4" : cut.applyCuts(xpi,cut4),
            "xL-4" : cut.applyCuts(xL,cut4),
            "y-4" : cut.applyCuts(y,cut4),
            "lumi-4" : cut.applyCuts(lumi,cut4),
        }        

        cut4Table = dict2df(cut4Dict)

        cut5Dict = {

            "Q2-5" : cut.applyCuts(Q2,cut5),
            "xbj-5" : cut.applyCuts(xbj,cut5),
            "xpi-5" : cut.applyCuts(xpi,cut5),
            "xL-5" : cut.applyCuts(xL,cut5),
            "y-5" : cut.applyCuts(y,cut5),
            "lumi-5" : cut.applyCuts(lumi,cut5),
        }        

        cut5Table = dict2df(cut5Dict)

        cut6Dict = {

            "Q2-6" : cut.applyCuts(Q2,cut6),
            "xbj-6" : cut.applyCuts(xbj,cut6),
            "xpi-6" : cut.applyCuts(xpi,cut6),
            "xL-6" : cut.applyCuts(xL,cut6),
            "y-6" : cut.applyCuts(y,cut6),
            "lumi-6" : cut.applyCuts(lumi,cut6),
        }        

        cut6Table = dict2df(cut6Dict)

    cut7Dict = {

        "Q2-7" : cut.applyCuts(Q2,cut7),
        "xbj-7" : cut.applyCuts(xbj,cut7),
        "xpi-7" : cut.applyCuts(xpi,cut7),
        "xL-7" : cut.applyCuts(xL,cut7),
        "y-7" : cut.applyCuts(y,cut7),
        "lumi-7" : cut.applyCuts(lumi,cut7),
    }

    cut7Table = dict2df(cut7Dict)
    
    cut15Dict = {

        "Q2-15" : cut.applyCuts(Q2,cut15),
        "xbj-15" : cut.applyCuts(xbj,cut15),
        "xpi-15" : cut.applyCuts(xpi,cut15),
        "xL-15" : cut.applyCuts(xL,cut15),
        "y-15" : cut.applyCuts(y,cut15),
        "lumi-15" : cut.applyCuts(lumi,cut15),
    }

    cut15Table = dict2df(cut15Dict)
    
    cut30Dict = {

        "Q2-30" : cut.applyCuts(Q2,cut30),
        "xbj-30" : cut.applyCuts(xbj,cut30),
        "xpi-30" : cut.applyCuts(xpi,cut30),
        "xL-30" : cut.applyCuts(xL,cut30),
        "y-30" : cut.applyCuts(y,cut30),
        "lumi-30" : cut.applyCuts(lumi,cut30),
    }

    cut30Table = dict2df(cut30Dict)
    
    cut60Dict = {

        "Q2-60" : cut.applyCuts(Q2,cut60),
        "xbj-60" : cut.applyCuts(xbj,cut60),
        "xpi-60" : cut.applyCuts(xpi,cut60),
        "xL-60" : cut.applyCuts(xL,cut60),
        "y-60" : cut.applyCuts(y,cut60),
        "lumi-60" : cut.applyCuts(lumi,cut60),
    }

    cut60Table = dict2df(cut60Dict)

    cut120Dict = {

        "Q2-120" : cut.applyCuts(Q2,cut120),
        "xbj-120" : cut.applyCuts(xbj,cut120),
        "xpi-120" : cut.applyCuts(xpi,cut120),
        "xL-120" : cut.applyCuts(xL,cut120),
        "y-120" : cut.applyCuts(y,cut120),
        "lumi-120" : cut.applyCuts(lumi,cut120),
    }

    cut120Table = dict2df(cut120Dict)

    cut240Dict = {

        "Q2-240" : cut.applyCuts(Q2,cut240),
        "xbj-240" : cut.applyCuts(xbj,cut240),
        "xpi-240" : cut.applyCuts(xpi,cut240),
        "xL-240" : cut.applyCuts(xL,cut240),
        "y-240" : cut.applyCuts(y,cut240),
        "lumi-240" : cut.applyCuts(lumi,cut240),
    }

    cut240Table = dict2df(cut240Dict)

    cut480Dict = {

        "Q2-480" : cut.applyCuts(Q2,cut480),
        "xbj-480" : cut.applyCuts(xbj,cut480),
        "xpi-480" : cut.applyCuts(xpi,cut480),
        "xL-480" : cut.applyCuts(xL,cut480),
        "y-480" : cut.applyCuts(y,cut480),
        "lumi-480" : cut.applyCuts(lumi,cut480),
    }

    cut480Table = dict2df(cut480Dict)

    cut1000Dict = {

        "Q2-1000" : cut.applyCuts(Q2,cut1000),
        "xbj-1000" : cut.applyCuts(xbj,cut1000),
        "xpi-1000" : cut.applyCuts(xpi,cut1000),
        "xL-1000" : cut.applyCuts(xL,cut1000),
        "y-1000" : cut.applyCuts(y,cut1000),
        "lumi-1000" : cut.applyCuts(lumi,cut1000),
    }

    cut1000Table = dict2df(cut1000Dict)

    if xLselection == "" or tselection == "":
        # Merge pandas df
        dataDict = {}
        for d in (cut1Table,cut2Table,cut3Table,cut4Table,cut5Table,cut6Table,cut7Table,cut15Table,cut30Table,cut60Table,cut120Table,cut240Table,cut480Table,cut1000Table):
            dataDict.update(d)
        data ={i : dataDict[i] for i in sorted(dataDict.keys())}        
        theory_table = dict2df(dataDict).sort_values(['Q2-1','Q2-2','Q2-3','Q2-4','Q2-5','Q2-6','Q2-7','Q2-15','Q2-30','Q2-60','Q2-120','Q2-240','Q2-480','Q2-1000'])
        print("Table created...\n",theory_table)
        print("\n\n",theory_table['xpi-1'].dropna().min())
        print("\n\n",theory_table['Q2-1'].dropna().min())        
    else:
        # Merge pandas df
        dataDict = {}
        for d in (cut7Table,cut15Table,cut30Table,cut60Table,cut120Table,cut240Table,cut480Table,cut1000Table):
            dataDict.update(d)
        data ={i : dataDict[i] for i in sorted(dataDict.keys())}
        theory_table = dict2df(dataDict).sort_values(['Q2-7','Q2-15','Q2-30','Q2-60','Q2-120','Q2-240','Q2-480','Q2-1000'])
        print("Table created...\n",theory_table)

    return theory_table
    
def main() :
    
    #dsigma_Plot()
    #fpivxpi_Plot()
    #phaseSpace_Plots()
    #plt.show()

    out_f = "OUTPUTS/theory_table_{}.csv".format(kinematics)
    theory_table = TheoryTable()
    theory_table.to_csv(out_f, index=False, header=True, mode='w+')
    
    
if __name__=='__main__': main()
