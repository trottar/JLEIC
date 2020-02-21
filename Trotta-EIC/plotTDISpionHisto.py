#! /usr/bin/python

#
# Description:This will read in the array data file that contains all the leave histogram information
# ================================================================
# Time-stamp: 2019-04-08 08:07:33 trottar
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

from __future__ import division
import logging

# Gets rid of matplot logging DEBUG messages
plt_logger = logging.getLogger('matplotlib')
plt_logger.setLevel(logging.WARNING)

# Suppresses unwanted numpy warning
import warnings
import numpy as np
warnings.simplefilter(action='ignore', category=FutureWarning)

import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
from matplotlib import interactive
from matplotlib import colors
import uproot as up
from sys import path
import time,math,sys
# np.set_printoptions(threshold=sys.maxsize)

# My class function
sys.path.insert(0,'/home/trottar/bin/python/root2py/')
from root2py import pyPlot, pyBranch

# rootName = "TDISpion_80k"
rootName_low = "/home/trottar/ResearchNP/JLEIC/Trotta-EIC/TDISpion_low.root"
rootName_mid = "/home/trottar/ResearchNP/JLEIC/Trotta-EIC/TDISpion_med.root"
rootName_high = "/home/trottar/ResearchNP/JLEIC/Trotta-EIC/TDISpion_high.root"
rootName = "/home/trottar/ResearchNP/JLEIC/Trotta-EIC/TDISpion.root"

# rootName = "chain_20k"
# rootName = "chain_low"

tree = up.open(rootName)["Evnts"]
tree_low = up.open(rootName_low)["Evnts"]
tree_mid = up.open(rootName_mid)["Evnts"]
tree_high = up.open(rootName_high)["Evnts"]

branch = pyBranch(tree)
branch_low = pyBranch(tree_low)
branch_mid = pyBranch(tree_mid)
branch_high = pyBranch(tree_high)

pdf = matplotlib.backends.backend_pdf.PdfPages("%s.pdf" % rootName)
pdf_q2 = matplotlib.backends.backend_pdf.PdfPages("Q2vxbj-xpi.pdf")
pdf_sigma = matplotlib.backends.backend_pdf.PdfPages("sigmavxbj-xpi.pdf")

Q2_low = branch_low.findBranch("invts","Q2")
Q2_mid = branch_mid.findBranch("invts","Q2")
Q2_high = branch_high.findBranch("invts","Q2")

TDIS_xbj_low = tree_low.array("TDIS_xbj")
TDIS_xbj_mid = tree_mid.array("TDIS_xbj")
TDIS_xbj_high = tree_high.array("TDIS_xbj")

xpi_low = tree_low.array("xpi")
xpi_mid = tree_mid.array("xpi")
xpi_high = tree_high.array("xpi")

sigma_dis_low = tree_low.array("sigma_dis")*(1e-6)
sigma_dis_mid = tree_mid.array("sigma_dis")*(1e-6)
sigma_dis_high = tree_high.array("sigma_dis")*(1e-6)

y_low = tree_low.array("y")
y_mid = tree_mid.array("y")
y_high = tree_high.array("y")

yplus_low = 1+((1-y_low)*(1-y_low))
yplus_mid = 1+((1-y_mid)*(1-y_mid))
yplus_high = 1+((1-y_high)*(1-y_high))

tot_sigma_low = (sigma_dis_low)*((TDIS_xbj_low*(Q2_low*Q2_low)*(137)*(137))/(2*math.pi*yplus_low))
tot_sigma_mid = (sigma_dis_mid)*((TDIS_xbj_mid*(Q2_mid*Q2_mid)*(137)*(137))/(2*math.pi*yplus_mid))
tot_sigma_high = (sigma_dis_high)*((TDIS_xbj_high*(Q2_high*Q2_high)*(137)*(137))/(2*math.pi*yplus_high))

# Define phyisics data
s_e = branch.findBranch("invts","s_e")
s_q = branch.findBranch("invts","s_q")
Q2 = branch.findBranch("invts","Q2")
xBj = branch.findBranch("invts","xBj")
t = -branch.findBranch("invts","tPrime")
y_D = branch.findBranch("invts","y_D")
nu = branch.findBranch("invts","nu")
TwoPdotk = branch.findBranch("invts","TwoPdotk")
TDIS_xbj = tree.array("TDIS_xbj")
sigma_dis = tree.array("sigma_dis")*(1e-5)
TDIS_y = tree.array("TDIS_y")
ppix_Lab = tree.array("ppix_Lab") # Long pion momentum
ppiy_Lab = tree.array("ppiy_Lab")
ppiz_Lab = tree.array("ppiz_Lab")
EpiE_Lab = tree.array("EpiE_Lab")
pprx_inc = tree.array("pprx_inc") # Long proton momentum
ppry_inc = tree.array("ppry_inc")
pprz_inc = tree.array("pprz_inc")
EprE_inc = tree.array("EprE_inc")
y = tree.array("y")
fpi = tree.array("fpi")
f2N = tree.array("f2N")
xpi = tree.array("xpi")
ypi = tree.array("ypi")
tpi = tree.array("tpi")
escat = tree.array("EScatRest")

# Q2_new = s_q/(xBj*y)
Q2_new = xBj/TwoPdotk
piQ2 = s_q/(xpi*ypi)
pNz_Lab = pprz_inc-ppiz_Lab # Long neutron momentum
ppiz_frac = ppiz_Lab/pprz_inc # Frac of proton momentum
pNz_frac = pNz_Lab/pprz_inc # Frac of proton momentum
yplus = 1+((1-y)*(1-y))
alt_sigma_dis = (1e1)*sigma_dis
tot_sigma = (sigma_dis)*((TDIS_xbj*(Q2*Q2)*(137)*(137))/(2*math.pi*yplus))
red_sig  = np.array([1.397,1.258,1.137,1.055,0.944,0.836,0.706,0.519])
red_y = np.array([0.657,0.416,0.263,0.163,0.103,0.066,0.033,0.008])
red_x = np.array([2.53e-4,4.0e-4,6.32e-4,1.02e-3,1.61e-3,2.52e-3,5.0e-3,2.10e-2])
red_Q2 = np.array([15,15,15,15,15,15,15,15])
red_yplus = 1+(1-red_y)*(1-red_y)
my_sigma = (red_sig)*((2*math.pi*red_yplus)/(red_x*(red_Q2*red_Q2)*(137)*(137)))

# arrPlot = arrPlot[(arrCut > low) & (arrCut < high)]
# To call item, cutDict.get(key,"Leaf name not found")
cutDict = {
    "xcut1" : ((1.03e-4 < TDIS_xbj) & (TDIS_xbj < 2.1e-2)),
    "xcut2" : ((2.53e-4 < TDIS_xbj) & (TDIS_xbj < 2.1e-2)),
    "xcut3" : ((8.0e-4 < TDIS_xbj) & (TDIS_xbj < 3.2e-2)),
    "xcut4" : ((1.3e-3 < TDIS_xbj) & (TDIS_xbj < 5.0e-2)),
    "xcut5" : ((2.1e-3 < TDIS_xbj) & (TDIS_xbj < 1.8e-1)),
    "xcut6" : ((5.0e-3 < TDIS_xbj) & (TDIS_xbj < 1.8e-1)),
    "Q2cut1" : ((6.499 < Q2) & (Q2 < 6.501)),
    "Q2cut2" : ((14.999 < Q2) & (Q2 < 15.001)),
    "Q2cut3" : ((34.999 < Q2) & (Q2 < 35.001)),
    "Q2cut4" : ((59.999 < Q2) & (Q2 < 60.001)),
    "Q2cut5" : ((119.999 < Q2) & (Q2 < 120.001)),
    "Q2cut6" : ((249.999 < Q2) & (Q2 < 250.001)),
    "Q2cutI" : ((6.499 < Q2) & (Q2 < 6.501)),
    "Q2cutII" : ((3.999 < Q2) & (Q2 < 4.001)),
    "Q2cutIII" : ((7.499 < Q2) & (Q2 < 7.501)),

}

lowDict = {
    "Q2low1" : ((6.499 < Q2_low) & (Q2_low < 6.501)),
    "Q2low2" : ((14.999 < Q2_low) & (Q2_low < 15.001)),
    "Q2low3" : ((34.999 < Q2_low) & (Q2_low < 35.001)),
    "Q2low4" : ((59.999 < Q2_low) & (Q2_low < 60.001)),
    "Q2low5" : ((119.999 < Q2_low) & (Q2_low < 120.001)),
    "Q2low6" : ((249.999 < Q2_low) & (Q2_low < 250.001)),

}

midDict = {
    "Q2mid1" : ((6.499 < Q2_mid) & (Q2_mid < 6.501)),
    "Q2mid2" : ((14.999 < Q2_mid) & (Q2_mid < 15.001)),
    "Q2mid3" : ((34.999 < Q2_mid) & (Q2_mid < 35.001)),
    "Q2mid4" : ((59.999 < Q2_mid) & (Q2_mid < 60.001)),
    "Q2mid5" : ((119.999 < Q2_mid) & (Q2_mid < 120.001)),
    "Q2mid6" : ((249.999 < Q2_mid) & (Q2_mid < 250.001)),

}

highDict = {
    "Q2high1" : ((6.499 < Q2_high) & (Q2_high < 6.501)),
    "Q2high2" : ((14.999 < Q2_high) & (Q2_high < 15.001)),
    "Q2high3" : ((34.999 < Q2_high) & (Q2_high < 35.001)),
    "Q2high4" : ((59.999 < Q2_high) & (Q2_high < 60.001)),
    "Q2high5" : ((119.999 < Q2_high) & (Q2_high < 120.001)),
    "Q2high6" : ((249.999 < Q2_high) & (Q2_high < 250.001)),

}

# Please find below the bins to make for Q2, xpi. I also include y here. For this, please run the simulation for 5 GeV electrons and 100 GeV protons. This should give a value for s= 4x 5 x 100 = 2000 GeV^2. Note that s is related to Q2 as s = Q2 x (x_B)/y. We will assume an implicit cut on y < 0.7 for now, which means Q2 x (xB) < 1400. We can relax that later. 

xarray   =[0.001,0.002,0.002,0.004,0.004,0.004,0.004,0.006,0.006,0.006,0.006,0.008,0.008,0.008,0.008,0.008]
Q2array  =[1.0,1.0,2.0,1.0,2.0,3.0,4.0,1.5,3.0,4.5,6.0,1.0,2.0,4.0,6.0,8.0]
yarray   =[0.8,0.4,0.8,0.2,0.4,0.6,0.8,0.2,0.4,0.6,0.8,0.1,0.2,0.4,0.6,0.8]

# xarray   =[0.01,0.01,0.01,0.01,0.01,0.02,0.02,0.02,0.02,0.02,0.04,0.04,0.04,0.04,0.04,0.06,0.06,0.06,0.06,0.06,0.08,0.08,0.08,0.08,0.08]
# Q2array  =[1.25,2.5,5.0,7.5,10.0,2.5,5.0,10.0,15.0,20.0,5.0,10.0,20.0,30.0,40.0,7.5,15.0,30.0,45.0,60.0,10.0,20.0,40.0,60.0,80.0]
# yarray   =[0.1,0.2,0.4,0.6,0.8,0.1,0.2,0.4,0.6,0.8,0.1,0.2,0.4,0.6,0.8,0.1,0.2,0.4,0.6,0.8,0.1,0.2,0.4,0.6,0.8]

# xarray   =[0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,0.6,0.6,0.6,0.6,0.8,0.8,0.8,0.8,0.8]
# Q2array  =[12.5,25.0,50.0,75.0,100.0,25.0,50.0,100.0,150.0,200.0,50.0,100.0,200.0,300.0,400.0,75.0,150.0,300.0,450.0,600.0,100.0,200.0,400.0,600.0,800.0]
# yarray   =[0.1,0.2,0.4,0.6,0.8,0.1,0.2,0.4,0.6,0.8,0.1,0.2,0.4,0.6,0.8,0.1,0.2,0.4,0.6,0.8,0.1,0.2,0.4,0.6,0.8]

cutBinDict = {}

i=0
for x in range(0,len(xarray)) :
    if i < (len(xarray)-1):
        xtmp = '{"xcut%ia" : ((%0.4f  < TDIS_xbj) & (TDIS_xbj < %0.4f))}' % (i,xarray[i]-0.0001,xarray[i]+0.0001)
        ytmp = '{"ycut%ia" : ((%0.1f < y) & (y < %0.1f))}' % (i,yarray[i]-0.1,yarray[i]+0.1)
        Q2tmp = '{"Q2cut%ia" : ((%0.1f < Q2) & (Q2 < %0.1f))}' % (i,Q2array[i]-0.1,Q2array[i]+0.1)
    else:
        xtmp = '{"xcut%ia" : ((%0.4f  < TDIS_xbj) & (TDIS_xbj < %0.4f))}' % (i,xarray[i]-0.0001,xarray[i]+0.0001)
        ytmp = '{"ycut%ia" : ((%0.1f < y) & (y < %0.1f))}' % (i,yarray[i]-0.1,yarray[i]+0.1)
        Q2tmp = '{"Q2cut%ia" : ((%0.1f < Q2) & (Q2 < %0.1f))}' % (i,Q2array[i]-0.1,Q2array[i]+0.1)
    cutBinDict.update(eval(xtmp))
    cutBinDict.update(eval(ytmp))
    cutBinDict.update(eval(Q2tmp))
    i+=1

c = pyPlot(cutDict)

clow = pyPlot(lowDict)
cmid = pyPlot(midDict)
chigh = pyPlot(highDict)

cbin = pyPlot(cutBinDict)

def uncern_Cut():

    low1 = ["Q2low1"]
    mid1 = ["Q2mid1"]
    high1 = ["Q2high1"]
    
    low2 = ["Q2low2"]
    mid2 = ["Q2mid2"]
    high2 = ["Q2high2"]
    
    low3 = ["Q2low3"]
    mid3 = ["Q2mid3"]
    high3 = ["Q2high3"]
    
    low4 = ["Q2low4"]
    mid4 = ["Q2mid4"]
    high4 = ["Q2high4"]
    
    low5 = ["Q2low5"]
    mid5 = ["Q2mid5"]
    high5 = ["Q2high5"]
    
    low6 = ["Q2low6"]
    mid6 = ["Q2mid6"]
    high6 = ["Q2high6"]


    return[low1,mid1,high1,low2,mid2,high2,low3,mid3,high3,low4,mid4,high4,low5,mid5,high5,low6,mid6,high6]
    
def sigmaDIS_Cut():
    
    cuts1 = ["xcut1", "Q2cut1"]
    cuts2 = ["xcut2", "Q2cut2"]
    cuts3 = ["xcut3", "Q2cut3"]
    cuts4 = ["xcut4", "Q2cut4"]
    cuts5 = ["xcut5", "Q2cut5"]
    cuts6 = ["xcut6", "Q2cut6"]
    cutsQ2 = ["Q2cutI"]
    cutsQ2a = ["Q2cutII"]
    cutsQ2b = ["Q2cutIII"]


    return[cuts1,cuts2,cuts3,cuts4,cuts5,cuts6,cutsQ2,cutsQ2a,cutsQ2b]

def sigmaBin_Cut():

    cuts1 = ["xcut1a", "ycut1a", "Q2cut1a"]
    cuts2 = ["xcut2a", "ycut2a", "Q2cut2a"]
    cuts3 = ["xcut3a", "ycut3a", "Q2cut3a"]
    cuts4 = ["xcut4a", "ycut4a", "Q2cut4a"]
    cuts5 = ["xcut5a", "ycut5a", "Q2cut5a"]
    cuts6 = ["xcut6a", "ycut6a", "Q2cut6a"]
    cuts7 = ["xcut7a", "ycut7a", "Q2cut7a"]
    cuts8 = ["xcut8a", "ycut8a", "Q2cut8a"]
    cuts9 = ["xcut9a", "ycut9a", "Q2cut9a"]
    cuts10 = ["xcut10a", "ycut10a", "Q2cut10a"]
    cuts11 = ["xcut11a", "ycut11a", "Q2cut11a"]
    cuts12 = ["xcut12a", "ycut12a", "Q2cut12a"]
    cuts13 = ["xcut13a", "ycut13a", "Q2cut13a"]
    cuts14 = ["xcut14a", "ycut14a", "Q2cut14a"]
    cuts15 = ["xcut15a", "ycut15a", "Q2cut15a"]
                                                    

    return[cuts1,cuts2,cuts3,cuts4,cuts5,cuts6,cuts7,cuts8,cuts9,cuts10,cuts11,cuts12,cuts13,cuts14,cuts15]

def densityPlot(x,y,title,xlabel,ylabel,binx,biny,xmin=None,xmax=None,ymin=None,ymax=None,cuts=None,figure=None,ax=None,layered=True):

    if cuts:
        xcut  = c.applyCuts(x,cuts)
        ycut = c.applyCuts(y,cuts)
    else:
        xcut = x
        ycut = y
    # if sub or figure:
    #     ax = figure.add_subplot(sub)
    if ax or figure:
        # ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
        print ""
    else:
        fig, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    if (xmin or xmax or ymin or ymax):
        # norm=colors.LogNorm() makes colorbar normed and logarithmic
        hist = ax.hist2d(xcut, ycut,bins=(p.setbin(x,binx,xmin,xmax),p.setbin(y,biny,ymin,ymax)))
    else:
        # norm=colors.LogNorm() makes colorbar normed and logarithmic
        hist = ax.hist2d(xcut, ycut,bins=(p.setbin(x,binx),p.setbin(y,biny)))
    if layered is True :    
        plt.colorbar(hist[3], ax=ax, spacing='proportional', label='Number of Events')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    
def sigmaDIS_Plot():

    [cuts1,cuts2,cuts3,cuts4,cuts5,cuts6,cutsQ2,cutsQ2a,cutsQ2b] = sigmaDIS_Cut()

    # print 'Q2 \n:', Q2
    # print 'TwoPdotk \n:', TwoPdotk
    # print 'TDIS_xbj \n:', TDIS_xbj, np.average(TDIS_xbj)
    # print 'xBj \n:', xBj, np.average(xBj)
    # print 'TDIS_xbj*TwoPdotk \n:', TDIS_xbj*TwoPdotk
    # print 'y \n:', y, np.average(y)
    # print 'y_D \n:', y_D, np.average(y_D)

    # print 'tot_sigma \n:', tot_sigma, np.average(tot_sigma)

    f = plt.figure(figsize=(11.69,8.27))
    f.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.7, wspace=0.35, hspace=0.45)
    ax = f.add_subplot(331)
    hsigma_dis = ax.hist(sigma_dis,bins=p.setbin(sigma_dis,200,0.,10.0),histtype='step', alpha=0.5, stacked=True, fill=True )
    plt.title("d$\sigma_{dis}$")
    ax = f.add_subplot(332)
    totsigplot = ax.hist(tot_sigma,bins=p.setbin(tot_sigma,200,0.,100.0),histtype='step', alpha=0.5, stacked=True, fill=True);
    plt.title('Reduced $\sigma_{dis}$',fontsize =16)
    ax = f.add_subplot(333)
    f2Nplot = ax.hist(f2N,bins=p.setbin(f2N,200,0.,100.0),histtype='step', alpha=0.5, stacked=True, fill=True);
    plt.title('$f^{2}_{N}$',fontsize =16)
    ax = f.add_subplot(334)
    hxbj = ax.hist(TDIS_xbj,bins=p.setbin(TDIS_xbj,200,0.,2.0),histtype='step', alpha=0.5, stacked=True, fill=True )
    plt.title("TDIS_xbj")
    ax = f.add_subplot(335)
    hy = ax.hist(y,bins=p.setbin(y,200,0.,2.0),histtype='step', alpha=0.5, stacked=True, fill=True )
    plt.title("y")
    ax = f.add_subplot(336)
    h1mmiss = ax.hist(yplus,bins=p.setbin(yplus,200,0.002,2.0),histtype='step', alpha=0.5, stacked=True, fill=True )
    plt.title("$Y_+$")

    f = plt.figure(figsize=(11.69,8.27))
    f.tight_layout()
    hQ2_sigma_dis = densityPlot(Q2, sigma_dis, 'd$\sigma_{dis}$ vs $Q^2$','$Q^2$','d$\sigma_{dis}$', 200, 200, 0., 350., 0., 1e-4, figure=f, sub=131)
    hxbj_sigma_dis = densityPlot(TDIS_xbj, sigma_dis, 'd$\sigma_{dis}$ vs TDIS_xbj','TDIS_xbj','d$\sigma_{dis}$', 200, 200, 0., 1., 0., 1e-4, figure=f, sub=132)
    hy_sigma_dis = densityPlot(y, sigma_dis, 'd$\sigma_{dis}$ vs y','y','d$\sigma_{dis}$', 200, 200, 0., 1., 0., 1e-4, figure=f, sub=133)
    
    hxbj_Q2 = densityPlot(TDIS_xbj, Q2, '$Q^{2}$ vs TDIS_xbj','TDIS_xbj','$Q^{2}$', 200, 200, 0., 1., 0., 350.)

    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    mysigscat = ax.scatter(red_x,my_sigma,label='Zeus')
    sigscat2 = ax.scatter(c.applyCuts(TDIS_xbj,cuts2),c.applyCuts(alt_sigma_dis,cuts2),label='$Q^2$=15.0 $GeV^2$')
    plt.title('d$\sigma_{DIS}$ vs TDIS_xbj', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.2,0.3), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('d$\sigma_{DIS}$ (fb)')
    plt.xlim(1e-4,1)
    plt.xscale('log')
    plt.yscale('log')
    
    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    sigscat1 = ax.scatter(c.applyCuts(TDIS_xbj,cuts1),c.applyCuts(sigma_dis,cuts1),label='$Q^2$=6.5 $GeV^2$')
    sigscat2 = ax.scatter(c.applyCuts(TDIS_xbj,cuts2),c.applyCuts(sigma_dis,cuts2),label='$Q^2$=15.0 $GeV^2$')
    sigscat3 = ax.scatter(c.applyCuts(TDIS_xbj,cuts3),c.applyCuts(sigma_dis,cuts3),label='$Q^2$=35.0 $GeV^2$')
    sigscat4 = ax.scatter(c.applyCuts(TDIS_xbj,cuts4),c.applyCuts(sigma_dis,cuts4),label='$Q^2$=60.0 $GeV^2$')
    sigscat5 = ax.scatter(c.applyCuts(TDIS_xbj,cuts5),c.applyCuts(sigma_dis,cuts5),label='$Q^2$=120.0 $GeV^2$')
    sigscat6 = ax.scatter(c.applyCuts(TDIS_xbj,cuts6),c.applyCuts(sigma_dis,cuts6),label='$Q^2$=250.0 $GeV^2$')
    plt.title('d$\sigma_{DIS}$ vs TDIS_xbj', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.2,0.3), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('d$\sigma_{DIS}$ (nb)')
    plt.xlim(1e-4,1)
    plt.xscale('log')
    # plt.yscale('log')
    print "sigma_dis Q^2=6.5:\n",sigscat1.get_offsets()
    print "sigma_dis Q^2=15.0:\n",sigscat2.get_offsets()
    print "sigma_dis Q^2=35.0:\n",sigscat3.get_offsets()
    print "sigma_dis Q^2=60.0:\n",sigscat4.get_offsets()
    print "sigma_dis Q^2=120.0:\n",sigscat5.get_offsets()
    print "sigma_dis Q^2=250.0:\n",sigscat6.get_offsets()

    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    f2Nscat1 = ax.scatter(c.applyCuts(TDIS_xbj,cuts1),c.applyCuts(f2N,cuts1),label='$Q^2$=6.5 $GeV^2$')
    f2Nscat2 = ax.scatter(c.applyCuts(TDIS_xbj,cuts2),c.applyCuts(f2N,cuts2),label='$Q^2$=15.0 $GeV^2$')
    f2Nscat3 = ax.scatter(c.applyCuts(TDIS_xbj,cuts3),c.applyCuts(f2N,cuts3),label='$Q^2$=35.0 $GeV^2$')
    f2Nscat4 = ax.scatter(c.applyCuts(TDIS_xbj,cuts4),c.applyCuts(f2N,cuts4),label='$Q^2$=60.0 $GeV^2$')
    f2Nscat5 = ax.scatter(c.applyCuts(TDIS_xbj,cuts5),c.applyCuts(f2N,cuts5),label='$Q^2$=120.0 $GeV^2$')
    f2Nscat6 = ax.scatter(c.applyCuts(TDIS_xbj,cuts6),c.applyCuts(f2N,cuts6),label='$Q^2$=250.0 $GeV^2$')
    plt.title('$F^{N}_{2}$ vs TDIS_xbj', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.2,0.3), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('$F^{N}_{2}$')
    plt.xlim(1e-4,1)
    plt.xscale('log')
    # plt.ylim(0.2,0.6)
    # plt.yscale('log')
    print "f2N Q^2=6.5:\n",f2Nscat1.get_offsets()
    print "f2N Q^2=15.0:\n",f2Nscat2.get_offsets()
    print "f2N Q^2=35.0:\n",f2Nscat3.get_offsets()
    print "f2N Q^2=60.0:\n",f2Nscat4.get_offsets()
    print "f2N Q^2=120.0:\n",f2Nscat5.get_offsets()
    print "f2N Q^2=250.0:\n",f2Nscat6.get_offsets()

    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    # totsigscat = ax.scatter(TDIS_xbj,tot_sigma,label='No cuts')
    totsigscat1 = ax.scatter(c.applyCuts(TDIS_xbj,cuts1),c.applyCuts(tot_sigma,cuts1),label='$Q^2$=6.5 $GeV^2$')
    totsigscat2 = ax.scatter(c.applyCuts(TDIS_xbj,cuts2),c.applyCuts(tot_sigma,cuts2),label='$Q^2$=15.0 $GeV^2$')
    totsigscat3 = ax.scatter(c.applyCuts(TDIS_xbj,cuts3),c.applyCuts(tot_sigma,cuts3),label='$Q^2$=35.0 $GeV^2$')
    totsigscat4 = ax.scatter(c.applyCuts(TDIS_xbj,cuts4),c.applyCuts(tot_sigma,cuts4),label='$Q^2$=60.0 $GeV^2$')
    totsigscat5 = ax.scatter(c.applyCuts(TDIS_xbj,cuts5),c.applyCuts(tot_sigma,cuts5),label='$Q^2$=120.0 $GeV^2$')
    totsigscat6 = ax.scatter(c.applyCuts(TDIS_xbj,cuts6),c.applyCuts(tot_sigma,cuts6),label='$Q^2$=250.0 $GeV^2$')
    plt.title('Reduced $\sigma_{DIS}$ vs TDIS_xbj', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.2,0.3), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('Reduced $\sigma_{DIS}$ ($nb^{-1}$)')
    plt.xlim(1e-4,1)
    plt.xscale('log')
    # plt.ylim(0.,15e-4)
    # plt.yscale('log')
    print "tot_sigma Q^2=6.5:\n",totsigscat1.get_offsets()
    print "tot_sigma Q^2=15.0:\n",totsigscat2.get_offsets()
    print "tot_sigma Q^2=35.0:\n",totsigscat3.get_offsets()
    print "tot_sigma Q^2=60.0:\n",totsigscat4.get_offsets()
    print "tot_sigma Q^2=120.0:\n",totsigscat5.get_offsets()
    print "tot_sigma Q^2=250.0:\n",totsigscat6.get_offsets()
    
    # f, ax = plt.subplots()
    # thist1 = ax.hist(c.applyCuts(t,cuts1),bins=p.setbin(t,100,0.,1.),histtype='step', alpha=0.5, stacked=True, fill=True,label='$Q^2$=6.5 $GeV^2$' )
    # thist2 = ax.hist(c.applyCuts(t,cuts2),bins=p.setbin(t,100,0.,1.),histtype='step', alpha=0.5, stacked=True, fill=True,label='$Q^2$=15.0 $GeV^2$' )
    # thist3 = ax.hist(c.applyCuts(t,cuts3),bins=p.setbin(t,100,0.,1.),histtype='step', alpha=0.5, stacked=True, fill=True,label='$Q^2$=35.0 $GeV^2$' )
    # plt.title("t cuts", fontsize =16)
    # leg = plt.legend(bbox_to_anchor=(0.6,0.5), loc="center right")
    # leg.get_frame().set_alpha(1.)
    
    for f in xrange(1, plt.figure().number):
        pdf.savefig(f)
    pdf.close()

def sigmaBin_Plot():

    [cuts1,cuts2,cuts3,cuts4,cuts5,cuts6,cuts7,cuts8,cuts9,cuts10,cuts11,cuts12,cuts13,cuts14,cuts15] = sigmaBin_Cut()

    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    sigbin1 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts1),cbin.applyCuts(sigma_dis,cuts1),label='cut1')
    sigbin2 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts2),cbin.applyCuts(sigma_dis,cuts2),label='cut2')
    sigbin3 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts3),cbin.applyCuts(sigma_dis,cuts3),label='cut3')
    sigbin4 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts4),cbin.applyCuts(sigma_dis,cuts4),label='cut4')
    sigbin5 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts5),cbin.applyCuts(sigma_dis,cuts5),label='cut5')
    sigbin6 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts6),cbin.applyCuts(sigma_dis,cuts6),label='cut6')
    sigbin7 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts7),cbin.applyCuts(sigma_dis,cuts7),label='cut7')
    sigbin8 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts8),cbin.applyCuts(sigma_dis,cuts8),label='cut8')
    sigbin9 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts9),cbin.applyCuts(sigma_dis,cuts9),label='cut9')
    sigbin10 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts10),cbin.applyCuts(sigma_dis,cuts10),label='cut10')
    sigbin11 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts11),cbin.applyCuts(sigma_dis,cuts11),label='cut11')
    sigbin12 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts12),cbin.applyCuts(sigma_dis,cuts12),label='cut12')
    sigbin13 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts13),cbin.applyCuts(sigma_dis,cuts13),label='cut13')
    sigbin14 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts14),cbin.applyCuts(sigma_dis,cuts14),label='cut14')
    sigbin15 = ax.scatter(cbin.applyCuts(TDIS_xbj,cuts15),cbin.applyCuts(sigma_dis,cuts15),label='cut15')
    plt.title('d$\sigma_{DIS}$ vs TDIS_xbj', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.2,0.3), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('d$\sigma_{DIS}$ (fb)')
    plt.xlim(1e-4,1)
    plt.xscale('log')
    plt.yscale('log')
    print "sigma_dis cut1:\n",sigbin1.get_offsets()
    print "sigma_dis cut2:\n",sigbin2.get_offsets()
    print "sigma_dis cut3:\n",sigbin3.get_offsets()
    print "sigma_dis cut4:\n",sigbin4.get_offsets()
    print "sigma_dis cut5:\n",sigbin5.get_offsets()
    print "sigma_dis cut6:\n",sigbin6.get_offsets()
    print "sigma_dis cut7:\n",sigbin7.get_offsets()
    print "sigma_dis cut8:\n",sigbin8.get_offsets()
    print "sigma_dis cut9:\n",sigbin9.get_offsets()
    print "sigma_dis cut10:\n",sigbin10.get_offsets()
    print "sigma_dis cut11:\n",sigbin11.get_offsets()
    print "sigma_dis cut12:\n",sigbin12.get_offsets()
    print "sigma_dis cut13:\n",sigbin13.get_offsets()
    print "sigma_dis cut14:\n",sigbin14.get_offsets()
    print "sigma_dis cut15:\n",sigbin15.get_offsets()
    
    for f in xrange(1, plt.figure().number):
        pdf.savefig(f)
    pdf.close()    

def neutronDist() :

    [cuts1,cuts2,cuts3,cuts4,cuts5,cuts6,cuts7,cuts8,cuts9,cuts10,cuts11,cuts12,cuts13,cuts14,cuts15] = sigmaBin_Cut()

    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    sigbin1 = ax.hist(cbin.applyCuts(f2N,cuts1),label='cut1',histtype='step', alpha=0.5, stacked=True, fill=True)
    sigbin2 = ax.hist(cbin.applyCuts(f2N,cuts2),label='cut2',histtype='step', alpha=0.5, stacked=True, fill=True)
    sigbin3 = ax.hist(cbin.applyCuts(f2N,cuts3),label='cut3',histtype='step', alpha=0.5, stacked=True, fill=True)
    sigbin4 = ax.hist(cbin.applyCuts(f2N,cuts4),label='cut4',histtype='step', alpha=0.5, stacked=True, fill=True)
    sigbin5 = ax.hist(cbin.applyCuts(f2N,cuts5),label='cut5',histtype='step', alpha=0.5, stacked=True, fill=True)
    sigbin6 = ax.hist(cbin.applyCuts(f2N,cuts6),label='cut6',histtype='step', alpha=0.5, stacked=True, fill=True)
    sigbin7 = ax.hist(cbin.applyCuts(f2N,cuts7),label='cut7',histtype='step', alpha=0.5, stacked=True, fill=True)
    sigbin8 = ax.hist(cbin.applyCuts(f2N,cuts8),label='cut8',histtype='step', alpha=0.5, stacked=True, fill=True)
    sigbin9 = ax.hist(cbin.applyCuts(f2N,cuts9),label='cut9',histtype='step', alpha=0.5, stacked=True, fill=True)
    sigbin10 = ax.hist(cbin.applyCuts(f2N,cuts10),label='cut10',histtype='step', alpha=0.5, stacked=True, fill=True)
    sigbin11 = ax.hist(cbin.applyCuts(f2N,cuts11),label='cut11',histtype='step', alpha=0.5, stacked=True, fill=True)
    sigbin12 = ax.hist(cbin.applyCuts(f2N,cuts12),label='cut12',histtype='step', alpha=0.5, stacked=True, fill=True)
    sigbin13 = ax.hist(cbin.applyCuts(f2N,cuts13),label='cut13',histtype='step', alpha=0.5, stacked=True, fill=True)
    sigbin14 = ax.hist(cbin.applyCuts(f2N,cuts14),label='cut14',histtype='step', alpha=0.5, stacked=True, fill=True)
    sigbin15 = ax.hist(cbin.applyCuts(f2N,cuts15),label='cut15',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.title('$f^{2}_{N}$ binned', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.2,0.3), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('$f^{2}_{N}$')
    plt.ylabel('Number of Events')
    
    # for f in xrange(1, plt.figure().number):
    #     pdf.savefig(f)
    # pdf.close()

def pionDist() :
    
    [cuts1,cuts2,cuts3,cuts4,cuts5,cuts6,cuts7,cuts8,cuts9,cuts10,cuts11,cuts12,cuts13,cuts14,cuts15] = sigmaBin_Cut()

    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    xpibin = ax.hist(xpi,bins=p.setbin(xpi,200,0.,1.0),histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.title('$x_{pi}$ binned', fontsize =16)
    plt.xlabel('$x_{pi}$')
    plt.ylabel('Number of Events')
    
    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    ypibin = ax.hist(ypi,bins=p.setbin(ypi,200,0.,1.0),histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.title('$y_{pi}$', fontsize =16)
    plt.xlabel('$y_{pi}$')
    plt.ylabel('Number of Events')
    
    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    tpibin = ax.hist(tpi,bins=p.setbin(tpi,200,-1.0,0.),histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.title('$t_{pi}$', fontsize =16)
    plt.xlabel('$t_{pi}$')
    plt.ylabel('Number of Events')
    
    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    fpibin = ax.hist(fpi,histtype='step', alpha=0.5, stacked=True, fill=True)
    # fpibin = ax.hist(fpi,bins=p.setbin(fpi,200,0.,1.0),histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.title('$f_{pi}$', fontsize =16)
    plt.xlabel('$f_{pi}$')
    plt.ylabel('Number of Events')
    
    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    xpibin1 = ax.hist(cbin.applyCuts(xpi,cuts1),label='cut1',histtype='step', alpha=0.5, stacked=True, fill=True)
    xpibin2 = ax.hist(cbin.applyCuts(xpi,cuts2),label='cut2',histtype='step', alpha=0.5, stacked=True, fill=True)
    xpibin3 = ax.hist(cbin.applyCuts(xpi,cuts3),label='cut3',histtype='step', alpha=0.5, stacked=True, fill=True)
    xpibin4 = ax.hist(cbin.applyCuts(xpi,cuts4),label='cut4',histtype='step', alpha=0.5, stacked=True, fill=True)
    xpibin5 = ax.hist(cbin.applyCuts(xpi,cuts5),label='cut5',histtype='step', alpha=0.5, stacked=True, fill=True)
    xpibin6 = ax.hist(cbin.applyCuts(xpi,cuts6),label='cut6',histtype='step', alpha=0.5, stacked=True, fill=True)
    xpibin7 = ax.hist(cbin.applyCuts(xpi,cuts7),label='cut7',histtype='step', alpha=0.5, stacked=True, fill=True)
    xpibin8 = ax.hist(cbin.applyCuts(xpi,cuts8),label='cut8',histtype='step', alpha=0.5, stacked=True, fill=True)
    xpibin9 = ax.hist(cbin.applyCuts(xpi,cuts9),label='cut9',histtype='step', alpha=0.5, stacked=True, fill=True)
    xpibin10 = ax.hist(cbin.applyCuts(xpi,cuts10),label='cut10',histtype='step', alpha=0.5, stacked=True, fill=True)
    xpibin11 = ax.hist(cbin.applyCuts(xpi,cuts11),label='cut11',histtype='step', alpha=0.5, stacked=True, fill=True)
    xpibin12 = ax.hist(cbin.applyCuts(xpi,cuts12),label='cut12',histtype='step', alpha=0.5, stacked=True, fill=True)
    xpibin13 = ax.hist(cbin.applyCuts(xpi,cuts13),label='cut13',histtype='step', alpha=0.5, stacked=True, fill=True)
    xpibin14 = ax.hist(cbin.applyCuts(xpi,cuts14),label='cut14',histtype='step', alpha=0.5, stacked=True, fill=True)
    xpibin15 = ax.hist(cbin.applyCuts(xpi,cuts15),label='cut15',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.title('$x_{pi}$ binned', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.2,0.3), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('$x_{pi}$')
    plt.ylabel('Number of Events')
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.05, 0.95, 'Num Evts cut 3 = %0.f' % (len(cbin.applyCuts(xpi,cuts3))), transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    
    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    ypibin1 = ax.hist(cbin.applyCuts(ypi,cuts1),label='cut1',histtype='step', alpha=0.5, stacked=True, fill=True)
    ypibin2 = ax.hist(cbin.applyCuts(ypi,cuts2),label='cut2',histtype='step', alpha=0.5, stacked=True, fill=True)
    ypibin3 = ax.hist(cbin.applyCuts(ypi,cuts3),label='cut3',histtype='step', alpha=0.5, stacked=True, fill=True)
    ypibin4 = ax.hist(cbin.applyCuts(ypi,cuts4),label='cut4',histtype='step', alpha=0.5, stacked=True, fill=True)
    ypibin5 = ax.hist(cbin.applyCuts(ypi,cuts5),label='cut5',histtype='step', alpha=0.5, stacked=True, fill=True)
    ypibin6 = ax.hist(cbin.applyCuts(ypi,cuts6),label='cut6',histtype='step', alpha=0.5, stacked=True, fill=True)
    ypibin7 = ax.hist(cbin.applyCuts(ypi,cuts7),label='cut7',histtype='step', alpha=0.5, stacked=True, fill=True)
    ypibin8 = ax.hist(cbin.applyCuts(ypi,cuts8),label='cut8',histtype='step', alpha=0.5, stacked=True, fill=True)
    ypibin9 = ax.hist(cbin.applyCuts(ypi,cuts9),label='cut9',histtype='step', alpha=0.5, stacked=True, fill=True)
    ypibin10 = ax.hist(cbin.applyCuts(ypi,cuts10),label='cut10',histtype='step', alpha=0.5, stacked=True, fill=True)
    ypibin11 = ax.hist(cbin.applyCuts(ypi,cuts11),label='cut11',histtype='step', alpha=0.5, stacked=True, fill=True)
    ypibin12 = ax.hist(cbin.applyCuts(ypi,cuts12),label='cut12',histtype='step', alpha=0.5, stacked=True, fill=True)
    ypibin13 = ax.hist(cbin.applyCuts(ypi,cuts13),label='cut13',histtype='step', alpha=0.5, stacked=True, fill=True)
    ypibin14 = ax.hist(cbin.applyCuts(ypi,cuts14),label='cut14',histtype='step', alpha=0.5, stacked=True, fill=True)
    ypibin15 = ax.hist(cbin.applyCuts(ypi,cuts15),label='cut15',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.title('$y_{pi}$ binned', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.2,0.3), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('$y_{pi}$')
    plt.ylabel('Number of Events')
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.05, 0.95, 'Num Evts cut 3 = %0.f' % (len(cbin.applyCuts(ypi,cuts3))), transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    
    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    tpibin1 = ax.hist(cbin.applyCuts(tpi,cuts1),label='cut1',histtype='step', alpha=0.5, stacked=True, fill=True)
    tpibin2 = ax.hist(cbin.applyCuts(tpi,cuts2),label='cut2',histtype='step', alpha=0.5, stacked=True, fill=True)
    tpibin3 = ax.hist(cbin.applyCuts(tpi,cuts3),label='cut3',histtype='step', alpha=0.5, stacked=True, fill=True)
    tpibin4 = ax.hist(cbin.applyCuts(tpi,cuts4),label='cut4',histtype='step', alpha=0.5, stacked=True, fill=True)
    tpibin5 = ax.hist(cbin.applyCuts(tpi,cuts5),label='cut5',histtype='step', alpha=0.5, stacked=True, fill=True)
    tpibin6 = ax.hist(cbin.applyCuts(tpi,cuts6),label='cut6',histtype='step', alpha=0.5, stacked=True, fill=True)
    tpibin7 = ax.hist(cbin.applyCuts(tpi,cuts7),label='cut7',histtype='step', alpha=0.5, stacked=True, fill=True)
    tpibin8 = ax.hist(cbin.applyCuts(tpi,cuts8),label='cut8',histtype='step', alpha=0.5, stacked=True, fill=True)
    tpibin9 = ax.hist(cbin.applyCuts(tpi,cuts9),label='cut9',histtype='step', alpha=0.5, stacked=True, fill=True)
    tpibin10 = ax.hist(cbin.applyCuts(tpi,cuts10),label='cut10',histtype='step', alpha=0.5, stacked=True, fill=True)
    tpibin11 = ax.hist(cbin.applyCuts(tpi,cuts11),label='cut11',histtype='step', alpha=0.5, stacked=True, fill=True)
    tpibin12 = ax.hist(cbin.applyCuts(tpi,cuts12),label='cut12',histtype='step', alpha=0.5, stacked=True, fill=True)
    tpibin13 = ax.hist(cbin.applyCuts(tpi,cuts13),label='cut13',histtype='step', alpha=0.5, stacked=True, fill=True)
    tpibin14 = ax.hist(cbin.applyCuts(tpi,cuts14),label='cut14',histtype='step', alpha=0.5, stacked=True, fill=True)
    tpibin15 = ax.hist(cbin.applyCuts(tpi,cuts15),label='cut15',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.title('$t_{pi}$ binned', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.2,0.3), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('$t_{pi}$')
    plt.ylabel('Number of Events')
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.05, 0.95, 'Num Evts cut 3 = %0.f' % (len(cbin.applyCuts(tpi,cuts3))), transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    
    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    fpibin1 = ax.hist(cbin.applyCuts(fpi,cuts1),label='cut1',histtype='step', alpha=0.5, stacked=True, fill=True)
    fpibin2 = ax.hist(cbin.applyCuts(fpi,cuts2),label='cut2',histtype='step', alpha=0.5, stacked=True, fill=True)
    fpibin3 = ax.hist(cbin.applyCuts(fpi,cuts3),label='cut3',histtype='step', alpha=0.5, stacked=True, fill=True)
    fpibin4 = ax.hist(cbin.applyCuts(fpi,cuts4),label='cut4',histtype='step', alpha=0.5, stacked=True, fill=True)
    fpibin5 = ax.hist(cbin.applyCuts(fpi,cuts5),label='cut5',histtype='step', alpha=0.5, stacked=True, fill=True)
    fpibin6 = ax.hist(cbin.applyCuts(fpi,cuts6),label='cut6',histtype='step', alpha=0.5, stacked=True, fill=True)
    fpibin7 = ax.hist(cbin.applyCuts(fpi,cuts7),label='cut7',histtype='step', alpha=0.5, stacked=True, fill=True)
    fpibin8 = ax.hist(cbin.applyCuts(fpi,cuts8),label='cut8',histtype='step', alpha=0.5, stacked=True, fill=True)
    fpibin9 = ax.hist(cbin.applyCuts(fpi,cuts9),label='cut9',histtype='step', alpha=0.5, stacked=True, fill=True)
    fpibin10 = ax.hist(cbin.applyCuts(fpi,cuts10),label='cut10',histtype='step', alpha=0.5, stacked=True, fill=True)
    fpibin11 = ax.hist(cbin.applyCuts(fpi,cuts11),label='cut11',histtype='step', alpha=0.5, stacked=True, fill=True)
    fpibin12 = ax.hist(cbin.applyCuts(fpi,cuts12),label='cut12',histtype='step', alpha=0.5, stacked=True, fill=True)
    fpibin13 = ax.hist(cbin.applyCuts(fpi,cuts13),label='cut13',histtype='step', alpha=0.5, stacked=True, fill=True)
    fpibin14 = ax.hist(cbin.applyCuts(fpi,cuts14),label='cut14',histtype='step', alpha=0.5, stacked=True, fill=True)
    fpibin15 = ax.hist(cbin.applyCuts(fpi,cuts15),label='cut15',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.title('$f_{pi}$ binned', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.2,0.3), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('$f_{pi}$')
    plt.ylabel('Number of Events')
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.05, 0.95, 'Num Evts cut 3 = %0.f' % (len(cbin.applyCuts(fpi,cuts3))), transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    
    # for f in xrange(1, plt.figure().number):
    #     pdf.savefig(f)
    # pdf.close() 

def momentumPlots():

    [cuts1,cuts2,cuts3,cuts4,cuts5,cuts6,cuts7,cuts8,cuts9,cuts10,cuts11,cuts12,cuts13,cuts14,cuts15] = sigmaBin_Cut()

    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    # pprx_incbin = ax.hist(pprx_inc,bins=p.setbin(pprx_inc,200,-150.,150.),histtype='step', alpha=0.5, stacked=True, fill=True)
    EprE_incbin = ax.hist(EprE_inc,bins=p.setbin(EprE_inc,200,-150.,150.),histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.title('$p_{pr}$', fontsize =16)
    plt.xlabel('$p_{pr}$')
    plt.ylabel('Number of Events')
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.05, 0.95, 'Num Evts = %0.f' % (len(pprz_inc)), transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    ppiz_fracbin = ax.hist(ppiz_frac,bins=p.setbin(ppiz_frac,200,0.,1.0),histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.title('$p_{\pi}$/$p_{pr}$', fontsize =16)
    plt.xlabel('$p_{\pi}$/$p_{pr}$')
    plt.ylabel('Number of Events')
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.75, 0.95, 'Num Evts = %0.f' % (len(ppiz_frac)), transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    pNz_fracbin = ax.hist(pNz_frac,bins=p.setbin(pNz_frac,200,0.,1.0),histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.title('$p_{N}$/$p_{pr}$', fontsize =16)
    plt.xlabel('$p_{N}$/$p_{pr}$')
    plt.ylabel('Number of Events')
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.05, 0.95, 'Num Evts = %0.f' % (len(pNz_frac)), transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)

def q2Plots():

    [cuts1,cuts2,cuts3,cuts4,cuts5,cuts6,cuts7,cuts8,cuts9,cuts10,cuts11,cuts12,cuts13,cuts14,cuts15] = sigmaBin_Cut()

    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    q2bin = ax.hist(Q2,bins=p.setbin(Q2,200,1.,8.0),histtype='step', alpha=0.5, stacked=True, fill=True,label='Q2')
    # q2_newbin = ax.hist(Q2_new,bins=p.setbin(Q2,200,1.,8.0),histtype='step', alpha=0.5, stacked=True, fill=True,label='Q2 new')
    piq2bin = ax.hist(piQ2,bins=p.setbin(Q2,200,1.,8.0),histtype='step', alpha=0.5, stacked=True, fill=True,label='$\pi$ Q2')
    plt.title('$Q^{2}$', fontsize =16)
    plt.xlabel('$Q^{2}$')
    plt.ylabel('Number of Events')
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax.text(0.05, 0.95, 'Num Evts = %0.f' % (len(pNz_frac)), transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
    # leg = plt.legend(bbox_to_anchor=(0.2,0.3), loc="center right")
    # leg.get_frame().set_alpha(1.)

def uncernCalc():

    [low1,mid1,high1,low2,mid2,high2,low3,mid3,high3,low4,mid4,high4,low5,mid5,high5,low6,mid6,high6] = uncern_Cut()
    
    # f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    # # xerr, yerr
    # # scat1 = 3ax.errorbar(c.applyCuts(TDIS_xbj,cutsQ2),c.applyCuts(Q2,cutsQ2),xerr=np.sqrt(c.applyCuts(TDIS_xbj,cutsQ2))/200,yerr=np.sqrt(c.applyCuts(Q2,cutsQ2))/200,fmt='o',label='$Q^2$=6.5 $GeV^2$')
    # scat1 = ax.errorbar(clow.applyCuts(TDIS_xbj_low,low1),clow.applyCuts(Q2_low,low1),xerr=np.sqrt(clow.applyCuts(TDIS_xbj_low,low1))/200,yerr=np.sqrt(clow.applyCuts(Q2_low,low1))/200,fmt='o',label='$Q^2$=6.5 $GeV^2$')
    # scat2 = ax.errorbar(cmid.applyCuts(TDIS_xbj_mid,mid2),cmid.applyCuts(Q2_mid,mid2),xerr=np.sqrt(cmid.applyCuts(TDIS_xbj_mid,mid2))/200,yerr=np.sqrt(cmid.applyCuts(Q2_mid,mid2))/200,fmt='o',label='$Q^2$=15.0 $GeV^2$')
    # scat3 = ax.errorbar(cmid.applyCuts(TDIS_xbj_mid,mid3),cmid.applyCuts(Q2_mid,mid3),xerr=np.sqrt(cmid.applyCuts(TDIS_xbj_mid,mid3))/200,yerr=np.sqrt(cmid.applyCuts(Q2_mid,mid3))/200,fmt='o',label='$Q^2$=35.0 $GeV^2$')
    # scat4 = ax.errorbar(cmid.applyCuts(TDIS_xbj_mid,mid4),cmid.applyCuts(Q2_mid,mid4),xerr=np.sqrt(cmid.applyCuts(TDIS_xbj_mid,mid4))/200,yerr=np.sqrt(cmid.applyCuts(Q2_mid,mid4))/200,fmt='o',label='$Q^2$=60.0 $GeV^2$')
    # scat5 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,high5),chigh.applyCuts(Q2_high,high5),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,high5))/200,yerr=np.sqrt(chigh.applyCuts(Q2_high,high5))/200,fmt='o',label='$Q^2$=120.0 $GeV^2$')
    # scat6 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,high6),chigh.applyCuts(Q2_high,high6),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,high6))/200,yerr=np.sqrt(chigh.applyCuts(Q2_high,high6))/200,fmt='o',label='$Q^2$=250.0 $GeV^2$')
    # plt.title('$Q^{2}$ vs TDIS_xbj', fontsize =16)
    # leg = plt.legend(bbox_to_anchor=(0.95,0.3), loc="center right")
    # leg.get_frame().set_alpha(1.)
    # plt.xlabel('TDIS_xbj')
    # plt.ylabel('$Q^{2}$')
    # plt.xlim(0.,1.)
    # # plt.ylim(0.,8.)
    # # print "uncern", np.sqrt(c.applyCuts(TDIS_xbj,cutsQ2a))/200


    f, ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27))
    # totsigscat = ax.scatter(TDIS_xbj,tot_sigma,label='No cuts')
    totsigscat1 = ax.errorbar(clow.applyCuts(TDIS_xbj_low,low1),clow.applyCuts(tot_sigma_low,low1),xerr=np.sqrt(clow.applyCuts(TDIS_xbj_low,low1))/200,fmt='o',label='$Q^2$=6.5 $GeV^2$')
    totsigscat2 = ax.errorbar(cmid.applyCuts(TDIS_xbj_mid,mid2),cmid.applyCuts(tot_sigma_mid,mid2),xerr=np.sqrt(cmid.applyCuts(TDIS_xbj_mid,mid2))/200,fmt='o',label='$Q^2$=15.0 $GeV^2$')
    totsigscat3 = ax.errorbar(cmid.applyCuts(TDIS_xbj_mid,mid3),cmid.applyCuts(tot_sigma_mid,mid3),xerr=np.sqrt(cmid.applyCuts(TDIS_xbj_mid,mid3))/200,fmt='o',label='$Q^2$=35.0 $GeV^2$')
    totsigscat4 = ax.errorbar(cmid.applyCuts(TDIS_xbj_mid,mid4),cmid.applyCuts(tot_sigma_mid,mid4),xerr=np.sqrt(cmid.applyCuts(TDIS_xbj_mid,mid4))/200,fmt='o',label='$Q^2$=60.0 $GeV^2$')
    totsigscat5 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,high5),chigh.applyCuts(tot_sigma_high,high5),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,high5))/200,fmt='o',label='$Q^2$=120.0 $GeV^2$')
    totsigscat6 = ax.errorbar(chigh.applyCuts(TDIS_xbj_high,high6),chigh.applyCuts(tot_sigma_high,high6),xerr=np.sqrt(chigh.applyCuts(TDIS_xbj_high,high6))/200,fmt='o',label='$Q^2$=250.0 $GeV^2$')
    plt.title('Reduced $\sigma_{DIS}$ vs TDIS_xbj', fontsize =16)
    leg = plt.legend(bbox_to_anchor=(0.95,0.3), loc="center right")
    leg.get_frame().set_alpha(1.)
    plt.xlabel('TDIS_xbj')
    plt.ylabel('Reduced $\sigma_{DIS}$ ($nb^{-1}$)')
    plt.xlim(0.,1.)
    # plt.xscale('log')
    # plt.ylim(0.,15e-4)
    # plt.yscale('log')
    # print "tot_sigma Q^2=6.5:\n",totsigscat1.get_offsets()
    # print "tot_sigma Q^2=15.0:\n",totsigscat2.get_offsets()
    # print "tot_sigma Q^2=35.0:\n",totsigscat3.get_offsets()
    # print "tot_sigma Q^2=60.0:\n",totsigscat4.get_offsets()
    # print "tot_sigma Q^2=120.0:\n",totsigscat5.get_offsets()
    # print "tot_sigma Q^2=250.0:\n",totsigscat6.get_offsets()

def phaseSpace():
    
    # Q2 vs xBj
    hxbj_Q2_low = densityPlot(TDIS_xbj_low, Q2_low, '$Q^{2}$[1-8 GeV] vs TDIS_xbj[0.001-0.008]','TDIS_xbj','$Q^{2}$', 20, 20, 0.001, 0.008, 1., 8.)
    plt.xscale('log')
    # plt.yscale('log')
    hxbj_Q2_mid = densityPlot(TDIS_xbj_mid, Q2_mid, '$Q^{2}$[10-80 GeV] vs TDIS_xbj[0.01-0.08]','TDIS_xbj','$Q^{2}$', 20, 20, 0.01, 0.08, 10., 80.)
    plt.xscale('log')
    # plt.yscale('log')
    hxbj_Q2_high = densityPlot(TDIS_xbj_high, Q2_high, '$Q^{2}$[100-800 GeV] vs TDIS_xbj[0.1-0.8]','TDIS_xbj','$Q^{2}$', 20, 20, 0.1, 0.8, 100., 800.)
    plt.xscale('log')
    # plt.yscale('log')

    # Q2 vs xBj chained
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    hxbj_Q2_low = densityPlot(TDIS_xbj_low, Q2_low, '','TDIS_xbj','$Q^{2}$', 20, 20, 0.001, 0.008, 1., 8., figure=f,ax=ax)
    hxbj_Q2_mid = densityPlot(TDIS_xbj_mid, Q2_mid, '','TDIS_xbj','$Q^{2}$', 20, 20,  0.01, 0.08, 10., 80., figure=f, ax=ax, layered=False)
    hxbj_Q2_high = densityPlot(TDIS_xbj_high, Q2_high, '$Q^{2}$ vs TDIS_xbj (all sets)','TDIS_xbj','$Q^{2}$', 20, 20, 0.1, 0.8, 100., 800., figure=f,ax=ax,layered=False)
    plt.xscale('log')
    plt.xlim(0.,1.)
    # plt.yscale('log')
    plt.ylim(0.,800.)

    # Q2 vs xpi
    hxpi_Q2_low = densityPlot(xpi_low, Q2_low, '$Q^{2}$[1-8 GeV] vs $x_{pi}$[0.001-0.008]','$x_{pi}$','$Q^{2}$', 20, 20, 0.001, 0.008, 1., 8.)
    plt.xscale('log')
    # plt.yscale('log')
    hxpi_Q2_mid = densityPlot(xpi_mid, Q2_mid, '$Q^{2}$[10-80 GeV] vs $x_{pi}$[0.01-0.08]','$x_{pi}$','$Q^{2}$', 20, 20, 0.01, 0.08, 10., 80.)
    plt.xscale('log')
    # plt.yscale('log')
    hxpi_Q2_high = densityPlot(xpi_high, Q2_high, '$Q^{2}$[100-800 GeV] vs $x_{pi}$[0.1-0.8]','$x_{pi}$','$Q^{2}$', 20, 20, 0.1, 0.8, 100., 800.)
    plt.xscale('log')
    # plt.yscale('log')

    # Q2 vs xpi chained
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    hxpi_Q2_low = densityPlot(xpi_low, Q2_low, '','$x_{pi}$','$Q^{2}$', 20, 20, 0.001, 0.008, 1., 8., figure=f,ax=ax)
    hxpi_Q2_mid = densityPlot(xpi_mid, Q2_mid, '','$x_{pi}$','$Q^{2}$', 20, 20,  0.01, 0.08, 10., 80., figure=f, ax=ax, layered=False)
    hxpi_Q2_high = densityPlot(xpi_high, Q2_high, '$Q^{2}$ vs $x_{pi}$ (all sets)','$x_{pi}$','$Q^{2}$', 20, 20, 0.1, 0.8, 100., 800., figure=f,ax=ax,layered=False)
    plt.xscale('log')
    plt.xlim(0.,1.)
    # plt.yscale('log')
    plt.ylim(0.,800.)

    
    for f in xrange(1, plt.figure().number):
        pdf_q2.savefig(f)
    pdf_q2.close()

def sigmavX():

    
    # Cross-section vs xBj
    hTDIS_xbj_tot_sigma_low = densityPlot(TDIS_xbj_low, tot_sigma_low, 'reduced $\sigma_{dis}$[1-8 GeV] vs TDIS_xbj[0.001-0.008]','TDIS_xbj','reduced $\sigma_{dis}$', 20, 20, 0.001, 0.008)
    plt.xscale('log')
    plt.yscale('log')
    hTDIS_xbj_tot_sigma_mid = densityPlot(TDIS_xbj_mid, tot_sigma_mid, 'reduced $\sigma_{dis}$[10-80 GeV] vs TDIS_xbj[0.01-0.08]','TDIS_xbj','reduced $\sigma_{dis}$', 20, 20, 0.01, 0.08)
    plt.xscale('log')
    plt.yscale('log')
    hTDIS_xbj_tot_sigma_high = densityPlot(TDIS_xbj_high, tot_sigma_high, 'reduced $\sigma_{dis}$[100-800 GeV] vs TDIS_xbj[0.1-0.8]','TDIS_xbj','reduced $\sigma_{dis}$', 20, 20, 0.1, 0.8)
    plt.xscale('log')
    plt.yscale('log')

    # Cross-section vs xBj chained
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    hTDIS_xbj_tot_sigma_low = densityPlot(TDIS_xbj_low, tot_sigma_low, '','TDIS_xbj','reduced $\sigma_{dis}$', 20, 20, 0.001, 0.008, figure=f,ax=ax)
    hTDIS_xbj_tot_sigma_mid = densityPlot(TDIS_xbj_mid, tot_sigma_mid, '','TDIS_xbj','reduced $\sigma_{dis}$', 20, 20, 0.01, 0.08, figure=f,ax=ax,layered=False)
    hTDIS_xbj_tot_sigma_high = densityPlot(TDIS_xbj_high, tot_sigma_high, 'reduced $\sigma_{dis}$ vs TDIS_xbj (all sets)','TDIS_xbj','reduced $\sigma_{dis}$', 20, 20, 0.1, 0.8, figure=f,ax=ax,layered=False)
    plt.xscale('log')
    plt.xlim(0.,1.)
    plt.yscale('log')
    plt.ylim(0.,0.15)
    
    # Cross-section vs xpi
    hxpi_tot_sigma_low = densityPlot(xpi_low, tot_sigma_low, 'reduced $\sigma_{dis}$[1-8 GeV] vs $x_{pi}$[0.001-.008]','$x_{pi}$','reduced $\sigma_{dis}$', 20, 20, 0.001, 0.008)
    plt.xscale('log')
    plt.yscale('log')
    hxpi_tot_sigma_mid = densityPlot(xpi_mid, tot_sigma_mid, 'reduced $\sigma_{dis}$[10-80 GeV] vs $x_{pi}$[0.01-0.08]','$x_{pi}$','reduced $\sigma_{dis}$', 20, 20, 0.01, 0.08)
    plt.xscale('log')
    plt.yscale('log')
    hxpi_tot_sigma_high = densityPlot(xpi_high, tot_sigma_high, 'reduced $\sigma_{dis}$[100-800 GeV] vs $x_{pi}$[0.1-0.8]','$x_{pi}$','reduced $\sigma_{dis}$', 20, 20, 0.1, 0.8)
    plt.xscale('log')
    plt.yscale('log')

    # Cross-section vs xpi chained
    f,ax = plt.subplots(tight_layout=True,figsize=(11.69,8.27));
    hxpi_tot_sigma_low = densityPlot(xpi_low, tot_sigma_low, '','$x_{pi}$','reduced $\sigma_{dis}$', 20, 20, 0.001, 0.008, figure=f,ax=ax)
    hxpi_tot_sigma_mid = densityPlot(xpi_mid, tot_sigma_mid, '','$x_{pi}$','reduced $\sigma_{dis}$', 20, 20, 0.01, 0.08, figure=f,ax=ax,layered=False)
    hxpi_tot_sigma_high = densityPlot(xpi_high, tot_sigma_high, 'reduced $\sigma_{dis}$ vs $x_{pi}$ (all sets)','$x_{pi}$','reduced $\sigma_{dis}$', 20, 20, 0.1, 0.8, figure=f,ax=ax,layered=False)
    plt.xscale('log')
    plt.xlim(0.,1.)
    plt.yscale('log')
    plt.ylim(0.,0.15)

    for f in xrange(1, plt.figure().number):
        pdf_sigma.savefig(f)
    pdf_sigma.close()
    
def main() :

    # q2Plots()
    # momentumPlots()
    # neutronDist()
    # pionDist()
    # sigmaBin_Plot()
    # sigmaDIS_Plot()
    # recreateLeaves()
    uncernCalc()
    # phaseSpace()
    # sigmavX()
    plt.show()
    
if __name__=='__main__': main()
