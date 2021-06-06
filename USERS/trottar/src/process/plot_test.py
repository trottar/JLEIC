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

import matplotlib.pyplot as plt
import matplotlib.colors as colors

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

fig = plt.figure(figsize=(17,12),facecolor='silver')

ax = fig.add_subplot(331)
#plt.scatter(data.t_qbin,data.fpi_qbin)
densityPlot(data.t_qbin,data.fpi_qbin, '$fpi$ vs $t$','$t$','$fpi$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(332)
#plt.scatter(TDIS_xbj_qbin,fpi_qbin)
densityPlot(data.TDIS_xbj_qbin,data.fpi_qbin, '$fpi$ vs $x$','$x$','$fpi$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(333)
densityPlot(data.TDIS_xbj_qbin,data.xL_qbin, '$xL$ vs $x$','$x$','$xL$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(334)
densityPlot(data.TDIS_xbj_qbin,data.t_qbin, '$t$ vs $x$','$x$','$t$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(335)
densityPlot(data.xL_qbin,data.t_qbin, '$t$ vs $xL$','$xL$','$t$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(336)
densityPlot(data.TDIS_xbj_xbin,data.Q2_xbin, '$Q^2$ vs $x$','$x$','$Q^{2}$', 200, 200, ax=ax, fig=fig)

ax = fig.add_subplot(337)
densityPlot(data.TDIS_xbj_qbin,data.Q2_qbin, '$Q^2$ vs $x$','$x$','$Q^{2}$', 200, 200, ax=ax, fig=fig)

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