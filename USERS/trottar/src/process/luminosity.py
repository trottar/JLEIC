#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2021-06-10 20:13:51 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

import numpy as np

def Lumi(cross_sec,xbinwidth,qbinwidth,tbinwidth,xLbinwidth):

    # Luminosity

    # all binned events
    sig_all = cross_sec

    nevt = [(100*((s)*(qbinwidth)*(xbinwidth)*(xLbinwidth))*1e6) for s in sig_all] # The 1e6 converts properly, integrated luminosiy: 100 fb^-1
    nevt  = np.asarray(nevt)
    #print("\nEvents expected running at 100 $fb^{-1}$: ", nevt)

    return nevt