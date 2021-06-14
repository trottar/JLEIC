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
    sig_all = cross_sec*1e5 # the 1e5 corrects the scaling up top
    evts_all = len(cross_sec)*[len(cross_sec)]

    lumi = [(e)/((s)*(qbinwidth)*(xbinwidth)*(tbinwidth)*(xLbinwidth)) for e,s in zip(evts_all,sig_all)]
    #print("---------------------------------\n")
    # print("\nLuminosity: ", lumi)

    nevt = [100/(l*1e-6) for l in lumi] # The 1e-6 converts properly, integrated luminosiy: 100 fb^-1
    nevt  = np.asarray(nevt)
    #print("\nEvents expected running at 100 $fb^{-1}$: ", nevt)

    return nevt