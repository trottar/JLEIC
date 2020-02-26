#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2020-02-25 16:25:06 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

from g4epy import Geant4Eic

g4e=Geant4Eic(detector='jleic',beamline='erhic')\
                        .source('../TDIS_lund_500k.dat')\
                        .output('jleic_mesonMC_500k')\
                        .beam_on(1)\
                        .vis()\
                        .run()
