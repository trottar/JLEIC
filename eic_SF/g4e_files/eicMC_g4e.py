#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2020-03-08 14:34:32 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

from g4epy import Geant4Eic
from pyjano.jana import Jana

g4e = Geant4Eic(detector='jleic', beamline='erhic')\
         .source('../OUTPUTS/pi_p_lund.dat')\
         .output('OUTPUTS/g4e_pi_p')\
         .beam_on(1)\
         .vis()\
         .run()

Jana().plugin('g4e_reader')\
      .plugin('meson_structure')\
      .plugin('jana', output='OUTPUTS/jana_pi_p.root')\
      .source('OUTPUTS/g4e_pi_p.root')\
      .run()
