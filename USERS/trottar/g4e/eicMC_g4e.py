#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2020-03-19 11:43:13 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

from g4epy import Geant4Eic
from pyjano.jana import Jana

# final_state="pi_p_18on275"
# final_state="pi_n_18on275"
# final_state="k_lambda_18on275"

# final_state="pi_p_10on275"
# final_state="pi_n_10on275"
# final_state="k_lambda_10on275"

final_state="pi_n_10on100"

g4e = Geant4Eic(detector='jleic', beamline='erhic')\
         .source('../OUTPUTS/%s_lund.dat' % final_state)\
         .output('OUTPUTS/g4e_%s' % final_state)\
         .beam_on(10000)
         # .vis()
g4e.run()

Jana().plugin('g4e_reader')\
      .plugin('meson_structure')\
      .plugin('jana', output='OUTPUTS/jana_%s.root' % final_state)\
      .source('OUTPUTS/g4e_%s.root' % final_state)\
      .run()
