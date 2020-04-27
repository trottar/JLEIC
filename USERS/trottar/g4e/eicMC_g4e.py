#! /usr/bin/python

#
# Description: Source escalate first....source ~/ResearchNP/gemc/escalate/escalate/escalate.csh
# ================================================================
# Time-stamp: "2020-04-26 21:28:19 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

from g4epy import Geant4Eic
from pyjano.jana import Jana, PluginFromSource

mstruct_general = PluginFromSource('/home/trottar/ResearchNP/gemc/my_plugins/meson_structure/mstruct_general')
# Name will be determined from folder name
# add name=<...> for custom name

# cmake clean
mstruct_general.builder.clean()
# change compiler version
mstruct_general.builder.config['cxx_standard'] = 11

# final_state="pi_p_18on275"
final_state="pi_n_18on275"
# final_state="k_lambda_18on275"

def runG4e(numEvts,g4e_flag=False):

    if(g4e_flag):
        g4e = Geant4Eic(detector='jleic', beamline='erhic')\
                                  .source('../OUTPUTS/%s_lund.dat' % final_state)\
                                  .output('OUTPUTS/g4e_%s' % final_state)\
                                  .beam_on(numEvts)\
        # g4e.vis()
        g4e.run()

    jana = Jana(nevents=numEvts, output='OUTPUTS/jana_%s.root' % final_state)

    # G4E reader here
    jana.plugin('g4e_reader') \
        .source('OUTPUTS/g4e_%s.root' % final_state)

    # Parameters:
    #     verbose   - Plugin output level. 0-almost nothing, 1-some, 2-everything
    #     smearing  - Particle smearing 0-true MC, 1-smearing, 2-reconstruction");
    #     e_beam_energy    -  Energy of colliding electron beam");
    #     ion_beam_energy  -  Energy of colliding ion beam");
    jana.plugin(mstruct_general, verbose=1)
    jana.run()

def main() :

    runG4e(10000)
    # runG4e(10000,True)
    
if __name__=='__main__': main()
