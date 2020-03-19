#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2020-03-19 12:49:35 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

import numpy as np
import uproot as up
import awkward
import matplotlib.pyplot as plt#
import matplotlib.backends.backend_pdf
import sys, math

sys.path.insert(0,'/home/trottar/bin/python/root2py/')
from root2py import pyPlot, pyBranch, pyBin

# kinematics = "pi_n_10on100"
kinematics = "pi_n_10on275"
# kinematics = "pi_n_18on275"

# pdf = matplotlib.backends.backend_pdf.PdfPages("g4e_pi_p.pdf")
pdf = matplotlib.backends.backend_pdf.PdfPages("g4e_%s_PIDcut.pdf") % kinematics

# rootName = "OUTPUTS/g4e_TDISpion_lund_p275_e10_10k.root" # Yulia root file
rootName = "OUTPUTS/g4e_%s.root" % kinematics

tree = up.open(rootName)["events"]
branch = pyBranch(tree)

gen_prt_charge = tree.array("gen_prt_charge")
gen_prt_id = tree.array("gen_prt_id")
gen_prt_dir_x = tree.array("gen_prt_dir_x")
gen_prt_dir_y = tree.array("gen_prt_dir_y")
gen_prt_dir_z = tree.array("gen_prt_dir_z")*(180/math.pi)
gen_prt_tot_mom = tree.array("gen_prt_tot_mom")

# Convert from awkward array to 1D
gen_prt_id = gen_prt_id.content
gen_prt_charge = gen_prt_charge.content
gen_prt_dir_z = gen_prt_dir_z.content
gen_prt_tot_mom = gen_prt_tot_mom.content

cutDict = {
    # "pion" : ((gen_prt_charge > 0.5)),
    "pion" : ((gen_prt_id > 0.5) & (gen_prt_id < 1.5)),
    # "neutron" : ((gen_prt_charge < 0.5) & (gen_prt_charge > -0.5)),
    "neutron" : ((gen_prt_id > 1.5) & (gen_prt_id < 2.5)),
    # "electron" : (gen_prt_charge < -0.5),
    "electron" : ((gen_prt_id < 0.5)),
}

b = pyBin()
c = pyPlot(cutDict)

picut = ["pion"]
ecut = ["electron"]
ncut = ["neutron"]

f = plt.figure(figsize=(11.69,8.27))
plt.style.use('default')

ax = f.add_subplot(211)
hist1 = ax.hist(c.applyCuts(gen_prt_tot_mom,picut),bins=b.setbin(gen_prt_tot_mom,200),label='pi cut',histtype='step', alpha=0.5, stacked=True, fill=True)
hist1a = ax.hist(c.applyCuts(gen_prt_tot_mom,ecut),bins=b.setbin(gen_prt_tot_mom,200),label='e cut',histtype='step', alpha=0.5, stacked=True, fill=True)
hist1b = ax.hist(c.applyCuts(gen_prt_tot_mom,ncut),bins=b.setbin(gen_prt_tot_mom,200),label='n cut',histtype='step', alpha=0.5, stacked=True, fill=True)
ax.legend(loc=1)
plt.title('gen_prt_tot_mom', fontsize =20)

ax = f.add_subplot(212)
hist2 = ax.hist(c.applyCuts(gen_prt_dir_z,picut),bins=b.setbin(gen_prt_dir_z,200),label='pi cut',histtype='step', alpha=0.5, stacked=True, fill=True)
hist2a = ax.hist(c.applyCuts(gen_prt_dir_z,ecut),bins=b.setbin(gen_prt_dir_z,200),label='e cut',histtype='step', alpha=0.5, stacked=True, fill=True)
hist2b = ax.hist(c.applyCuts(gen_prt_dir_z,ncut),bins=b.setbin(gen_prt_dir_z,200),label='n cut',histtype='step', alpha=0.5, stacked=True, fill=True)
ax.legend(loc=2)
plt.title('gen_prt_dir_z', fontsize =20)

phaseSpace = c.densityPlot(c.applyCuts(gen_prt_tot_mom,picut), c.applyCuts(gen_prt_dir_z,picut), 'dir_z vs tot_mom','tot_mom','dir_z', 200, 200,  b)
plt.ylim(-180.,180.)
plt.xlim(0.,600.)
plt.xlabel('tot_mom (GeV)')
plt.ylabel('Theta (deg)')
plt.title('pi cut', fontsize =20)
phaseSpace = c.densityPlot(c.applyCuts(gen_prt_tot_mom,ecut), c.applyCuts(gen_prt_dir_z,ecut), 'dir_z vs tot_mom','tot_mom','dir_z', 200, 200,  b)
plt.ylim(-180.,180.)
# plt.xlim(0.,50.)
plt.xlabel('tot_mom (GeV)')
plt.ylabel('Theta (deg)')
plt.title('e cut', fontsize =20)
phaseSpace = c.densityPlot(c.applyCuts(gen_prt_tot_mom,ncut), c.applyCuts(gen_prt_dir_z,ncut), 'dir_z vs tot_mom','tot_mom','dir_z', 200, 200,  b)
plt.ylim(-180.,180.)
# plt.xlim(0.,50.)
plt.xlabel('tot_mom (GeV)')
plt.ylabel('Theta (deg)')
plt.title('n cut', fontsize =20)
phaseSpace = c.densityPlot(gen_prt_tot_mom, gen_prt_dir_z, 'dir_z vs tot_mom','tot_mom','dir_z', 200, 200,  b, 0., 50., -180., 180.)
plt.ylim(-180.,180.)
# plt.xlim(0.,50.)
plt.xlabel('tot_mom (GeV)')
plt.ylabel('Theta (deg)')
plt.title('no cut', fontsize =20)

for f in xrange(1, plt.figure().number):
    pdf.savefig(f)
pdf.close()

plt.show()
