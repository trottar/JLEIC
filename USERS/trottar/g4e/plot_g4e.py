#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2020-05-03 17:33:10 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

import numpy as np
import uproot as up
import awkward
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import sys, math

sys.path.insert(0,'/home/trottar/bin/python/root2py/')
import root2py as r2p

kinematics = "pi_n_18on275"
# kinematics = "k_lambda_18on275"
# kinematics = "lambda_cone275_output"

pdf = matplotlib.backends.backend_pdf.PdfPages("g4e_%s.pdf" % kinematics) 

# rootName = "OUTPUTS/g4e_TDISpion_lund_p275_e10_10k.root" # Yulia root file
# rootName = "OUTPUTS/g4e_%s.root" % kinematics
rootName = "OUTPUTS/g4e_%s_lund_new_output.root" % kinematics # Yulia root file

tree = up.open(rootName)["events"]
branch = r2p.pyBranch(tree)

gen_prt_charge = tree.array("gen_prt_charge")
gen_prt_id = tree.array("gen_prt_id")
gen_prt_dir_x = np.arccos(tree.array("gen_prt_dir_x"))*(180/math.pi)
gen_prt_dir_y = np.arccos(tree.array("gen_prt_dir_y"))*(180/math.pi)
gen_prt_dir_z = np.arccos(tree.array("gen_prt_dir_z"))*(180/math.pi)
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

c = r2p.pyPlot(cutDict)

picut = ["pion"]
ecut = ["electron"]
ncut = ["neutron"]

f = plt.figure(figsize=(11.69,8.27))
plt.style.use('default')

ax = f.add_subplot(211)
hist1 = ax.hist(c.applyCuts(gen_prt_tot_mom,picut),bins=c.setbin(gen_prt_tot_mom,200),label='pi cut',histtype='step', alpha=0.5, stacked=True, fill=True)
hist1a = ax.hist(c.applyCuts(gen_prt_tot_mom,ecut),bins=c.setbin(gen_prt_tot_mom,200),label='e cut',histtype='step', alpha=0.5, stacked=True, fill=True)
hist1b = ax.hist(c.applyCuts(gen_prt_tot_mom,ncut),bins=c.setbin(gen_prt_tot_mom,200),label='n cut',histtype='step', alpha=0.5, stacked=True, fill=True)
ax.legend(loc=1)
plt.title('gen_prt_tot_mom', fontsize =20)

ax = f.add_subplot(212)
hist2 = ax.hist(c.applyCuts(gen_prt_dir_z,picut),bins=c.setbin(gen_prt_dir_z,200),label='pi cut',histtype='step', alpha=0.5, stacked=True, fill=True)
hist2a = ax.hist(c.applyCuts(gen_prt_dir_z,ecut),bins=c.setbin(gen_prt_dir_z,200),label='e cut',histtype='step', alpha=0.5, stacked=True, fill=True)
hist2b = ax.hist(c.applyCuts(gen_prt_dir_z,ncut),bins=c.setbin(gen_prt_dir_z,200),label='n cut',histtype='step', alpha=0.5, stacked=True, fill=True)
ax.legend(loc=2)
plt.title('gen_prt_dir_z', fontsize =20)

phaseSpace = c.densityPlot(c.applyCuts(gen_prt_tot_mom,picut), c.applyCuts(gen_prt_dir_z,picut), 'dir_z vs tot_mom','tot_mom','dir_z', 200, 200,  c)
# plt.ylim(0.,60.)
# plt.xlim(0.,600.)
plt.xlabel('tot_mom (GeV)', fontsize =20)
plt.ylabel('Theta (deg)', fontsize =20)
plt.title('pi cut', fontsize =20)
phaseSpace = c.densityPlot(c.applyCuts(gen_prt_tot_mom,ecut), c.applyCuts(gen_prt_dir_z,ecut), 'dir_z vs tot_mom','tot_mom','dir_z', 200, 200,  c)
# plt.ylim(0.,180.)
# plt.xlim(0.,50.)
plt.xlabel('tot_mom (GeV)', fontsize =20)
plt.ylabel('Theta (deg)', fontsize =20)
plt.title('e cut', fontsize =20)
phaseSpace = c.densityPlot(c.applyCuts(gen_prt_tot_mom,ncut), c.applyCuts(gen_prt_dir_z,ncut), 'dir_z vs tot_mom','tot_mom','dir_z', 200, 200,  c)
# plt.ylim(0.,180.)
# plt.xlim(0.,500.)
plt.xlabel('tot_mom (GeV)', fontsize =20)
plt.ylabel('Theta (deg)', fontsize =20)
plt.title('n cut', fontsize =20)
phaseSpace = c.densityPlot(gen_prt_tot_mom, gen_prt_dir_z, 'dir_z vs tot_mom','tot_mom','dir_z', 200, 200,  c)
# plt.ylim(-30.,180.)
# plt.xlim(0.,50.)
plt.xlabel('tot_mom (GeV)', fontsize =20)
plt.ylabel('Theta (deg)', fontsize =20)
plt.title('no cut', fontsize =20)

for f in range(1, plt.figure().number):
    pdf.savefig(f)
pdf.close()

plt.show()
