#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2020-04-07 11:47:20 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import ROOT, math, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from ROOT import Math, TLorentzVector, TFile
from array import array

sys.path.insert(0,'/home/trottar/bin/python/root2py/')
from root2py import pyPlot, pyBin

b = pyBin()
c = pyPlot(None)

# kinematics = "pi_n_18on275"
# kinematics = "pi_n_5on100"
# kinematics = "k_lambda_5on100"
kinematics = "k_lambda_18on275"
# kinematics = "pi_p_18on275"

pdf = matplotlib.backends.backend_pdf.PdfPages("MC_%s_PIDcut.pdf" % kinematics) 

myfile = TFile('../OUTPUTS/%s.root' % kinematics)
mytree = myfile.Get("Evnts")

print(mytree)

scat_electron_theta = np.array([])
pion_theta = np.array([])
neutron_theta = np.array([])

scat_electron_mom = np.array([])
pion_mom = np.array([])
neutron_mom = np.array([])

for entryNum in range(0 , mytree.GetEntries()):
    mytree.GetEntry(entryNum)
    lepton0 = TLorentzVector()
    lepton1 = TLorentzVector()
    scat_electron = getattr(mytree ,"e_Scat.")
    pion = getattr(mytree ,"pi.")
    neutron = getattr(mytree ,"p1_Sp.")

    scat_electron_theta = np.append(scat_electron_theta,scat_electron.Theta()*(180/math.pi))
    pion_theta = np.append(pion_theta,pion.Theta()*(180/math.pi))
    neutron_theta = np.append(neutron_theta,neutron.Theta()*(180/math.pi))

    scat_electron_mom = np.append(scat_electron_mom,scat_electron.E())
    pion_mom = np.append(pion_mom,pion.E())
    neutron_mom = np.append(neutron_mom,neutron.E())

f = plt.figure(figsize=(11.69,8.27))
plt.style.use('default')

ax = f.add_subplot(211)
hist1 = ax.hist(pion_theta,bins=b.setbin(pion_theta,200),label='pi cut',histtype='step', alpha=0.5, stacked=True, fill=True)
hist1a = ax.hist(scat_electron_theta,bins=b.setbin(scat_electron_theta,200),label='e cut',histtype='step', alpha=0.5, stacked=True, fill=True)
hist1b = ax.hist(neutron_theta,bins=b.setbin(neutron_theta,200),label='n cut',histtype='step', alpha=0.5, stacked=True, fill=True)
ax.legend(loc=1)
plt.title('theta', fontsize =20)

ax = f.add_subplot(212)
hist1 = ax.hist(pion_mom,bins=b.setbin(pion_mom,200),label='pi cut',histtype='step', alpha=0.5, stacked=True, fill=True)
hist1a = ax.hist(scat_electron_mom,bins=b.setbin(scat_electron_mom,200),label='e cut',histtype='step', alpha=0.5, stacked=True, fill=True)
hist1b = ax.hist(neutron_mom,bins=b.setbin(neutron_mom,200),label='n cut',histtype='step', alpha=0.5, stacked=True, fill=True)
ax.legend(loc=1)
plt.title('tot_mom', fontsize =20)

phaseSpace = c.densityPlot(scat_electron_mom, scat_electron_theta, 'dir_z vs tot_mom','tot_mom','dir_z', 200, 200,  b)
# plt.ylim(-180.,180.)
# plt.xlim(0.,50.)
plt.xlabel('tot_mom (GeV)', fontsize =20)
plt.ylabel('Theta (deg)', fontsize =20)
plt.title('e cut', fontsize =20)
phaseSpace = c.densityPlot(pion_mom, pion_theta, 'dir_z vs tot_mom','tot_mom','dir_z', 200, 200,  b)
# plt.ylim(0.,60.)
# plt.xlim(0.,600.)
plt.xlabel('tot_mom (GeV)', fontsize =20)
plt.ylabel('Theta (deg)', fontsize =20)
plt.title('pi cut', fontsize =20)
phaseSpace = c.densityPlot(neutron_mom, neutron_theta, 'dir_z vs tot_mom','tot_mom','dir_z', 200, 200,  b)
# plt.ylim(-180.,180.)
plt.xlim(0.,1000.)
plt.xlabel('tot_mom (GeV)', fontsize =20)
plt.ylabel('Theta (deg)', fontsize =20)
plt.title('$\Lambda$ cut', fontsize =20)

for f in range(1, plt.figure().number):
    pdf.savefig(f)
pdf.close()

plt.show()
