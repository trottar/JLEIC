#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2020-09-14 11:03:13 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import ROOT, math, sys
import numpy as np
import uproot as up
import matplotlib.pyplot as plt
from ROOT import Math, TLorentzVector, TFile
from array import array

sys.path.insert(0,'/home/trottar/bin/python/')
import root2py as r2p

kinematics = sys.argv[1]

# kinematics = "pi_n_18on275"
# kinematics = "pi_n_10on100"
# kinematics = "pi_n_5on100"
# kinematics = "pi_n_5on41"
# kinematics = "k_lambda_18on275"

num = 0

rootName = '../OUTPUTS/%s.root' % kinematics
tree = up.open(rootName)["Evnts"]
branch = r2p.pyBranch(tree)

Q2 = branch.findBranch("invts","Q2")
x = branch.findBranch("invts","xBj")
# t = -branch.findBranch("invts","tPrime")
t = -branch.findBranch("invts","tSpectator")
TDIS_znq = tree.array("TDIS_znq")
# TDIS_Mx2 = tree.array("TDIS_Mx2")
EeE_Lab = tree.array("EeE_Lab")

if "pi" in kinematics:
    xpi = tree.array("xpi")
    EnE_Lab = tree.array("EnE_Lab")
    EpiE_Lab = tree.array("EpiE_Lab")
    ppiz_Lab = tree.array("ppiz_Lab")
if "k" in kinematics:
    xpi = tree.array("xk")
    EnE_Lab = tree.array("ElambE_Lab")
    EpiE_Lab = tree.array("EkE_Lab")
    ppiz_Lab = tree.array("pkz_Lab")

c = r2p.pyPlot(None)

myfile = TFile(rootName)
mytree = myfile.Get("Evnts")

print(mytree)

scat_electron_theta = np.array([])
pion_theta = np.array([])
neutron_theta = np.array([])

scat_electron_mom = np.array([])
pion_mom = np.array([])
neutron_mom = np.array([])

Q2_cut = np.array([])
x_cut = np.array([])
t_cut = np.array([])
TDIS_znq_cut = np.array([])

for entryNum in range(0,mytree.GetEntries()):
    mytree.GetEntry(entryNum)
    # lepton0 = TLorentzVector()
    # lepton1 = TLorentzVector()
    # virtual = TLorentzVector()
    scat_electron = getattr(mytree ,"e_Scat.")
    if "pi" in kinematics:
        pion = getattr(mytree ,"pi.")
    if "k" in kinematics:
        pion = getattr(mytree ,"k.")
    neutron = getattr(mytree ,"p1_Sp.")
    
    # if Q2[entryNum] > num:
    # if t[entryNum]<-0. and t[entryNum]>-1.0:
    
    # if EpiE_Lab[entryNum]+EpiE_Lab[entryNum] < 275:
    if ppiz_Lab[entryNum] == pion.Z():
        # if pion_theta[entryNum] < 85.0:
        
        pion_theta = np.append(pion_theta,pion.Theta()*(180/math.pi))
        # pion_mom = np.append(pion_mom,pion.E())
        pion_mom = np.append(pion_mom,EpiE_Lab[entryNum])
    
        # if EnE_Lab[entryNum] == neutron.E():
        neutron_theta = np.append(neutron_theta,neutron.Theta()*(180/math.pi))
        # neutron_mom = np.append(neutron_mom,neutron.E())
        neutron_mom = np.append(neutron_mom,EnE_Lab[entryNum])
    
        scat_electron_theta = np.append(scat_electron_theta,scat_electron.Theta()*(180/math.pi))
        # scat_electron_mom = np.append(scat_electron_mom,scat_electron.E())
        scat_electron_mom = np.append(scat_electron_mom,EeE_Lab[entryNum])
        
        Q2_cut = np.append(Q2_cut,Q2[entryNum])
        x_cut = np.append(x_cut,x[entryNum])
        t_cut = np.append(t_cut,t[entryNum])
        TDIS_znq_cut = np.append(TDIS_znq_cut,TDIS_znq[entryNum])
    
def plot_physics():
    
    f = plt.figure(figsize=(11.69,8.27))
    plt.style.use('default')
    plt.rcParams.update({'font.size': 15})
    
    ax = f.add_subplot(211)
    hist1 = ax.hist(pion_theta,bins=c.setbin(pion_theta,200),label='pi cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    hist1a = ax.hist(scat_electron_theta,bins=c.setbin(scat_electron_theta,200),label='e cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    hist1b = ax.hist(neutron_theta,bins=c.setbin(neutron_theta,200),label='n cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    ax.legend(loc=1)
    plt.title('theta', fontsize =20)

    ax = f.add_subplot(212)
    hist1 = ax.hist(pion_mom,bins=c.setbin(pion_mom,200),label='pi cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    hist1a = ax.hist(scat_electron_mom,bins=c.setbin(scat_electron_mom,200),label='e cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    hist1b = ax.hist(neutron_mom,bins=c.setbin(neutron_mom,200),label='n cut',histtype='step', alpha=0.5, stacked=True, fill=True)
    plt.title('tot_mom', fontsize =20)

    f.savefig('hist.png')
    
    xQ2 = c.densityPlot(x_cut,Q2_cut, '$Q^2$ vs $x_{Bj}$','$x_{Bj}$','$Q^2$', 200, 200,  c, 0.001, 1.0, 0., 100.0)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.001,1.)
    plt.ylim(10.,100.)
    plt.xlabel('$x_{Bj}$', fontsize =20)
    plt.ylabel('$Q^2$ ($GeV^2$)', fontsize =20)
    plt.title('$Q^2$ vs $x_{Bj}$', fontsize =20)

    xQ2[1].savefig('xQ2.png')

    pimomt = c.densityPlot(pion_mom,t_cut, 't vs pi tot_mom','tot_mom','t', 200, 200,  c)
    # plt.ylim(-180.,180.)
    # plt.xlim(0.,50.)
    plt.xlabel('tot_mom (GeV)', fontsize =20)
    plt.ylabel('t ($GeV^2$)', fontsize =20)
    plt.title('t vs pi tot_mom', fontsize =20)

    pimomt[1].savefig('pimomt.png')
    
    pithetat = c.densityPlot(pion_theta,t_cut, 't vs pi $\Theta$','$\Theta$','t', 200, 200,  c)
    # plt.ylim(-180.,180.)
    # plt.xlim(0.,50.)
    plt.xlabel('$\Theta$', fontsize =20)
    plt.ylabel('t ($GeV^2$)', fontsize =20)
    plt.title('t vs pi $\Theta$', fontsize =20)

    pithetat[1].savefig('pithetat.png')
    
    nmomt = c.densityPlot(neutron_mom,t_cut, 't vs n tot_mom','tot_mom','t', 200, 200,  c)
    # plt.ylim(-180.,180.)
    # plt.xlim(0.,50.)
    plt.xlabel('tot_mom (GeV)', fontsize =20)
    plt.ylabel('t ($GeV^2$)', fontsize =20)
    plt.title('t vs n tot_mom', fontsize =20)

    nmomt[1].savefig('nmomt.png')
    
    nthetat = c.densityPlot(neutron_theta,t_cut, 't vs n $\Theta$','$\Theta$','t', 200, 200,  c)
    # plt.ylim(-180.,180.)
    # plt.xlim(0.,50.)
    plt.xlabel('$\Theta$', fontsize =20)
    plt.ylabel('t ($GeV^2$)', fontsize =20)
    plt.title('t vs n $\Theta$', fontsize =20)

    nthetat[1].savefig('nthetat.png')
    
    phaseSpace = c.densityPlot(scat_electron_mom, scat_electron_theta, 'dir_z vs tot_mom','tot_mom','dir_z', 200, 200,  c)
    # plt.ylim(-180.,180.)
    # plt.xlim(0.,50.)
    plt.xlabel('tot_mom (GeV)', fontsize =20)
    plt.ylabel('Theta (deg)', fontsize =20)
    plt.title('e cut', fontsize =20)

    phaseSpace[1].savefig('ephaseSpace.png')
    
    phaseSpace = c.densityPlot(pion_mom, pion_theta, 'dir_z vs tot_mom','tot_mom','dir_z', 200, 200,  c)
    # plt.ylim(0.0,360.0)
    plt.xlim(0.,15.)
    plt.xlabel('tot_mom (GeV)', fontsize =20)
    plt.ylabel('Theta (deg)', fontsize =20)
    plt.title('pi cut', fontsize =20)

    phaseSpace[1].savefig('piphaseSpace.png')
    
    phaseSpace = c.densityPlot(neutron_mom, neutron_theta, 'dir_z vs tot_mom','tot_mom','dir_z', 200, 200,  c)
    # plt.ylim(0.,180.)
    # plt.xlim(200.,500.)
    plt.xlabel('tot_mom (GeV)', fontsize =20)
    plt.ylabel('Theta (deg)', fontsize =20)
    # plt.title('$\Lambda$ cut', fontsize =20)
    plt.title('n cut', fontsize =20)

    phaseSpace[1].savefig('nphaseSpace.png')
    
    # polar = c.polarPlot(pion_theta,pion_mom, 'P vs pi $\Theta$','$\Theta$','P')

    # polar.savefig('pipolar.png')
    
    polar = c.polarPlot(scat_electron_theta,scat_electron_mom, 'P vs scat e $\Theta$','$\Theta$','P [GeV]',0.,180.,9)

    polar.savefig('epolar.png')
        
    polar = c.polarPlot(neutron_theta,neutron_mom, 'P vs n $\Theta$','$\Theta$','P [GeV]',0.0,5.0,4)

    polar.savefig('npolar.png')
    
    # phaseSpace = c.densityPlot(pion_mom, neutron_mom, 'n vs pi','pi','n', 200, 200,  b,0,3.,0.,10.)
    # # plt.ylim(-180.,180.)
    # # plt.xlim(200.,500.)
    # plt.xlabel('pi tot_mom (GeV)', fontsize =20)
    # plt.ylabel('n tot_mom (GeV)', fontsize =20)
    # # plt.title('$\Lambda$ cut', fontsize =20)
    # plt.title('n vs pi', fontsize =20)

    # phaseSpace = c.densityPlot(pion_theta, neutron_theta, 'n $\Theta$ vs pi $\Theta$','pi $\Theta$','n $\Theta$', 200, 200,  c)
    # # plt.ylim(-180.,180.)
    # # plt.xlim(200.,500.)
    # plt.xlabel('pi $\Theta$', fontsize =20)
    # plt.ylabel('n $\Theta$', fontsize =20)
    # # plt.title('$\Lambda$ cut', fontsize =20)
    # plt.title('n $\Theta$ vs pi $\Theta$', fontsize =20)

meta  = up.open(rootName)["Meta"]
metaBranch = r2p.pyBranch(meta)

p2_pt = metaBranch.findBranch("P2","p2_pt")
p2_z = metaBranch.findBranch("P2","p2_z")
phi_p2 = metaBranch.findBranch("P2","phi_p2")
theta_p2 = metaBranch.findBranch("P2","theta_p2")
Pz_p2 = metaBranch.findBranch("P2","Pz_p2")
P_p2 = metaBranch.findBranch("P2","P_p2")
E_p2 = metaBranch.findBranch("P2","E_p2")
Px_p2 = metaBranch.findBranch("P2","Px_p2")
Py_p2 = metaBranch.findBranch("P2","Py_p2")

def plot_meta():
    f = plt.figure(figsize=(11.69,8.27))
    plt.style.use('default')

    ax = f.add_subplot(211)
    ax.hist(Pz_p2,bins=c.setbin(Pz_p2,200),label='Pz_p2',histtype='step', alpha=0.5, stacked=True, fill=True)
    ax.hist(P_p2,bins=c.setbin(P_p2,200),label='P_p2',histtype='step', alpha=0.5, stacked=True, fill=True)
    ax.legend(loc=1)
    plt.title('P(z)_p2', fontsize =20)

    c.densityPlot(P_p2,theta_p2, 'theta_p2 vs P_p2','P_p2','theta_p2', 200, 200,  c)
    # plt.ylim(-180.,180.)
    # plt.xlim(0.,50.)
    plt.xlabel('P_p2', fontsize =20)
    plt.ylabel('theta_p2', fontsize =20)
    plt.title('theta_p2 vs P_p2', fontsize =20)

    c.densityPlot(Pz_p2,theta_p2, 'theta_p2 vs Pz_p2','Pz_p2','theta_p2', 200, 200,  c)
    # plt.ylim(-180.,180.)
    # plt.xlim(0.,50.)
    plt.xlabel('Pz_p2', fontsize =20)
    plt.ylabel('theta_p2', fontsize =20)
    plt.title('theta_p2 vs Pz_p2', fontsize =20)

    c.densityPlot(P_p2,pion_theta, 'pion_theta vs P_p2','P_p2','pion_theta', 200, 200,  c)
    # plt.ylim(-180.,180.)
    # plt.xlim(0.,50.)
    plt.xlabel('P_p2', fontsize =20)
    plt.ylabel('pion_theta', fontsize =20)
    plt.title('pion_theta vs P_p2', fontsize =20)

    c.densityPlot(P_p2,pion_mom, 'pion_mom vs P_p2','P_p2','pion_mom', 200, 200,  c)
    # plt.ylim(-180.,180.)
    # plt.xlim(0.,50.)
    plt.xlabel('P_p2', fontsize =20)
    plt.ylabel('pion_mom', fontsize =20)
    plt.title('pion_mom vs P_p2', fontsize =20)

    c.densityPlot(P_p2,neutron_theta, 'neutron_theta vs P_p2','P_p2','neutron_theta', 200, 200,  c)
    # plt.ylim(0,5000.)
    # plt.xlim(0.,5.)
    plt.xlabel('P_p2', fontsize =20)
    plt.ylabel('neutron_theta', fontsize =20)
    plt.title('neutron_theta vs P_p2', fontsize =20)
    
    c.densityPlot(P_p2,neutron_mom, 'neutron_mom vs P_p2','P_p2','neutron_mom', 200, 200,  b,0,5,0,5000)
    plt.ylim(0,5000.)
    plt.xlim(0.,5.)
    plt.xlabel('P_p2', fontsize =20)
    plt.ylabel('neutron_mom', fontsize =20)
    plt.title('neutron_mom vs P_p2', fontsize =20)

        
def main() :

    # plot_meta()
    plot_physics()
    plt.show()
    
if __name__=='__main__': main()
