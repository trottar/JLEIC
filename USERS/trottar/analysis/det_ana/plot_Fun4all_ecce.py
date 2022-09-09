#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2022-05-20 12:41:03 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import MaxNLocator
import uproot as up
import numpy as np
import pandas as pd
import sys
import ltsep as lt

numEvts = sys.argv[1]
IP = sys.argv[2]

#e_list = ["5on41","5on100","10on100","10on135","18on275"]
e_list = ["5on41","10on100"]

rootName = [None]*len(e_list)
tdata = [None]*len(e_list)
print("Energy settings...")
for i,e in enumerate(e_list):
    rootName[i]="./INPUTS/{0}_{1}_{2}.root".format(e,numEvts,IP) # Note: this is relative to the bash script NOT this python script!
    print(rootName[i])
    tdata[i] = up.open(rootName[i])    

def plot_mom(scat_eTruth_flag=True,nTruth_flag=True,n_flag=True):
    '''
    Plot momentum
    '''
    
    '''
    scattered eTruth
    '''

    if scat_eTruth_flag==True:
    
        scat_eTruth = [None]*len(e_list)
        epxpy = [None]*len(e_list)
        epxpz = [None]*len(e_list)
        epypz = [None]*len(e_list)
        for i,e in enumerate(e_list):
            scat_eTruth[i] = tdata[i]['Scattered_Electron_Truth_Info']

            epxpy[i] = scat_eTruth[i]['eTruth_pxpy'].to_hist()
            epxpz[i] = scat_eTruth[i]['eTruth_pxpz'].to_hist()
            epypz[i] = scat_eTruth[i]['eTruth_pypz'].to_hist()

        f = plt.figure(figsize=(11.69,8.27))
        
        for i,e in enumerate(e_list):
            
            ax = f.add_subplot(3,len(e_list),i+1)
            epxpy[i].plot(cbar=False,cmin=1,ax=ax)
            ax.set_ylim(-10,10)
            ax.set_xlim(-10,10)
            ax.set_yticks(range(-5,5+1,10))
            ax.set_xticks(range(-5,5+1,10))
            if i==0 :
                plt.ylabel(r"e'$\frac{\Delta p_{y}}{Truth p_{y}}$",fontsize=20)
            else:
                plt.ylabel("",fontsize=20)
                ax.yaxis.set_major_formatter(plt.NullFormatter())
            plt.xlabel(r"e'$\frac{\Delta p_{x}}{Truth p_{x}}$",fontsize=20)
            plt.annotate("{}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

            ax = f.add_subplot(3,len(e_list),i+3)
            epxpz[i].plot(cbar=False,cmin=1,ax=ax)
            ax.set_ylim(-10,10)
            ax.set_xlim(-10,10)
            ax.set_yticks(range(-5,5+1,10))
            ax.set_xticks(range(-5,5+1,10))
            if i==0 :
                plt.ylabel(r"e'$\frac{\Delta p_{z}}{Truth p_{z}}$",fontsize=20)
            else:
                plt.ylabel("",fontsize=20)
                ax.yaxis.set_major_formatter(plt.NullFormatter())
            plt.xlabel(r"e'$\frac{\Delta p_{x}}{Truth p_{x}}$",fontsize=20)

            ax = f.add_subplot(3,len(e_list),i+5)
            epypz[i].plot(cbar=False,cmin=1,ax=ax)
            ax.set_ylim(-10,10)
            ax.set_xlim(-10,10)
            ax.set_yticks(range(-5,5+1,10))
            ax.set_xticks(range(-5,5+1,10))
            if i==0 :
                plt.ylabel(r"e'$\frac{\Delta p_{z}}{Truth p_{z}}$",fontsize=20)
            else:
                plt.ylabel("",fontsize=20)
                ax.yaxis.set_major_formatter(plt.NullFormatter())
            plt.xlabel(r"e'$\frac{\Delta p_{y}}{Truth p_{y}}$",fontsize=20)            
                        
        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.subplots_adjust(wspace=0.0)
        plt.savefig('OUTPUTS/MOM/scat_eTruth_{}.png'.format(IP))
        plt.close(f)
        
    '''
    nTruth
    '''

    if nTruth_flag==True:
    
        scat_nTruth = [None]*len(e_list)
        npxpy = [None]*len(e_list)
        npxpz = [None]*len(e_list)
        npypz = [None]*len(e_list)
        for i,e in enumerate(e_list):
            scat_nTruth[i] = tdata[i]['Neutron_Truth_Info']

            npxpy[i] = scat_nTruth[i]['nTruth_pxpy'].to_hist()
            npxpz[i] = scat_nTruth[i]['nTruth_pxpz'].to_hist()
            npypz[i] = scat_nTruth[i]['nTruth_pypz'].to_hist()

        f = plt.figure(figsize=(11.69,8.27))
        
        for i,e in enumerate(e_list):
            
            ax = f.add_subplot(3,len(e_list),i+1)
            npxpy[i].plot(cbar=False,cmin=1,ax=ax)
            ax.set_ylim(-10,10)
            ax.set_xlim(-10,10)
            ax.set_yticks(range(-5,5+1,10))
            ax.set_xticks(range(-5,5+1,10))
            if i==0 :
                plt.ylabel(r"n$\frac{\Delta p_{y}}{Truth p_{y}}$",fontsize=20)
            else:
                plt.ylabel("",fontsize=20)
                ax.yaxis.set_major_formatter(plt.NullFormatter())
            plt.xlabel(r"n$\frac{\Delta p_{x}}{Truth p_{x}}$",fontsize=20)
            plt.annotate("{}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

            ax = f.add_subplot(3,len(e_list),i+3)
            npxpz[i].plot(cbar=False,cmin=1,ax=ax)
            ax.set_ylim(-10,10)
            ax.set_xlim(-10,10)
            ax.set_yticks(range(-5,5+1,10))
            ax.set_xticks(range(-5,5+1,10))
            if i==0 :
                plt.ylabel(r"n$\frac{\Delta p_{z}}{Truth p_{z}}$",fontsize=20)
            else:
                plt.ylabel("",fontsize=20)
                ax.yaxis.set_major_formatter(plt.NullFormatter())
            plt.xlabel(r"n$\frac{\Delta p_{x}}{Truth p_{x}}$",fontsize=20)

            ax = f.add_subplot(3,len(e_list),i+5)
            npypz[i].plot(cbar=False,cmin=1,ax=ax)
            ax.set_ylim(-10,10)
            ax.set_xlim(-10,10)
            ax.set_yticks(range(-5,5+1,10))
            ax.set_xticks(range(-5,5+1,10))
            if i==0 :
                plt.ylabel(r"n$\frac{\Delta p_{z}}{Truth p_{z}}$",fontsize=20)
            else:
                plt.ylabel("",fontsize=20)
                ax.yaxis.set_major_formatter(plt.NullFormatter())
            plt.xlabel(r"n$\frac{\Delta p_{y}}{Truth p_{y}}$",fontsize=20)            
            
        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.subplots_adjust(wspace=0.0)
        plt.savefig('OUTPUTS/MOM/scat_nTruth_{}.png'.format(IP))
        plt.close(f)
        
        '''
        Neutron kinematics
        '''

    if n_flag==True:        

        nInfo = [None]*len(e_list)
        n_px = [None]*len(e_list)
        n_py = [None]*len(e_list)
        n_pz = [None]*len(e_list)
        n_p = [None]*len(e_list)
        n_E = [None]*len(e_list)
        n_Theta = [None]*len(e_list)
        n_Phi = [None]*len(e_list)
        nTrack_ThetaPhi = [None]*len(e_list)
        nTrack_pTheta = [None]*len(e_list)

        for i,e in enumerate(e_list):
            nInfo[i] = tdata[i]['Neutron_Info']
            n_px[i] = nInfo[i]['n_px'].to_hist()
            n_py[i] = nInfo[i]['n_py'].to_hist()
            n_pz[i] = nInfo[i]['n_pz'].to_hist()
            n_p[i] = nInfo[i]['n_p'].to_hist()
            n_E[i] = nInfo[i]['n_E'].to_hist()
            n_Theta[i] = nInfo[i]['n_Theta'].to_hist()
            n_Phi[i] = nInfo[i]['n_Phi'].to_hist()
            nTrack_ThetaPhi[i] = nInfo[i]['nTrack_ThetaPhi'].to_hist()
            nTrack_pTheta[i] = nInfo[i]['nTrack_pTheta'].to_hist()

        f = plt.figure(figsize=(11.69,8.27))            
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            n_px[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('n_px',fontsize=20)
            plt.ylabel('counts',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/MOM/n_px_{}.png'.format(IP))
        plt.clf()
        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            n_py[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('n_py',fontsize=20)
            plt.ylabel('counts',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/MOM/n_py_{}.png'.format(IP))        
        plt.clf()
        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            n_pz[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('n_pz',fontsize=20)
            plt.ylabel('counts',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/MOM/n_pz_{}.png'.format(IP))        
        plt.clf()
        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            n_p[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('n_p',fontsize=20)
            plt.ylabel('counts',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/MOM/n_p_{}.png'.format(IP))
        plt.clf()
        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            n_E[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('n_E',fontsize=20)
            plt.ylabel('counts',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/MOM/n_E_{}.png'.format(IP))
        plt.clf()
        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            n_Theta[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('$\Theta$ [deg]',fontsize=20)
            plt.ylabel('counts',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/MOM/n_Theta_{}.png'.format(IP))
        plt.clf()
        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            n_Phi[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('$\Phi$',fontsize=20)
            plt.ylabel('counts',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/MOM/n_Phi_{}.png'.format(IP))
        plt.clf()
        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            nTrack_ThetaPhi[i].plot(cmin=5,ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('$\Theta$ [deg]',fontsize=20)
            plt.ylabel('$\Phi$ [deg]',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/MOM/nTrack_ThetaPhi_{}.png'.format(IP))
        plt.clf()
        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            nTrack_pTheta[i].plot(cmin=5,ax=ax)
            ax.set_ylim(0,120)
            #ax.set_xlim(0,5)
            ax.set_yticks(range(0,120+1,10))
            #ax.set_xticks(range(0,5+1,10))
            plt.xlabel('$\Theta$ [deg]',fontsize=20)
            plt.ylabel('P [GeV/c]',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))
            #plt.gca().set_aspect('equal')
            #plt.axis('scaled')
            #ax.set_aspect(abs(500)/abs(5000))
            ax.set_aspect('auto')

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/MOM/nTrack_pTheta_{}.png'.format(IP))
        plt.close(f)
        
def plot_kin(k1d_flag=True,k1dT_flag=True,k2_flag=True,k3_flag=True):

    if k1d_flag==True:        

        Q2_Dist = [None]*len(e_list)
        W_Dist = [None]*len(e_list)
        t_Dist = [None]*len(e_list)
        xb_Dist = [None]*len(e_list)
        xi_Dist = [None]*len(e_list)
        Delta_t = [None]*len(e_list)

        for i,e in enumerate(e_list):
            Q2_Dist[i] = tdata[i]['Kinematics_Info/Q2_Dist'].to_hist()
            W_Dist[i] = tdata[i]['Kinematics_Info/W_Dist'].to_hist()
            t_Dist[i] = tdata[i]['Kinematics_Info/t_Dist'].to_hist()
            xb_Dist[i] = tdata[i]['Kinematics_Info/xb_Dist'].to_hist()
            xi_Dist[i] = tdata[i]['Kinematics_Info/xi_Dist'].to_hist()
            Delta_t[i] = tdata[i]['Kinematics_Info/Delta_t'].to_hist()

        f = plt.figure(figsize=(11.69,8.27))
        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            Q2_Dist[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('Q2_Dist',fontsize=20)
            plt.ylabel('counts',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/KIN/Q2_Dist_{}.png'.format(IP))
        plt.clf()

        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            Delta_t[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('$\Delta$t',fontsize=20)
            plt.ylabel('counts',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/KIN/Delta_t_{}.png'.format(IP))
        plt.clf()

        
        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            W_Dist[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('W_Dist',fontsize=20)
            plt.ylabel('counts',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/KIN/W_Dist_{}.png'.format(IP))
        plt.clf()

        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            t_Dist[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('t_Dist',fontsize=20)
            plt.ylabel('counts',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/KIN/t_Dist_{}.png'.format(IP))
        plt.clf()

        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            xb_Dist[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('xb_Dist',fontsize=20)
            plt.ylabel('counts',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/KIN/xb_Dist_{}.png'.format(IP))
        plt.clf()
        
        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            xi_Dist[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('xi_Dist',fontsize=20)
            plt.ylabel('counts',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/KIN/xi_Dist_{}.png'.format(IP))
        plt.close(f)

    if k1dT_flag==True:        

        Q2Truth_Dist = [None]*len(e_list)
        WTruth_Dist = [None]*len(e_list)
        tTruth_Dist = [None]*len(e_list)
        xbTruth_Dist = [None]*len(e_list)
        xiTruth_Dist = [None]*len(e_list)
        t_tTruth = [None]*len(e_list)
        Q2_zdc_eff = [None]*len(e_list)

        for i,e in enumerate(e_list):
            Q2Truth_Dist[i] = tdata[i]['Kinematics_Truth_Info/Q2Truth_Dist'].to_hist()
            WTruth_Dist[i] = tdata[i]['Kinematics_Truth_Info/WTruth_Dist'].to_hist()
            tTruth_Dist[i] = tdata[i]['Kinematics_Truth_Info/tTruth_Dist'].to_hist()
            xbTruth_Dist[i] = tdata[i]['Kinematics_Truth_Info/xbTruth_Dist'].to_hist()
            xiTruth_Dist[i] = tdata[i]['Kinematics_Truth_Info/xiTruth_Dist'].to_hist()
            t_tTruth[i] = tdata[i]['Kinematics_Analysis/t_tTruth'].to_hist()
            Q2_zdc_eff[i] = tdata[i]['Kinematics_Analysis/Q2_zdc_eff'].to_hist()

        f = plt.figure(figsize=(11.69,8.27))

        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            t_tTruth[i].plot(cmin=5,ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('t',fontsize=20)
            plt.ylabel('tTruth',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))
            plt.axis('scaled')

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/KIN/t_tTruth_{}.png'.format(IP))
        plt.clf()
        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            Q2_zdc_eff[i].plot(cmin=5,ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('zdc eff',fontsize=20)
            plt.ylabel('$Q^{2}$',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/KIN/Q2_zdc_eff_{}.png'.format(IP))
        plt.clf()
        
        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            Q2Truth_Dist[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('Q2Truth_Dist',fontsize=20)
            plt.ylabel('counts',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/KIN/Q2Truth_Dist_{}.png'.format(IP))
        plt.clf()

        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            WTruth_Dist[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('WTruth_Dist',fontsize=20)
            plt.ylabel('counts',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/KIN/WTruth_Dist_{}.png'.format(IP))
        plt.clf()

        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            tTruth_Dist[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('tTruth_Dist',fontsize=20)
            plt.ylabel('counts',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/KIN/tTruth_Dist_{}.png'.format(IP))
        plt.clf()

        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            xbTruth_Dist[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('xbTruth_Dist',fontsize=20)
            plt.ylabel('counts',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/KIN/xbTruth_Dist_{}.png'.format(IP))
        plt.clf()
        
        
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            xiTruth_Dist[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('xiTruth_Dist',fontsize=20)
            plt.ylabel('counts',fontsize=20)
            plt.annotate("{0}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/KIN/xiTruth_Dist_{}.png'.format(IP))
        plt.close(f)
        
    if k2_flag==True:

        arrQ2 = [7,15,30,60]
        t_Q2 = [[None for _ in range(len(arrQ2))] for _ in range(len(e_list))]
        delta_t_t_Q2 = [[None for _ in range(len(arrQ2))] for _ in range(len(e_list))]
        zdc_eff_Q2 = [[None for _ in range(len(arrQ2))] for _ in range(len(e_list))]

        for i,e in enumerate(e_list):
            for j,val in enumerate(arrQ2):
                t_Q2[i][j] = tdata[i]['Kinematics_Analysis/t_Q2_{}'.format(val)].to_hist()
                delta_t_t_Q2[i][j] = tdata[i]['Kinematics_Analysis/delta_t_t_Q2_{}'.format(val)].to_hist()
                zdc_eff_Q2[i][j] = tdata[i]['Kinematics_Analysis/zdc_eff_Q2_{}'.format(val)].to_hist()
            
        t_Q2_Dict = [None]*len(e_list)
        df = [None]*len(e_list)

        f = plt.figure(figsize=(11.69,8.27))
        for i,e in enumerate(e_list):
            t_Q2_Dict[i] = {}
            for j,val in enumerate(arrQ2):
                t_Q2_Dict[i].update({"{0} < $Q^2$ < {1}".format(val-5,val+5) : tdata[i]['Kinematics_Analysis/t_Q2_{}'.format(val)].values()})
            #t_Q2_Dict[i] = {k : t_Q2_Dict[i][k] for k in sorted(t_Q2_Dict[i].keys())}
            df[i] = pd.DataFrame(t_Q2_Dict[i])
            print("-------------------------------------\n{}".format(e))
            print(df[i])

            ax = f.add_subplot(int("12{}".format(i+1)),projection='3d')

            #determine the number of columns
            ncol = df[i].shape[1]

            #define the yticks, i.e., the column numbers
            yticks = np.arange(ncol)

            #we create evenly spaced bins between the minimum and maximum of the entire dataframe
            xbins = np.linspace(0, 0.3, 201)
            
            #and calculate the center and widths of the bars
            xcenter = np.convolve(xbins, np.ones(2), "valid")/2
            xwidth = np.diff(xbins)

            #calculate now the histogram and plot it for each column
            for ytick in yticks:

                #extract the current column from your df[i] by its number
                col =  df[i].iloc[:, ytick].tolist()

                ax.bar(left=xcenter, height=col, width=xwidth, zs=ytick, zdir="y", alpha=0.75, color='blue')

            ax.set_xlabel("-t",fontsize=20)
            ax.set_zlabel("counts",fontsize=20)

            ylabel = list(df[i])
            ax.set_yticks(yticks)
            ax.set_yticklabels(ylabel)
            plt.setp(ax.yaxis.get_majorticklabels(), rotation=-15, ha="left", rotation_mode="anchor")
            plt.annotate("{}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))
            
        plt.savefig('OUTPUTS/KIN/t_Q2_3d_{}.png'.format(IP))
        plt.clf()
        plt.close(f)
        
        f = plt.figure(figsize=(11.69,8.27))
        print('\nt_Q2_{}'.format(IP))
        for i,e in enumerate(e_list):
        
            numfig = (len(arrQ2)//2)+(len(arrQ2)%2>0)
            for j in range(1,(numfig+1)):
                ax = f.add_subplot(int("{0}{1}{2}".format(numfig,2,j)))
                t_Q2[i][j-1].plot(ax=ax)
                if j==1 :
                    plt.ylabel("Counts",fontsize=20)
                else:
                    plt.ylabel("",fontsize=20)
                    ax.yaxis.set_major_formatter(plt.NullFormatter())
                plt.xlabel("",fontsize=20)
                ax.xaxis.set_major_formatter(plt.NullFormatter())
                plt.annotate("{0} < $Q^2$ < {1}".format(arrQ2[j-1]-5,arrQ2[j-1]+5), xy=(0.50, 0.90), xycoords='axes fraction',fontsize=20)
                
                ax = f.add_subplot(int("{0}{1}{2}".format(numfig,2,j+numfig)))
                t_Q2[i][j+numfig-1].plot(ax=ax)
                if j==1 :
                    plt.ylabel("Counts",fontsize=20)
                else:
                    plt.ylabel("",fontsize=20)
                    ax.yaxis.set_major_formatter(plt.NullFormatter())
                plt.xlabel("t",fontsize=20)
                plt.annotate("{0} < $Q^2$ < {1}".format(arrQ2[j+numfig-1]-5,arrQ2[j+numfig-1]+5), xy=(0.50, 0.90), xycoords='axes fraction',fontsize=20)
                
            plt.subplots_adjust(hspace=0.0,wspace=0.0)
            plt.savefig('OUTPUTS/KIN/t_Q2_{0}_{1}.png'.format(e,IP))
            plt.clf()
            lt.Misc.progressBar(i,len(e_list)-1)

        print('\ndelta_t_t_Q2_{}'.format(IP))            
        for i,e in enumerate(e_list):

            numfig = (len(arrQ2)//2)+(len(arrQ2)%2>0)
            for j in range(1,(numfig+1)):
                ax = f.add_subplot(int("{0}{1}{2}".format(numfig,2,j)))
                delta_t_t_Q2[i][j-1].plot(cmin=1,cbar=False,ax=ax)
                if j==1 :
                    plt.ylabel("-t",fontsize=20)
                    ax.xaxis.set_major_locator(MaxNLocator(prune='upper'))
                    plt.annotate("{0}".format(e), xy=(0.10, 0.90), xycoords='axes fraction',fontsize=24)
                    plt.annotate("{0} < $Q^2$ < {1}".format(arrQ2[j-1]-5,arrQ2[j-1]+5), xy=(0.50, 0.90), xycoords='axes fraction',fontsize=20)
                else:
                    plt.ylabel("",fontsize=20)
                    ax.yaxis.set_major_formatter(plt.NullFormatter())
                    ax.xaxis.set_major_locator(MaxNLocator(prune='lower'))
                    plt.annotate("{0} < $Q^2$ < {1}".format(arrQ2[j-1]-5,arrQ2[j-1]+5), xy=(0.50, 0.90), xycoords='axes fraction',fontsize=20)
                plt.xlabel("",fontsize=20)
                ax.xaxis.set_major_formatter(plt.NullFormatter())
                
                ax = f.add_subplot(int("{0}{1}{2}".format(numfig,2,j+numfig)))
                delta_t_t_Q2[i][j+numfig-1].plot(cmin=1,cbar=False,ax=ax)
                if j==1 :
                    plt.ylabel("-t",fontsize=20)
                    ax.xaxis.set_major_locator(MaxNLocator(prune='upper')) 
                else:
                    plt.ylabel("",fontsize=20)
                    ax.yaxis.set_major_formatter(plt.NullFormatter())
                    ax.xaxis.set_major_locator(MaxNLocator(prune='lower')) 
                plt.xlabel(r"$\Delta$t",fontsize=20)
                plt.annotate("{0} < $Q^2$ < {1}".format(arrQ2[j+numfig-1]-5,arrQ2[j+numfig-1]+5), xy=(0.50, 0.90), xycoords='axes fraction',fontsize=20)
            plt.subplots_adjust(hspace=0.0,wspace=0.0)
            plt.savefig('OUTPUTS/KIN/delta_t_t_Q2_{0}_{1}.png'.format(e,IP))
            plt.clf()
            lt.Misc.progressBar(i,len(e_list)-1)

        print('\nzdc_eff_Q2_{}'.format(IP))            
        for i,e in enumerate(e_list):

            numfig = (len(arrQ2)//2)+(len(arrQ2)%2>0)
            for j in range(1,(numfig+1)):
                ax = f.add_subplot(int("{0}{1}{2}".format(numfig,2,j)))
                zdc_eff_Q2[i][j-1].plot(cmin=1,ax=ax)
                if j==1 :
                    plt.ylabel("ZDC Y",fontsize=20)
                else:
                    plt.ylabel("",fontsize=20)
                    ax.yaxis.set_major_formatter(plt.NullFormatter())
                plt.xlabel("",fontsize=20)
                ax.xaxis.set_major_formatter(plt.NullFormatter())
                plt.annotate("{0} < $Q^2$ < {1}".format(arrQ2[j-1]-5,arrQ2[j-1]+5), xy=(0.50, 0.90), xycoords='axes fraction',fontsize=20)
                
                ax = f.add_subplot(int("{0}{1}{2}".format(numfig,2,j+numfig)))
                zdc_eff_Q2[i][j+numfig-1].plot(cmin=1,ax=ax)
                if j==1 :
                    plt.ylabel("ZDC Y",fontsize=20)
                else:
                    plt.ylabel("",fontsize=20)
                    ax.yaxis.set_major_formatter(plt.NullFormatter())
                plt.xlabel(r"ZDC X",fontsize=20)
                plt.annotate("{0} < $Q^2$ < {1}".format(arrQ2[j+numfig-1]-5,arrQ2[j+numfig-1]+5), xy=(0.50, 0.90), xycoords='axes fraction',fontsize=20)
            plt.subplots_adjust(hspace=0.1,wspace=0.1)
            plt.savefig('OUTPUTS/KIN/zdc_eff_Q2_{0}_{1}.png'.format(e,IP))
            plt.clf()
            lt.Misc.progressBar(i,len(e_list)-1)
        plt.close(f)
        
    if k3_flag==True:
        
        ZDC_delta_t = [None]*len(e_list)

        for i,e in enumerate(e_list):
            ZDC_delta_t[i] = tdata[i]['Kinematics_Analysis/ZDC_delta_t'].to_hist()
            #ZDC_delta_t[i] = tdata[i]['Kinematics_Analysis/ZDC_delta_t'].to_numpy()

        ZDC_X = [None]*len(e_list)
        ZDC_Y = [None]*len(e_list)
        for i,e in enumerate(e_list):
            ZDC_X[i] = tdata[i]['ZDC/ZDC_X'].values()
            ZDC_Y[i] = tdata[i]['ZDC/ZDC_Y'].values()
        
        delta_t = [None]*len(e_list)
        for i,e in enumerate(e_list):
            delta_t[i] = tdata[i]['Kinematics_Info/Delta_t'].values()

        f = plt.figure(figsize=(11.69,8.27))
        print('\nzdc_x_t_accept_{}'.format(IP))
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12{}".format(i+1)),projection='3d')

            #define the colormap 
            my_cmap = plt.cm.viridis
            
            x = tdata[i]['ZDC/ZDC_X'].to_numpy()[0]
            y = tdata[i]['ZDC/ZDC_Y'].to_numpy()[0]
            t = tdata[i]['Kinematics_Info/t_Dist'].to_numpy()[0]
            dx = tdata[i]['ZDC/ZDC_X'].to_numpy()[1]
            dy = tdata[i]['ZDC/ZDC_Y'].to_numpy()[1]
            dt = tdata[i]['Kinematics_Info/t_Dist'].to_numpy()[1]
            
            #we create evenly spaced bins between the minimum and maximum of the entire dataframe
            tbins = dt
            
            #and calculate the center and widths of the bars
            tcenter = np.convolve(tbins, np.ones(2), "valid")/2
            twidth = np.diff(tbins)

            j=0
            print("{}".format(e))
            #calculate now the histogram and plot it for each column
            for xi,di in zip(x,dx):
                lt.Misc.progressBar(j,len(x)-1)
                accept_map = xi/max(x)
                ax.bar(left=tcenter, height=t, width=twidth, zs=di, zdir="y", color=my_cmap(accept_map), alpha=accept_map)
                j+=1

            plt.annotate("{}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))
            ax.set_xlabel("-t",fontsize=20)
            ax.set_ylabel("x",fontsize=20)
            ax.set_zlabel("counts",fontsize=20)
        plt.savefig('OUTPUTS/KIN/zdc_x_t_accept_{}.png'.format(IP))
        plt.clf()
        
        #f = plt.figure(figsize=(11.69,8.27))
        print('\nQ2_t_accept_{}'.format(IP))
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12{}".format(i+1)),projection='3d')

            #define the colormap 
            my_cmap = plt.cm.viridis
            
            x = tdata[i]['ZDC/ZDC_X'].to_numpy()[0]
            y = tdata[i]['ZDC/ZDC_Y'].to_numpy()[0]
            q = tdata[i]['Kinematics_Info/Q2_Dist'].to_numpy()[0]
            t = tdata[i]['Kinematics_Info/t_Dist'].to_numpy()[0]
            delta_t = tdata[i]['Kinematics_Info/Delta_t'].to_numpy()[0]
            dx = tdata[i]['ZDC/ZDC_X'].to_numpy()[1]
            dq = tdata[i]['Kinematics_Info/Q2_Dist'].to_numpy()[1]
            dy = tdata[i]['ZDC/ZDC_Y'].to_numpy()[1]
            dt = tdata[i]['Kinematics_Info/t_Dist'].to_numpy()[1]
            d_delta_t = tdata[i]['Kinematics_Info/Delta_t'].to_numpy()[1]
            #we create evenly spaced bins between the minimum and maximum of the entire dataframe
            tbins = dt
            
            #and calculate the center and widths of the bars
            tcenter = np.convolve(tbins, np.ones(2), "valid")/2
            twidth = np.diff(tbins)
            
            j=0
            print("{}".format(e))
            maxDelBin = []
            #calculate now the histogram and plot it for each column
            for del_ti,qi,di in zip(delta_t,q,dq):
                lt.Misc.progressBar(j,len(delta_t)-1)
                x_max_index = np.argmax(x)
                #accept_map = qi/(q[x_max_index]*10)
                accept_map = del_ti/max(delta_t)
                if del_ti > 0 and d_delta_t[j] < 0.02:
                    maxDelBin.append(d_delta_t[j])
                ax.bar(left=tcenter, height=t, width=twidth, zs=di, zdir="y", color=my_cmap(accept_map), alpha=accept_map)
                j+=1
            plt.annotate("{}".format(e),fontsize=20,xycoords='axes fraction',xy=(0.50, 0.90))
            ax.set_xlabel("-t",fontsize=20)
            ax.set_ylabel("Q2",fontsize=20)
            ax.set_zlabel("counts",fontsize=20)
        plt.savefig('OUTPUTS/KIN/Q2_t_accept_{}.png'.format(IP))
        plt.clf()
        plt.close(f)
        print("-"*50,"\n")
            
        
def plot_det(zdc_flag=True,rp_flag=True,offmom_flag=True,b0_flag=True):
    '''
    Plot detectors
    '''
    '''
    ZDC
    '''

    if zdc_flag==True:
        ZDC_XY = [None]*len(e_list)
        for i,e in enumerate(e_list):
            ZDC_XY[i] = tdata[i]['ZDC/ZDC_XY_l'].to_hist()

        f = plt.figure(figsize=(11.69,8.27))
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            # cmin is the minimum value to include in histogram (cmax can be used for a max limit), here I am excluding 0 values
            print(tdata[i]['ZDC_XY'])
            ZDC_XY[i].plot(cmin=5,ax=ax)
            ax.set_ylim(-50,50)
            ax.set_xlim(-50,50)
            ax.set_yticks(range(-50,50+1,10))
            ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('X [cm]',fontsize=20)
            plt.ylabel('Y [cm]',fontsize=20)
            zdc_counts = ZDC_XY[i].counts()
            zdc_counts = sum(zdc_counts[~np.isnan(zdc_counts)])
            zdc_accp = 100*(zdc_counts/float(numEvts))
            if zdc_accp > 100:
                zdc_accp = 100
            plt.annotate("{0}, {1:3.3f}% acceptance".format(e,zdc_accp),fontsize=17,xycoords='axes fraction',xy=(0.05, 0.95))
            plt.axis('scaled')

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/DET/zdc_all_{}.png'.format(IP))
        plt.clf()
        
    '''
    RP
    '''

    if rp_flag==True:
        RP_XY = [None]*len(e_list)
        for i,e in enumerate(e_list):
            RP_XY[i] = tdata[i]['RP/RP_XY_l'].to_hist()

        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            RP_XY[i].plot(cmin=5,ax=ax)
            ax.set_ylim(-50,50)
            ax.set_xlim(-50,50)
            ax.set_yticks(range(-50,50+1,10))
            ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('X [cm]',fontsize=20)
            plt.ylabel('Y [cm]',fontsize=20)
            rp_counts = RP_XY[i].counts()
            rp_counts = sum(rp_counts[~np.isnan(rp_counts)])
            rp_accp = 100*(rp_counts/float(numEvts))
            plt.annotate("{0}, {1:3.3f}% occupancy".format(e,rp_accp),fontsize=17,xycoords='axes fraction',xy=(0.05, 0.95))
            plt.axis('scaled')

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/DET/rp_all_{}.png'.format(IP))
        plt.clf()
        
    '''
    OFFMOM
    '''

    if offmom_flag==True:
        OFFMOM_XY = [None]*len(e_list)
        for i,e in enumerate(e_list):
            OFFMOM_XY[i] = tdata[i]['OFFMOM_XY'].to_hist()

        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            OFFMOM_XY[i].plot(cmin=5,ax=ax)
            ax.set_ylim(-50,50)
            ax.set_xlim(-50,50)
            ax.set_yticks(range(-50,50+1,10))
            ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('X [cm]',fontsize=20)
            plt.ylabel('Y [cm]',fontsize=20)
            offmom_counts = OFFMOM_XY[i].counts()
            offmom_counts = sum(offmom_counts[~np.isnan(offmom_counts)])
            offmom_accp = 100*(offmom_counts/float(numEvts))
            plt.annotate("{0}, {1:3.3f}% occupancy".format(e,offmom_accp),fontsize=17,xycoords='axes fraction',xy=(0.05, 0.95))
            plt.axis('scaled')

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/DET/offmom_all_{}.png'.format(IP))
        plt.clf()

    '''
    B0
    '''

    if b0_flag==True:
        B0_XY = [None]*len(e_list)
        for i,e in enumerate(e_list):
            B0_XY[i] = tdata[i]['B0/B0_XY_l'].to_hist()

        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("12%i" % (i+1)))
            B0_XY[i].plot(cmin=5,ax=ax)
            ax.set_ylim(-25,25)
            ax.set_xlim(-25,25)
            ax.set_yticks(range(-25,25+1,10))
            ax.set_xticks(range(-25,25+1,10))
            plt.xlabel('X [cm]',fontsize=20)
            plt.ylabel('Y [cm]',fontsize=20)
            b0_counts = B0_XY[i].counts()
            b0_counts = sum(b0_counts[~np.isnan(b0_counts)])
            b0_accp = 100*(b0_counts/float(numEvts))
            plt.annotate("{0}, {1:3.3f}% occupancy".format(e,b0_accp),fontsize=17,xycoords='axes fraction',xy=(0.05, 0.95))
            plt.axis('scaled')

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/DET/b0_all_{}.png'.format(IP))
        plt.close(f)


def main() :
    plot_mom()
    #plot_mom(scat_eTruth_flag=False,nTruth_flag=False)
    plot_kin()
    #plot_kin(k3_flag=False)
    plot_det()
    #plot_det(rp_flag=False,offmom_flag=False,b0_flag=False)
    plt.show()

if __name__=='__main__': main()    

