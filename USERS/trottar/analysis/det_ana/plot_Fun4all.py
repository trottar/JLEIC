#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2021-11-08 09:38:05 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.ticker import MaxNLocator
#import mplhep as hep
#hep.style.use(hep.style.ROOT)
import uproot as up
import numpy as np
import sys

numEvts = sys.argv[1]
IP = sys.argv[2]

#e_list = ["5on41","5on100","10on100","10on135","18on275"]
e_list = ["5on41","5on100","10on100","18on275"]

rootName = [None]*len(e_list)
tdata = [None]*len(e_list)
print("Energy settings...")
for i,e in enumerate(e_list):
    rootName[i]="./INPUTS/%s_%s.root" % (e,numEvts) # Note: this is relative to the bash script NOT this python script!
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

        for i,e in enumerate(e_list):
            f = plt.figure(figsize=(11.69,8.27))
            f.suptitle("{}".format(IP))
            ax = f.add_subplot(221)
            epxpy[i].plot(cbar=False,cmap=plt.cm.BuPu,ax=ax)
            ax.set_ylim(-10,10)
            ax.set_xlim(-10,10)
            ax.set_yticks(range(-10,10+1,10))
            ax.set_xticks(range(-10,10+1,10))
            ax.xaxis.set_major_formatter(plt.NullFormatter())
            ax.text(0.25, 0.85, r"e'$\frac{\Delta p_{x}}{Truth p_{x}}$ vs $\frac{\Delta p_{y}}{Truth p_{y}}$", transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
            plt.title(r"%s" % (e))

            ax = f.add_subplot(222)
            epxpz[i].plot(cmap=plt.cm.BuPu,ax=ax)
            ax.set_ylim(-10,10)
            ax.set_xlim(-10,10)
            ax.set_yticks(range(-10,10+1,10))
            ax.set_xticks(range(-10,10+1,10))
            ax.yaxis.set_major_formatter(plt.NullFormatter())
            ax.text(0.25, 0.85, r"e'$\frac{\Delta p_{x}}{Truth p_{x}}$ vs $\frac{\Delta p_{z}}{Truth p_{z}}$", transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')

            ax = f.add_subplot(223)
            epypz[i].plot(cbar=False,cmap=plt.cm.BuPu,ax=ax)
            ax.set_ylim(-10,10)
            ax.set_xlim(-10,10)
            ax.set_yticks(range(-10,10+1,10))
            ax.set_xticks(range(-10,10+1,10))
            ax.text(0.25, 0.85, r"e'$\frac{\Delta p_{y}}{Truth p_{y}}$ vs $\frac{\Delta p_{z}}{Truth p_{z}}$", transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.subplots_adjust(hspace=0.0,wspace=0.0)
        plt.savefig('OUTPUTS/MOM/scat_eTruth{}.png'.format(IP))

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

        for i,e in enumerate(e_list):
            f = plt.figure(figsize=(11.69,8.27))
            f.suptitle("{}".format(IP))
            ax = f.add_subplot(221)
            npxpy[i].plot(cbar=False,cmap=plt.cm.BuPu,ax=ax)
            ax.set_ylim(-10,10)
            ax.set_xlim(-10,10)
            ax.set_yticks(range(-10,10+1,10))
            ax.set_xticks(range(-10,10+1,10))
            ax.xaxis.set_major_formatter(plt.NullFormatter())
            ax.text(0.25, 0.85, r"n'$\frac{\Delta p_{x}}{Truth p_{x}}$ vs $\frac{\Delta p_{y}}{Truth p_{y}}$", transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')
            plt.title(r"%s" % (e))

            ax = f.add_subplot(222)
            npxpz[i].plot(cmap=plt.cm.BuPu,ax=ax)
            ax.set_ylim(-10,10)
            ax.set_xlim(-10,10)
            ax.set_yticks(range(-10,10+1,10))
            ax.set_xticks(range(-10,10+1,10))
            ax.yaxis.set_major_formatter(plt.NullFormatter())
            ax.text(0.25, 0.85, r"n'$\frac{\Delta p_{x}}{Truth p_{x}}$ vs $\frac{\Delta p_{z}}{Truth p_{z}}$", transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')

            ax = f.add_subplot(223)
            npypz[i].plot(cbar=False,cmap=plt.cm.BuPu,ax=ax)
            ax.set_ylim(-10,10)
            ax.set_xlim(-10,10)
            ax.set_yticks(range(-10,10+1,10))
            ax.set_xticks(range(-10,10+1,10))
            ax.text(0.25, 0.85, r"n'$\frac{\Delta p_{y}}{Truth p_{y}}$ vs $\frac{\Delta p_{z}}{Truth p_{z}}$", transform=ax.transAxes, fontsize=15, verticalalignment='top', horizontalalignment='left')

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.subplots_adjust(hspace=0.0,wspace=0.0)
        plt.savefig('OUTPUTS/MOM/nTruth{}.png'.format(IP))

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
        f.suptitle("{}".format(IP))
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("22%i" % (i+1)))
            n_px[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('n_px')
            plt.ylabel('counts')
            plt.title("IP6: n_px ({0})".format(e),fontsize=10)

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/MOM/n_px{}.png'.format(IP))

        f = plt.figure(figsize=(11.69,8.27))
        f.suptitle("{}".format(IP))
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("22%i" % (i+1)))
            n_py[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('n_py')
            plt.ylabel('counts')
            plt.title("IP6: n_py ({0})".format(e),fontsize=10)

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/MOM/n_py{}.png'.format(IP))        

        f = plt.figure(figsize=(11.69,8.27))
        f.suptitle("{}".format(IP))
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("22%i" % (i+1)))
            n_pz[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('n_pz')
            plt.ylabel('counts')
            plt.title("IP6: n_pz ({0})".format(e),fontsize=10)

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/MOM/n_pz{}.png'.format(IP))        

        f = plt.figure(figsize=(11.69,8.27))
        f.suptitle("{}".format(IP))
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("22%i" % (i+1)))
            n_p[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('n_p')
            plt.ylabel('counts')
            plt.title("IP6: n_p ({0})".format(e),fontsize=10)

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/MOM/n_p{}.png'.format(IP))

        f = plt.figure(figsize=(11.69,8.27))
        f.suptitle("{}".format(IP))
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("22%i" % (i+1)))
            n_E[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('n_E')
            plt.ylabel('counts')
            plt.title("IP6: n_E ({0})".format(e),fontsize=10)

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/MOM/n_E{}.png'.format(IP))

        f = plt.figure(figsize=(11.69,8.27))
        f.suptitle("{}".format(IP))
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("22%i" % (i+1)))
            n_Theta[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('$\Theta$')
            plt.ylabel('counts')
            plt.title("IP6: n_Theta ({0})".format(e),fontsize=10)

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/MOM/n_Theta{}.png'.format(IP))

        f = plt.figure(figsize=(11.69,8.27))
        f.suptitle("{}".format(IP))
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("22%i" % (i+1)))
            n_Phi[i].plot(ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('$\Phi$')
            plt.ylabel('counts')
            plt.title("IP6: n_Phi ({0})".format(e),fontsize=10)

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/MOM/n_Phi{}.png'.format(IP))

        f = plt.figure(figsize=(11.69,8.27))
        f.suptitle("{}".format(IP))
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("22%i" % (i+1)))
            nTrack_ThetaPhi[i].plot(cmin=5,ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('$\Theta$')
            plt.ylabel('$\Phi$')
            plt.title("IP6: nTrack_ThetaPhi ({0})".format(e),fontsize=10)

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/MOM/nTrack_ThetaPhi{}.png'.format(IP))

        f = plt.figure(figsize=(11.69,8.27))
        f.suptitle("{}".format(IP))
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("22%i" % (i+1)))
            nTrack_pTheta[i].plot(cmin=5,ax=ax)
            #ax.set_ylim(-50,50)
            #ax.set_xlim(-50,50)
            #ax.set_yticks(range(-50,50+1,10))
            #ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('$\Theta$')
            plt.ylabel('P')
            plt.title("IP6: nTrack_pTheta ({0})".format(e),fontsize=10)        

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/MOM/nTrack_pTheta{}.png'.format(IP))
    
        


        
        
        
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
            ZDC_XY[i] = tdata[i]['ZDC_XY'].to_hist()

        f = plt.figure(figsize=(11.69,8.27))
        f.suptitle("{}".format(IP))
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("22%i" % (i+1)))
            # cmin is the minimum value to include in histogram (cmax can be used for a max limit), here I am excluding 0 values
            print(tdata[i]['ZDC_XY'])
            ZDC_XY[i].plot(cmin=5,ax=ax)
            ax.set_ylim(-50,50)
            ax.set_xlim(-50,50)
            ax.set_yticks(range(-50,50+1,10))
            ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('X')
            plt.ylabel('Y')
            zdc_counts = ZDC_XY[i].counts()
            zdc_counts = sum(zdc_counts[~np.isnan(zdc_counts)])
            zdc_accp = 100*(zdc_counts/float(numEvts))
            if zdc_accp > 100:
                zdc_accp = 100
            plt.title("IP6: ZDC XY ({0}, {1:3.3f}% acceptance)".format(e,zdc_accp),fontsize=10)

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/DET/zdc_all_{}.png'.format(IP))

    '''
    RP
    '''

    if rp_flag==True:
        RP_XY = [None]*len(e_list)
        for i,e in enumerate(e_list):
            RP_XY[i] = tdata[i]['RP_XY'].to_hist()

        f = plt.figure(figsize=(11.69,8.27))
        f.suptitle("{}".format(IP))
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("22%i" % (i+1)))
            RP_XY[i].plot(cmin=5,ax=ax)
            ax.set_ylim(-50,50)
            ax.set_xlim(-50,50)
            ax.set_yticks(range(-50,50+1,10))
            ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('X')
            plt.ylabel('Y')
            rp_counts = RP_XY[i].counts()
            rp_counts = sum(rp_counts[~np.isnan(rp_counts)])
            rp_accp = 100*(rp_counts/float(numEvts))
            plt.title("IP6: RP XY ({0}, {1:3.3f}% acceptance)".format(e,rp_accp),fontsize=10)

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/DET/rp_all_{}.png'.format(IP))
        
    '''
    OFFMOM
    '''

    if offmom_flag==True:
        OFFMOM_XY = [None]*len(e_list)
        for i,e in enumerate(e_list):
            OFFMOM_XY[i] = tdata[i]['OFFMOM_XY'].to_hist()

        f = plt.figure(figsize=(11.69,8.27))
        f.suptitle("{}".format(IP))
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("22%i" % (i+1)))
            OFFMOM_XY[i].plot(cmin=5,ax=ax)
            ax.set_ylim(-50,50)
            ax.set_xlim(-50,50)
            ax.set_yticks(range(-50,50+1,10))
            ax.set_xticks(range(-50,50+1,10))
            plt.xlabel('X')
            plt.ylabel('Y')
            offmom_counts = OFFMOM_XY[i].counts()
            offmom_counts = sum(offmom_counts[~np.isnan(offmom_counts)])
            offmom_accp = 100*(offmom_counts/float(numEvts))
            plt.title("IP6: OFFMOM XY ({0}, {1:3.3f}% acceptance)".format(e,offmom_accp),fontsize=10)

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/DET/offmom_all_{}.png'.format(IP))        

    '''
    B0
    '''

    if b0_flag==True:
        B0_XY = [None]*len(e_list)
        for i,e in enumerate(e_list):
            B0_XY[i] = tdata[i]['B0_XY'].to_hist()

        f = plt.figure(figsize=(11.69,8.27))
        f.suptitle("{}".format(IP))
        for i,e in enumerate(e_list):
            ax = f.add_subplot(int("22%i" % (i+1)))
            B0_XY[i].plot(cmin=5,ax=ax)
            ax.set_ylim(-25,25)
            ax.set_xlim(-50,0)
            ax.set_yticks(range(-25,25+1,10))
            ax.set_xticks(range(-50,0+1,10))
            plt.xlabel('X')
            plt.ylabel('Y')            
            b0_counts = B0_XY[i].counts()
            b0_counts = sum(b0_counts[~np.isnan(b0_counts)])
            b0_accp = 100*(b0_counts/float(numEvts))
            plt.title("IP6: B0 XY ({0}, {1:3.3f}% acceptance)".format(e,b0_accp),fontsize=10)

        plt.tight_layout(rect=[0.0,0.03,1,0.95])
        plt.savefig('OUTPUTS/DET/b0_all_{}.png'.format(IP))


def main() :
    #plot_mom()
    #plot_mom(scat_eTruth_flag=False,nTruth_flag=False,n_flag=True)
    plot_det()
    #plot_det(zdc_flag=True,rp_flag=False,offmom_flag=False,b0_flag=False)
    plt.show()

if __name__=='__main__': main()    
