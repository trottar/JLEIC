#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2021-02-11 16:47:03 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

import uproot as up
import scipy as sci
import scipy.optimize as opt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import ticker
from collections import namedtuple
from scipy.interpolate import griddata
from sys import path
import time,math,sys,itertools

sys.path.insert(0,'/home/trottar/bin/python/')
import root2py as r2p

kinematics = sys.argv[1]

# kinematics="pi_n_18on275"
# kinematics="pi_n_10on100"
# kinematics="pi_n_5on100"
# kinematics="pi_n_5on41"
# kinematics="k_lambda_5on100"
# kinematics="k_lambda_18on275"

rootName="/home/trottar/ResearchNP/JLEIC/USERS/trottar/OUTPUTS/%s.root" % kinematics

# Patrick's interpolate grid
points,values=np.load('points.npy'),np.load('values.npy')

tree = up.open(rootName)["Evnts"]
branch = r2p.pyBranch(tree)

# Define phyisics data
s_e = branch.findBranch("invts","s_e")
s_q = branch.findBranch("invts","s_q")
Q2 = branch.findBranch("invts","Q2")
xBj = branch.findBranch("invts","xBj")
t = -branch.findBranch("invts","tSpectator")
tPrime = -branch.findBranch("invts","tPrime")
y_D = branch.findBranch("invts","y_D")
nu = branch.findBranch("invts","nu")
TwoPdotk = branch.findBranch("invts","TwoPdotk")
# xBj = tree.array("TDIS_xbj")
TDIS_xbj = tree.array("TDIS_xbj")
sigma_dis = tree.array("sigma_dis")*(1e-5) # used to compare to HERA cross-section
TDIS_y = tree.array("TDIS_y")
ppix_Lab = tree.array("ppix_Lab") # Long pion momentum
ppiy_Lab = tree.array("ppiy_Lab")
ppiz_Lab = tree.array("ppiz_Lab")
EpiE_Lab = tree.array("EpiE_Lab")
pprx_inc = tree.array("pprx_inc") # Long proton momentum
ppry_inc = tree.array("ppry_inc")
pprz_inc = tree.array("pprz_inc")
EprE_inc = tree.array("EprE_inc")
EnE_Lab = tree.array("EnE_Lab")
# pprx_inc = tree.array("pnx_Lab") # Long proton momentum
# ppry_inc = tree.array("pny_Lab")
# pprz_inc = tree.array("pnz_Lab")
# EprE_inc = tree.array("EnE_Lab")
y = tree.array("TDIS_y")
fpi = tree.array("fpi")
f2N = tree.array("f2N")
xpi = tree.array("xpi")
# xpi = TDIS_xbj
ypi = tree.array("ypi")
tpi = tree.array("tpi")
escat = tree.array("EScatRest")
pprz_inc = tree.array("pprz_inc")
ppiz_Lab = tree.array("ppiz_Lab")

binx  = 0.001
# binx  = 0.01
# binx  = 0.1
binQ2 = 10.0
# binQ2 = 3.0
# binQ2 = 1.0

print("\nx binning:",binx)
print("Q^2 binning:",binQ2,"\n")

# Create cut dictionary
cutDict = {}

# More coarse, set bin size
xpiarray = np.arange(binx/2,1.0,binx).tolist()

for i,x in enumerate(xpiarray):
    xpitmp = '{"xpicut%i" : ((%0.5f <= xpi) & (xpi <= %0.5f))}' % (i,xpiarray[i]-binx/2,xpiarray[i]+binx/2) 
    # (i,xpiarray[i]-binx/2,xpiarray[i]+binx/2) # no binning
    # (i,xpiarray[i]-binx/20,xpiarray[i]+binx/20) # for proper binning
    print('{"xpicut%i" : ((%0.5f <= xpi) & (xpi <= %0.5f))}' % (i,xpiarray[i]-binx/2,xpiarray[i]+binx/2))
    cutDict.update(eval(xpitmp))
c = r2p.pyPlot(cutDict)

xpicut = []
for i,evt in enumerate(xpiarray):
    xpicut.append("xpicut%s" % i)
    tmp  = "xpicut_%s = [\"xpicut%s\"]" % (i,i)
    exec(tmp)
    
s_e = c.applyCuts(s_e,xpicut,either=True)
s_q = c.applyCuts(s_q,xpicut,either=True)
Q2 = c.applyCuts(Q2,xpicut,either=True)
xBj = c.applyCuts(xBj,xpicut,either=True)
t = c.applyCuts(t,xpicut,either=True)
tPrime = c.applyCuts(tPrime,xpicut,either=True)
y_D = c.applyCuts(y_D,xpicut,either=True)
nu = c.applyCuts(nu,xpicut,either=True)
TwoPdotk = c.applyCuts(TwoPdotk,xpicut,either=True)
TDIS_xbj = c.applyCuts(TDIS_xbj,xpicut,either=True)
sigma_dis = c.applyCuts(sigma_dis,xpicut,either=True)
TDIS_y = c.applyCuts(TDIS_y,xpicut,either=True)
ppix_Lab = c.applyCuts(ppix_Lab,xpicut,either=True)
ppiy_Lab = c.applyCuts(ppiy_Lab,xpicut,either=True)
ppiz_Lab = c.applyCuts(ppiz_Lab,xpicut,either=True)
EpiE_Lab = c.applyCuts(EpiE_Lab,xpicut,either=True)
pprx_inc = c.applyCuts(pprx_inc,xpicut,either=True)
ppry_inc = c.applyCuts(ppry_inc,xpicut,either=True)
pprz_inc = c.applyCuts(pprz_inc,xpicut,either=True)
EprE_inc = c.applyCuts(EprE_inc,xpicut,either=True)
EnE_Lab = c.applyCuts(EnE_Lab,xpicut,either=True)
y = c.applyCuts(y,xpicut,either=True)
fpi = c.applyCuts(fpi,xpicut,either=True)
f2N = c.applyCuts(f2N,xpicut,either=True)
ypi = c.applyCuts(ypi,xpicut,either=True)
tpi = c.applyCuts(tpi,xpicut,either=True)
escat = c.applyCuts(escat,xpicut,either=True)
# pprz_inc = c.applyCuts(pprz_inc,xpicut,either=True)
# ppiz_Lab = c.applyCuts(ppiz_Lab,xpicut,either=True)

xpi = c.applyCuts(xpi,xpicut,either=True)

Q2array = np.arange(0.0,500.0,binQ2).tolist()
# Q2array = np.arange(1.0,3.0,1.0).tolist() 
for i,x in enumerate(Q2array) :
    Q2tmp = '{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,Q2array[i]-binQ2/2,Q2array[i]+binQ2/2)
    # (i,Q2array[i]-binQ2/20,Q2array[i]+binQ2/20) # for proper binning
    print('{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,Q2array[i]-binQ2/2,Q2array[i]+binQ2/2))
    cutDict.update(eval(Q2tmp))
c = r2p.pyPlot(cutDict)

Q2cut = []
for i,evt in enumerate(Q2array):
    Q2cut.append("Q2cut%s" % i)

s_e = c.applyCuts(s_e,Q2cut,either=True)
s_q = c.applyCuts(s_q,Q2cut,either=True)
xBj = c.applyCuts(xBj,Q2cut,either=True)
t = c.applyCuts(t,Q2cut,either=True)
tPrime = c.applyCuts(tPrime,Q2cut,either=True)
y_D = c.applyCuts(y_D,Q2cut,either=True)
nu = c.applyCuts(nu,Q2cut,either=True)
TwoPdotk = c.applyCuts(TwoPdotk,Q2cut,either=True)
TDIS_xbj = c.applyCuts(TDIS_xbj,Q2cut,either=True)
sigma_dis = c.applyCuts(sigma_dis,Q2cut,either=True)
TDIS_y = c.applyCuts(TDIS_y,Q2cut,either=True)
ppix_Lab = c.applyCuts(ppix_Lab,Q2cut,either=True)
ppiy_Lab = c.applyCuts(ppiy_Lab,Q2cut,either=True)
ppiz_Lab = c.applyCuts(ppiz_Lab,Q2cut,either=True)
EpiE_Lab = c.applyCuts(EpiE_Lab,Q2cut,either=True)
pprx_inc = c.applyCuts(pprx_inc,Q2cut,either=True)
ppry_inc = c.applyCuts(ppry_inc,Q2cut,either=True)
pprz_inc = c.applyCuts(pprz_inc,Q2cut,either=True)
EprE_inc = c.applyCuts(EprE_inc,Q2cut,either=True)
EnE_Lab = c.applyCuts(EnE_Lab,Q2cut,either=True)
y = c.applyCuts(y,Q2cut,either=True)
fpi = c.applyCuts(fpi,Q2cut,either=True)
f2N = c.applyCuts(f2N,Q2cut,either=True)
xpi = c.applyCuts(xpi,Q2cut,either=True)
ypi = c.applyCuts(ypi,Q2cut,either=True)
tpi = c.applyCuts(tpi,Q2cut,either=True)
escat = c.applyCuts(escat,Q2cut,either=True)
# pprz_inc = c.applyCuts(pprz_inc,Q2cut,either=True)
# ppiz_Lab = c.applyCuts(ppiz_Lab,Q2cut,either=True)

Q2 = c.applyCuts(Q2,Q2cut,either=True)

sigma_tdis = sigma_dis*(fpi/f2N)
# xL = EnE_Lab/EprE_inc # Frac of proton momentum
xL = 1-(xBj/xpi) # Frac of proton momentum

Q2array1 = np.arange(0.0,10,1.0).tolist()
Q2array10 = np.arange(10,100.0,10.0).tolist()
Q2array100 = np.arange(100,600.0,100.0).tolist()
# Q2binarray = Q2array1 + Q2array10 + Q2array100 # xbin 0.1
# Q2binarray = Q2array1 + Q2array10 # xbin 0.001
# Q2binarray = np.arange(1.0,3.0,1.0).tolist() #

Q2binarray = np.arange(40.0,80.0,10.0).tolist() 
for i,x in enumerate(Q2binarray) :
    Q2tmp = '{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,Q2binarray[i]-binQ2/2,Q2binarray[i]+binQ2/2)
    print('{"Q2cut%i" : ((%0.1f <= Q2) & (Q2 <= %0.1f))}' % (i,Q2binarray[i]-binQ2/2,Q2binarray[i]+binQ2/2))
    cutDict.update(eval(Q2tmp))

xarray = np.arange(0.05,1.0,0.1).tolist() # 0.1 binning x=0.1-1
# xarray = np.arange(0.001,0.1,0.001).tolist() # 0.001 binning x=0.001-0.1
# xarray = np.arange(0.05,0.1,0.1).tolist() 
for i,x in enumerate(xarray):
    xtmp = '{"xcut%i" : ((%0.4f <= xpi) & (xpi <= %0.4f))}' % (i,xarray[i]-0.005,xarray[i]+0.005)
    print('{"xcut%i" : ((%0.4f <= xpi) & (xpi <= %0.4f))}' % (i,xarray[i]-0.005,xarray[i]+0.005))
    cutDict.update(eval(xtmp))

xcut = []
for i,evt in enumerate(xarray):
    xcut.append("xcut%s" % i)
    tmp  = "xcut_%s = [\"xcut%s\"]" % (i,i)
    exec(tmp)

    
ytmp = '{"ycut" : ((0.01 <= y) & (y <= 0.95))}'
ttmp = '{"tcut" : ((-1.00 <= t) & (t <= 0.00))}'
cutDict.update(eval(ttmp))
cutDict.update(eval(ytmp))
c = r2p.pyPlot(cutDict)

ycut1 = ["ycut"]
tcut1 = ["tcut"]

for i,evt in enumerate(Q2binarray):
    for j,nevt in enumerate(xarray):
        tmp = "cutx%s_q%s = [\"Q2cut%s\",\"xcut%s\"]" % (j,int(evt),i,j)
        exec(tmp)

def F2pi(xpi, Q2):
    points,values=np.load('xpiQ2.npy'),np.load('F2pi.npy')
    F2pi=lambda xpi,Q2: griddata(points,values,(np.log10(xpi),np.log10(Q2)))
    return F2pi(xpi,Q2)
        
# Calculate cross-section using Patrick's interpolate grid
def ds_dxdQ2dxLdt(x, Q2, xL,t):
    
    sigma=lambda x,Q2,xL,t: griddata(points,values,(np.log(x),np.log(Q2),xL,t))/Q2/x**2.5*0.389379372e12

    return sigma(x,Q2,xL,t)

def Lumi(evts, x, Q2, xL, t):

    # Luminosity
    tot_sigma = ds_dxdQ2dxLdt(x, Q2, xL, t)
    
    # all binned events
    sig_all = tot_sigma 
    evts_all = evts

    bint = 0.1
    binxL = 0.1

    lumi = evts_all/(sig_all*binQ2*binx)

    nevt = 100/(lumi)

    return nevt
    # return tot_sigma

def generateTable():

    # print(type(ycut1))        
    # print(ycut1)
    # print(c.applyCuts(y,ycut1))
    # print(type(cutx0_q1))        
    # print(cutx0_q1)
    # # print(c.applyCuts(y,["Q2cut0","xcut0"]))
    # print(c.applyCuts(y,["xcut0"]))
    # print(c.applyCuts(y,cutx0_q1))
    
    Q2Dict = {}
    xBjDict = {}
    xpiDict = {}
    xLDict = {}
    tDict = {}
    yDict = {}
    yieldDict = {}
    theoryDict = {}
    
    for i,evt in enumerate(Q2binarray):
        for j,nevt in enumerate(xarray):
            try:
                # Q2
                Q2_tmp = '{"Q2avg_cutx%s_q%s" : (c.applyCuts(Q2,cutx%s_q%s)[0])}' % (j,int(evt),j,int(evt))
                # print('{"Q2avg_cutx%s_q%s" : (c.applyCuts(Q2,cutx%s_q%s)[0])}' % (j,int(evt),j,int(evt)))
                Q2Dict.update(eval(Q2_tmp))        
            except IndexError as e:
                print(e)
                continue
                # print(sys.exc_type)
            # xBj
            xBj_tmp = '{"xBjavg_cutx%s_q%s" : (c.applyCuts(xBj,cutx%s_q%s)[0])}' % (j,int(evt),j,int(evt))
            # print('{"xBjavg_cutx%s_q%s" : (c.applyCuts(xBj,cutx%s_q%s)[0])}' % (j,int(evt),j,int(evt)))
            xBjDict.update(eval(xBj_tmp))
            # xpi
            xpi_tmp = '{"xpiavg_cutx%s_q%s" : (c.applyCuts(xpi,cutx%s_q%s)[0])}' % (j,int(evt),j,int(evt))
            # print('{"xpiavg_cutx%s_q%s" : (c.applyCuts(xpi,cutx%s_q%s)[0])}' % (j,int(evt),j,int(evt)))
            xpiDict.update(eval(xpi_tmp))
            # xL
            xL_tmp = '{"xLavg_cutx%s_q%s" : (c.applyCuts(xL,cutx%s_q%s)[0])}' % (j,int(evt),j,int(evt))
            # print('{"xLavg_cutx%s_q%s" : (c.applyCuts(xL,cutx%s_q%s)[0])}' % (j,int(evt),j,int(evt)))
            xLDict.update(eval(xL_tmp))
            if xLDict["xLavg_cutx%s_q%s" % (j,int(evt))] < 0.8:
                print("xL = %s" % xLDict["xLavg_cutx%s_q%s" % (j,int(evt))])
                continue            
            # t
            t_tmp = '{"tavg_cutx%s_q%s" : (c.applyCuts(t,cutx%s_q%s)[0])}' % (j,int(evt),j,int(evt))
            # print('{"tavg_cutx%s_q%s" : (c.applyCuts(t,cutx%s_q%s)[0])}' % (j,int(evt),j,int(evt)))
            tDict.update(eval(t_tmp))
            # y
            y_tmp = '{"yavg_cutx%s_q%s" : (c.applyCuts(y,cutx%s_q%s)[0])}' % (j,int(evt),j,int(evt))
            # print('{"yavg_cutx%s_q%s" : (c.applyCuts(y,cutx%s_q%s)[0])}' % (j,int(evt),j,int(evt)))
            yDict.update(eval(y_tmp))            
            binevts = eval("c.applyCuts(t,cutx%s_q%s)" % (j,int(evt)))
            # # check lumi
            lumi = Lumi(len(binevts), xBjDict["xBjavg_cutx%s_q%s" % (j,int(evt))], Q2Dict["Q2avg_cutx%s_q%s" % (j,int(evt))], xLDict["xLavg_cutx%s_q%s" % (j,int(evt))], tDict["tavg_cutx%s_q%s" % (j,int(evt))])              
            # if math.isnan(lumi):
            #     print("\nNan value found")
            #     continue        
            # yield
            # yield_tmp = '{"yieldavg_cutx%s_q%s" : len(binevts)}' % (j,int(evt))
            yield_tmp = '{"yieldavg_cutx%s_q%s" : lumi}' % (j,int(evt))
            yieldDict.update(eval(yield_tmp))
            # # print('{"yieldavg_cutx%s_q%s" : (c.applyCuts(lumi,cutx%s_q%s)[0])}' % (j,int(evt),j,int(evt)))
            # print("\n¬¬¬¬¬¬¬¬¬¬¬¬ end of loop: ", j)                
            # print("\n",xBjDict["xBjavg_cutx%s_q%s" % (j,int(evt))],Q2Dict["Q2avg_cutx%s_q%s" % (j,int(evt))],xLDict["xLavg_cutx%s_q%s" % (j,int(evt))]) 
            
            theoryDict["cutx%s_q%s" % (j,int(evt))] ={
                "Q2"    : Q2Dict["Q2avg_cutx%s_q%s" % (j,int(evt))],
                "x"     : xBjDict["xBjavg_cutx%s_q%s" % (j,int(evt))],
                "xpi"   : xpiDict["xpiavg_cutx%s_q%s" % (j,int(evt))],
                "xL"    : xLDict["xLavg_cutx%s_q%s" % (j,int(evt))],
                "t"     : tDict["tavg_cutx%s_q%s" % (j,int(evt))],
                "y"     : yDict["yavg_cutx%s_q%s" % (j,int(evt))],
                "events per bin" : yieldDict["yieldavg_cutx%s_q%s" % (j,int(evt))]
            }
            
        c.progressBar(i,len(Q2binarray)-1,70)            

    print(theoryDict)
        
    tmp = []
    frames = []

    # for key,val in theoryDict.items():
    #     tmp.append(key)
    #     # frames.append(pd.DataFrame.from_dict(val))
    #     frames.append(pd.DataFrame.from_dict(val,orient='index'))            
    # theory_df_table = pd.concat(frames,keys=tmp)
    # theory_df_table = theory_df_table.T
        
    theory_df_table = pd.DataFrame(data=theoryDict)
    theory_df_table = theory_df_table.T
    print(theory_df_table)
    return theory_df_table



def main() :

    df_table = generateTable()

    #render dataframe as html
    html_table = df_table.to_html()
    csv_table = df_table.to_csv()

    #write html to file
    text_file = open("theory_table.html", "w")
    text_file.write(html_table)
    text_file.close()

    text_file = open("theory_table.csv", "w")
    text_file.write(csv_table)
    text_file.close()
    
if __name__=='__main__': main()

