#!/usr/bin/env python

import time,sys,os,argparse,atexit,subprocess,math

from collections import namedtuple

import matplotlib.pyplot as plt

from array import array

def getValues() :

    binNum = array('i')          
    yieldVal = array('d')
    NEvts = array('i')

    filename = "Yield_Data.dat"
    
    f = open(filename)
    
    for line in f:
        data = line.split()
        #print(str(data))
        if "#" not in line :
            binNum.append(int(data[0]))    
            #print("Run %s" % (binNum))
            yieldVal.append(float(data[1]))
            NEvts.append(int(data[2]))
 
    f.close()

    return[binNum,yieldVal,NEvts]


def main() :

    [binNum,yieldVal,NEvts] = getValues()

    foutname = 'plot_yield.png'

    yieldPlot = plt.figure()

    #plt.subplot(1,2,1)    
    plt.grid(zorder=1)
    #plt.xlim(0,60)
    #plt.ylim(0.98,1.02)
    #plt.plot([0,60], [1,1], 'r-',zorder=2)
    #plt.errorbar(current,yieldValRel_HMS,yerr=uncerEvts_HMS,color='black',linestyle='None',zorder=3)
    plt.scatter(binNum,yieldVal,color='blue',zorder=4)
    plt.ylabel('Yield [events/sec]', fontsize=16)
    plt.title('Yield vs. Bin for %s events' % (NEvts[0]), fontsize =16)
    plt.xlabel('Bin Number', fontsize =16)
    
    plt.tight_layout()
    plt.show()
    yieldPlot.savefig(foutname)
    
if __name__=='__main__': main()
