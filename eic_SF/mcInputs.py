#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2020-02-26 13:18:18 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

import sys

mcinputs="TDISMC_constants.input"
f    = open(mcinputs)
fout = open('tmp','wb')

def main():
    for line in f:
        data = line.split('=')
        fout.write(data[1])
    f.close()
    return data
    
if __name__=='__main__': main()
