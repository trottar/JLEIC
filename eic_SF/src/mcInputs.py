#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2020-02-26 14:43:58 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

import sys

mcinputs="inputs/kinematics.input"
f    = open(mcinputs)
fout = open('tmp','wb')

def main():
    for line in f:
        data = line.split('=')
        fout.write(bytes(data[1], 'UTF-8'))
    f.close()
    return data
    
if __name__=='__main__': main()
