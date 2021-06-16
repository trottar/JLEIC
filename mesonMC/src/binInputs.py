#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2021-06-14 20:03:46 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

import sys

bininputs="inputs/bins.input"
f    = open(bininputs)
fout = open('tmp2','wb')

def main():
    for line in f:
        if "#" in line:
            continue
        else:
            data = line.split('=')
            fout.write(bytes(data[1], 'UTF-8'))
    f.close()
    fout.close()
    return data
    
if __name__=='__main__': main()
