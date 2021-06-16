#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2021-06-05 23:40:57 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#

import uproot as up
import sys

kinematics = sys.argv[1]
rootName="./OUTPUTS/%s.root" % kinematics # Note: this is relative to the bash script NOT this python script!

tName = "Process"
print("Reading tree {} in file {}".format(tName, rootName))
tdata = up.open(rootName)[tName]

bnames = tdata.keys()
print ("Tree has the following branches:")
print ("  [{}]".format(', '.join(bnames)))

tdict = {}
print("Tree has the following data in it's branches:")
for i,val in enumerate(bnames):
    print("{0} = {1}".format(val,tdata[val].array()))
    tdict[val] = tdata[val].array()

#print("Tree dictionary of all branches:")
#print("{}".format(tdict)

def findKey(key='invts'):
    """
    findKey(key='invts')
    
    """
    if key=='invts':
        # [key for key in tdict][0] finds the first key in the dictionary
        return list(tdict[[k for k in tdict][0]])
    else:
        # [k for k in tdict if k==key] finds the key matching the input
        return list(tdict[[k for k in tdict if k==key][0]])