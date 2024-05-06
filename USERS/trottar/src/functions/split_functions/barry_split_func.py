#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-05-06 11:05:15 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trottar.iii@gmail.com>
#
# Copyright (c) trottar
#
import os,sys
import threading
import copy
import numpy as np
import pandas as pd
import time
from mpmath import fp

from scipy.integrate import quad,fixed_quad
from scipy.special import jv, zeta
from scipy import interpolate
from scipy.optimize import leastsq,least_squares

from tools.tools import load, save, checkdir,lprint
from tools.config import conf, load_config, load_config2
from qcdlib.aux import AUX
from qcdlib.eweak import EWEAK
from qcdlib.mellin import MELLIN
from qcdlib.alphaS import ALPHAS
from fitlib.resman import RESMAN
from fitlib.parman import PARMAN
from fitlib.maxlike import MAXLIKE
from analysis.corelib import core

'''
I have provided some parameters for the regulators (written as $\Lambda$ in the note and JAM 2018 PRL). In my code, they are written as “par” in get_regulator. I list below the mean values of my fits to the smallest kT and largest xL bins.

 

IMF exp: 1.532

Cov exp: 0.774

Cov mon: 0.532

Regge: 0.977

Pauli-Villars: 0.551
'''

class PP2LAMBDA_SPLITTING:
    
    def __init__(self):
        self.aux = conf['aux']
        self.mK2  = self.aux.Mk**2
        self.mN2  = self.aux.M2
        self.ML   = 1.115683 #--Lambda mass particle from PDG
        
        self.mN   = self.aux.M
        self.mK   = self.aux.Mk
        self.Mbar = self.ML+self.mN
        self.DelM = self.ML-self.mN

        self.model = conf['model']

    def get_regulator(self,kT2,xL,par):
        
        mK2,mN2,ML=self.mK2,self.mN2,self.ML
        
        sKL = (kT2+mK2)/(1-xL)+(kT2+ML**2)/xL
        t = -kT2/xL - (1-xL)/xL*(ML**2-xL*mN2)
        DKY = -(kT2+(1-xL)*ML**2+xL*mK2-(1-xL)*xL*mN2)/xL

        if self.model=='IMF exp':
            reg = np.exp(2*(mN2-sKL)/par**2)
        elif self.model=='cov exp':
            reg = np.exp(2*(t-mK2)/par**2)
        elif self.model=='cov mon':
            reg = ((par**2-mK2)/(par**2-t))**2
        elif self.model=='cov dip':
            reg = ((par**2-mK2)/(par**2-t))**4
        elif self.model=='Regge': 
            reg = (1-xL)**(-2*t) * np.exp(2*(t-mK2)/par**2)
        elif self.model=='Pauli-Villars':
            reg = (1 - (t-mK2)**2/(t-par**2)**2)
        return reg


    def get_fL(self,kT2,xL):

        #--prefactors

        psdc = 0.093 #--pseudoscalar decay constant

        D = 0.85
        F = 0.41
        CKL2 = ((D+3*F)/2.0/3**0.5)**2 #--K^+ Lambda coefficient C_{KY}

        prefact=CKL2*self.Mbar**2/(4*np.pi*psdc)**2
        DKY = -(kT2+(1-xL)*self.ML**2+xL*self.mK2-(1-xL)*xL*self.mN2)/xL
        fKpL = prefact*(1-xL)*(kT2+(self.mN*(1-xL)+self.DelM)**2)/(xL**2*DKY**2)

        return fKpL
    
    def get_theory(self,par,xL,kT):
    
        sig_tot = 19.9 #--mb

        return self.get_fL(kT**2,xL)/np.pi*xL*sig_tot * self.get_regulator(kT**2,xL,par)
    
    
class PP2N_SPLITTING:
    
    def __init__(self):
        self.aux  = conf['aux']
        self.mPi2 = self.aux.Mpi2        
        self.mN   = self.aux.M
        self.mN2  = self.aux.M2

        self.model = conf['model']

    def get_regulator(self,kT2,xL,par):
        
        mpi2,mN2=self.mPi2,self.mN2
        
        spiN = (kT2+mpi2)/(1-xL)+(kT2+mN2)/xL
        t = -(kT2 + (1-xL)**2 * mN2)/xL
        DpiN = -(kT2 + (1-xL)**2*mN2 + xL*mpi2)/xL

        if self.model=='IMF exp':
            reg = np.exp(2*(mN2-spiN)/par**2)
        elif self.model=='cov exp':
            reg = np.exp(2*(t-mpi2)/par**2)
        elif self.model=='cov mon':
            reg = ((par**2-mpi2)/(par**2-t))**2
        elif self.model=='cov dip':
            reg = ((par**2-mpi2)/(par**2-t))**4
        elif self.model=='Regge': 
            reg = (1-xL)**(-2*t) * np.exp(2*(t-mpi2)/par**2)
        elif self.model=='Pauli-Villars':
            reg = (1 - (t-mpi2)**2/(t-par**2)**2)
        return reg


    def get_fN(self,kT2,xL):
                
        prefact=(13.7)/(4*np.pi) #--comes from gA**2 * M**2 / fpi**2 = 13.7 * 4*pi
        DpiN = -(kT2 + (1-xL)**2*self.mN2 + xL*self.mPi2)/xL
        fppin = prefact*(1-xL)*(kT2 + self.mN2*(1-xL)**2)/(xL**2*DpiN**2)

        return 2 * fppin #--factor 2 for the charged pion
    
    def get_theory(self,par,xL,kT):
    
        sig_tot = 23.8 #--mb

        return self.get_fN(kT**2,xL)/np.pi*xL*sig_tot * self.get_regulator(kT**2,xL,par)
    
