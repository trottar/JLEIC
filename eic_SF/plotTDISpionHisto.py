#!/usr/bin/env python

#import ROOT import *
from ROOT import TCanvas, TPad, TFile, TPaveLabel, TPaveText, TTreeReader, TTreeReaderValue
from ROOT import gROOT
from rootpy.interactive import wait
import numpy as np

outputpdf = 'OUTPUT/plot_TDISpionHisto.pdf'

cypi = TCanvas( 'cypi', 'hypi', 600, 400)
cse = TCanvas( 'cse', 'hse', 600, 400)
cQ2vsX = TCanvas( 'cQ2vsX', 'hQ2vsX', 600, 400)
ctS = TCanvas( 'ctS', 'htS', 600, 400)
ctSY = TCanvas( 'ctSY', 'htSY', 600, 400)
ctSZ = TCanvas( 'ctSZ', 'htSZ', 600, 400)
cnu = TCanvas( 'cnu', 'hnu', 600, 400)
calphaS = TCanvas( 'calphaS', 'halphaS', 600, 400)
cxpi = TCanvas( 'cxpi', 'hxpi', 600, 400)
cypi = TCanvas( 'cypi', 'hypi', 600, 400)
ctpi = TCanvas( 'ctpi', 'htpi', 600, 400)
cptp2 = TCanvas( 'cptp2', 'hptp2', 600, 400)
czp2 = TCanvas( 'czp2', 'hzp2', 600, 400)
csigmaDIS = TCanvas( 'csigmaDIS', 'hsigmaDIS', 600, 400)
csigmaTDIS = TCanvas( 'csigmaTDIS', 'hsigmaTDIS', 600, 400)

f = TFile.Open("TDISpion.root","read")
if not f or f.IsZombie() :
	exit
f.ls()

#tree = f.Get("Evnts")
#r = TTreeReader("Evnts", f)
#xpi = TTreeReaderValue(float)(r, "xpi")
#while r.Next():
#    print(xpi)

cypi.cd()
cypi.SetLogy()
hypi = gROOT.FindObject( 'hypi' )
hypi.SetFillColor( 45 )
hypi.DrawCopy()
label1 = TPaveLabel( -3.5, 700, -1, 800, 'Default option' )
label1.SetFillColor( 42 )
label1.Draw()
cypi.Print(outputpdf+"(")

cse.cd()
hse = gROOT.FindObject( 'hse' )
hse.SetFillColor( 45 )
hse.DrawCopy()
cse.Print(outputpdf)

cQ2vsX.cd()
hQ2vsX = gROOT.FindObject( 'hQ2vsX' )
hQ2vsX.SetFillColor( 45 )
hQ2vsX.DrawCopy()
cQ2vsX.Print(outputpdf)

ctS.cd()
htS = gROOT.FindObject( 'htS' )
htS.SetFillColor( 45 )
htS.DrawCopy()
ctS.Print(outputpdf)

ctSY.cd()
htSY = gROOT.FindObject( 'htSY' )
htSY.SetFillColor( 45 )
htSY.DrawCopy()
ctSY.Print(outputpdf)

ctSZ.cd()
htSZ = gROOT.FindObject( 'htSZ' )
htSZ.SetFillColor( 45 )
htSZ.DrawCopy()
ctSZ.Print(outputpdf)

cnu.cd()
hnu = gROOT.FindObject( 'hnu' )
hnu.SetFillColor( 45 )
hnu.DrawCopy()
cnu.Print(outputpdf)

calphaS.cd()
halphaS = gROOT.FindObject( 'halphaS' )
halphaS.SetFillColor( 45 )
halphaS.DrawCopy()
calphaS.Print(outputpdf)

cxpi.cd()
hxpi = gROOT.FindObject( 'hxpi' )
hxpi.SetFillColor( 45 )
hxpi.DrawCopy()
cxpi.Print(outputpdf)

cypi.cd()
hypi = gROOT.FindObject( 'hypi' )
hypi.SetFillColor( 45 )
hypi.DrawCopy()
cypi.Print(outputpdf)

ctpi.cd()
htpi = gROOT.FindObject( 'htpi' )
htpi.SetFillColor( 45 )
htpi.DrawCopy()
ctpi.Print(outputpdf)

cptp2.cd()
hptp2 = gROOT.FindObject( 'hptp2' )
hptp2.SetFillColor( 45 )
hptp2.DrawCopy()
cptp2.Print(outputpdf)

czp2.cd()
hzp2 = gROOT.FindObject( 'hzp2' )
hzp2.SetFillColor( 45 )
hzp2.DrawCopy()
czp2.Print(outputpdf)

csigmaDIS.cd()
hsigmaDIS = gROOT.FindObject( 'hsigmaDIS' )
hsigmaDIS.SetFillColor( 45 )
hsigmaDIS.DrawCopy()
csigmaDIS.Print(outputpdf)

csigmaTDIS.cd()
hsigmaTDIS = gROOT.FindObject( 'hsigmaTDIS' )
hsigmaTDIS.SetFillColor( 45 )
hsigmaTDIS.DrawCopy()
csigmaTDIS.Print(outputpdf+")")

#if not gROOT.IsBatch():
#    input("press [Enter] to continue ")

