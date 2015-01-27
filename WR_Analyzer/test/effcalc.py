#!/usr/bin/env python

import os,sys
from itertools import product
from math import sqrt
import JChaves.Tools.Samples as Samples
sys.argv.append('-b')
import ROOT

signals = ['800','1000','2000','3000','4000','5000']
anas = ['ana','ana2','ana3']
xsignals = []
signal_events = 50000
ttbar_events = Samples.ttbar_samples[0].size
dyjets_events = (Samples.dyjets_samples[0].size+Samples.dyjets_samples[1].size+Samples.dyjets_samples[2].size+Samples.dyjets_samples[3].size) # FIX this number
luminosity = 10000.0 # 10fb-1 in 2015?

for s,t in product(signals,Samples.signal_samples):
    if '_'+s in t.name:
        xsignals.append([s,t.xsection*luminosity/signal_events])

fs = []
fs.append([ROOT.TFile('out_ttbar.root'),Samples.ttbar_samples[0].xsection*luminosity/ttbar_events,'ttbar'])
fs.append([ROOT.TFile('out_dyjets.root'),Samples.dyjets_samples[0].xsection*luminosity/dyjets_events,'dyjets'])
for s in xsignals:
    fs.append([ROOT.TFile('out_'+s[0]+'.root'),s[1],'signal_'+s[0]])    

sb1 = 0
sb2 = 0
sb3 = 0
for a in anas:
    print a
    ds = []
    for f in fs:
        ds.append([f[0].GetDirectory(a),f[1],f[2]])

    h1 = ROOT.TH1F("h","",100,0,1)
    h1 = ds[0][0].Get('N_pv')
    h2 = ROOT.TH1F("h","",100,0,1)
    h2 = ds[1][0].Get('N_pv')
    b = h1.Integral()*ds[0][1] + h2.Integral()*ds[1][1]
    for d in ds[2:]:
        h = ROOT.TH1F("h","",100,0,1)
        h = d[0].Get('N_pv')
        s = h.Integral()*d[1]
        print d[2],s/sqrt(b)
        
