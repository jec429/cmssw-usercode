#!/usr/bin/env python

import os,sys
from itertools import product
from JChaves.Tools.ROOTTools import differentiate_stat_box
import JChaves.Tools.Samples as Samples
sys.argv.append('-b')
import ROOT

anas = ['ana','ana2','ana3','ana4']
logs = [True,False]
signals = ['800','1000','2000','3000','4000','5000']
signals = ['2000']
xsignals = []
signal_events = 50000
ttbar_events = Samples.ttbar_samples[0].size
dyjets_events = (Samples.dyjets_samples[0].size+Samples.dyjets_samples[1].size+Samples.dyjets_samples[2].size+Samples.dyjets_samples[3].size) # FIX this number
luminosity = 10000.0 # 10fb-1 in 2015?
hists=['M_ee','M_jj','M_eeqq','M_eejj','M_N1_ee','M_N2_ee','M_N1_mumu','M_N2_mumu','mMET','jet_pt','jet1_pt','jet2_pt','jet1_eta','jet2_eta','jet1_phi','jet2_phi','e1_pt','e2_pt','e1_eta','e2_eta','e1_phi','e2_phi','deta_ee','dphi_ee','deta_mumu','dphi_mumu','deta_jets','dphi_jets','N_jets','M_mumu','M_emu','M_emujj','M_mumujj','mu1_pt','mu2_pt','mu1_eta','mu2_eta','mu1_phi','mu2_phi']
edges=[[0,3000],[0,3000],[0,6000],[0,6000],[0,5000],[0,5000],[0,5000],[0,5000],[0,1000],[0,2000],[0,2000],[0,2000],[-3,3],[-3,3],[-3.15,3.15],[-3.15,3.15],[0,2000],[0,2000],[-3,3],[-3,3],[-3.15,3.15],[-3.15,3.15],[-3,3],[-3.15,3.15],[-3,3],[-3.15,3.15],[-3,3],[-3.15,3.15],[0,40],[0,3000],[0,3000],[0,6000],[0,6000],[0,2000],[0,2000],[-3,3],[-3,3],[-3.15,3.15],[-3.15,3.15]]


for s,t in product(signals,Samples.signal_samples):
    if '_'+s in t.name:
        xsignals.append([s,t.xsection*luminosity/signal_events])

fs = []
fs.append([ROOT.TFile('plots/rootfiles/54_10_02_24_15_out_EXO_ttbar.root'),Samples.ttbar_samples[0].xsection*luminosity/ttbar_events,'ttbar'])
fs.append([ROOT.TFile('plots/rootfiles/54_10_02_24_15_out_EXO_dyjets.root'),Samples.dyjets_samples[0].xsection*luminosity/dyjets_events,'dyjets'])
fs.append([ROOT.TFile('plots/rootfiles/29_08_02_24_15_out_EXO_EE.root'),7.07E-02*luminosity/signal_events,'ee_2000'])
fs.append([ROOT.TFile('plots/rootfiles/29_08_02_24_15_out_EXO_MUMU.root'),6.22E-02*luminosity/signal_events,'mumu_2000'])
#fs.append([ROOT.TFile('plots/rootfiles/43_08_02_23_15_genhistos_EE.root'),7.07E-02*luminosity/signal_events,'ee_2000'])
#fs.append([ROOT.TFile('plots/rootfiles/22_08_02_23_15_genhistos_MUMU.root'),6.22E-02*luminosity/signal_events,'mumu_2000'])
#for s in xsignals:
#    fs.append([ROOT.TFile('out_'+s[0]+'.root'),s[1],'signal_'+s[0]])    
    
for a,log in product(anas,logs):
    os.system('mkdir -p plots/eemumu/'+a)
    ds = []
    for f in fs:
        ds.append([f[0].GetDirectory(a),f[1],f[2]])        
    for h,e in zip(hists,edges):
        hs = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        ROOT.gStyle.SetOptStat(1111111)
        c1 = ROOT.TCanvas("c1","Canvas",600,600)
        pad1 = ROOT.TPad("pad1","This is pad1",0.02,0.52,0.48,0.98);
        pad2 = ROOT.TPad("pad2","This is pad2",0.52,0.52,0.98,0.98);
        pad3 = ROOT.TPad("pad3","This is pad3",0.02,0.02,0.48,0.48);
        pad4 = ROOT.TPad("pad4","This is pad4",0.52,0.02,0.98,0.48);
        th = ROOT.THStack("th","stacked histograms");
        for i,d in enumerate(ds):
            hs[i] = ROOT.TH1F("h"+str(i+1),"",100,e[0],e[1])
            hs[i] = d[0].Get(h)
            hs[i].SetLineColor(i+1)
            hs[i].Scale(d[1])
            hs[i].SetName(d[2])
            hs[i].SetTitle(h)
        #pad1.SetFillColor(ROOT.kWhite)
        #pad2.SetFillColor(ROOT.kWhite)
        #pad3.SetFillColor(ROOT.kWhite)
        #pad4.SetFillColor(ROOT.kWhite)
        pad1.Draw()
        pad2.Draw()
        pad3.Draw()
        pad4.Draw()
        pad1.cd()
        if log:
            pad1.SetLogy()
            pad2.SetLogy()
            h = h + '_log'
        
        histos = sorted(hs[:len(ds)], key=lambda x:x.GetMaximum())
        #print histos
        #print histos[:-1]
        histos[-1].Draw()
        histos[-1].SetMinimum(0.1)
        
        for histo in histos[:-1]:      
            histo.Draw("sames")
            c1.Update()
        th.Add(hs[0])
        th.Add(hs[1]) 
        pad2.cd()
        th.Draw()
        th.SetMinimum(0.1)
        for histo in hs[2:len(ds)]:
            histo.Draw("same")            

        pad3.cd()
        hs[0].DrawNormalized()
        for histo in hs[1:len(ds)]:
            histo.DrawNormalized("same")
        pad4.cd()

        for i,hi in enumerate(hs[:len(ds)]):
            differentiate_stat_box(hi,i)
                  
        c1.Print('plots/eemumu/'+a+'/histos_'+h+'_scaled.pdf')
        c1.Print('plots/eemumu/'+a+'/histos_'+h+'_scaled.png')
        c1=0
        th=0
