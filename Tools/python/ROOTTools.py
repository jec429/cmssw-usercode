#!/usr/bin/env python

import os,sys
sys.argv.append('-b')
import ROOT

def differentiate_stat_box(hist, movement=1, new_color=None, new_size=None, color_from_hist=True):
    """Move hist's stat box and change its line/text color. If
    movement is just an int, that number specifies how many units to
    move the box downward. If it is a 2-tuple of ints (m,n), the stat
    box will be moved to the left m units and down n units. A unit is
    the width or height of the stat box.

    Call TCanvas::Update first (and use TH1::Draw('sames') if
    appropriate) or else the stat box will not exist."""

    s = hist.FindObject('stats')

    if color_from_hist:
        new_color = hist.GetLineColor()

    if new_color is not None:
        s.SetTextColor(new_color)
        s.SetLineColor(new_color)

    if type(movement) == int:
        movement = (0,movement)
    m,n = movement
    
    x1,x2 = s.GetX1NDC(), s.GetX2NDC()
    y1,y2 = s.GetY1NDC(), s.GetY2NDC()

    if new_size is not None:
        x1 = x2 - new_size[0]
        y1 = y2 - new_size[1]

    s.SetX1NDC(x1 - (x2-x1)*m)
    s.SetX2NDC(x2 - (x2-x1)*m)
    s.SetY1NDC(y1 - (y2-y1)*n)
    s.SetY2NDC(y2 - (y2-y1)*n)

def print_histos(root_file,formats):
    """Print all the histograms from a ROOT file
    in each of the formats."""
    
    f = ROOT.TFile(root_file)
    os.system('mkdir -p plots/print_histos/'+root_file.split('/')[-1][:-5])
    for x in f.GetListOfKeys():
        d = f.GetDirectory(x.GetName())
        for y in d.GetListOfKeys():
            c1 = ROOT.TCanvas("c1","Canvas",600,600)        
            a = ''+y.GetName()
            h = d.Get(a)
            h.Draw()
            for i in formats:
                c1.Print('plots/print_histos/'+root_file.split('/')[-1][:-5]+'/'+ a +'.'+i)
            c1 = 0

