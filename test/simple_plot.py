import numpy as np
import ROOT as R
import os
import shutil
from argparse import ArgumentParser

www = '/eos/user/r/rlopezru/www/SignalModelStudies/SimHits_Analysis/'
indexfilepath = '/eos/user/r/rlopezru/index.php'
if not os.path.exists(www): 
    print(f'Creating directory {www}')
    os.makedirs(www); 
shutil.copyfile(indexfilepath, www+'index.php')

vardict = {}
vardict['firstHit_R']        = {'varexp': 'sqrt(simTrack_firstHit_x**2+simTrack_firstHit_y**2)', 'title': 'R of first hit', 'name': 'R (cm)', 'xmin': 0, 'xmax': 50, 'nbins': 25}
vardict['firstHit_x']        = {'varexp': 'simTrack_firstHit_x', 'title': 'X of first hit', 'name': 'x (cm)', 'xmin': -50, 'xmax': 50, 'nbins': 25}
vardict['firstHit_y']        = {'varexp': 'simTrack_firstHit_y', 'title': 'Y of first hit', 'name': 'y (cm)', 'xmin': -50, 'xmax': 50, 'nbins': 25}
vardict['firstHit_z']        = {'varexp': 'simTrack_firstHit_z', 'title': 'Z of first hit', 'name': 'z (cm)', 'xmin': -100, 'xmax': 100, 'nbins': 25}
vardict['firstHit_tof']      = {'varexp': 'simTrack_firstHit_tof', 'title': 'Time of flight of first hit', 'name': 't (ns)', 'xmin': -1, 'xmax': 24, 'nbins': 25}
vardict['simHit_R']        = {'varexp': 'sqrt(simHit_globalPosition_x**2+simHit_globalPosition_y**2)', 'title': 'R of sim hit', 'name': 'R (cm)', 'xmin': 0, 'xmax': 100, 'nbins': 50}
vardict['simHit_x']        = {'varexp': 'simHit_globalPosition_x', 'title': 'X of sim hit', 'name': 'x (cm)', 'xmin': -150, 'xmax': 150, 'nbins': 75}
vardict['simHit_y']        = {'varexp': 'simHit_globalPosition_y', 'title': 'Y of sim hit', 'name': 'y (cm)', 'xmin': -150, 'xmax': 150, 'nbins': 75}
vardict['simHit_z']        = {'varexp': 'simHit_globalPosition_z', 'title': 'Z of sim hit', 'name': 'z (cm)', 'xmin': -100, 'xmax': 100, 'nbins': 50}
vardict['simHit_tof']      = {'varexp': 'simHit_tof', 'title': 'Time of flight of sim hit', 'name': 't (ns)', 'xmin': -1, 'xmax': 24, 'nbins': 25}

###################################################################################
## Parsing
parser = ArgumentParser()
parser.add_argument('--infile', dest='inputfile', type=str, help='Input .root NTuple')
parser.add_argument('--var'   , dest='var_'     , type=str, help='Variable(s) to plot', nargs='*', choices=vardict.keys())
parser.add_argument('--sel'   , dest='sel_'     , type=str, help='Selection')
parser.add_argument('--imgname', dest='imgname',  type=str, help='Name of output plot image', default='plot')
args = parser.parse_args()
inputfile = args.inputfile
var_ = list(args.var_)
sel_ = args.sel_
imgname = args.imgname
###################################################################################

if len(var_) == 1: plotType = '1D'
if len(var_) == 2: plotType = '2D'
if len(var_) == 3: plotType = '3D'
treeName = 'Events'

def plot1D(h, option='HIST', outputfilename='plot'):
    c = R.TCanvas("c","c")
    c.cd()
    h.SetLineWidth(2)
    h.Draw(option)
    
    latex = R.TLatex()
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(R.kBlack);
    latex.SetTextFont(42);
    latex.SetTextAlign(11);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.15, 0.9, r"#bf{CMS} #it{Private work} #sqrt{s} = 13.6 TeV")
    c.SaveAs(www+outputfilename+".png")
    c.SaveAs(www+outputfilename+".pdf")

def plot2D(h, option='COLZ', outputfilename='plot'):
    c = R.TCanvas("c","c")
    c.cd()
    h.Draw(option)
    
    latex = R.TLatex()
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(R.kBlack);
    latex.SetTextFont(42);
    latex.SetTextAlign(11);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.15, 0.9, r"#bf{CMS} #it{Private work} #sqrt{s} = 13.6 TeV")
    c.SaveAs(www+outputfilename+".png")
    c.SaveAs(www+outputfilename+".pdf")

if __name__=='__main__':
    R.gROOT.ProcessLine('.L ./tdrstyle2D.C')
    R.gROOT.SetBatch(1)
    R.setTDRStyle()

    chain = R.TChain(treeName)
    chain.Add(inputfile)
    print(f"Entries: {chain.GetEntries()}")

    #f = R.TFile.Open(inputfile)
    #tree = f.Get("Events")

    #from tdrStyle import tdrGrid, setTDRStyle
    #setTDRStyle()

    if plotType == '1D':
        title          = f"{vardict[var_[0]]['title']};{vardict[var_[0]]['name']};N events"
        nbinsx = vardict[var_[0]]['nbins']
        xmin   = vardict[var_[0]]['xmin'] 
        xmax   = vardict[var_[0]]['xmax'] 
        hToDraw = R.TH1F("hToDraw", title, nbinsx, xmin, xmax)
        varDraw = f"{vardict[var_[0]]['varexp']}>>+hToDraw"
        selDraw = sel_
        optDraw = ''
        print(varDraw, selDraw, optDraw)
        chain.Draw(varDraw, selDraw, optDraw)
        hToDraw.SetTitle(hToDraw.GetTitle() + " [{0}]".format(hToDraw.GetEntries()))
        plot1D(hToDraw, option='HIST', outputfilename=imgname)

    if plotType == '2D':
        title          = f"{vardict[var_[0]]['title']} vs {vardict[var_[1]]['title']};{vardict[var_[0]]['name']};{vardict[var_[1]]['name']}"
        nbinsx, nbinsy = vardict[var_[0]]['nbins'], vardict[var_[1]]['nbins']
        xmin, ymin     = vardict[var_[0]]['xmin'] , vardict[var_[1]]['xmin']
        xmax, ymax     = vardict[var_[0]]['xmax'] , vardict[var_[1]]['xmax']
        hToDraw = R.TH2F("hToDraw", title, nbinsx, xmin, xmax, nbinsy, ymin, ymax)
        varDraw = f"{vardict[var_[1]]['varexp']}:{vardict[var_[0]]['varexp']}>>+hToDraw"
        selDraw = sel_
        optDraw = ''
        print(varDraw, selDraw, optDraw)
        chain.Draw(varDraw, selDraw, optDraw)
        hToDraw.SetTitle(hToDraw.GetTitle() + " [{0}]".format(hToDraw.GetEntries()))
        plot2D(hToDraw, option='COLZ', outputfilename=imgname)
