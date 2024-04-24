import numpy as np
import ROOT as R
import os
import shutil
from argparse import ArgumentParser

www = '/eos/user/r/rlopezru/www/MuonReco/TkAlIssue/'
indexfilepath = '/eos/user/r/rlopezru/index.php'
if not os.path.exists(www): 
    print(f'Creating directory {www}')
    os.makedirs(www); 
shutil.copyfile(indexfilepath, www+'index.php')

vardict = {}
vardict['mu_phi']        = {'name': '#phi', 'xmin': -3.2, 'xmax': 3.2, 'nbins': 100}
vardict['mu_eta']        = {'name': '#eta', 'xmin': -3., 'xmax': 3., 'nbins': 100}
vardict['mu_iso03_emEt'] = {'name': 'Isolation EmEt', 'xmin': 0., 'xmax': 10, 'nbins': 100}
vardict['mu_iso03_emEt/mu_pt'] = {'name': 'Isolation EmEt / pT', 'xmin': 0., 'xmax': 0.2, 'nbins': 100}
vardict['mu_iso03_hadEt'] = {'name': 'Isolation hadEt', 'xmin': 0., 'xmax': 0.2, 'nbins': 50}
vardict['mu_iso03_hadEt/mu_pt'] = {'name': 'Isolation hadEt / pT', 'xmin': 0., 'xmax': 0.2, 'nbins': 50}
vardict['mu_iso03_sumPt'] = {'name': 'Isolation sumPt', 'xmin': 0., 'xmax': 0.2, 'nbins': 50}
vardict['mu_iso03_sumPt/mu_pt'] = {'name': 'Isolation sumPt / pT', 'xmin': 0., 'xmax': 0.2, 'nbins': 50}

###################################################################################
## Parsing
parser = ArgumentParser()
parser.add_argument('--infile', dest='inputfile', type=str, help='Input .root NTuple')
parser.add_argument('--var'   , dest='var_'     , type=str, help='Variable(s) to plot', nargs='*', choices=vardict.keys())
parser.add_argument('--imgname', dest='imgname',  type=str, help='Name of output plot image', default='plot')
args = parser.parse_args()
inputfile = args.inputfile
var_ = list(args.var_)
imgname = args.imgname
###################################################################################

if len(var_) == 1: plotType = '1D'
if len(var_) == 2: plotType = '2D'
if len(var_) == 3: plotType = '3D'
treeName = 'Events'

def plot1D(h, option='HIST', outputfilename='plot'):
    c = R.TCanvas("c","c")
    c.cd()
    h.GetXaxis().SetTitle(h.GetXaxis().GetTitle())
    h.SetTitle(h.GetTitle())
    h.GetYaxis().SetTitle("N events")
    h.SetLineWidth(2)
    h.Draw(option)
    
    '''latex = R.TLatex()
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(R.kBlack);
    latex.SetTextFont(42);
    latex.SetTextAlign(11);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.15, 0.94, r"#bf{CMS} #it{Simulation} #sqrt{s} = 13.6 TeV")'''
    c.SaveAs(www+outputfilename+".png")
    c.SaveAs(www+outputfilename+".pdf")

def plot2D(h, option='COLZ', outputfilename='plot'):
    c = R.TCanvas("c","c")
    c.cd()
    h.GetXaxis().SetTitle(h.GetXaxis().GetTitle())
    h.GetYaxis().SetTitle(h.GetYaxis().GetTitle())
    h.SetTitle(h.GetTitle())
    h.Draw(option)
    
    '''latex = R.TLatex()
    latex.SetNDC();
    latex.SetTextAngle(0);
    latex.SetTextColor(R.kBlack);
    latex.SetTextFont(42);
    latex.SetTextAlign(11);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.15, 0.94, r"#bf{CMS} #it{Simulation} #sqrt{s} = 13.6 TeV")'''
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
        title          = f"{vardict[var_[0]]['name']}"
        nbinsx = vardict[var_[0]]['nbins']
        xmin   = vardict[var_[0]]['xmin'] 
        xmax   = vardict[var_[0]]['xmax'] 
        hToDraw = R.TH1F("hToDraw", title, nbinsx, xmin, xmax)
        varDraw = f'{var_[0]}>>+hToDraw'
        selDraw = '1'
        optDraw = ''
        print(varDraw, selDraw, optDraw)
        chain.Draw(varDraw, selDraw, optDraw)
        hToDraw.SetTitle(hToDraw.GetTitle() + " [{0}]".format(hToDraw.GetEntries()))
        plot1D(hToDraw, option='HIST', outputfilename=imgname)

    if plotType == '2D':
        title          = f"ZMM 14TeV: {vardict[var_[0]]['name']} vs {vardict[var_[1]]['name']};{vardict[var_[0]]['name']};{vardict[var_[1]]['name']}"
        nbinsx, nbinsy = vardict[var_[0]]['nbins'], vardict[var_[1]]['nbins']
        xmin, ymin     = vardict[var_[0]]['xmin'] , vardict[var_[1]]['xmin']
        xmax, ymax     = vardict[var_[0]]['xmax'] , vardict[var_[1]]['xmax']
        hToDraw = R.TH2F("hToDraw", title, nbinsx, xmin, xmax, nbinsy, ymin, ymax)
        varDraw = f'{var_[1]}:{var_[0]}>>+hToDraw'
        selDraw = '1'
        optDraw = ''
        print(varDraw, selDraw, optDraw)
        chain.Draw(varDraw, selDraw, optDraw)
        hToDraw.SetTitle(hToDraw.GetTitle() + " [{0}]".format(hToDraw.GetEntries()))
        plot2D(hToDraw, option='COLZ', outputfilename=imgname)
