# Standard Ntuplizer

This code serves as a template for new Ntuplizers to work with CMSSW. Basic instructions for installation and standard modifications can be found below.
It presents an example where a ROOT tree if filled with plain Ntuples made from pat::Muon variables read from MiniAOD. It is configured to read Cosmic data from the NoBPTX dataset.

The Ntuplizer is an EDAnalyzer. More information about this class and its structure can be found in https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookWriteFrameworkModule.

## How to install

This branch is to be used in release CMSSW_12_4_11_patch3. Commands to setup the analyzer are:

```
cmsrel CMSSW_12_4_11_patch3

cd CMSSW_12_4_11_patch3/src

cmsenv

mkdir Analysis

cd Analysis

git clone git@github.com:24LopezR/Muon-Ntuplizer.git -b GenSimInfoReader

scram b -j 8
```

## How to run

In file test/HTo2LongLivedTo4mu_1600_50_100000_PREMIXRAW_runNtuplizer.py, modify input file name.

```
cd Muon-Ntuplizer/test

cmsRun HTo2LongLivedTo4mu_1600_50_100000_PREMIXRAW_runNtuplizer.py &> logHTo2LL_1600_50_100000.log
```
