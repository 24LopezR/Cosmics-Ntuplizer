# Standard Ntuplizer

This code serves as a template for new Ntuplizers to work with CMSSW. Basic instructions for installation and standard modifications can be found below.
It presents an example where a ROOT tree if filled with plain Ntuples made from pat::Muon variables read from MiniAOD. It is configured to read Cosmic data from the NoBPTX dataset.

The Ntuplizer is an EDAnalyzer. More information about this class and its structure can be found in https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookWriteFrameworkModule.

## How to install

This branch is to be used in release CMSSW_14_1_0_pre2 or later. Commands to setup the analyzer are:

```
cmsrel CMSSW_14_1_0_pre2

cd CMSSW_14_1_0_pre2/src

cmsenv

mkdir Analysis

cd Analysis

git clone git@github.com:24LopezR/Muon-Ntuplizer.git -b TkAlIssue_IsoEm

scram b -j 8
```
