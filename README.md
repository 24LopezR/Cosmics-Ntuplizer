# Standard Ntuplizer

This code serves as a template for new Ntuplizers to work with CMSSW. Basic instructions for installation and standard modifications can be found below.
It presents an example where a ROOT tree if filled with plain Ntuples made from pat::Muon variables read from MiniAOD. It is configured to read Cosmic data from the NoBPTX dataset.

The Ntuplizer is an EDAnalyzer. More information about this class and its structure can be found in https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookWriteFrameworkModule.

## How to install together with my debugging changes in CMSSW

This branch is to be used in release CMSSW_14_1_0_pre3. Commands to setup the analyzer are:

```
cmsrel CMSSW_14_1_0_pre3

cd CMSSW_14_1_0_pre3/src

cmsenv

git-cms-init

git checkout -b [new-branch]

git-cms-addpkg Configuration/Generator SimG4Core/Generators SimG4Core/Application SimG4Core/Notification SimG4CMS/Tracker

git cms-merge-topic 24LopezR:SimDisplacedSUSY

mkdir Analysis

cd Analysis

git clone git@github.com:24LopezR/Muon-Ntuplizer.git -b GenSimInfoReader

cd ..

scram b -j 8
```

## How to run over SMuon samples

Generate some events (2 for this cfg)

```
cd /afs/cern.ch/user/r/rlopezru/public/ForSlava/SMuonToMuGravitino-M_500_ctau_100000mm_TuneCP5_13p6TeV_pythia8_cff.py ./SMuonToMuGravitino-M_500_ctau_100000mm_TuneCP5_13p6TeV_pythia8_cff.py

cmsRun SMuonToMuGravitino-M_500_ctau_100000mm_TuneCP5_13p6TeV_pythia8_cff.py &> log_GENSIM.log
```

Analyze GENSIM output file to print SimHits and SimTracks.
In file test/SMuonToMuGravitino_500_100000mm_GENSIM_runNtuplizer.py, you may need to modify input file name.

```
cd Analysis/Muon-Ntuplizer/test

cmsRun SMuonToMuGravitino_500_100000mm_GENSIM_runNtuplizer.py &> log_Analyzer.log
```
