# Standard Ntuplizer

This code serves as a template for new Ntuplizers to work with CMSSW. Basic instructions for installation and standard modifications can be found below.
It presents an example where a ROOT tree if filled with plain Ntuples made from pat::Muon variables read from MiniAOD. It is configured to read Cosmic data from the NoBPTX dataset.

The Ntuplizer is an EDAnalyzer. More information about this class and its structure can be found in https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookWriteFrameworkModule.

## How to install

Recommended release for this analyzer is CMSSW_12_4_0 or later. Commands to setup the analyzer are:

```
cmsrel CMSSW_12_4_0

cd CMSSW_12_4_0/src

cmsenv

mkdir Analysis

cd Analysis

git clone git@github.com:CeliaFernandez/standard-Ntuplizer.git

scram b -j 8
```


## Ntuplizer structure

<p> The analyzer consists of three folders: </p> 
<ul>
  <li> plugins: which contains the plugins (EDAnalyzer's) where the analyzers are defined in `.cc` files. These are the main code.</li>
  <li> python : which contains cfi files to setup the sequences that run with the plugins contained in `plugins/`. A sequence is an specific configuration of the parameters that run with one of the plugins defined in `plugins`. One single plugin may have different sequences defined in the same or multiple files.</li> 
  <li> test: which contains cfg files to run the sequences defined in the `python/` folder.</li>
</ul>

## How to run

This example runs with a file of the 2023 NoBPTX dataset that may need to be accessed throught xrootd. Make sure that you have a valid proxy before running and do at least once:

```
voms-proxy-init --voms cms
```

Then you can run the Ntuplizer with the setup configuration through the cfg file:

```
cmsRun test/runNtuplizer_cfg.py
```


## Quick start: How to modify the analyzer

In this section (to be completed) there are several examples of how modify the existing analyzer.

### How to add new variables

### How to read a new collection

