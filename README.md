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
  <li> <strong>plugins/</strong>: which contains the plugins (EDAnalyzer's) where the analyzers are defined in .cc files. These are the main code.</li>
  <li> <strong>python/</strong>: which contains cfi files to setup the sequences that run with the plugins contained in plugins/. A sequence is an specific configuration of the parameters that run with one of the plugins defined in plugins. One single plugin may have different sequences defined in the same or multiple files.</li> 
  <li> <strong>test/</strong>: which contains cfg files to run the sequences defined in the python/ folder.</li>
</ul>

### EDAnalyzer plugin

(to be completed)

### Configuration cfi files and parameters

Parameters are values that are defined "per sequence" and serve to configure how the code should run. For example, if we want to run the same EDAnalyzer for both data and Monte Carlo we may need to know if the generation variables can be accesed or not as if we try to access them in data we may likely get an error. This could be done via parameters.

The parameter values are defined in a cfi file, whose structure is as follows e.g. python/ntuples_cfi.py:
https://github.com/CeliaFernandez/standard-Ntuplizer/blob/ac4da89ab55d8d01a4b70a8a585680ea5a8961d7/python/ntuples_cfi.py#L1-L12
where ```ntuples``` is the name of the sequence (instance of EDAnalyzer) and ```'ntuplizer'``` matches the name of the plugin we want to run.

Each parameter as a variable that is declared in the EDAnalyzer constructor as a private variable that can be used when the code is running. For example, to indicate if we are running on data samples we have can define a bool variable ```isData```:
https://github.com/CeliaFernandez/standard-Ntuplizer/blob/5e3b77f976d88d9c812b7a5cff1a32b70b0cfe25/plugins/ntuplizer.cc#L73
and we also define a ```isData``` parameter in the cfi file ntuples_cfi.py:
https://github.com/CeliaFernandez/standard-Ntuplizer/blob/ac4da89ab55d8d01a4b70a8a585680ea5a8961d7/python/ntuples_cfi.py#L5

The ```isData``` variable is initiated with the value set in the cfi file. The values defined there can be accessed in the constructor with ```iConfig``` variable:
https://github.com/CeliaFernandez/standard-Ntuplizer/blob/5e3b77f976d88d9c812b7a5cff1a32b70b0cfe25/plugins/ntuplizer.cc#L115

To access ```iConfig``` in other parts of the code is useful to define a ```edm::ParameterSet``` variable, which in our case is called ```parameters```and it is declared in the class definition as
https://github.com/CeliaFernandez/standard-Ntuplizer/blob/5e3b77f976d88d9c812b7a5cff1a32b70b0cfe25/plugins/ntuplizer.cc#L49
and initiated in the constructor as a copy of ```iConfig```:
https://github.com/CeliaFernandez/standard-Ntuplizer/blob/5e3b77f976d88d9c812b7a5cff1a32b70b0cfe25/plugins/ntuplizer.cc#L119

Then we can assign the correct value to ```isData``` before the analyzer runs in ```beginJob()``` function like:
https://github.com/CeliaFernandez/standard-Ntuplizer/blob/5e3b77f976d88d9c812b7a5cff1a32b70b0cfe25/plugins/ntuplizer.cc#L147


### Configuration cfg files

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

### How to add a parameter


### How to add new variables of an existing collection

1) We first need to declare a new variable that will act as a container for the value we want to store e.g. the number of displacedGlobalMuon tracks ```ndgl```. It is defined in the constructor of the EDAnalyzer as a private variable (although it could be also a global variable):
https://github.com/CeliaFernandez/standard-Ntuplizer/blob/8656711d7fa7d640a9ec160daa955738d283720e/plugins/ntuplizer.cc#L79

2) We then need to link this variable's address ```&ndlg``` to the TTree branch. This is done at the beginning, where the TTree is created in ```beginJob()```:
https://github.com/CeliaFernandez/standard-Ntuplizer/blob/8656711d7fa7d640a9ec160daa955738d283720e/plugins/ntuplizer.cc#L147

3) This variable will be saved inside the TTree once the Fill() command is executed:
https://github.com/CeliaFernandez/standard-Ntuplizer/blob/8656711d7fa7d640a9ec160daa955738d283720e/plugins/ntuplizer.cc#L244
So the value of this variable should be assigned before that like:
https://github.com/CeliaFernandez/standard-Ntuplizer/blob/8656711d7fa7d640a9ec160daa955738d283720e/plugins/ntuplizer.cc#L210
https://github.com/CeliaFernandez/standard-Ntuplizer/blob/8656711d7fa7d640a9ec160daa955738d283720e/plugins/ntuplizer.cc#L216

4) It is possible to save an array of values. In this case we must define a container array with a long enough length (as it is declared when the analyzer is defined and bound to be always the same for every event). For example, for the pt of the stored displacedGlobalMuon tracks we define a default array container of 200 entries (no event has more than 200 displaced global muons):
https://github.com/CeliaFernandez/standard-Ntuplizer/blob/df839c1226bc0d8acb6d3d5408cd9dad6daa94a6/plugins/ntuplizer.cc#L80
And then set the array length for the branch. If the array length is dependent of the event, we can make the array branch length dependent of another TTree variable. In this case we have as much pt measurements as number of displacedGlobalMuon tracks i.e. ```ndgl```:
https://github.com/CeliaFernandez/standard-Ntuplizer/blob/df839c1226bc0d8acb6d3d5408cd9dad6daa94a6/plugins/ntuplizer.cc#L148
Since the array is itself a contained of adresses, it is not needed to include the ```&```. The pt values are filled per displacedGlobalMuon track in a loop:
https://github.com/CeliaFernandez/standard-Ntuplizer/blob/df839c1226bc0d8acb6d3d5408cd9dad6daa94a6/plugins/ntuplizer.cc#L213

### How to read a new collection

