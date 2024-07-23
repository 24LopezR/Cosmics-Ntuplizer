import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("demo")
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.Services_cff')

# Debug printout and summary.
process.load("FWCore.MessageService.MessageLogger_cfi")

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  # Set up multi-threaded run. Must be consistent with config.JobType.numCores in crab_cfg.py.
  numberOfThreads=cms.untracked.uint32(8)
)

from Configuration.AlCa.GlobalTag import GlobalTag

# Select number of events to be processed
nEvents = -1
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nEvents) )

# Read events
#listOfFiles = [f'/store/user/rlopezru/RelValZMM_14/GEN-SIM-RECO/240710_154415/0000/output_{i}.root' for i in range(1,20)]
listOfFiles = ['file:/eos/user/r/rlopezru/MuonReco/TKAlIssue/RECO/FinalTest_ZMM_14TeV_GEN-SIM-RECO_target_v1.root']
#listOfFiles = ['file:/afs/cern.ch/user/r/rlopezru/private/MuonPOG/TkAlIssue/CMSSW_14_1_0_pre3/src/output.root']
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( listOfFiles ),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(0)
  )
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2023_realistic_v14')

## Define the process to run 
## 
process.load("Analysis.Muon-Ntuplizer.Muon_RECO_RelValZMM_14_ntuples_cfi")

process.ntuples.nameOfOutput = '/eos/user/r/rlopezru/MuonReco/TKAlIssue/NTuples_FinalTest_RelValZMM_14_target_v2.root'
#process.ntuples.nameOfOutput = '2ndFix_test.root'

process.p = cms.EndPath(process.ntuples)
