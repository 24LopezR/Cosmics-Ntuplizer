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

from Configuration.AlCa.GlobalTag import GlobalTag

# Select number of events to be processed
nEvents = 2000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nEvents) )

# Read events
#listOfFiles = ['file:/eos/user/r/rlopezru/DisplacedDimuons/Stop_500_100000mm_WF_test_LimitTau0/output_GS.root']
listOfFiles = ['file:/eos/user/r/rlopezru/DisplacedDimuons/Stop_500_100000mm_WF_test_PR40820/output_GS.root']
'''listOfFiles = ['/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_163.root',
'/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_8.root',
'/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_7.root',
'/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_6.root',
'/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_9.root',
'/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_1.root',
'/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_18.root',
'/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_37.root',
'/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_146.root',
'/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_190.root',
'/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_118.root',
'/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_47.root',
'/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_193.root',
'/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_22.root',
'/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_17.root',
'/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_103.root',
'/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_10.root',
'/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_104.root',
'/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_45.root',
'/store/user/rlopezru/StopToMuB_M_500_100000mm_13p6TeV_2022MC/GenSim/240115_093806/0000/output_79.root']'''
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( listOfFiles ),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(0)
  )
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2023_realistic_v14')

## Define the process to run 
## 
process.load("Analysis.Muon-Ntuplizer.Muon_HTo2LongLivedTo4Mu_PREMIXRAW_ntuples_cfi")
process.p = cms.EndPath(process.ntuples)
