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
nEvents = 100
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nEvents) )

# Read events
listOfFiles = [
'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/0f73e103-273c-416c-8fdd-e322aa8de063.root',
'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/1abd3dd2-bc65-4028-995e-9ad505672d08.root',
'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/1ffbad78-d003-4a61-869d-46ac2087e80c.root',
'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/4085c749-1d0b-4264-9a46-6ed926475a6b.root',
'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/560602b8-bd74-4228-bd95-a49099cc031e.root',
'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/59292bf8-70e0-4982-a4af-6f27370dfba8.root',
'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/70b39b48-f2b8-4584-b268-e464aaba3531.root',
'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/734e7003-ac98-462c-87c8-66dedc28fe8a.root',
'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/a807ed3e-e844-41b3-b6e1-ea16e6ae9632.root',
'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/af7f0f95-05f8-42cf-9792-0f8653c651eb.root',
'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/b0a8ead8-0fcd-48d5-822f-a323925b101a.root',
'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/b579e360-bb6a-4707-a4c9-9b01f2ce4e7d.root',
'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/b996f92a-73fd-4e45-9e37-20c7e86c2061.root',
'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/b9a33790-ba83-4715-b28f-b6c95e19c7ca.root',
'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/bc8feccd-4e0b-402d-8741-edcb0e7e571b.root',
'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/bc90c644-bb51-4488-b1d0-b42b2ec7a991.root',
'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/d3fecec7-889d-442a-9249-fdf2eadbd24b.root',
'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/d99c2acd-545b-4095-80a8-8f57e77b4551.root',
'/store/relval/CMSSW_14_0_0_pre2/RelValSingleMuFlatPt2To100/GEN-SIM-RECO/133X_mcRun4_realistic_v1_STD_2026D98_noPU_RV229-v1/2580000/db45efef-81c1-4da8-99d8-6e491edf0237.root']
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( listOfFiles ),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(0)
  )
process.GlobalTag = GlobalTag(process.GlobalTag, '130X_mcRun3_2023_realistic_v14')

## Define the process to run 
## 
process.load("Muon_RECO_ntuples_cfi")

process.p = cms.EndPath(process.ntuples)
