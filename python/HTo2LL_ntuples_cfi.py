import FWCore.ParameterSet.Config as cms

ntuples = cms.EDAnalyzer('ntuplizer',
    nameOfOutput = cms.string('HTo2LongLivedTo2mu2jets_Ntuples_CMSSW_12_4_11_patch3_nsegments2.root'),
    isData = cms.bool(False),
    EventInfo = cms.InputTag("generator"),
    RunInfo = cms.InputTag("generator"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    displacedGlobalCollection = cms.InputTag("displacedGlobalMuons"),
    displacedStandAloneCollection = cms.InputTag("displacedStandAloneMuons"),
    displacedMuonCollection = cms.InputTag("slimmedDisplacedMuons"),

    #prescales  = cms.InputTag("patTrigger"),
    bits       = cms.InputTag("TriggerResults","","HLT"),
    #objects    = cms.InputTag("slimmedPatTrigger")
)
