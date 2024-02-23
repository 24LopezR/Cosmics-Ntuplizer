import FWCore.ParameterSet.Config as cms

ntuples = cms.EDAnalyzer('ntuplizer',
    nameOfOutput = cms.string('NTuples_RelValSingleMuFlatPt2To100_v0.root'),
    isData                        = cms.bool(True),
    EventInfo                     = cms.InputTag("generator"),
    RunInfo                       = cms.InputTag("generator"),
    BeamSpot                      = cms.InputTag("offlineBeamSpot"),
    muonCollection                = cms.InputTag("muons"),
    genParticleCollection         = cms.InputTag("prunedGenParticles"),
    PrimaryVertexCollection       = cms.InputTag("offlinePrimaryVertices"),

    #prescales  = cms.InputTag("patTrigger"),
    bits       = cms.InputTag("TriggerResults","","HLT"),
    #objects    = cms.InputTag("slimmedPatTrigger")
)
