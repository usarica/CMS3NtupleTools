import FWCore.ParameterSet.Config as cms

patMuonMaker = cms.EDFilter("PATMuonMaker",
    # pat muon collection
    patMuonsInputTag  = cms.InputTag("selectedPatMuons"),
    recoMuonsInputTag = cms.InputTag("muons")                            
)


