import FWCore.ParameterSet.Config as cms

patMuonMaker = cms.EDProducer("PATMuonMaker",
	aliasPrefix = cms.untracked.string("mus_pat"),
    # pat muon collection
    patMuonsInputTag  = cms.InputTag("selectedPatMuons"),
    recoMuonsInputTag = cms.InputTag("muons")                            
)


