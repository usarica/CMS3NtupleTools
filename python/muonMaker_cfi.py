import FWCore.ParameterSet.Config as cms

muonMaker = cms.EDFilter("MuonMaker",
	aliasPrefix = cms.untracked.string("mus"),
                         muonsInputTag    = cms.InputTag("muons"        ),                         
                         beamSpotInputTag = cms.InputTag("beamSpotMaker"),
                         tevMuonsName     = cms.string("tevMuons")
)


