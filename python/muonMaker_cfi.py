import FWCore.ParameterSet.Config as cms

muonMaker = cms.EDProducer("MuonMaker",
	aliasPrefix = cms.untracked.string("mus"),
                         muonsInputTag    = cms.InputTag("muons"        ),                         
                         beamSpotInputTag = cms.InputTag("beamSpotMaker"),
                         pfCandsInputTag = cms.InputTag("particleFlow"),
                         vtxInputTag = cms.InputTag("offlinePrimaryVertices"),
                         tevMuonsName     = cms.string("tevMuons")
)


