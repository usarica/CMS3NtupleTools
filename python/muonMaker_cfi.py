import FWCore.ParameterSet.Config as cms

muonMaker = cms.EDFilter("MuonMaker",
                         muonsInputTag    = cms.InputTag("muons"        ),                         
                         beamSpotInputTag = cms.InputTag("beamSpotMaker"),
                         tevMuonsName     = cms.string("tevMuons")
)


