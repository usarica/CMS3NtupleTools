import FWCore.ParameterSet.Config as cms

muonMaker = cms.EDFilter("MuonMaker",
    # muon collection
    muonsInputTag = cms.InputTag("selectedLayer1Muons"),
    #track collection
    tracksInputTag = cms.InputTag("generalTracks"),
    #beamSpot (from CMS2)
    beamSpotInputTag = cms.InputTag("beamSpotMaker")
)


