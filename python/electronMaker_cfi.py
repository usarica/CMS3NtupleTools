import FWCore.ParameterSet.Config as cms

electronMaker = cms.EDFilter("ElectronMaker",
    # jet collection
    electronsInputTag = cms.InputTag("allLayer1Electrons"),
    #track collection
    tracksInputTag = cms.InputTag("generalTracks"),
    # MC particles
    genParticlesInputTag = cms.InputTag("genParticles"),
    beamSpotInputTag = cms.InputTag("beamSpotMaker")

)


