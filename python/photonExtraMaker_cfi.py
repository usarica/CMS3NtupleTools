import FWCore.ParameterSet.Config as cms

photonExtraMaker = cms.EDProducer("PhotonExtraMaker",
                             aliasPrefix = cms.untracked.string("photons"),
                             minEt       = cms.double(10.), #gev, min to keep
                             photonsInputTag = cms.InputTag("slimmedPhotons"),
                             )

