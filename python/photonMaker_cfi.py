import FWCore.ParameterSet.Config as cms

photonMaker = cms.EDProducer("PhotonMaker",
                             aliasPrefix = cms.untracked.string("photons"),
                             minEt       = cms.double(10.), #gev, min to keep
                             # Photon collection
                             # photonsInputTag = cms.InputTag("photons"),
                             photonsInputTag = cms.InputTag("slimmedPhotons"),
                             ecalRecHitsInputTag_EE = cms.InputTag("reducedEgamma","reducedEERecHits","PAT"),
                             ecalRecHitsInputTag_EB = cms.InputTag("reducedEgamma","reducedEBRecHits","PAT"),
                             # cms2scsseeddetidInputTag = cms.InputTag("scMaker"),
                             )

