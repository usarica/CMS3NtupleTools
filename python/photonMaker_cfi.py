import FWCore.ParameterSet.Config as cms
import CMS3.NtupleMaker.configProcessName as configProcessName

photonMaker = cms.EDProducer("PhotonMaker",
                             aliasPrefix = cms.untracked.string("photons"),
                             minEt       = cms.double(10.), #gev, min to keep
                             # Photon collection
                             # photonsInputTag = cms.InputTag("photons"),
                             photonsInputTag = cms.InputTag("slimmedPhotons"),
                             ecalRecHitsInputTag_EE = cms.InputTag("reducedEgamma","reducedEERecHits",configProcessName.name),
                             ecalRecHitsInputTag_EB = cms.InputTag("reducedEgamma","reducedEBRecHits",configProcessName.name),
                             # cms2scsseeddetidInputTag = cms.InputTag("scMaker"),

                             ebReducedRecHitCollectionTag = cms.InputTag("reducedEgamma:reducedEBRecHits"),
                             eeReducedRecHitCollectionTag = cms.InputTag("reducedEgamma:reducedEERecHits"),
                             esReducedRecHitCollectionTag = cms.InputTag("reducedEgamma:reducedESRecHits"),

                             )

