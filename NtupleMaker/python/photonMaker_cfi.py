import FWCore.ParameterSet.Config as cms

photonMaker = cms.EDProducer(
   "PhotonMaker",

   aliasPrefix = cms.untracked.string("photons"),

   year = cms.int32(-1), # Must be overriden by main_pset

   MVACuts = cms.VPSet(),

   photonsInputTag = cms.InputTag("slimmedPhotons"),

   pfcandsInputTag = cms.InputTag("packedPFCandidates"),

   rhoInputTag = cms.InputTag("fixedGridRhoFastjetAll"),

   EBHitsInputTag = cms.InputTag("reducedEgamma:reducedEBRecHits"),
   EEHitsInputTag = cms.InputTag("reducedEgamma:reducedEERecHits"),

   )

