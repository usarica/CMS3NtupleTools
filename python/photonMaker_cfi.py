import FWCore.ParameterSet.Config as cms

photonMaker = cms.EDProducer(
   "PhotonMaker",

   aliasPrefix = cms.untracked.string("photons"),

   year = cms.int32(-1), # Must be overriden by main_pset

   MVACuts = cms.VPSet(),

   photonsInputTag = cms.InputTag("slimmedPhotons"),

   rhoInputTag = cms.InputTag("fixedGridRhoFastjetAll"),

   )

