import FWCore.ParameterSet.Config as cms


electronMaker = cms.EDProducer(
   "ElectronMaker",

   aliasPrefix = cms.untracked.string("electrons"),
   year = cms.int32(-1), # Must be overriden by main_pset

   MVACuts = cms.VPSet(),

   # Electron collection
   electronsInputTag   = cms.InputTag("slimmedElectrons"),

   # Beamspot and vertex
   beamSpotInputTag  = cms.InputTag("offlineBeamSpot"),
   vtxInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),

   # Track collections
   trksInputTag      = cms.InputTag("generalTracks"),
   gsftracksInputTag = cms.InputTag("electronGsfTracks"),

   # Reco conversions
   recoConversionInputTag = cms.InputTag("reducedEgamma:reducedConversions"),

   # Rec. hits
   EBHitsInputTag = cms.InputTag("reducedEgamma:reducedEBRecHits"),
   EEHitsInputTag = cms.InputTag("reducedEgamma:reducedEERecHits"),
   ESHitsInputTag = cms.InputTag("reducedEgamma:reducedESRecHits"),

   # Event rhos
   rhoInputTag = cms.InputTag("fixedGridRhoFastjetAll"),
   rhoCaloInputTag = cms.InputTag("fixedGridRhoFastjetAllCalo"),

   )

