import FWCore.ParameterSet.Config as cms


electronMaker = cms.EDProducer(
   "ElectronMaker",

   aliasPrefix = cms.untracked.string("electrons"),
   year = cms.int32(-1), # Must be overriden by main_pset

   MVACuts = cms.VPSet(),

   # Electron collection
   electronsInputTag   = cms.InputTag("slimmedElectrons"),

   # Beamspot and vertex
   beamSpotInputTag  = cms.InputTag("beamSpotMaker","evtbsp4"),
   vtxInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),

   # Track collections
   trksInputTag      = cms.InputTag("generalTracks"),
   gsftracksInputTag = cms.InputTag("electronGsfTracks"),

   # reco conversions
   recoConversionInputTag = cms.InputTag("reducedEgamma:reducedConversions"),

   ebReducedRecHitCollectionTag = cms.InputTag("reducedEgamma:reducedEBRecHits"),
   eeReducedRecHitCollectionTag = cms.InputTag("reducedEgamma:reducedEERecHits"),
   esReducedRecHitCollectionTag = cms.InputTag("reducedEgamma:reducedESRecHits"),

   # Event rho
   rhoInputTag = cms.InputTag("fixedGridRhoFastjetAll"),

   )

