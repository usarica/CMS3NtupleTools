import FWCore.ParameterSet.Config as cms

cms3ntuple = cms.EDAnalyzer(
   "CMS3Ntuplizer",

   year = cms.int32(-1), # Must be overriden by main_pset
   treename = cms.untracked.string("Events"),

   isMC = cms.bool(False),
   is80x = cms.bool(False),

   prefiringWeightsTag = cms.untracked.string(""),

   electronSrc = cms.InputTag("electronMaker"),
   photonSrc = cms.InputTag("photonMaker"),
   muonSrc = cms.InputTag("muonMaker"),
   ak4jetSrc = cms.InputTag("pfJetMaker"),
   ak8jetSrc = cms.InputTag("subJetMaker"),
   isotrackSrc = cms.InputTag("isoTrackMaker"),
   pfcandSrc = cms.InputTag("packedPFCandidates"),

   pfmetSrc = cms.InputTag("pfmetMaker"),
   puppimetSrc = cms.InputTag("pfmetpuppiMaker"),

   vtxSrc = cms.InputTag("offlineSlimmedPrimaryVertices"),

   rhoSrc = cms.InputTag("fixedGridRhoFastjetAll"),
   triggerInfoSrc = cms.InputTag("hltMaker"),
   puInfoSrc = cms.InputTag("slimmedAddPileupInfo"),
   metFilterInfoSrc = cms.InputTag("metFilterMaker"),

   genInfoSrc = cms.InputTag("genMaker"),

   keepGenParticles = cms.untracked.string("reducedfinalstates"),
   prunedGenParticlesSrc = cms.InputTag("prunedGenParticles"),
   packedGenParticlesSrc = cms.InputTag("packedGenParticles"),

   keepGenJets = cms.bool(True),
   genAK4JetsSrc  = cms.InputTag("slimmedGenJets"),
   genAK8JetsSrc  = cms.InputTag("slimmedGenJetsAK8"),

   minNmuons = cms.int32(-1),
   minNelectrons = cms.int32(-1),
   minNleptons = cms.int32(-1),
   minNphotons = cms.int32(-1),
   minNak4jets = cms.int32(-1),
   minNak8jets = cms.int32(-1),

   )

