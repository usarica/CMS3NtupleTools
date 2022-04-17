import FWCore.ParameterSet.Config as cms


genMaker = cms.EDProducer(
   "GenMaker",

   aliasprefix = cms.untracked.string("genps"),
   year = cms.int32(-1),
   recoMode = cms.untracked.string("Run2_preUL"),

   xsec = cms.double(-1),
   BR = cms.double(-1),

   LHEInputTag = cms.InputTag("externalLHEProducer"),
   genEvtInfoInputTag  = cms.InputTag("generator"),
   prunedGenParticlesInputTag  = cms.InputTag("prunedGenParticles"),
   packedGenParticlesInputTag  = cms.InputTag("packedGenParticles"),

   genJetsInputTag = cms.InputTag("slimmedGenJets"),

   genMETInputTag = cms.InputTag("slimmedMETs"), # For genMET

   ntuplePackedGenParticles    = cms.bool(False), # default is False

   superMH = cms.double(125), # Higgs mass used in MELA SuperMELA

   candVVmode = cms.untracked.string("none"), # Has to correspond to MELAEvent::nCandidateVVModes
   decayVVmode = cms.int32(-1), # -1 means any decay mode
   doHiggsKinematics = cms.bool(False),
   lheMElist = cms.vstring(),

   kfactors = cms.VPSet()

   )

