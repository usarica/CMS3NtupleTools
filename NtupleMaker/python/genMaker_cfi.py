import FWCore.ParameterSet.Config as cms


genMaker = cms.EDProducer(
   "GenMaker",

   aliasprefix = cms.untracked.string("genps"),
   year = cms.int32(-1),

   xsec = cms.double(-1),
   BR = cms.double(-1),

   LHEInputTag = cms.InputTag("externalLHEProducer"),
   genEvtInfoInputTag  = cms.InputTag("generator"),
   prunedGenParticlesInputTag  = cms.InputTag("prunedGenParticles"),
   packedGenParticlesInputTag  = cms.InputTag("packedGenParticles"),

   genMETInputTag = cms.InputTag("slimmedMETs"), # For genMET

   ntuplePackedGenParticles    = cms.bool(False), # default is False

   sqrts = cms.int32(13), # CoM in TeV
   superMH = cms.double(125), # Higgs mass used in MELA SuperMELA

   candVVmode = cms.int32(7), # Has to correspond to MELAEvent::nCandidateVVModes
   decayVVmode = cms.int32(-1), # -1 means any decay mode
   doHiggsKinematics = cms.bool(False),
   lheMElist = cms.vstring(),

   )

