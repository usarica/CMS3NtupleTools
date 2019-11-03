import FWCore.ParameterSet.Config as cms


genMaker = cms.EDProducer(
   "GenMaker",

   aliasprefix = cms.untracked.string("genps"),
   year = cms.int32(-1),

   LHEInputTag = cms.InputTag("externalLHEProducer"),
   genEvtInfoInputTag  = cms.InputTag("generator" ),
   prunedGenParticlesInputTag  = cms.InputTag("prunedGenParticles" ),
   packedGenParticlesInputTag  = cms.InputTag("packedGenParticles" ), # Assign Status "1111" to these to avoid duplication. Only save p4, ID, status

   genMETInputTag = cms.InputTag("slimmedMETs"), # For genMET

   ntuplePackedGenParticles    = cms.bool(False), # default is False

   candVVmode = cms.int32(7), # Has to correspond to MELAEvent::nCandidateVVModes
   decayVVmode = cms.int32(-1), # -1 means any decay mode
   doHiggsKinematics = cms.bool(False),
   lheMElist = cms.vstring(),

   )

