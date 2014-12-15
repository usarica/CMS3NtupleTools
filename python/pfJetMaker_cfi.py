import FWCore.ParameterSet.Config as cms

pfJetMaker = cms.EDProducer("PFJetMaker",
  pfJetsInputTag                   = cms.InputTag("slimmedJets"),
  pfCandidatesTag                  = cms.InputTag("packedPFCandidates"),
  pfJetPtCut                       = cms.double(5.),
  L1File          = cms.string("CMS3/NtupleMaker/data/JEC_phys14/PHYS14_V1_MC_L1FastJet_AK4PFchs.txt"),
  L2File          = cms.string("CMS3/NtupleMaker/data/JEC_phys14/PHYS14_V1_MC_L2Relative_AK4PFchs.txt"),
  L3File          = cms.string("CMS3/NtupleMaker/data/JEC_phys14/PHYS14_V1_MC_L3Absolute_AK4PFchs.txt")

  #PFJetCorrectorL2L3               = cms.string("ak5PFCHSL2L3"),
)


pfchsJetMaker = cms.EDProducer("PFJetMaker",
  pfJetsInputTag                   = cms.InputTag("prunedUncorrectedCMS2Jets", "pfchsjet"),
  pfCandidatesTag                  = cms.InputTag("particleFlow"),
  pfJetPtCut                       = cms.double(5.),
  PFJetCorrectorL2L3               = cms.string("ak5PFCHSL2L3"),
  PFJetCorrectorL1L2L3             = cms.string("ak5PFCHSL1L2L3"),
  PFJetCorrectorL1FastL2L3         = cms.string("ak5PFCHSL1FastL2L3"),
  PFJetCorrectorL1Fast             = cms.string("ak5PFCHSL1Fastjet"),
  PFJetCorrectorL1FastL2L3residual = cms.string("ak5PFCHSL1FastL2L3Residual"),
  AliasPrefix                      = cms.string("pfchsjets")
                               
)
