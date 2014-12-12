import FWCore.ParameterSet.Config as cms

pfJetMaker = cms.EDProducer("PFJetMaker",
  pfJetsInputTag                   = cms.InputTag("slimmedJets"),
  pfCandidatesTag                  = cms.InputTag("packedPFCandidates"),
  pfJetPtCut                       = cms.double(5.),
  L1File          = cms.string("CMS2/NtupleMaker/data/JEC_phys14/PHYS14_V1_MC_L1FastJet_AK4PFchs.txt"),
  L2File          = cms.string("CMS2/NtupleMaker/data/JEC_phys14/PHYS14_V1_MC_L2Relative_AK4PFchs.txt"),
  L3File          = cms.string("CMS2/NtupleMaker/data/JEC_phys14/PHYS14_V1_MC_L3Absolute_AK4PFchs.txt")

  #PFJetCorrectorL2L3               = cms.string("ak5PFCHSL2L3"),
)


