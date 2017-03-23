import FWCore.ParameterSet.Config as cms
import CMS3.NtupleMaker.configProcessName as configProcessName

pfJetMaker = cms.EDProducer("PFJetMaker",
  aliasPrefix = cms.untracked.string("pfjets"),
  pfJetsInputTag                   = cms.InputTag("slimmedJets","",configProcessName.name),
  pfJetPtCut                       = cms.double(5.),
  pfCandidatesTag     = cms.InputTag("packedPFCandidates","",configProcessName.name),
)

pfJetPUPPIMaker = cms.EDProducer("PFJetMaker",
  aliasPrefix = cms.untracked.string("pfjets_puppi"),
  pfJetsInputTag                   = cms.InputTag("slimmedJetsPuppi","",configProcessName.name),
  pfJetPtCut                       = cms.double(20.), #miniAOD doesn't go lower than 20
  pfCandidatesTag     = cms.InputTag("packedPFCandidates","",configProcessName.name),
)
