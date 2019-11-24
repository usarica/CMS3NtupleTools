import FWCore.ParameterSet.Config as cms

ca12subJetMaker = cms.EDProducer("CA12SubJetMaker",
  pfJetsInputTag                   = cms.InputTag("selectedPatJetsCA10PFCHS","","CMS3"),
  #pfJetsInputTag                   = cms.InputTag("selectedPatJetsAK4PF","","CMS3"),
  pfCandidatesTag                  = cms.InputTag("packedPFCandidates"),
  pfJetPtCut                       = cms.double(5.),
)


