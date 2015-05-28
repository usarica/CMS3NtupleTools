import FWCore.ParameterSet.Config as cms

ca12subJetMaker = cms.EDProducer("CA12SubJetMaker",
  srcJet                           = cms.InputTag('selectedPatJetsCA10PFCHS',"","CMS3"),
  pfJetsInputTag                   = cms.InputTag("selectedPatJetsCA10PFCHS","","CMS3"),
  pfCandidatesTag                  = cms.InputTag("packedPFCandidates"),
  pfJetPtCut                       = cms.double(5.),
)


