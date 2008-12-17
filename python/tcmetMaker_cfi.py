import FWCore.ParameterSet.Config as cms

tcmetMaker = cms.EDFilter("TCMETMaker",
     metInputTag = cms.InputTag("tcMet")
)


