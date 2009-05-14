import FWCore.ParameterSet.Config as cms

tcmetMaker = cms.EDFilter("TCMETMaker",
                          muon_tag_     = cms.InputTag("muons"),
                          tcmet_tag_    = cms.InputTag("tcMet"),
                          tcmet_vm_tag_ = cms.InputTag("muonTCMETValueMapProducer", "muCorrData")
)


