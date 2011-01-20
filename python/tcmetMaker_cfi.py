import FWCore.ParameterSet.Config as cms

tcmetMaker = cms.EDProducer("TCMETMaker",
                            aliasPrefix = cms.untracked.string("evt"),
                            muon_tag_     = cms.InputTag("muons"),
                            tcmet_tag_    = cms.InputTag("tcMet"),
                            pftcmet_tag_  = cms.InputTag("tcMetWithPFclusters"),
                            tcmet_vm_tag_ = cms.InputTag("muonTCMETValueMapProducer", "muCorrData")
)


