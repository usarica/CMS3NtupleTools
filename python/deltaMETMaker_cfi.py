import FWCore.ParameterSet.Config as cms

deltaMETMaker = cms.EDProducer("DeltaMETMaker",
                             aliasPrefix          = cms.untracked.string(""),
                             cms2_metInputTag_    = cms.InputTag("metMaker", "evtmet"   ),
                             cms2_metphiInputTag_ = cms.InputTag("metMaker", "evtmetPhi"),
                             cms2_sumetInputTag_  = cms.InputTag("metMaker", "evtsumet" ),
                             metInputTag_ = cms.InputTag("")
)                                                              
