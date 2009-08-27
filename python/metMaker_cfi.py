import FWCore.ParameterSet.Config as cms

metMaker = cms.EDFilter("METMaker",
                        met_tag_               = cms.InputTag("met"                  ),               
                        metHO_tag_             = cms.InputTag("metHO"                ),             
                        metNoHF_tag_           = cms.InputTag("metNoHF"              ),           
                        metNoHFHO_tag_         = cms.InputTag("metNoHFHO"            ),                        
                        metOpt_tag_            = cms.InputTag("metOpt"               ),            
                        metOptHO_tag_          = cms.InputTag("metOptHO"             ),          
                        metOptNoHF_tag_        = cms.InputTag("metOptNoHF"           ),        
                        metOptNoHFHO_tag_      = cms.InputTag("metOptNoHFHO"         ),     
                        corMetGlobalMuons_tag_ = cms.InputTag("corMetGlobalMuons"    ),
                        MuonJEScorMET_tag_     = cms.InputTag("metMuonJESCorSC5CMS2" ),
                        muon_tag_              = cms.InputTag("muons"                ),
                        muon_vm_tag_           = cms.InputTag("muonMETValueMapProducer", "muCorrData")
)                                                              


