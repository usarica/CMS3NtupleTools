import FWCore.ParameterSet.Config as cms

pftauExtraMaker = cms.EDProducer("PFTauExtraMaker",
                            aliasPrefix = cms.untracked.string("taus_pf"),
                            pftausInputTag = cms.InputTag("slimmedTaus"),
)




