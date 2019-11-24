import FWCore.ParameterSet.Config as cms

pftauMaker = cms.EDProducer("PFTauMaker",
                            aliasPrefix = cms.untracked.string("taus_pf"),
                            pftausInputTag = cms.InputTag("slimmedTaus"),
)




