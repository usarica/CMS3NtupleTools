import FWCore.ParameterSet.Config as cms

jptMaker = cms.EDFilter("JPTMaker",
	aliasPrefix = cms.untracked.string("jpts"),
                        jptInputTag       = cms.InputTag("JetPlusTrackZSPCorJetAntiKt5")
)
