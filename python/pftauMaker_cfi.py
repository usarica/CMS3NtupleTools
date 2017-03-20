import FWCore.ParameterSet.Config as cms

pftauMaker = cms.EDProducer("PFTauMaker",
                            aliasPrefix = cms.untracked.string("taus_pf"),
                            # cms2PFJetsTag = cms.InputTag("prunedUncorrectedCMS2Jets", "pfjet"),
                            # referencePFJetsTag = cms.InputTag("ak5PFJets"),
                            # particleFlowTag = cms.InputTag("particleFlow"),

                            # PFTau collection
                            # pftausInputTag = cms.InputTag("hpsPFTauProducer"),
                            pftausInputTag = cms.InputTag("slimmedTaus"),
                            do_full = cms.untracked.bool(False)
)




