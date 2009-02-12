import FWCore.ParameterSet.Config as cms

jptMaker = cms.EDFilter("JPTMaker",
    # jpt collection
    jptInputTag = cms.InputTag("JetPlusTrackZSPCorJetIcone5"),

    # ic5jet collection
#    ic5jetInputTag = cms.InputTag("iterativeCone5CaloJets")
    ic5jetInputTag = cms.InputTag("sisCone5CaloJets")
)


