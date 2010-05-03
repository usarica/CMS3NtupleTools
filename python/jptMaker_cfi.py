import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.DefaultJEC_cff import *

jptMaker = cms.EDFilter("JPTMaker",
	aliasPrefix = cms.untracked.string("jpts"),
        jptInputTag = cms.InputTag("JetPlusTrackZSPCorJetAntiKt5"),
        minUncorPt  = cms.double("5."),
        JPTCorrectorL2L3 = cms.string("ak5JPTL2L3")                        
)
