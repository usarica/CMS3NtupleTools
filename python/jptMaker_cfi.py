import FWCore.ParameterSet.Config as cms

from JetMETCorrections.Configuration.DefaultJEC_cff import *

jptMaker = cms.EDProducer("JPTMaker",
	aliasPrefix = cms.untracked.string("jpts"),
        jptInputTag = cms.InputTag("prunedUncorrectedCMS2Jets", "jpt"),
        JPTCorrectorL2L3 = cms.string("ak5JPTL2L3")                        
)
