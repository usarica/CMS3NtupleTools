import FWCore.ParameterSet.Config as cms

patMETMaker = cms.EDFilter("PATMETMaker",
	aliasPrefix = cms.untracked.string("met_pat"),
                           patMETsInputTag = cms.InputTag("patMETs")
)


