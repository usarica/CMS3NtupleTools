import FWCore.ParameterSet.Config as cms

triggerFilter = cms.EDFilter("TriggerFilter",
	
# 	8E29
	menu = cms.string("HLT8E29"),
#	1E31
#        menu = cms.string("HLT"),   
	filterExpressions = cms.vstring("HLT_Ele10_LW_L1R",
					"HLT_Mu9"),
	hlt = cms.bool(True)
)
