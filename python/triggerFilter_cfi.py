import FWCore.ParameterSet.Config as cms

triggerFilter = cms.EDFilter("TriggerFilter",
                           filterExpressions = cms.vstring(	"HLT_Ele15_SW_L1R",
								"HLT_Mu9"),
                           hlt = cms.bool(True)
)
