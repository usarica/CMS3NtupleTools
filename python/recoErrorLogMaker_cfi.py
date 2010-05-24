import FWCore.ParameterSet.Config as cms

recoErrorLogMaker = cms.EDFilter("RECOErrorLogMaker",
                              errorSummaryCollInputTag = cms.InputTag("logErrorHarvester"),
                              #supported severity is warning, error
                              #recommended is error to keep the event size down
                              minSeverity             = cms.string("error"), 
)
