import FWCore.ParameterSet.Config as cms

weightinfoMaker = cms.EDProducer( "WeightInfoMaker",
                                  LHEEventInputTag = cms.string("source"),                               
                                  aliasPrefix = cms.untracked.string("weightinfo")
                                  )


