import FWCore.ParameterSet.Config as cms

weightMaker = cms.EDProducer( "WeightMaker",
                              LHEEventInputTag = cms.string("source"),                               
                              aliasPrefix = cms.untracked.string("weight")
                              )


