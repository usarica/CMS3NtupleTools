import FWCore.ParameterSet.Config as cms

PFClusterMaker = cms.EDProducer("PFClusterMaker",
                                aliasPrefix = cms.untracked.string("pfc"),
                                PFClustersECAL = cms.InputTag("particleFlowClusterECAL"),
                                PFClustersHCAL = cms.InputTag("particleFlowClusterHCAL"),
                                PFClustersHFEM = cms.InputTag("particleFlowClusterHFEM"),
                                PFClustersHFHAD = cms.InputTag("particleFlowClusterHFHAD")
                                )

