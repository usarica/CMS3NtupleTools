import FWCore.ParameterSet.Config as cms

scMaker = cms.EDFilter("SCMaker",
    # sc collection for EE and EB
    scInputTag_EE = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
    scInputTag_EB = cms.InputTag("correctedHybridSuperClusters")
)

