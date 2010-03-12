import FWCore.ParameterSet.Config as cms

bcMaker = cms.EDFilter("BCMaker",
	aliasPrefix = cms.untracked.string("bcs"),
    # sc collection for EE and EB

#    scInputTag_EE = cms.InputTag("multi5x5BasicClusters", "multi5x5EndcapBasicClusters"),
    scInputTag_EE = cms.InputTag("islandBasicClusters", "islandEndcapBasicClusters"),

#    scInputTag_EE = cms.InputTag("multi5x5BasicClusters", "multi5x5EndcapBasicClusters", "testRecoEcal"),
#    scInputTag_EB = cms.InputTag("multi5x5BasicClusters", "multi5x5BarrelBasicClusters", "testRecoEcal"),
#    scInputTag_EB = cms.InputTag("hybridSuperClusters", "hybridBarrelBasicClusters"),
    scInputTag_EB = cms.InputTag("islandBasicClusters", "islandBarrelBasicClusters"),

    ecalRecHitsInputTag_EE = cms.InputTag("ecalRecHitFit","EcalRecHitsEE"),
    ecalRecHitsInputTag_EB = cms.InputTag("ecalRecHitFit","EcalRecHitsEB"),
    primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices"),
    scEtMin = cms.double(0.0)

)

