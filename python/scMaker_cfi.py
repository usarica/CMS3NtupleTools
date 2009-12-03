import FWCore.ParameterSet.Config as cms

scMaker = cms.EDFilter("SCMaker",
    # sc collection for EE and EB
    scInputTag_EE = cms.InputTag("correctedMulti5x5SuperClustersWithPreshower"),
    scInputTag_EB = cms.InputTag("correctedHybridSuperClusters"),
    hcalRecHitsInputTag_HBHE = cms.InputTag("hbhereco"),
    ecalRecHitsInputTag_EE = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
    ecalRecHitsInputTag_EB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices"),
    electronsInputTag = cms.InputTag("gsfElectrons"),
    scEtMin = cms.double(10.0),
    debug = cms.bool(False),
    MCTruthCollection = cms.InputTag("generator")

)

