import FWCore.ParameterSet.Config as cms
from RecoMuon.MuonIsolationProducers.caloExtractorByAssociatorBlocks_cff import *
from RecoMuon.MuonIsolationProducers.trackExtractorBlocks_cff import *
from RecoMuon.MuonIsolationProducers.jetExtractorBlock_cff import *

randomConeIsoMaker = cms.EDProducer("RandomConeIsoMaker",
	aliasPrefix = cms.untracked.string("ran"),
    primaryVertexInputTag = cms.InputTag("offlinePrimaryVertices"),
    ecalBarrelRecHitProducer = cms.InputTag("ecalRecHit"),
    ecalBarrelRecHitCollection = cms.InputTag("EcalRecHitsEB"),
    ecalEndcapRecHitProducer = cms.InputTag("ecalRecHit"),
    ecalEndcapRecHitCollection = cms.InputTag("EcalRecHitsEE"),
    towerProducer = cms.InputTag("towerMaker"),
    BeamspotProducer = cms.InputTag("offlineBeamSpot"),
    srProducerEE = cms.InputTag("ecalDigis"),
    srProducerEB = cms.InputTag("ecalDigis"),
    trackProducer = cms.InputTag("generalTracks"),

    CaloExtractorPSet = cms.PSet(
    MIsoCaloExtractorByAssociatorTowersBlock
    ),
    TrackExtractorPSet = cms.PSet(
    MIsoTrackExtractorBlock   
    ),
    JetExtractorPSet = cms.PSet(
    MIsoJetExtractorBlock
    ),
   
 
    egammaPSet = cms.PSet(
   
    useNumCrystals = cms.bool(True),
    intRadiusBarrel = cms.double(3.0),
    intRadiusEndcap = cms.double(3.0),
    jurassicWidth = cms.double(1.5),    #dEta strip width
    extRadius = cms.double(0.3),
    etMinBarrel = cms.double(0.0),
    eMinBarrel = cms.double(0.08),
    etMinEndcap = cms.double(0.100),
    eMinEndcap = cms.double(0.0),
    
    useIsolEt = cms.bool(True),
    tryBoth   = cms.bool(True),
    subtract  = cms.bool(False),
    vetoClustered  = cms.bool(False),

    intRadius = cms.double(0.15), # to be orthogonal with the H/E ID cut
    #extRadius = cms.double(0.3),
    etMin = cms.double(0.0),
    Depth = cms.int32(-1),

    intRadius_trk  = cms.double(0.04),
    extRadius_trk  = cms.double(0.3),
    ptMin_trk      = cms.double(0.7),
    maxVtxDist_trk = cms.double(0.2),
    maxVtxDistXY_trk= cms.double(9999.0)
    
    ),

#
# for matching and propagation
#
TrackAssociatorParameters = cms.PSet(
        muonMaxDistanceSigmaX = cms.double(0.0),
        muonMaxDistanceSigmaY = cms.double(0.0),
        CSCSegmentCollectionLabel = cms.InputTag("cscSegments"),
        dRHcal = cms.double(9999.0),
        dREcal = cms.double(9999.0),
        CaloTowerCollectionLabel = cms.InputTag("towerMaker"),
        useEcal = cms.bool(True),
        dREcalPreselection = cms.double(0.05),
        HORecHitCollectionLabel = cms.InputTag("horeco"),
        dRMuon = cms.double(9999.0),
	trajectoryUncertaintyTolerance = cms.double(-1.0),
        propagateAllDirections = cms.bool(True),
        muonMaxDistanceX = cms.double(5.0),
        muonMaxDistanceY = cms.double(5.0),
        useHO = cms.bool(True),
        accountForTrajectoryChangeCalo = cms.bool(False),
        DTRecSegment4DCollectionLabel = cms.InputTag("dt4DSegments"),
        EERecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEE"),
        dRHcalPreselection = cms.double(0.2),
        useMuon = cms.bool(False),
        useCalo = cms.bool(False),
        EBRecHitCollectionLabel = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
        dRMuonPreselection = cms.double(0.2),
	usePreshower = cms.bool(False),
	dRPreshowerPreselection = cms.double(0.2),
        truthMatch = cms.bool(False),
        HBHERecHitCollectionLabel = cms.InputTag("hbhereco"),
        useHcal = cms.bool(True)
)
                                 
                                  
)

randomConeIsoMaker.CaloExtractorPSet.TrackAssociatorParameters.dREcal = 0.5
randomConeIsoMaker.CaloExtractorPSet.TrackAssociatorParameters.dRHcal = 0.5
randomConeIsoMaker.CaloExtractorPSet.TrackAssociatorParameters.dREcalPreselection = 0.5
randomConeIsoMaker.CaloExtractorPSet.TrackAssociatorParameters.dRHcalPreselection = 0.5

