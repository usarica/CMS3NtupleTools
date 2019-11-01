import FWCore.ParameterSet.Config as cms


electronMaker = cms.EDProducer(
    "ElectronMaker",
    aliasPrefix = cms.untracked.string("els"),
    year = cms.int32(-1), # Must be overriden by main_pset

    # Electron collection
    electronsInputTag   = cms.InputTag("slimmedElectrons"),
    #electronVIDFall17V2NoIsoMvaValueMap   = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2Values"),
    #electronVIDFall17V2NoIsoMvaCatMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV2Categories"),
    #electronVIDFall17V2IsoMvaValueMap   = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV2Values"),
    #electronVIDFall17V2IsoMvaCatMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV2Categories"),

    miniIsoChgValueMap     = cms.InputTag("isoForEle:miniIsoChg"),
    miniIsoAllValueMap     = cms.InputTag("isoForEle:miniIsoAll"),

    # Beamspot
    beamSpotInputTag  = cms.InputTag("beamSpotMaker","evtbsp4"),
    #beamSpotTag = cms.InputTag("offlineBeamSpot"),
    # reco Track collection
    trksInputTag      = cms.InputTag("generalTracks"),
    gsftracksInputTag = cms.InputTag("electronGsfTracks"),
    # pfCandidate and Vertex collection
    pfJetsInputTag = cms.InputTag("slimmedJets"),
    #pfCandsInputTag = cms.InputTag("packedPFCandidates"),
    vtxInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),

    bFieldInputTag = cms.InputTag("eventMaker", "evtbField"),

    # isolations from external
    # pfIsoCharged03InputTag = cms.InputTag("elPFIsoValueCharged03PFIdPFIso"),
    # pfIsoGamma03InputTag   = cms.InputTag("elPFIsoValueGamma03PFIdPFIso"),
    # pfIsoNeutral03InputTag = cms.InputTag("elPFIsoValueNeutral03PFIdPFIso"),
    # pfIsoCharged04InputTag = cms.InputTag("elPFIsoValueCharged04PFIdPFIso"),
    # pfIsoGamma04InputTag   = cms.InputTag("elPFIsoValueGamma04PFIdPFIso"),
    # pfIsoNeutral04InputTag = cms.InputTag("elPFIsoValueNeutral04PFIdPFIso"),

    # reco conversions
    recoConversionInputTag = cms.InputTag("reducedEgamma:reducedConversions"),
    # egamma ID
    eidLHTag = cms.InputTag("egammaIDLikelihood"),
    # cms2scsseeddetidInputTag = cms.InputTag("scMaker"),

    #conversion stuff
    minAbsDist  = cms.double(0.02),
    minAbsDcot  = cms.double(0.02),
    minSharedFractionOfHits = cms.double(0.45),

    rhoInputTag = cms.InputTag("fastJetMaker", "evtrho"),

    ebReducedRecHitCollectionTag = cms.InputTag("reducedEgamma:reducedEBRecHits"),
    eeReducedRecHitCollectionTag = cms.InputTag("reducedEgamma:reducedEERecHits"),
    esReducedRecHitCollectionTag = cms.InputTag("reducedEgamma:reducedESRecHits"),

)

