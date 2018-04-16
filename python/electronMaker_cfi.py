import FWCore.ParameterSet.Config as cms
import CMS3.NtupleMaker.configProcessName as configProcessName


electronMaker = cms.EDProducer(   
    "ElectronMaker",
    aliasPrefix = cms.untracked.string("els"),
    # Electron collection
    electronsInputTag   = cms.InputTag("slimmedElectrons"),
    electronVetoIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
    electronLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    electronTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
    electronHEEPIdMap                = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV60"),
    electronVIDNonTrigMvaWP80IdMap   = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
    electronVIDNonTrigMvaWP90IdMap   = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
    electronVIDTrigMvaWP80IdMap      = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"),
    electronVIDTrigMvaWP90IdMap      = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
    electronVIDNonTrigMvaValueMap    = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
    electronVIDTrigMvaValueMap       = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
    electronVIDNonTrigMvaCatMap      = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
    electronVIDTrigMvaCatMap         = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
    electronVIDSpring16GPMvaValueMap    = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
    electronVIDSpring16GPMvaCatMap      = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),
    electronVIDSpring16HZZMvaValueMap   = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Values"),
    electronVIDSpring16HZZMvaCatMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Categories"),
    electronVIDFall17NoIsoMvaValueMap   = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
    electronVIDFall17NoIsoMvaCatMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17NoIsoV1Categories"),
    electronVIDFall17IsoMvaValueMap   = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Values"),
    electronVIDFall17IsoMvaCatMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Categories"),

    miniIsoChgValueMap     = cms.InputTag("isoForEle:miniIsoChg"),
    miniIsoAllValueMap     = cms.InputTag("isoForEle:miniIsoAll"),
    
    # Beamspot
    beamSpotInputTag  = cms.InputTag("beamSpotMaker","evtbsp4"),
    # reco Track collection
    trksInputTag      = cms.InputTag("generalTracks"),
    gsftracksInputTag = cms.InputTag("electronGsfTracks"),
    # pfCandidate and Vertex collection
    pfJetsInputTag = cms.InputTag("slimmedJets","",configProcessName.name),
    pfCandsInputTag = cms.InputTag("packedPFCandidates"),
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
    cms2scsseeddetidInputTag = cms.InputTag("scMaker"),
    #conversion stuff    
    minAbsDist  = cms.double(0.02),        
    minAbsDcot  = cms.double(0.02),
    minSharedFractionOfHits = cms.double(0.45),
    rhoInputTag = cms.InputTag("fastJetMaker", "evtrho"),
    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    
    ebReducedRecHitCollectionTag = cms.InputTag("reducedEgamma:reducedEBRecHits"),
    eeReducedRecHitCollectionTag = cms.InputTag("reducedEgamma:reducedEERecHits"),
    esReducedRecHitCollectionTag = cms.InputTag("reducedEgamma:reducedESRecHits"),

)

