import FWCore.ParameterSet.Config as cms


electronMaker = cms.EDProducer(   
    "ElectronMaker",
    aliasPrefix = cms.untracked.string("els"),
    # Electron collection
    electronsInputTag   = cms.InputTag("slimmedElectrons"),
    electronVetoIdMap   = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto"),
    electronLooseIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-loose"),
    electronMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-medium"),
    electronTightIdMap  = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-tight"),
    # Beamspot
    beamSpotInputTag  = cms.InputTag("beamSpotMaker"),
    # reco Track collection
    trksInputTag      = cms.InputTag("generalTracks"),
    gsftracksInputTag = cms.InputTag("electronGsfTracks"),
    # pfCandidate and Vertex collection
    pfCandsInputTag = cms.InputTag("packedPFCandidates"),
    vtxInputTag = cms.InputTag("offlineSlimmedPrimaryVertices"),

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

