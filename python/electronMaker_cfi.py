import FWCore.ParameterSet.Config as cms

electronMaker = cms.EDProducer(
	"ElectronMaker",
	aliasPrefix = cms.untracked.string("els"),
    # Electron collection
    electronsInputTag = cms.InputTag("gsfElectrons"),
    # Beamspot
    beamSpotInputTag  = cms.InputTag("beamSpotMaker"),
    # reco Track collection
    trksInputTag      = cms.InputTag("generalTracks"),
    gsftracksInputTag = cms.InputTag("electronGsfTracks"),
    # pfCandidate and Vertex collection
    pfCandsInputTag = cms.InputTag("particleFlow"),
    vtxInputTag = cms.InputTag("offlinePrimaryVerticesDA"),
    # egamma ID
    eidLHTag = cms.InputTag("egammaIDLikelihood"),
    cms2scsseeddetidInputTag = cms.InputTag("scMaker"),
    #conversion stuff    
    minAbsDist  = cms.double(0.02),        
    minAbsDcot  = cms.double(0.02),
    minSharedFractionOfHits = cms.double(0.45)    
)

