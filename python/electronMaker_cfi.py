import FWCore.ParameterSet.Config as cms

electronMaker = cms.EDFilter(
	"ElectronMaker",
	aliasPrefix = cms.untracked.string("els"),
    # Electron collection
    electronsInputTag = cms.InputTag("gsfElectrons"),
    # Beamspot
    beamSpotInputTag  = cms.InputTag("beamSpotMaker"),
    # reco Track collection
    trksInputTag      = cms.InputTag("generalTracks"),
    # egamma ID
    eidRobustLooseTag = cms.InputTag("eidRobustLooseCMS2"),
    eidRobustTightTag = cms.InputTag("eidRobustTightCMS2"),
    eidRobustHighEnergyTag = cms.InputTag("eidRobustHighEnergyCMS2"),
    eidLooseTag = cms.InputTag("eidLooseCMS2"),
    eidTightTag = cms.InputTag("eidTightCMS2"),
    cms2scsseeddetidInputTag = cms.InputTag("scMaker"),
    #conversion stuff    
    minAbsDist  = cms.double(0.02),        
    minAbsDcot  = cms.double(0.02),
    minSharedFractionOfHits = cms.double(0.45)    
)

