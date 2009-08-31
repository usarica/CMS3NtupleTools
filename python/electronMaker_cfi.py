import FWCore.ParameterSet.Config as cms

electronMaker = cms.EDFilter(
	"ElectronMaker",
    # Electron collection
    electronsInputTag = cms.InputTag("gsfElectrons"),
    # Beamspot
    beamSpotInputTag = cms.InputTag("beamSpotMaker"),
    # egamma ID
    eidRobustLooseTag = cms.InputTag("eidRobustLooseCMS2"),
    eidRobustTightTag = cms.InputTag("eidRobustTightCMS2"),
    eidRobustHighEnergyTag = cms.InputTag("eidRobustHighEnergyCMS2"),
    eidLooseTag = cms.InputTag("eidLooseCMS2"),
    eidTightTag = cms.InputTag("eidTightCMS2")

)

