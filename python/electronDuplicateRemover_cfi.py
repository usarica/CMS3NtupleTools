import FWCore.ParameterSet.Config as cms

uniqueElectrons = cms.EDFilter("ElectronDuplicateRemover",
    # Electron collection
    electronsInputTag = cms.InputTag("pixelMatchGsfElectrons")
)

