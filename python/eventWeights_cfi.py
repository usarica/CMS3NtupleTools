import FWCore.ParameterSet.Config as cms

eventWeights = cms.EDFilter("EventWeights",
    # number of events in DiElectron skim ntuple
    topDilepton2Electron_nevts = cms.untracked.double(1.0),
    # flag for DiMuon skim
    IsTopDileptonMuonX = cms.untracked.bool(True),
    # CSA07 module label
    csa07InputTag = cms.InputTag("csa07InfoMaker"),
    # otherwise: provide xsec (in pb) and number of events in the sample
    IsNotSoup = cms.untracked.bool(False),
    # flag for DiElectron skim
    IsTopDilepton2Electron = cms.untracked.bool(False),
    # event module label
    eventInputTag = cms.InputTag("eventMaker"),
    # number of events in DiMuon skim ntuple
    topDileptonMuonX_nevts = cms.untracked.double(25000.0),
    # number of events in DiMuon skim data file
    Chowder_topDileptonMuonX_nevts = cms.untracked.double(725295.0),
    # number of events in the sample
    Nevts = cms.untracked.int32(1000),
    # Xsec in pb
    Xsec = cms.untracked.double(1000.0),
    # number of events in DiElectron skim data file
    Chowder_topDilepton2Electron_nevts = cms.untracked.double(639013.0)
)


