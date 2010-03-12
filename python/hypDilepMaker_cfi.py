import FWCore.ParameterSet.Config as cms

hypDilepMaker = cms.EDFilter("HypDilepMaker",
	aliasPrefix = cms.untracked.string("hyp"),
    #are we using pat jets in the JetMaker?
    #need to know if the jets were corrected
    usingCorrectedJets = cms.bool(True),
    # met collection
    tcmetInputTag = cms.InputTag("tcmetMaker"),
    metInputTag = cms.InputTag("metMaker"),
    # muon to gen particle association
    muToGenInputTag = cms.InputTag("muToGenAssMaker"),
    # electrons collection
    electronsInputTag = cms.InputTag("electronMaker"),
    hypJetMaxEtaCut = cms.double(2.4),
    #jet collection
    jetsInputTag = cms.InputTag("jptMaker"),
    hypJetMinPtCut = cms.double(30.0), ##this is a corrected pt cut!

    #tight lepton pt cut
    TightLepton_PtCut = cms.double(20.0),
    #tracks collection
    trksInputTag = cms.InputTag("trackMaker"),
    #loose lepton pt cut
    LooseLepton_PtCut = cms.double(10.0),
    # muons collection
    muonsInputTag = cms.InputTag("muonMaker")
)


