import FWCore.ParameterSet.Config as cms

hypDilepMaker = cms.EDFilter("HypDilepMaker",
    #are we using pat jets in the JetMaker?
    #matters for JetCorrection
    usingPATJets = cms.bool(True),
    #Candidate to Generator association maker
    candToGenAssTag = cms.InputTag("candToGenAssMaker"),
    # met collection
    metInputTag = cms.InputTag("metMaker"),
    # muon to gen particle association
    muToGenInputTag = cms.InputTag("muToGenAssMaker"),
    # electrons collection
    electronsInputTag = cms.InputTag("electronMaker"),
    hypJetMaxEtaCut = cms.double(3.0),
    #jet collection
    jetsInputTag = cms.InputTag("jetMaker"),
    hypJetMinPtCut = cms.double(15.0), ##this is an uncorrected pt cut!

    #tight lepton pt cut
    TightLepton_PtCut = cms.double(20.0),
    #tracks collection
    trksInputTag = cms.InputTag("trackMaker"),
    #cuts on the hyp jets
    hypJetMinEtaCut = cms.double(-3.0),
    #patjets collection
    patJetsInputTag = cms.InputTag("patJetMaker"),
    #loose lepton pt cut
    LooseLepton_PtCut = cms.double(5.0),
    # muons collection
    muonsInputTag = cms.InputTag("muonMaker")
)


