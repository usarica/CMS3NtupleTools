import FWCore.ParameterSet.Config as cms

hypTrilepMaker = cms.EDFilter("HypTrilepMaker",
    # met collection
    metInputTag = cms.InputTag("metMaker"),
    hypJetMinPtCut = cms.double(15.0), ##this is an uncorrected pt cut!

    # electrons collection
    electronsInputTag = cms.InputTag("electronMaker"),
    hypJetMaxEtaCut = cms.double(3.0),
    #jet collection
    jetsInputTag = cms.InputTag("jetMaker"),
    # muons collection
    muonsInputTag = cms.InputTag("muonMaker"),
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
    #are we using pat jets in the JetMaker?
    #matters for JetCorrection, pat jets are corrected, hyp jets cut is on uncorrected jets
    usingPATJets = cms.bool(True)
)


