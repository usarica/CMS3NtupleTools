import FWCore.ParameterSet.Config as cms

hypDilepMaker = cms.EDProducer("HypDilepMaker",
    aliasPrefix       = cms.untracked.string("hyp"),

    metInputTag       = cms.InputTag("pfmetMaker"),
    # muon to gen particle association
    muToGenInputTag   = cms.InputTag("muToGenAssMaker"),
    # electrons collection
    electronsInputTag = cms.InputTag("electronMaker"),
    hypJetMaxEtaCut   = cms.double(5.0),
    #jet collection
    jetsInputTag      = cms.InputTag("pfJetMaker","pfjetsp4"),
    hypJetMinPtCut    = cms.double(30.0), ##this is a corrected pt cut!

    #tight lepton pt cut
    TightLepton_PtCut = cms.double(20.0),
    #loose lepton pt cut
    LooseLepton_PtCut = cms.double(10.0),
    # muons collection
    muonsInputTag     = cms.InputTag("muonMaker"),
    useSTAMuons       = cms.bool(False)
)


