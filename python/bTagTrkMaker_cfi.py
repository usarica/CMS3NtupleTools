import FWCore.ParameterSet.Config as cms

bTagTrkMaker = cms.EDFilter(
        "BTagMaker",
        AliasPrefix                        = cms.string("trkjets"),
        uncorRecoJetsTag                   = cms.InputTag("ak5TrackJets"),
        combinedSecondaryVertexBJetTags    = cms.InputTag("CMS2TrkCombinedSecondaryVertexBJetTags"),
	combinedSecondaryVertexMVABJetTags = cms.InputTag("CMS2TrkCombinedSecondaryVertexMVABJetTags"),
        jetBProbabilityBJetTags            = cms.InputTag("CMS2TrkJetBProbabilityBJetTags"),
	jetProbabilityBJetTags             = cms.InputTag("CMS2TrkJetProbabilityBJetTags"),
	simpleSecondaryVertexBJetTags      = cms.InputTag("CMS2TrkSimpleSecondaryVertexBJetTags"),
        softElectronByIP3dBJetTags         = cms.InputTag("CMS2TrkSoftElectronByIP3dBJetTags"),
        softElectronByPtBJetTags           = cms.InputTag("CMS2TrkSoftElectronByPtBJetTags"),
	softMuonBJetTags                   = cms.InputTag("CMS2TrkSoftMuonBJetTags"),
        softMuonByIP3dBJetTags             = cms.InputTag("CMS2TrkSoftMuonByIP3dBJetTags"),
        softMuonByPtBJetTags               = cms.InputTag("CMS2TrkSoftMuonByPtBJetTags"),
	softMuonNoIPBJetTags               = cms.InputTag("CMS2TrkSoftMuonNoIPBJetTags"),
	trackCountingHighEffBJetTags       = cms.InputTag("CMS2TrkTrackCountingHighEffBJetTags"),
	trackCountingHighPurBJetTags       = cms.InputTag("CMS2TrkTrackCountingHighPurBJetTags")
)

