import FWCore.ParameterSet.Config as cms

bTagMaker = cms.EDFilter(
	"BTagMaker",
        AliasPrefix                        = cms.string("jets"),
	uncorRecoJetsTag                   = cms.InputTag("prunedUncorrectedCMS2Jets"),
	combinedSecondaryVertexBJetTags    = cms.InputTag("CMS2CombinedSecondaryVertexBJetTags"),
	combinedSecondaryVertexMVABJetTags = cms.InputTag("CMS2CombinedSecondaryVertexMVABJetTags"),
        jetBProbabilityBJetTags            = cms.InputTag("CMS2JetBProbabilityBJetTags"),
	jetProbabilityBJetTags             = cms.InputTag("CMS2JetProbabilityBJetTags"),
	simpleSecondaryVertexBJetTags      = cms.InputTag("CMS2SimpleSecondaryVertexBJetTags"),
        softElectronByIP3dBJetTags         = cms.InputTag("CMS2SoftElectronByIP3dBJetTags"),
        softElectronByPtBJetTags           = cms.InputTag("CMS2SoftElectronByPtBJetTags"),
	softMuonBJetTags                   = cms.InputTag("CMS2SoftMuonBJetTags"),
        softMuonByIP3dBJetTags             = cms.InputTag("CMS2SoftMuonByIP3dBJetTags"),
        softMuonByPtBJetTags               = cms.InputTag("CMS2SoftMuonByPtBJetTags"),
	softMuonNoIPBJetTags               = cms.InputTag("CMS2SoftMuonNoIPBJetTags"),
	trackCountingHighEffBJetTags       = cms.InputTag("CMS2TrackCountingHighEffBJetTags"),
	trackCountingHighPurBJetTags       = cms.InputTag("CMS2TrackCountingHighPurBJetTags")
)

