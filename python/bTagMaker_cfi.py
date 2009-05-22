import FWCore.ParameterSet.Config as cms

bTagMaker = cms.EDFilter(
	"BTagMaker",
	uncorRecoJetsTag                   = cms.InputTag("prunedUncorrectedCMS2Jets"),
	combinedSecondaryVertexBJetTags    = cms.InputTag("CMS2CombinedSecondaryVertexBJetTags"),
	combinedSecondaryVertexMVABJetTags = cms.InputTag("CMS2CombinedSecondaryVertexMVABJetTags"),
	impactParameterMVABJetTags         = cms.InputTag("CMS2ImpactParameterMVABJetTags"),
	jetBProbabilityBJetTags            = cms.InputTag("CMS2JetBProbabilityBJetTags"),
	jetProbabilityBJetTags             = cms.InputTag("CMS2JetProbabilityBJetTags"),
	simpleSecondaryVertexBJetTags      = cms.InputTag("CMS2SimpleSecondaryVertexBJetTags"),
	softElectronBJetTags               = cms.InputTag("CMS2SoftElectronBJetTags"),
	softMuonBJetTags                   = cms.InputTag("CMS2SoftMuonBJetTags"),
	softMuonNoIPBJetTags               = cms.InputTag("CMS2SoftMuonNoIPBJetTags"),
	trackCountingHighEffBJetTags       = cms.InputTag("CMS2TrackCountingHighEffBJetTags"),
	trackCountingHighPurBJetTags       = cms.InputTag("CMS2TrackCountingHighPurBJetTags"),
)

