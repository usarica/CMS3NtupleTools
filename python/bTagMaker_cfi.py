import FWCore.ParameterSet.Config as cms

bTagMaker = cms.EDFilter(
	"BTagMaker",
        AlgoSuffix                         = cms.string("Calo"),
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


bTagTrkMaker = cms.EDFilter(
        "BTagMaker",
        AlgoSuffix                         = cms.string("Trk"),
        uncorRecoJetsTag                   = cms.InputTag("SISCone5TrkJets"),
        combinedSecondaryVertexBJetTags    = cms.InputTag("CMS2TrkCombinedSecondaryVertexBJetTags"),
        combinedSecondaryVertexMVABJetTags = cms.InputTag("CMS2TrkCombinedSecondaryVertexMVABJetTags"),
        impactParameterMVABJetTags         = cms.InputTag("CMS2TrkImpactParameterMVABJetTags"),
        jetBProbabilityBJetTags            = cms.InputTag("CMS2TrkJetBProbabilityBJetTags"),
        jetProbabilityBJetTags             = cms.InputTag("CMS2TrkJetProbabilityBJetTags"),
        simpleSecondaryVertexBJetTags      = cms.InputTag("CMS2TrkSimpleSecondaryVertexBJetTags"),
        softElectronBJetTags               = cms.InputTag("CMS2TrkSoftElectronBJetTags"),
        softMuonBJetTags                   = cms.InputTag("CMS2TrkSoftMuonBJetTags"),
        softMuonNoIPBJetTags               = cms.InputTag("CMS2TrkSoftMuonNoIPBJetTags"),
        trackCountingHighEffBJetTags       = cms.InputTag("CMS2TrkTrackCountingHighEffBJetTags"),
        trackCountingHighPurBJetTags       = cms.InputTag("CMS2TrkTrackCountingHighPurBJetTags"),
)

