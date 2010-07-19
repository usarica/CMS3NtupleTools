import FWCore.ParameterSet.Config as cms

bTagJPTJetMaker = cms.EDProducer(
        "BTagMaker",
        AliasPrefix                           = cms.string("jpts"),
	cms2CaloJetsTag                       = cms.InputTag("prunedUncorrectedCMS2Jets", "jpt"),
        referenceCaloJetsTag                  = cms.InputTag("JetPlusTrackZSPCorJetAntiKt5"),
        combinedSecondaryVertexBJetTags       = cms.InputTag("CMS2JPTCombinedSecondaryVertexBJetTags"),
	combinedSecondaryVertexMVABJetTags    = cms.InputTag("CMS2JPTCombinedSecondaryVertexMVABJetTags"),
#        ghostTrackBJetTags                    = cms.InputTag("CMS2ghostTrackBJetTags"),
        jetBProbabilityBJetTags               = cms.InputTag("CMS2JPTJetBProbabilityBJetTags"),
	jetProbabilityBJetTags                = cms.InputTag("CMS2JPTJetProbabilityBJetTags"),
        simpleSecondaryVertexHighEffBJetTags  = cms.InputTag("CMS2JPTSimpleSecondaryVertexHighEffBJetTags"),
        simpleSecondaryVertexHighPurBJetTags  = cms.InputTag("CMS2JPTSimpleSecondaryVertexHighPurBJetTags"),  
        softElectronByIP3dBJetTags            = cms.InputTag("CMS2JPTSoftElectronByIP3dBJetTags"),
        softElectronByPtBJetTags              = cms.InputTag("CMS2JPTSoftElectronByPtBJetTags"),
	softMuonBJetTags                      = cms.InputTag("CMS2JPTSoftMuonBJetTags"),
        softMuonByIP3dBJetTags                = cms.InputTag("CMS2JPTSoftMuonByIP3dBJetTags"),
        softMuonByPtBJetTags                  = cms.InputTag("CMS2JPTSoftMuonByPtBJetTags"),
	trackCountingHighEffBJetTags          = cms.InputTag("CMS2JPTTrackCountingHighEffBJetTags"),
	trackCountingHighPurBJetTags          = cms.InputTag("CMS2JPTTrackCountingHighPurBJetTags")
)

