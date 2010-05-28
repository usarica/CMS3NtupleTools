import FWCore.ParameterSet.Config as cms

bTagPFJetMaker = cms.EDFilter(
        "BTagMaker",
        AliasPrefix                           = cms.string("pfjets"),
	cms2CaloJetsTag                       = cms.InputTag("prunedUncorrectedCMS2Jets","pfjet"),        
        referenceCaloJetsTag                  = cms.InputTag("ak5PFJets"),
        combinedSecondaryVertexBJetTags       = cms.InputTag("CMS2PFCombinedSecondaryVertexBJetTags"),
	combinedSecondaryVertexMVABJetTags    = cms.InputTag("CMS2PFCombinedSecondaryVertexMVABJetTags"),
#        ghostTrackBJetTags                    = cms.InputTag("CMS2ghostTrackBJetTags"),
        jetBProbabilityBJetTags               = cms.InputTag("CMS2PFJetBProbabilityBJetTags"),
	jetProbabilityBJetTags                = cms.InputTag("CMS2PFJetProbabilityBJetTags"),
        simpleSecondaryVertexHighEffBJetTags  = cms.InputTag("CMS2PFSimpleSecondaryVertexHighEffBJetTags"),
        simpleSecondaryVertexHighPurBJetTags  = cms.InputTag("CMS2PFSimpleSecondaryVertexHighPurBJetTags"),  
        softElectronByIP3dBJetTags            = cms.InputTag("CMS2PFSoftElectronByIP3dBJetTags"),
        softElectronByPtBJetTags              = cms.InputTag("CMS2PFSoftElectronByPtBJetTags"),
	softMuonBJetTags                      = cms.InputTag("CMS2PFSoftMuonBJetTags"),
        softMuonByIP3dBJetTags                = cms.InputTag("CMS2PFSoftMuonByIP3dBJetTags"),
        softMuonByPtBJetTags                  = cms.InputTag("CMS2PFSoftMuonByPtBJetTags"),
	trackCountingHighEffBJetTags          = cms.InputTag("CMS2PFTrackCountingHighEffBJetTags"),
	trackCountingHighPurBJetTags          = cms.InputTag("CMS2PFTrackCountingHighPurBJetTags")
)

