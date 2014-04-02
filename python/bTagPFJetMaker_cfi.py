import FWCore.ParameterSet.Config as cms

bTagPFJetMaker = cms.EDProducer(
        "BTagMaker",
        AliasPrefix                           = cms.string("pfjets"),
	cms2CaloJetsTag                       = cms.InputTag("prunedUncorrectedCMS2Jets","pfjet"),        
        referenceCaloJetsTag                  = cms.InputTag("ak5PFJets"),
        combinedSecondaryVertexBJetTags       = cms.InputTag("CMS2PFCombinedSecondaryVertexBJetTags"),
	combinedSecondaryVertexMVABJetTags    = cms.InputTag("CMS2PFCombinedSecondaryVertexMVABJetTags"),
#        ghostTrackBJetTags                    = cms.InputTag("CMS2ghostTrackBJetTags"),
        jetBProbabilityBJetTags               = cms.InputTag("CMS2PFJetBProbabilityBJetTags"),
	jetProbabilityBJetTags                = cms.InputTag("CMS2PFJetProbabilityBJetTags"),
        simpleSecondaryVertexBJetTags         = cms.InputTag("CMS2PFSimpleSecondaryVertexBJetTags"),
        simpleSecondaryVertexHighEffBJetTags  = cms.InputTag("CMS2PFSimpleSecondaryVertexHighEffBJetTags"),
        simpleSecondaryVertexHighPurBJetTags  = cms.InputTag("CMS2PFSimpleSecondaryVertexHighPurBJetTags"),  
        softElectronTags                      = cms.InputTag("CMS2PFSoftElectronBJetTags"),
        softElectronByIP3dBJetTags            = cms.InputTag("CMS2PFSoftElectronByIP3dBJetTags"),
        softElectronByPtBJetTags              = cms.InputTag("CMS2PFSoftElectronByPtBJetTags"),
	softMuonBJetTags                      = cms.InputTag("CMS2PFSoftMuonBJetTags"),
        softMuonByIP3dBJetTags                = cms.InputTag("CMS2PFSoftMuonByIP3dBJetTags"),
        softMuonByPtBJetTags                  = cms.InputTag("CMS2PFSoftMuonByPtBJetTags"),
	trackCountingHighEffBJetTags          = cms.InputTag("CMS2PFTrackCountingHighEffBJetTags"),
	trackCountingHighPurBJetTags          = cms.InputTag("CMS2PFTrackCountingHighPurBJetTags")
)


bTagPFCHSJetMaker = cms.EDProducer(
        "BTagMaker",
        AliasPrefix                           = cms.string("pfchsjets"),
	cms2CaloJetsTag                       = cms.InputTag("prunedUncorrectedCMS2Jets","pfchsjet"),        
        referenceCaloJetsTag                  = cms.InputTag("ak5PFJetsCHS"),
        combinedSecondaryVertexBJetTags       = cms.InputTag("CMS2PFCHSCombinedSecondaryVertexBJetTags"),
	combinedSecondaryVertexMVABJetTags    = cms.InputTag("CMS2PFCHSCombinedSecondaryVertexMVABJetTags"),
#        ghostTrackBJetTags                    = cms.InputTag("CMS2ghostTrackBJetTags"),
        jetBProbabilityBJetTags               = cms.InputTag("CMS2PFCHSJetBProbabilityBJetTags"),
	jetProbabilityBJetTags                = cms.InputTag("CMS2PFCHSJetProbabilityBJetTags"),
        simpleSecondaryVertexBJetTags         = cms.InputTag("CMS2PFCHSSimpleSecondaryVertexBJetTags"),
        simpleSecondaryVertexHighEffBJetTags  = cms.InputTag("CMS2PFCHSSimpleSecondaryVertexHighEffBJetTags"),
        simpleSecondaryVertexHighPurBJetTags  = cms.InputTag("CMS2PFCHSSimpleSecondaryVertexHighPurBJetTags"),  
        softElectronTags                      = cms.InputTag("CMS2PFCHSSoftElectronBJetTags"),
        softElectronByIP3dBJetTags            = cms.InputTag("CMS2PFCHSSoftElectronByIP3dBJetTags"),
        softElectronByPtBJetTags              = cms.InputTag("CMS2PFCHSSoftElectronByPtBJetTags"),
	softMuonBJetTags                      = cms.InputTag("CMS2PFCHSSoftMuonBJetTags"),
        softMuonByIP3dBJetTags                = cms.InputTag("CMS2PFCHSSoftMuonByIP3dBJetTags"),
        softMuonByPtBJetTags                  = cms.InputTag("CMS2PFCHSSoftMuonByPtBJetTags"),
	trackCountingHighEffBJetTags          = cms.InputTag("CMS2PFCHSTrackCountingHighEffBJetTags"),
	trackCountingHighPurBJetTags          = cms.InputTag("CMS2PFCHSTrackCountingHighPurBJetTags")
)

