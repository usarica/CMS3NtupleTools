import FWCore.ParameterSet.Config as cms

bTagMaker = cms.EDProducer(
        "BTagMaker",
        AliasPrefix                           = cms.string("jets"),
	cms2CaloJetsTag                       = cms.InputTag("prunedUncorrectedCMS2Jets", "calojet"),        
        referenceCaloJetsTag                  = cms.InputTag("ak5CaloJets"),
        combinedSecondaryVertexBJetTags       = cms.InputTag("CMS2CombinedSecondaryVertexBJetTags"),
	combinedSecondaryVertexMVABJetTags    = cms.InputTag("CMS2CombinedSecondaryVertexMVABJetTags"),
#        ghostTrackBJetTags                    = cms.InputTag("CMS2ghostTrackBJetTags"),
        jetBProbabilityBJetTags               = cms.InputTag("CMS2JetBProbabilityBJetTags"),
	jetProbabilityBJetTags                = cms.InputTag("CMS2JetProbabilityBJetTags"),
        simpleSecondaryVertexHighEffBJetTags  = cms.InputTag("CMS2SimpleSecondaryVertexHighEffBJetTags"),
        simpleSecondaryVertexHighPurBJetTags  = cms.InputTag("CMS2SimpleSecondaryVertexHighPurBJetTags"),  
        softElectronByIP3dBJetTags            = cms.InputTag("CMS2SoftElectronByIP3dBJetTags"),
        softElectronByPtBJetTags              = cms.InputTag("CMS2SoftElectronByPtBJetTags"),
	softMuonBJetTags                      = cms.InputTag("CMS2SoftMuonBJetTags"),
        softMuonByIP3dBJetTags                = cms.InputTag("CMS2SoftMuonByIP3dBJetTags"),
        softMuonByPtBJetTags                  = cms.InputTag("CMS2SoftMuonByPtBJetTags"),
	trackCountingHighEffBJetTags          = cms.InputTag("CMS2TrackCountingHighEffBJetTags"),
	trackCountingHighPurBJetTags          = cms.InputTag("CMS2TrackCountingHighPurBJetTags")
)

