import FWCore.ParameterSet.Config as cms

# sequence comes from:
# https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookBTagging#DifferentJets


# b-tagging general configuration
from RecoJets.JetAssociationProducers.ic5PFJetTracksAssociatorAtVertex_cfi import *
from RecoBTag.Configuration.RecoBTag_cff import *
from RecoBTag.SoftLepton.softElectronCandProducer_cfi import *



#The first step is to clone the definition of the association between jets and tracks. This is where most of the customization will take place. As an example, here generalTracks and sisCone5CaloJets are used - to try different ones, change them to what you need throughout the configuration (check later the soft lepton part, too):
# create a CMS2 jets and tracks association
CMS2PFJetTracksAssociatorAtVertex = ic5PFJetTracksAssociatorAtVertex.clone()
CMS2PFJetTracksAssociatorAtVertex.jets = "ak5PFJets"
CMS2PFJetTracksAssociatorAtVertex.tracks = "generalTracks"
#The one needs to clone the b-tag producers and instruct them to use this CMS2 collection. First the impact parameter-based b-tag:
# impact parameter b-tag
CMS2PFImpactParameterTagInfos = impactParameterTagInfos.clone()
CMS2PFImpactParameterTagInfos.jetTracks = "CMS2PFJetTracksAssociatorAtVertex"
CMS2PFTrackCountingHighEffBJetTags = trackCountingHighEffBJetTags.clone()
CMS2PFTrackCountingHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFImpactParameterTagInfos") )
CMS2PFTrackCountingHighPurBJetTags = trackCountingHighPurBJetTags.clone()
CMS2PFTrackCountingHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFImpactParameterTagInfos") )
CMS2PFJetProbabilityBJetTags = jetProbabilityBJetTags.clone()
CMS2PFJetProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFImpactParameterTagInfos") )
CMS2PFJetBProbabilityBJetTags = jetBProbabilityBJetTags.clone()
CMS2PFJetBProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFImpactParameterTagInfos") )
#Then the secondary vertex-based b-tag. Note that these producers inherit the jets and tracks they use from the impact parameter modules: # secondary vertex b-tag
CMS2PFSecondaryVertexTagInfos = secondaryVertexTagInfos.clone()
CMS2PFSecondaryVertexTagInfos.trackIPTagInfos = cms.InputTag("CMS2PFImpactParameterTagInfos")
#CMS2PFSimpleSecondaryVertexBJetTags = simpleSecondaryVertexBJetTags.clone()
#CMS2PFSimpleSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFSecondaryVertexTagInfos") )
CMS2PFSimpleSecondaryVertexHighEffBJetTags = simpleSecondaryVertexHighEffBJetTags.clone()
CMS2PFSimpleSecondaryVertexHighEffBJetTags.tagInfos = cms.VInputTag(cms.InputTag("CMS2PFSecondaryVertexTagInfos"))
CMS2PFSimpleSecondaryVertexHighPurBJetTags = simpleSecondaryVertexHighPurBJetTags.clone()
CMS2PFSimpleSecondaryVertexHighPurBJetTags.tagInfos = cms.VInputTag(cms.InputTag("CMS2PFSecondaryVertexTagInfos"))
CMS2PFCombinedSecondaryVertexBJetTags = combinedSecondaryVertexBJetTags.clone()
CMS2PFCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFImpactParameterTagInfos"), cms.InputTag("CMS2PFSecondaryVertexTagInfos") )
CMS2PFCombinedSecondaryVertexMVABJetTags = combinedSecondaryVertexMVABJetTags.clone()
CMS2PFCombinedSecondaryVertexMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFImpactParameterTagInfos"), cms.InputTag("CMS2PFSecondaryVertexTagInfos") )
# add some ghost b-tagging stuff
#CMS2PFghostVertexTagInfos = ghostTrackVertexTagInfos.clone()
#CMS2PFghostVertexTagInfos.trackIPTagInfos = "CMS2PFImpactParameterTagInfos"
#CMS2PFghostTrackBJetTags = ghostTrackBJetTags.clone()

#CMS2PFghostTrackBJetTags.tagInfos = cms.VInputTag(cms.InputTag("CMS2PFImpactParameterTagInfos"),
#                                                cms.InputTag("CMS2PFghostVertexTagInfos"))
#And the soft lepton b-tag. These producers will accept as input either the raw jets, or the association collection:
# soft electron b-tag
CMS2PFSoftElectronTagInfos = softElectronTagInfos.clone()
CMS2PFSoftElectronTagInfos.jets = "ak5PFJets"
#CMS2PFSoftElectronBJetTags = softElectronBJetTags.clone()
#CMS2PFSoftElectronBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFSoftElectronTagInfos") )
CMS2PFSoftElectronByIP3dBJetTags = softElectronByIP3dBJetTags.clone()
CMS2PFSoftElectronByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFSoftElectronTagInfos") )
CMS2PFSoftElectronByPtBJetTags = softElectronByPtBJetTags.clone()
CMS2PFSoftElectronByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFSoftElectronTagInfos") )


# soft muon b-tag
CMS2PFSoftMuonTagInfos = softMuonTagInfos.clone()
CMS2PFSoftMuonTagInfos.jets = "ak5PFJets"
CMS2PFSoftMuonBJetTags = softMuonBJetTags.clone()
CMS2PFSoftMuonBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFSoftMuonTagInfos") )
CMS2PFSoftMuonByIP3dBJetTags = softMuonByIP3dBJetTags.clone()
CMS2PFSoftMuonByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFSoftMuonTagInfos") )
CMS2PFSoftMuonByPtBJetTags = softMuonByPtBJetTags.clone()
CMS2PFSoftMuonByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFSoftMuonTagInfos") )
#Finally, there needs to be a CMS2 path running all these modules
# prepare a path running the CMS2 modules
CMS2PFJetTracksAssociator = cms.Sequence(
    CMS2PFJetTracksAssociatorAtVertex
)

CMS2PFJetBtaggingIP = cms.Sequence(
    CMS2PFImpactParameterTagInfos * (
        CMS2PFTrackCountingHighEffBJetTags +
        CMS2PFTrackCountingHighPurBJetTags +
        CMS2PFJetProbabilityBJetTags +
        CMS2PFJetBProbabilityBJetTags
    )
)

CMS2PFJetBtaggingSV = cms.Sequence(
    CMS2PFImpactParameterTagInfos *
    CMS2PFSecondaryVertexTagInfos * (
#        CMS2PFSimpleSecondaryVertexBJetTags +
        CMS2PFSimpleSecondaryVertexHighEffBJetTags +
        CMS2PFSimpleSecondaryVertexHighPurBJetTags +
        CMS2PFCombinedSecondaryVertexBJetTags +
        CMS2PFCombinedSecondaryVertexMVABJetTags
    )
)

#CMS2PFJetghostBTagging = cms.Sequence(
#    CMS2PFghostVertexTagInfos *
#    CMS2PFghostTrackBJetTags
#)

CMS2PFJetBtaggingEle = cms.Sequence(
    #btagSoftElectrons *
    softElectronCands *
    CMS2PFSoftElectronTagInfos * (
    CMS2PFSoftElectronByIP3dBJetTags + 
#    CMS2PFSoftElectronBJetTags
    CMS2PFSoftElectronByPtBJetTags
    )
)

CMS2PFJetBtaggingMu = cms.Sequence(
    CMS2PFSoftMuonTagInfos * (
        CMS2PFSoftMuonBJetTags +
        CMS2PFSoftMuonByIP3dBJetTags +
        CMS2PFSoftMuonByPtBJetTags
    )
)

CMS2PFJetBtagging = cms.Sequence(
    CMS2PFJetBtaggingIP +
    CMS2PFJetBtaggingSV +
#    CMS2PFJetghostBTagging +
    CMS2PFJetBtaggingEle +
    CMS2PFJetBtaggingMu
)

CMS2PFBtagging = cms.Sequence(
    CMS2PFJetTracksAssociator *
    CMS2PFJetBtagging
)
