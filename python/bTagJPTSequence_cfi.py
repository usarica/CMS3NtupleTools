import FWCore.ParameterSet.Config as cms

# sequence comes from:
# https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookBTagging#DifferentJets


# b-tagging general configuration
from RecoJets.JetAssociationProducers.ic5PFJetTracksAssociatorAtVertex_cfi import *
from RecoBTag.Configuration.RecoBTag_cff import *
from RecoBTag.SoftLepton.softElectronCandProducer_cfi import *



#The first step is to clone the definition of the association between jets and tracks. This is where most of the customization will take place. As an example, here generalTracks and sisCone5CaloJets are used - to try different ones, change them to what you need throughout the configuration (check later the soft lepton part, too):
# create a CMS2 jets and tracks association
CMS2JPTJetTracksAssociatorAtVertex = ic5PFJetTracksAssociatorAtVertex.clone()
CMS2JPTJetTracksAssociatorAtVertex.jets = "JetPlusTrackZSPCorJetAntiKt5"
CMS2JPTJetTracksAssociatorAtVertex.tracks = "generalTracks"
#The one needs to clone the b-tag producers and instruct them to use this CMS2 collection. First the impact parameter-based b-tag:
# impact parameter b-tag
CMS2JPTImpactParameterTagInfos = impactParameterTagInfos.clone()
CMS2JPTImpactParameterTagInfos.jetTracks = "CMS2JPTJetTracksAssociatorAtVertex"
CMS2JPTTrackCountingHighEffBJetTags = trackCountingHighEffBJetTags.clone()
CMS2JPTTrackCountingHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2JPTImpactParameterTagInfos") )
CMS2JPTTrackCountingHighPurBJetTags = trackCountingHighPurBJetTags.clone()
CMS2JPTTrackCountingHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2JPTImpactParameterTagInfos") )
CMS2JPTJetProbabilityBJetTags = jetProbabilityBJetTags.clone()
CMS2JPTJetProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2JPTImpactParameterTagInfos") )
CMS2JPTJetBProbabilityBJetTags = jetBProbabilityBJetTags.clone()
CMS2JPTJetBProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2JPTImpactParameterTagInfos") )
#Then the secondary vertex-based b-tag. Note that these producers inherit the jets and tracks they use from the impact parameter modules: # secondary vertex b-tag
CMS2JPTSecondaryVertexTagInfos = secondaryVertexTagInfos.clone()
CMS2JPTSecondaryVertexTagInfos.trackIPTagInfos = cms.InputTag("CMS2JPTImpactParameterTagInfos")
#CMS2JPTSimpleSecondaryVertexBJetTags = simpleSecondaryVertexBJetTags.clone()
#CMS2JPTSimpleSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2JPTSecondaryVertexTagInfos") )
CMS2JPTSimpleSecondaryVertexHighEffBJetTags = simpleSecondaryVertexHighEffBJetTags.clone()
CMS2JPTSimpleSecondaryVertexHighEffBJetTags.tagInfos = cms.VInputTag(cms.InputTag("CMS2JPTSecondaryVertexTagInfos"))
CMS2JPTSimpleSecondaryVertexHighPurBJetTags = simpleSecondaryVertexHighPurBJetTags.clone()
CMS2JPTSimpleSecondaryVertexHighPurBJetTags.tagInfos = cms.VInputTag(cms.InputTag("CMS2JPTSecondaryVertexTagInfos"))
CMS2JPTCombinedSecondaryVertexBJetTags = combinedSecondaryVertexBJetTags.clone()
CMS2JPTCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2JPTImpactParameterTagInfos"), cms.InputTag("CMS2JPTSecondaryVertexTagInfos") )
CMS2JPTCombinedSecondaryVertexMVABJetTags = combinedSecondaryVertexMVABJetTags.clone()
CMS2JPTCombinedSecondaryVertexMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2JPTImpactParameterTagInfos"), cms.InputTag("CMS2JPTSecondaryVertexTagInfos") )
# add some ghost b-tagging stuff
#CMS2JPTghostVertexTagInfos = ghostTrackVertexTagInfos.clone()
#CMS2JPTghostVertexTagInfos.trackIPTagInfos = "CMS2JPTImpactParameterTagInfos"
#CMS2JPTghostTrackBJetTags = ghostTrackBJetTags.clone()

#CMS2JPTghostTrackBJetTags.tagInfos = cms.VInputTag(cms.InputTag("CMS2JPTImpactParameterTagInfos"),
#                                                cms.InputTag("CMS2JPTghostVertexTagInfos"))
#And the soft lepton b-tag. These producers will accept as input either the raw jets, or the association collection:
# soft electron b-tag
CMS2JPTSoftElectronTagInfos = softElectronTagInfos.clone()
CMS2JPTSoftElectronTagInfos.jets = "JetPlusTrackZSPCorJetAntiKt5"
#CMS2JPTSoftElectronBJetTags = softElectronBJetTags.clone()
#CMS2JPTSoftElectronBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2JPTSoftElectronTagInfos") )
CMS2JPTSoftElectronByIP3dBJetTags = softElectronByIP3dBJetTags.clone()
CMS2JPTSoftElectronByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2JPTSoftElectronTagInfos") )
CMS2JPTSoftElectronByPtBJetTags = softElectronByPtBJetTags.clone()
CMS2JPTSoftElectronByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2JPTSoftElectronTagInfos") )


# soft muon b-tag
CMS2JPTSoftMuonTagInfos = softMuonTagInfos.clone()
CMS2JPTSoftMuonTagInfos.jets = "JetPlusTrackZSPCorJetAntiKt5"
CMS2JPTSoftMuonBJetTags = softMuonBJetTags.clone()
CMS2JPTSoftMuonBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2JPTSoftMuonTagInfos") )
CMS2JPTSoftMuonByIP3dBJetTags = softMuonByIP3dBJetTags.clone()
CMS2JPTSoftMuonByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2JPTSoftMuonTagInfos") )
CMS2JPTSoftMuonByPtBJetTags = softMuonByPtBJetTags.clone()
CMS2JPTSoftMuonByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2JPTSoftMuonTagInfos") )
#Finally, there needs to be a CMS2 path running all these modules
# prepare a path running the CMS2 modules
CMS2JPTJetTracksAssociator = cms.Sequence(
    CMS2JPTJetTracksAssociatorAtVertex
)

CMS2JPTJetBtaggingIP = cms.Sequence(
    CMS2JPTImpactParameterTagInfos * (
        CMS2JPTTrackCountingHighEffBJetTags +
        CMS2JPTTrackCountingHighPurBJetTags +
        CMS2JPTJetProbabilityBJetTags +
        CMS2JPTJetBProbabilityBJetTags
    )
)

CMS2JPTJetBtaggingSV = cms.Sequence(
    CMS2JPTImpactParameterTagInfos *
    CMS2JPTSecondaryVertexTagInfos * (
#        CMS2JPTSimpleSecondaryVertexBJetTags +
        CMS2JPTSimpleSecondaryVertexHighEffBJetTags +
        CMS2JPTSimpleSecondaryVertexHighPurBJetTags +
        CMS2JPTCombinedSecondaryVertexBJetTags +
        CMS2JPTCombinedSecondaryVertexMVABJetTags
    )
)

#CMS2JPTJetghostBTagging = cms.Sequence(
#    CMS2JPTghostVertexTagInfos *
#    CMS2JPTghostTrackBJetTags
#)

CMS2JPTJetBtaggingEle = cms.Sequence(
    #btagSoftElectrons *
    softElectronCands *
    CMS2JPTSoftElectronTagInfos * (
    CMS2JPTSoftElectronByIP3dBJetTags + 
#    CMS2JPTSoftElectronBJetTags
    CMS2JPTSoftElectronByPtBJetTags
    )
)

CMS2JPTJetBtaggingMu = cms.Sequence(
    CMS2JPTSoftMuonTagInfos * (
        CMS2JPTSoftMuonBJetTags +
        CMS2JPTSoftMuonByIP3dBJetTags +
        CMS2JPTSoftMuonByPtBJetTags
    )
)

CMS2JPTJetBtagging = cms.Sequence(
    CMS2JPTJetBtaggingIP +
    CMS2JPTJetBtaggingSV +
#    CMS2JPTJetghostBTagging +
    CMS2JPTJetBtaggingEle +
    CMS2JPTJetBtaggingMu
)

CMS2JPTBtagging = cms.Sequence(
    CMS2JPTJetTracksAssociator *
    CMS2JPTJetBtagging
)
