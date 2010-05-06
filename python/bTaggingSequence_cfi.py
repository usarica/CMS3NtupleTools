import FWCore.ParameterSet.Config as cms

# sequence comes from:
# https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookBTagging#DifferentJets


# b-tagging general configuration
from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import *
from RecoBTag.Configuration.RecoBTag_cff import *
from RecoBTag.SoftLepton.softElectronCandProducer_cfi import *



#The first step is to clone the definition of the association between jets and tracks. This is where most of the customization will take place. As an example, here generalTracks and sisCone5CaloJets are used - to try different ones, change them to what you need throughout the configuration (check later the soft lepton part, too):
# create a CMS2 jets and tracks association
CMS2JetTracksAssociatorAtVertex = ic5JetTracksAssociatorAtVertex.clone()
CMS2JetTracksAssociatorAtVertex.jets = "ak5CaloJets"
CMS2JetTracksAssociatorAtVertex.tracks = "generalTracks"
#The one needs to clone the b-tag producers and instruct them to use this CMS2 collection. First the impact parameter-based b-tag:
# impact parameter b-tag
CMS2ImpactParameterTagInfos = impactParameterTagInfos.clone()
CMS2ImpactParameterTagInfos.jetTracks = "CMS2JetTracksAssociatorAtVertex"
CMS2TrackCountingHighEffBJetTags = trackCountingHighEffBJetTags.clone()
CMS2TrackCountingHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2ImpactParameterTagInfos") )
CMS2TrackCountingHighPurBJetTags = trackCountingHighPurBJetTags.clone()
CMS2TrackCountingHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2ImpactParameterTagInfos") )
CMS2JetProbabilityBJetTags = jetProbabilityBJetTags.clone()
CMS2JetProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2ImpactParameterTagInfos") )
CMS2JetBProbabilityBJetTags = jetBProbabilityBJetTags.clone()
CMS2JetBProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2ImpactParameterTagInfos") )
#Then the secondary vertex-based b-tag. Note that these producers inherit the jets and tracks they use from the impact parameter modules: # secondary vertex b-tag
CMS2SecondaryVertexTagInfos = secondaryVertexTagInfos.clone()
CMS2SecondaryVertexTagInfos.trackIPTagInfos = cms.InputTag("CMS2ImpactParameterTagInfos")
#CMS2SimpleSecondaryVertexBJetTags = simpleSecondaryVertexBJetTags.clone()
#CMS2SimpleSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2SecondaryVertexTagInfos") )
CMS2SimpleSecondaryVertexHighEffBJetTags = simpleSecondaryVertexHighEffBJetTags.clone()
CMS2SimpleSecondaryVertexHighEffBJetTags.tagInfos = cms.VInputTag(cms.InputTag("CMS2SecondaryVertexTagInfos"))
CMS2SimpleSecondaryVertexHighPurBJetTags = simpleSecondaryVertexHighPurBJetTags.clone()
CMS2SimpleSecondaryVertexHighPurBJetTags.tagInfos = cms.VInputTag(cms.InputTag("CMS2SecondaryVertexTagInfos"))
CMS2CombinedSecondaryVertexBJetTags = combinedSecondaryVertexBJetTags.clone()
CMS2CombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2ImpactParameterTagInfos"), cms.InputTag("CMS2SecondaryVertexTagInfos") )
CMS2CombinedSecondaryVertexMVABJetTags = combinedSecondaryVertexMVABJetTags.clone()
CMS2CombinedSecondaryVertexMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2ImpactParameterTagInfos"), cms.InputTag("CMS2SecondaryVertexTagInfos") )
# add some ghost b-tagging stuff
#CMS2ghostVertexTagInfos = ghostTrackVertexTagInfos.clone()
#CMS2ghostVertexTagInfos.trackIPTagInfos = "CMS2ImpactParameterTagInfos"
#CMS2ghostTrackBJetTags = ghostTrackBJetTags.clone()
#CMS2ghostTrackBJetTags.tagInfos = cms.VInputTag(cms.InputTag("CMS2ImpactParameterTagInfos"),
#                                                cms.InputTag("CMS2ghostVertexTagInfos"))
#And the soft lepton b-tag. These producers will accept as input either the raw jets, or the association collection:
# soft electron b-tag
CMS2SoftElectronTagInfos = softElectronTagInfos.clone()
CMS2SoftElectronTagInfos.jets = "ak5TrackJets"
#CMS2SoftElectronBJetTags = softElectronBJetTags.clone()
#CMS2SoftElectronBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2SoftElectronTagInfos") )
CMS2SoftElectronByIP3dBJetTags = softElectronByIP3dBJetTags.clone()
CMS2SoftElectronByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2SoftElectronTagInfos") )
CMS2SoftElectronByPtBJetTags = softElectronByPtBJetTags.clone()
CMS2SoftElectronByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2SoftElectronTagInfos") )


# soft muon b-tag
CMS2SoftMuonTagInfos = softMuonTagInfos.clone()
CMS2SoftMuonTagInfos.jets = "ak5TrackJets"
CMS2SoftMuonBJetTags = softMuonBJetTags.clone()
CMS2SoftMuonBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2SoftMuonTagInfos") )
CMS2SoftMuonByIP3dBJetTags = softMuonByIP3dBJetTags.clone()
CMS2SoftMuonByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2SoftMuonTagInfos") )
CMS2SoftMuonByPtBJetTags = softMuonByPtBJetTags.clone()
CMS2SoftMuonByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2SoftMuonTagInfos") )
#Finally, there needs to be a CMS2 path running all these modules
# prepare a path running the CMS2 modules
CMS2JetTracksAssociator = cms.Sequence(
    CMS2JetTracksAssociatorAtVertex
)

CMS2JetBtaggingIP = cms.Sequence(
    CMS2ImpactParameterTagInfos * (
        CMS2TrackCountingHighEffBJetTags +
        CMS2TrackCountingHighPurBJetTags +
        CMS2JetProbabilityBJetTags +
        CMS2JetBProbabilityBJetTags
    )
)

CMS2JetBtaggingSV = cms.Sequence(
    CMS2ImpactParameterTagInfos *
    CMS2SecondaryVertexTagInfos * (
#        CMS2SimpleSecondaryVertexBJetTags +
        CMS2SimpleSecondaryVertexHighEffBJetTags +
        CMS2SimpleSecondaryVertexHighPurBJetTags +
        CMS2CombinedSecondaryVertexBJetTags +
        CMS2CombinedSecondaryVertexMVABJetTags
    )
)

#CMS2JetghostBTagging = cms.Sequence(
#    CMS2ghostVertexTagInfos *
#    CMS2ghostTrackBJetTags
#)

CMS2JetBtaggingEle = cms.Sequence(
    #btagSoftElectrons *
    softElectronCands *
    CMS2SoftElectronTagInfos * (
    CMS2SoftElectronByIP3dBJetTags + 
#    CMS2SoftElectronBJetTags
    CMS2SoftElectronByPtBJetTags
    )
)

CMS2JetBtaggingMu = cms.Sequence(
    CMS2SoftMuonTagInfos * (
        CMS2SoftMuonBJetTags +
        CMS2SoftMuonByIP3dBJetTags +
        CMS2SoftMuonByPtBJetTags
    )
)

CMS2JetBtagging = cms.Sequence(
    CMS2JetBtaggingIP +
    CMS2JetBtaggingSV +
#    CMS2JetghostBTagging +
    CMS2JetBtaggingEle +
    CMS2JetBtaggingMu
)

CMS2Btagging = cms.Sequence(
    CMS2JetTracksAssociator *
    CMS2JetBtagging
)
