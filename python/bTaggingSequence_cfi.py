import FWCore.ParameterSet.Config as cms

# sequence comes from:
# https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookBTagging#DifferentJets

# b-tagging general configuration
from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import *
from RecoBTag.Configuration.RecoBTag_cff import *

# create a CMS2 jets and tracks association
CMS2JetTracksAssociatorAtVertex = ic5JetTracksAssociatorAtVertex.clone()
CMS2JetTracksAssociatorAtVertex.jets = "prunedUncorrectedCMS2Jets"
CMS2JetTracksAssociatorAtVertex.tracks = "generalTracks"

#The one needs to clone the b-tag producers and instruct them to use this new collection. First the impact parameter-based b-tag:

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
CMS2ImpactParameterMVABJetTags = impactParameterMVABJetTags.clone()
CMS2ImpactParameterMVABJetTags.tagInfos = cms.VInputTag(cms.InputTag("CMS2ImpactParameterTagInfos"))


#Then the secondary vertex-based b-tag. Note that these producers inherit the jets and tracks they use from the impact parameter modules:

# secondary vertex b-tag
CMS2SecondaryVertexTagInfos = secondaryVertexTagInfos.clone()
CMS2SecondaryVertexTagInfos.trackIPTagInfos = "CMS2ImpactParameterTagInfos"
CMS2SimpleSecondaryVertexBJetTags = simpleSecondaryVertexBJetTags.clone()
CMS2SimpleSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2SecondaryVertexTagInfos") )
CMS2CombinedSecondaryVertexBJetTags = combinedSecondaryVertexBJetTags.clone()
CMS2CombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2ImpactParameterTagInfos"), cms.InputTag("CMS2SecondaryVertexTagInfos") )
CMS2CombinedSecondaryVertexMVABJetTags = combinedSecondaryVertexMVABJetTags.clone()
CMS2CombinedSecondaryVertexMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2ImpactParameterTagInfos"), cms.InputTag("CMS2SecondaryVertexTagInfos") )



# soft electron b-tag
CMS2SoftElectronTagInfos = softElectronTagInfos.clone()
CMS2SoftElectronTagInfos.jets = "prunedUncorrectedCMS2Jets"
CMS2SoftElectronBJetTags = softElectronBJetTags.clone()
CMS2SoftElectronBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2SoftElectronTagInfos") )

# soft muon b-tag
CMS2SoftMuonTagInfos = softMuonTagInfos.clone()
CMS2SoftMuonTagInfos.jets = "prunedUncorrectedCMS2Jets"
CMS2SoftMuonBJetTags = softMuonBJetTags.clone()
CMS2SoftMuonBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2SoftMuonTagInfos") )
CMS2SoftMuonByIP3dBJetTags = softMuonByIP3dBJetTags.clone()
CMS2SoftMuonByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2SoftMuonTagInfos") )
CMS2SoftMuonByPtBJetTags = softMuonByPtBJetTags.clone()
CMS2SoftMuonByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2SoftMuonTagInfos") )
CMS2SoftMuonNoIPBJetTags = softMuonNoIPBJetTags.clone()
CMS2SoftMuonNoIPBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2SoftMuonTagInfos") )

CMS2JetTracksAssociator = cms.Sequence(
    CMS2JetTracksAssociatorAtVertex
)

CMS2JetBtaggingIP = cms.Sequence(
    CMS2ImpactParameterTagInfos * (
        CMS2TrackCountingHighEffBJetTags +
        CMS2TrackCountingHighPurBJetTags +
        CMS2JetProbabilityBJetTags +
        CMS2JetBProbabilityBJetTags +
        CMS2ImpactParameterMVABJetTags
    )
)

CMS2JetBtaggingSV = cms.Sequence(
    CMS2ImpactParameterTagInfos *
    CMS2SecondaryVertexTagInfos * (
        CMS2SimpleSecondaryVertexBJetTags +
        CMS2CombinedSecondaryVertexBJetTags +
        CMS2CombinedSecondaryVertexMVABJetTags
    )
)

CMS2JetBtaggingEle = cms.Sequence(
    btagSoftElectrons *
    CMS2SoftElectronTagInfos *
    CMS2SoftElectronBJetTags
)

CMS2JetBtaggingMu = cms.Sequence(
    CMS2SoftMuonTagInfos * (
        CMS2SoftMuonBJetTags +
        CMS2SoftMuonByIP3dBJetTags +
        CMS2SoftMuonByPtBJetTags +
        CMS2SoftMuonNoIPBJetTags
    )
)

CMS2JetBtagging = cms.Sequence(
    CMS2JetBtaggingIP +
    CMS2JetBtaggingSV +
    CMS2JetBtaggingEle +
    CMS2JetBtaggingMu
)

CMS2Btagging = cms.Sequence(
    CMS2JetTracksAssociator *
    CMS2JetBtagging
)
