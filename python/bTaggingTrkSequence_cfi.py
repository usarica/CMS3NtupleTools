import FWCore.ParameterSet.Config as cms

# sequence comes from:
# https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookBTagging#DifferentJets

# b-tagging general configuration
from RecoJets.JetAssociationProducers.ic5JetTracksAssociatorAtVertex_cfi import *
from RecoBTag.Configuration.RecoBTag_cff import *

# create a CMS2 jets and tracks association
CMS2TrkJetTracksAssociatorAtVertex = ic5JetTracksAssociatorAtVertex.clone()
CMS2TrkJetTracksAssociatorAtVertex.jets = "SISCone5TrkJets"
CMS2TrkJetTracksAssociatorAtVertex.tracks = "generalTracks"

#The one needs to clone the b-tag producers and instruct them to use this new collection. First the impact parameter-based b-tag:

# impact parameter b-tag
CMS2TrkImpactParameterTagInfos = impactParameterTagInfos.clone()
CMS2TrkImpactParameterTagInfos.jetTracks = "CMS2TrkJetTracksAssociatorAtVertex"
CMS2TrkTrackCountingHighEffBJetTags = trackCountingHighEffBJetTags.clone()
CMS2TrkTrackCountingHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2TrkImpactParameterTagInfos") )
CMS2TrkTrackCountingHighPurBJetTags = trackCountingHighPurBJetTags.clone()
CMS2TrkTrackCountingHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2TrkImpactParameterTagInfos") )
CMS2TrkJetProbabilityBJetTags = jetProbabilityBJetTags.clone()
CMS2TrkJetProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2TrkImpactParameterTagInfos") )
CMS2TrkJetBProbabilityBJetTags = jetBProbabilityBJetTags.clone()
CMS2TrkJetBProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2TrkImpactParameterTagInfos") )
#removed for 312
#CMS2TrkImpactParameterMVABJetTags = impactParameterMVABJetTags.clone()
#CMS2TrkImpactParameterMVABJetTags.tagInfos = cms.VInputTag(cms.InputTag("CMS2TrkImpactParameterTagInfos"))


#Then the secondary vertex-based b-tag. Note that these producers inherit the jets and tracks they use from the impact parameter modules:

# secondary vertex b-tag
CMS2TrkSecondaryVertexTagInfos = secondaryVertexTagInfos.clone()
CMS2TrkSecondaryVertexTagInfos.trackIPTagInfos = "CMS2TrkImpactParameterTagInfos"
CMS2TrkSimpleSecondaryVertexBJetTags = simpleSecondaryVertexBJetTags.clone()
CMS2TrkSimpleSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2TrkSecondaryVertexTagInfos") )
CMS2TrkCombinedSecondaryVertexBJetTags = combinedSecondaryVertexBJetTags.clone()
CMS2TrkCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2TrkImpactParameterTagInfos"), cms.InputTag("CMS2TrkSecondaryVertexTagInfos") )
CMS2TrkCombinedSecondaryVertexMVABJetTags = combinedSecondaryVertexMVABJetTags.clone()
CMS2TrkCombinedSecondaryVertexMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2TrkImpactParameterTagInfos"), cms.InputTag("CMS2TrkSecondaryVertexTagInfos") )



# soft electron b-tag
CMS2TrkSoftElectronTagInfos = softElectronTagInfos.clone()
CMS2TrkSoftElectronTagInfos.jets = "SISCone5TrkJets"
#removed for 312
#CMS2TrkSoftElectronBJetTags = softElectronBJetTags.clone()
#CMS2TrkSoftElectronBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2TrkSoftElectronTagInfos") )
###########################################new
CMS2TrkSoftElectronByPtBJetTags = softElectronByPtBJetTags.clone()
CMS2TrkSoftElectronByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2TrkSoftElectronTagInfos") )
CMS2TrkSoftElectronByIP3dBJetTags = softElectronByIP3dBJetTags.clone()
CMS2TrkSoftElectronByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2TrkSoftElectronTagInfos") )

# soft muon b-tag
CMS2TrkSoftMuonTagInfos = softMuonTagInfos.clone()
CMS2TrkSoftMuonTagInfos.jets = "SISCone5TrkJets"
CMS2TrkSoftMuonBJetTags = softMuonBJetTags.clone()
CMS2TrkSoftMuonBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2TrkSoftMuonTagInfos") )
CMS2TrkSoftMuonByIP3dBJetTags = softMuonByIP3dBJetTags.clone()
CMS2TrkSoftMuonByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2TrkSoftMuonTagInfos") )
CMS2TrkSoftMuonByPtBJetTags = softMuonByPtBJetTags.clone()
CMS2TrkSoftMuonByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2TrkSoftMuonTagInfos") )
CMS2TrkSoftMuonNoIPBJetTags = softMuonNoIPBJetTags.clone()
CMS2TrkSoftMuonNoIPBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2TrkSoftMuonTagInfos") )

CMS2TrkJetTracksAssociator = cms.Sequence(
    CMS2TrkJetTracksAssociatorAtVertex
)

CMS2TrkJetBtaggingIP = cms.Sequence(
    CMS2TrkImpactParameterTagInfos * (
        CMS2TrkTrackCountingHighEffBJetTags +
        CMS2TrkTrackCountingHighPurBJetTags +
        CMS2TrkJetProbabilityBJetTags +
        CMS2TrkJetBProbabilityBJetTags
        #CMS2TrkImpactParameterMVABJetTags
    )
)

CMS2TrkJetBtaggingSV = cms.Sequence(
    CMS2TrkImpactParameterTagInfos *
    CMS2TrkSecondaryVertexTagInfos * (
        CMS2TrkSimpleSecondaryVertexBJetTags +
        CMS2TrkCombinedSecondaryVertexBJetTags +
        CMS2TrkCombinedSecondaryVertexMVABJetTags
    )
)

CMS2TrkJetBtaggingEle = cms.Sequence(
    #btagSoftElectrons *
    CMS2TrkSoftElectronTagInfos * (
    #CMS2TrkSoftElectronBJetTags
    CMS2TrkSoftElectronByIP3dBJetTags +
    CMS2TrkSoftElectronByPtBJetTags
    )
)

CMS2TrkJetBtaggingMu = cms.Sequence(
    CMS2TrkSoftMuonTagInfos * (
        CMS2TrkSoftMuonBJetTags +
        CMS2TrkSoftMuonByIP3dBJetTags +
        CMS2TrkSoftMuonByPtBJetTags +
        CMS2TrkSoftMuonNoIPBJetTags #this isn't even in the official reco sequence for 312, but it runs. WTF?
    )
)

CMS2TrkJetBtagging = cms.Sequence(
    CMS2TrkJetBtaggingIP +
    CMS2TrkJetBtaggingSV +
    CMS2TrkJetBtaggingEle +
    CMS2TrkJetBtaggingMu
)

CMS2TrkBtagging = cms.Sequence(
    CMS2TrkJetTracksAssociator *
    CMS2TrkJetBtagging
)
