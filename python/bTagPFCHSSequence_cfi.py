import FWCore.ParameterSet.Config as cms

# sequence comes from:
# https://twiki.cern.ch/twiki/bin/view/CMS/WorkBookBTagging#DifferentJets


# b-tagging general configuration
from RecoJets.JetAssociationProducers.ic5PFJetTracksAssociatorAtVertex_cfi import *
from RecoBTag.Configuration.RecoBTag_cff import *
#from RecoBTag.SoftLepton.softElectronCandProducer_cfi import *



#The first step is to clone the definition of the association between jets and tracks. This is where most of the customization will take place. As an example, here generalTracks and sisCone5CaloJets are used - to try different ones, change them to what you need throughout the configuration (check later the soft lepton part, too):
# create a CMS2 jets and tracks association
CMS2PFCHSJetTracksAssociatorAtVertex = ic5PFJetTracksAssociatorAtVertex.clone()
CMS2PFCHSJetTracksAssociatorAtVertex.jets = "ak5PFJetsCHS"
CMS2PFCHSJetTracksAssociatorAtVertex.tracks = "generalTracks"
#The one needs to clone the b-tag producers and instruct them to use this CMS2 collection. First the impact parameter-based b-tag:
# impact parameter b-tag
CMS2PFCHSImpactParameterTagInfos = impactParameterTagInfos.clone()
CMS2PFCHSImpactParameterTagInfos.jetTracks = "CMS2PFCHSJetTracksAssociatorAtVertex"
CMS2PFCHSTrackCountingHighEffBJetTags = trackCountingHighEffBJetTags.clone()
CMS2PFCHSTrackCountingHighEffBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFCHSImpactParameterTagInfos") )
CMS2PFCHSTrackCountingHighPurBJetTags = trackCountingHighPurBJetTags.clone()
CMS2PFCHSTrackCountingHighPurBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFCHSImpactParameterTagInfos") )
CMS2PFCHSJetProbabilityBJetTags = jetProbabilityBJetTags.clone()
CMS2PFCHSJetProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFCHSImpactParameterTagInfos") )
CMS2PFCHSJetBProbabilityBJetTags = jetBProbabilityBJetTags.clone()
CMS2PFCHSJetBProbabilityBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFCHSImpactParameterTagInfos") )

#Then the secondary vertex-based b-tag. Note that these producers inherit the jets and tracks they use from the impact parameter modules: # secondary vertex b-tag
CMS2PFCHSSecondaryVertexTagInfos = secondaryVertexTagInfos.clone()
CMS2PFCHSSecondaryVertexTagInfos.trackIPTagInfos = cms.InputTag("CMS2PFCHSImpactParameterTagInfos")

CMS2PFCHSSimpleSecondaryVertexBJetTags = simpleSecondaryVertexBJetTags.clone()
CMS2PFCHSSimpleSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFCHSSecondaryVertexTagInfos") )

CMS2PFCHSSimpleSecondaryVertexHighEffBJetTags = simpleSecondaryVertexHighEffBJetTags.clone()
CMS2PFCHSSimpleSecondaryVertexHighEffBJetTags.tagInfos = cms.VInputTag(cms.InputTag("CMS2PFCHSSecondaryVertexTagInfos"))

CMS2PFCHSSimpleSecondaryVertexHighPurBJetTags = simpleSecondaryVertexHighPurBJetTags.clone()
CMS2PFCHSSimpleSecondaryVertexHighPurBJetTags.tagInfos = cms.VInputTag(cms.InputTag("CMS2PFCHSSecondaryVertexTagInfos"))

CMS2PFCHSCombinedSecondaryVertexBJetTags = combinedSecondaryVertexBJetTags.clone()
CMS2PFCHSCombinedSecondaryVertexBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFCHSImpactParameterTagInfos"), cms.InputTag("CMS2PFCHSSecondaryVertexTagInfos") )

CMS2PFCHSCombinedSecondaryVertexMVABJetTags = combinedSecondaryVertexMVABJetTags.clone()
CMS2PFCHSCombinedSecondaryVertexMVABJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFCHSImpactParameterTagInfos"), cms.InputTag("CMS2PFCHSSecondaryVertexTagInfos") )

# add some ghost b-tagging stuff
#CMS2PFghostVertexTagInfos = ghostTrackVertexTagInfos.clone()
#CMS2PFghostVertexTagInfos.trackIPTagInfos = "CMS2PFImpactParameterTagInfos"
#CMS2PFghostTrackBJetTags = ghostTrackBJetTags.clone()

#CMS2PFghostTrackBJetTags.tagInfos = cms.VInputTag(cms.InputTag("CMS2PFImpactParameterTagInfos"),
#                                                cms.InputTag("CMS2PFghostVertexTagInfos"))
#And the soft lepton b-tag. These producers will accept as input either the raw jets, or the association collection:
# soft electron b-tag
CMS2PFCHSSoftElectronTagInfos = softPFElectronsTagInfos.clone()
CMS2PFCHSSoftElectronTagInfos.jets = "ak5PFJetsCHS"
CMS2PFCHSSoftElectronBJetTags = softPFElectronBJetTags.clone()
CMS2PFCHSSoftElectronBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFCHSSoftElectronTagInfos") )
CMS2PFCHSSoftElectronByIP3dBJetTags = softPFElectronByIP3dBJetTags.clone()
CMS2PFCHSSoftElectronByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFCHSSoftElectronTagInfos") )
CMS2PFCHSSoftElectronByPtBJetTags = softPFElectronByPtBJetTags.clone()
CMS2PFCHSSoftElectronByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFCHSSoftElectronTagInfos") )


# soft muon b-tag
CMS2PFCHSSoftMuonTagInfos = softPFMuonsTagInfos.clone()
CMS2PFCHSSoftMuonTagInfos.jets = "ak5PFJetsCHS"
CMS2PFCHSSoftMuonBJetTags = softPFMuonBJetTags.clone()
CMS2PFCHSSoftMuonBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFCHSSoftMuonTagInfos") )
CMS2PFCHSSoftMuonByIP3dBJetTags = softPFMuonByIP3dBJetTags.clone()
CMS2PFCHSSoftMuonByIP3dBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFCHSSoftMuonTagInfos") )
CMS2PFCHSSoftMuonByPtBJetTags = softPFMuonByPtBJetTags.clone()
CMS2PFCHSSoftMuonByPtBJetTags.tagInfos = cms.VInputTag( cms.InputTag("CMS2PFCHSSoftMuonTagInfos") )
#Finally, there needs to be a CMS2 path running all these modules
# prepare a path running the CMS2 modules
CMS2PFCHSJetTracksAssociator = cms.Sequence(
    CMS2PFCHSJetTracksAssociatorAtVertex
)

CMS2PFCHSJetBtaggingIP = cms.Sequence(
    CMS2PFCHSImpactParameterTagInfos * (
        CMS2PFCHSTrackCountingHighEffBJetTags +
        CMS2PFCHSTrackCountingHighPurBJetTags +
        CMS2PFCHSJetProbabilityBJetTags +
        CMS2PFCHSJetBProbabilityBJetTags
    )
)

CMS2PFCHSJetBtaggingSV = cms.Sequence(
    CMS2PFCHSImpactParameterTagInfos *
    CMS2PFCHSSecondaryVertexTagInfos * (
        CMS2PFCHSSimpleSecondaryVertexBJetTags +
        CMS2PFCHSSimpleSecondaryVertexHighEffBJetTags +
        CMS2PFCHSSimpleSecondaryVertexHighPurBJetTags +
        CMS2PFCHSCombinedSecondaryVertexBJetTags +
        CMS2PFCHSCombinedSecondaryVertexMVABJetTags
    )
)

#CMS2PFCHSJetghostBTagging = cms.Sequence(
#    CMS2PFCHSghostVertexTagInfos *
#    CMS2PFCHSghostTrackBJetTags
#)

CMS2PFCHSJetBtaggingEle = cms.Sequence(
    #btagSoftElectrons *
#    softElectronCands *
    CMS2PFCHSSoftElectronTagInfos * (
    CMS2PFCHSSoftElectronByIP3dBJetTags + 
    CMS2PFCHSSoftElectronBJetTags +
    CMS2PFCHSSoftElectronByPtBJetTags
    )
)

CMS2PFCHSJetBtaggingMu = cms.Sequence(
     CMS2PFCHSSoftMuonTagInfos * (
         CMS2PFCHSSoftMuonBJetTags +
         CMS2PFCHSSoftMuonByIP3dBJetTags +
         CMS2PFCHSSoftMuonByPtBJetTags
     )
 )

CMS2PFCHSJetBtagging = cms.Sequence(
    CMS2PFCHSJetBtaggingIP +
    CMS2PFCHSJetBtaggingSV +
#    CMS2PFCHSJetghostBTagging +
    CMS2PFCHSJetBtaggingEle +
    CMS2PFCHSJetBtaggingMu
)

CMS2PFCHSBtagging = cms.Sequence(
    CMS2PFCHSJetTracksAssociator *
    CMS2PFCHSJetBtagging
)
