import FWCore.ParameterSet.Config as cms

from CMS2.NtupleMaker.pfElectronSequence_cff import *
from CMS2.NtupleMaker.pfMuonSequence_cff import *
from CMS2.NtupleMaker.pfJetMaker_cfi import *
#from CMS2.NtupleMaker.pfmetMaker_cfi import *
from CMS2.NtupleMaker.pfmetSequence_cff import *
from CMS2.NtupleMaker.pftauMaker_cfi import *
from CMS2.NtupleMaker.pfMuToMuAssMaker_cfi import *
from CMS2.NtupleMaker.muToPFMuAssMaker_cfi import *
from CMS2.NtupleMaker.pfElToElAssMaker_cfi import *
from CMS2.NtupleMaker.elToPFElAssMaker_cfi import *
from CMS2.NtupleMaker.bTagPFSequence_cfi import *
from CMS2.NtupleMaker.bTagPFJetMaker_cfi import *
from CMS2.NtupleMaker.pfCandidateMaker_cfi import *
from CMS2.NtupleMaker.trackIsolationMaker_cfi import *
from CMS2.NtupleMaker.trkMetSequence_cff   import *
from CMS2.NtupleMaker.mvaJetIdMaker_cfi import *


from CommonTools.ParticleFlow.TopProjectors.pfNoMuon_cfi import *
from CommonTools.ParticleFlow.Isolation.pfIsolatedMuons_cfi import *
from CommonTools.ParticleFlow.pfNoPileUp_cff  import *
from CommonTools.ParticleFlow.pfNoPileUpIso_cff  import *
from CommonTools.ParticleFlow.ParticleSelectors.pfAllMuons_cfi import *
from CommonTools.ParticleFlow.ParticleSelectors.pfAllElectrons_cfi import *
from CommonTools.ParticleFlow.ParticleSelectors.pfAllNeutralHadrons_cfi  import *
from CommonTools.ParticleFlow.ParticleSelectors.pfAllChargedHadrons_cfi import *
from CommonTools.ParticleFlow.ParticleSelectors.pfAllPhotons_cfi import *
from CommonTools.ParticleFlow.pfMuons_cff import *
from JetMETCorrections.Configuration.DefaultJEC_cff import *

#VERTEX_SEL=("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2")
#        
#goodPrimaryVertices = cms.EDFilter("VertexSelector",
#                                           src = cms.InputTag("offlinePrimaryVertices"),
#                                           cut = cms.string(VERTEX_SEL),
#                                           filter = cms.bool(True),
#                                           )
#pfPileUp.Vertices = cms.InputTag("goodPrimaryVertices")
#pfPileUp.PFCandidates = cms.InputTag("particleFlow")
#pfNoPileUp.bottomCollection = cms.InputTag("particleFlow")

CMS2pfIsolatedMuons = pfIsolatedMuons.clone()
CMS2pfIsolatedMuons.src = cms.InputTag("pfAllMuons")
CMS2pfIsolatedMuons.isolationValueMapsCharged = cms.VInputTag(
  cms.InputTag("CMS2isoValMuonWithCharged03")
)
CMS2pfIsolatedMuons.isolationValueMapsNeutral = cms.VInputTag(
  cms.InputTag("CMS2isoValMuonWithNeutral03"),
  cms.InputTag("CMS2isoValMuonWithPhotons03")
)
#CMS2pfIsolatedMuons.isolationValueMaps = cms.VInputTag(
#  cms.InputTag("CMS2isoValMuonWithCharged03"),
#  cms.InputTag("CMS2isoValMuonWithNeutral03"),
#  cms.InputTag("CMS2isoValMuonWithPhotons03")
#)

CMS2pfNoMuon = pfNoMuon.clone()
CMS2pfNoMuon.topCollection = cms.InputTag("CMS2pfIsolatedMuons")

CMS2pfAllElectrons = pfAllElectrons.clone()
CMS2pfAllElectrons.src = cms.InputTag("CMS2pfNoMuon")

#cms2PFSequence = cms.Sequence(pfJetMaker + pfmetMaker + pftauMaker + CMS2PFBtagging + bTagPFJetMaker + pfNoPileUpSequence + pfAllNeutralHadrons + pfAllChargedHadrons + pfAllPhotons + pfAllMuons + 
##CMS2pfMuonIsolationSequence + CMS2pfIsolatedMuons + 
#CMS2pfNoMuon + CMS2pfAllElectrons + pfMuonMaker + pfMuToMuAssMaker + muToPFMuAssMaker + 
##CMS2pfElectronIsolationSequence + 
#pfElectronMaker + pfElToElAssMaker + elToPFElAssMaker + pfCandidateMaker + trkMetSequence)

   
#pfAllNeutralHadrons.src = cms.InputTag("pfNoPileUp")
#pfAllChargedHadrons.src = cms.InputTag("pfNoPileUp")
#pfAllPhotons.src        = cms.InputTag("pfNoPileUp")

from CMS2.NtupleMaker.muonIsolationMaker_cfi import *
from CMS2.NtupleMaker.electronIsolationMaker_cfi import *
pfNoPileUpClones = cms.EDProducer("PFCandidateFromFwdPtrProducer", src = cms.InputTag("pfNoPileUp") )
muonIsolationMaker.pfNoPileUpInputTag_ = cms.InputTag("pfNoPileUpClones")
electronIsolationMaker.pfNoPileUpInputTag_ = cms.InputTag("pfNoPileUpClones")

cms2PFNoTauSequence = cms.Sequence( 
  pfJetMaker + 
#  pfmetMaker + 
  CMS2pfMetSequence +
  CMS2PFBtagging + 
  bTagPFJetMaker + 
#  goodPrimaryVertices +
  pfNoPileUpIsoSequence +
  pfNoPileUpSequence +
  pfNoPileUpClones + 
  pfAllNeutralHadrons + 
  pfAllChargedHadrons + 
  pfAllPhotons + 
  pfAllMuons +
  CMS2pfMuonIsolationSequence +
  CMS2pfIsolatedMuons +
  CMS2pfNoMuon + 
  CMS2pfAllElectrons + 
  pfMuonMaker + 
  pfMuToMuAssMaker + 
  muToPFMuAssMaker + 
  CMS2pfElectronIsolationSequence +
  pfElectronMaker + 
  pfElToElAssMaker + 
  elToPFElAssMaker + 
  pfCandidateMaker + 
  trackIsolationMaker + 
  trkMetSequence +
  ak5PFJetsL1FastL2L3 +
  ak5PFJetsL1FastL2L3Residual +
  mvaJetIdMaker +
  mvaJetIdMakerFull5x +
  mvaJetIdMakerFull53x +
  muonIsolationMaker +
  electronIsolationMaker
)
