import FWCore.ParameterSet.Config as cms

#from CMS3.NtupleMaker.pfElectronSequence_cff import *
#from CMS3.NtupleMaker.pfMuonSequence_cff import *
from CMS3.NtupleMaker.pfJetMaker_cfi import *
from CMS3.NtupleMaker.pfmetMaker_cfi import *
#from CMS3.NtupleMaker.pfmetSequence_cff import *
from CMS3.NtupleMaker.pftauMaker_cfi import *
#from CMS3.NtupleMaker.pfMuToMuAssMaker_cfi import *
#from CMS3.NtupleMaker.muToPFMuAssMaker_cfi import *
#from CMS3.NtupleMaker.pfElToElAssMaker_cfi import *
#from CMS3.NtupleMaker.elToPFElAssMaker_cfi import *
#from CMS3.NtupleMaker.bTagPFSequence_cfi import *
#from CMS3.NtupleMaker.bTagPFJetMaker_cfi import *
from CMS3.NtupleMaker.pfCandidateMaker_cfi import *
#from CMS3.NtupleMaker.trackIsolationMaker_cfi import *
#from CMS3.NtupleMaker.trkMetSequence_cff   import *
#from CMS3.NtupleMaker.mvaJetIdMaker_cfi import *
#
#
#from CommonTools.ParticleFlow.TopProjectors.pfNoMuon_cfi import *
#from CommonTools.ParticleFlow.Isolation.pfIsolatedMuons_cfi import *
#from CommonTools.ParticleFlow.pfNoPileUp_cff  import *
#from CommonTools.ParticleFlow.pfNoPileUpIso_cff  import *
#from CommonTools.ParticleFlow.ParticleSelectors.pfAllMuons_cfi import *
#from CommonTools.ParticleFlow.ParticleSelectors.pfAllElectrons_cfi import *
#from CommonTools.ParticleFlow.ParticleSelectors.pfAllNeutralHadrons_cfi  import *
#from CommonTools.ParticleFlow.ParticleSelectors.pfAllChargedHadrons_cfi import *
#from CommonTools.ParticleFlow.ParticleSelectors.pfAllPhotons_cfi import *
#from CommonTools.ParticleFlow.pfMuons_cff import *
#from JetMETCorrections.Configuration.DefaultJEC_cff import *
#
#CMS2pfIsolatedMuons = pfIsolatedMuons.clone()
#CMS2pfIsolatedMuons.src = cms.InputTag("pfAllMuons")
#CMS2pfIsolatedMuons.isolationValueMapsCharged = cms.VInputTag(
#  cms.InputTag("CMS2isoValMuonWithCharged03")
#)
#CMS2pfIsolatedMuons.isolationValueMapsNeutral = cms.VInputTag(
#  cms.InputTag("CMS2isoValMuonWithNeutral03"),
#  cms.InputTag("CMS2isoValMuonWithPhotons03")
#)
#
#CMS2pfNoMuon = pfNoMuon.clone()
#CMS2pfNoMuon.topCollection = cms.InputTag("CMS2pfIsolatedMuons")
#
#CMS2pfAllElectrons = pfAllElectrons.clone()
#CMS2pfAllElectrons.src = cms.InputTag("CMS2pfNoMuon")
#
#from CMS3.NtupleMaker.muonIsolationMaker_cfi import *
#from CMS3.NtupleMaker.electronIsolationMaker_cfi import *
#pfNoPileUpClones = cms.EDProducer("PFCandidateFromFwdPtrProducer", src = cms.InputTag("pfNoPileUp") )
#muonIsolationMaker.pfNoPileUpInputTag_ = cms.InputTag("pfNoPileUpClones")
#electronIsolationMaker.pfNoPileUpInputTag_ = cms.InputTag("pfNoPileUpClones")
#
#cms2PFNoTauSequence = cms.Sequence( pfCandidateMaker )
