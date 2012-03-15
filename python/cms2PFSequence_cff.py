import FWCore.ParameterSet.Config as cms

from CMS2.NtupleMaker.pfElectronSequence_cff import *
from CMS2.NtupleMaker.pfMuonSequence_cff import *
from CMS2.NtupleMaker.pfJetMaker_cfi import *
from CMS2.NtupleMaker.pfmetMaker_cfi import *
from CMS2.NtupleMaker.pftauMaker_cfi import *
from CMS2.NtupleMaker.pfMuToMuAssMaker_cfi import *
from CMS2.NtupleMaker.muToPFMuAssMaker_cfi import *
from CMS2.NtupleMaker.pfElToElAssMaker_cfi import *
from CMS2.NtupleMaker.elToPFElAssMaker_cfi import *
from CMS2.NtupleMaker.bTagPFSequence_cfi import *
from CMS2.NtupleMaker.bTagPFJetMaker_cfi import *
from CMS2.NtupleMaker.pfCandidateMaker_cfi import *
from CMS2.NtupleMaker.trkMetSequence_cff   import *

from CommonTools.ParticleFlow.ParticleSelectors.pfAllMuons_cfi import *
from CommonTools.ParticleFlow.ParticleSelectors.pfAllElectrons_cfi import *
from CommonTools.ParticleFlow.TopProjectors.pfNoMuon_cfi import *
from CommonTools.ParticleFlow.Isolation.pfIsolatedMuons_cfi import *
from CommonTools.ParticleFlow.pfNoPileUp_cff  import *
from CommonTools.ParticleFlow.ParticleSelectors.pfAllNeutralHadrons_cfi  import *
from CommonTools.ParticleFlow.ParticleSelectors.pfAllChargedHadrons_cfi import *
from CommonTools.ParticleFlow.ParticleSelectors.pfAllPhotons_cfi import *
from CommonTools.ParticleFlow.ParticleSelectors.pfAllMuons_cfi import *
from CommonTools.ParticleFlow.ParticleSelectors.pfAllElectrons_cfi import *
from CommonTools.ParticleFlow.pfMuons_cff import *

CMS2pfIsolatedMuons = pfIsolatedMuons.clone()
CMS2pfIsolatedMuons.src = cms.InputTag("pfAllMuons")
CMS2pfIsolatedMuons.isolationValueMaps = cms.VInputTag(cms.InputTag("CMS2isoValMuonWithCharged"),
                                                       cms.InputTag("CMS2isoValMuonWithNeutral"),
                                                       cms.InputTag("CMS2isoValMuonWithPhotons"))

CMS2pfNoMuon = pfNoMuon.clone()
CMS2pfNoMuon.topCollection = cms.InputTag("CMS2pfIsolatedMuons")

CMS2pfAllElectrons = pfAllElectrons.clone()
CMS2pfAllElectrons.src = cms.InputTag("CMS2pfNoMuon")

cms2PFSequence = cms.Sequence(pfJetMaker + pfmetMaker + pftauMaker + CMS2PFBtagging + bTagPFJetMaker + pfNoPileUpSequence + pfAllNeutralHadrons + pfAllChargedHadrons + pfAllPhotons + pfAllMuons + 
#CMS2pfMuonIsolationSequence + CMS2pfIsolatedMuons + 
CMS2pfNoMuon + CMS2pfAllElectrons + pfMuonMaker + pfMuToMuAssMaker + muToPFMuAssMaker + CMS2pfElectronIsolationSequence + pfElectronMaker + pfElToElAssMaker + elToPFElAssMaker + pfCandidateMaker + trkMetSequence)

cms2PFNoTauSequence = cms.Sequence(pfJetMaker + pfmetMaker + CMS2PFBtagging + bTagPFJetMaker + pfNoPileUpSequence + pfAllNeutralHadrons + pfAllChargedHadrons + pfAllPhotons + pfAllMuons + 
#CMS2pfMuonIsolationSequence + CMS2pfIsolatedMuons + 
CMS2pfNoMuon + CMS2pfAllElectrons + pfMuonMaker + pfMuToMuAssMaker + muToPFMuAssMaker + CMS2pfElectronIsolationSequence + pfElectronMaker + pfElToElAssMaker + elToPFElAssMaker + pfCandidateMaker + trkMetSequence)
