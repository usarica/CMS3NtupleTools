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

from PhysicsTools.PFCandProducer.ParticleSelectors.pfAllMuons_cfi import *
from PhysicsTools.PFCandProducer.ParticleSelectors.pfAllElectrons_cfi import *
from PhysicsTools.PFCandProducer.TopProjectors.pfNoMuon_cfi import *
from PhysicsTools.PFCandProducer.Isolation.pfIsolatedMuons_cfi import *
from PhysicsTools.PFCandProducer.pfNoPileUp_cff  import *
from PhysicsTools.PFCandProducer.ParticleSelectors.pfAllNeutralHadrons_cfi  import *
from PhysicsTools.PFCandProducer.ParticleSelectors.pfAllChargedHadrons_cfi import *
from PhysicsTools.PFCandProducer.ParticleSelectors.pfAllPhotons_cfi import *
from PhysicsTools.PFCandProducer.ParticleSelectors.pfAllMuons_cfi import *
from PhysicsTools.PFCandProducer.ParticleSelectors.pfAllElectrons_cfi import *
from PhysicsTools.PFCandProducer.pfMuons_cff import *

CMS2pfIsolatedMuons = pfIsolatedMuons.clone()
CMS2pfIsolatedMuons.src = cms.InputTag("pfAllMuons")
CMS2pfIsolatedMuons.isolationValueMaps = cms.VInputTag(cms.InputTag("CMS2isoValMuonWithCharged"),
                                                       cms.InputTag("CMS2isoValMuonWithNeutral"),
                                                       cms.InputTag("CMS2isoValMuonWithPhotons"))

CMS2pfNoMuon = pfNoMuon.clone()
CMS2pfNoMuon.topCollection = cms.InputTag("CMS2pfIsolatedMuons")

CMS2pfAllElectrons = pfAllElectrons.clone()
CMS2pfAllElectrons.src = cms.InputTag("CMS2pfNoMuon")

cms2PFSequence = cms.Sequence(pfJetMaker + pfmetMaker + pftauMaker + CMS2PFBtagging + bTagPFJetMaker + pfNoPileUpSequence + pfAllNeutralHadrons + pfAllChargedHadrons + pfAllPhotons + pfAllMuons + CMS2pfMuonIsolationSequence + CMS2pfIsolatedMuons + CMS2pfNoMuon + CMS2pfAllElectrons + pfMuonMaker + CMS2pfElectronIsolationSequence + pfElectronMaker + pfElToElAssMaker + elToPFElAssMaker + pfCandidateMaker )

cms2PFNoTauSequence = cms.Sequence(pfJetMaker + pfmetMaker + CMS2PFBtagging + bTagPFJetMaker + pfNoPileUpSequence + pfAllNeutralHadrons + pfAllChargedHadrons + pfAllPhotons + pfAllMuons + CMS2pfMuonIsolationSequence + CMS2pfIsolatedMuons + CMS2pfNoMuon + CMS2pfAllElectrons + pfMuonMaker + pfMuToMuAssMaker + muToPFMuAssMaker + CMS2pfElectronIsolationSequence + pfElectronMaker + pfElToElAssMaker + elToPFElAssMaker + pfCandidateMaker )
