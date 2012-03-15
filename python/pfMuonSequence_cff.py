import FWCore.ParameterSet.Config as cms

from CommonTools.ParticleFlow.Isolation.tools_cfi import *
#from CommonTools.ParticleFlow.Isolation.pfMuonIsolationFromDeposits_cff import *
from CommonTools.ParticleFlow.Isolation.isoValMuonWithCharged_cfi import *
from CommonTools.ParticleFlow.Isolation.isoValMuonWithNeutral_cfi import *
from CommonTools.ParticleFlow.Isolation.isoValMuonWithPhotons_cfi import *
from CMS2.NtupleMaker.pfMuonMaker_cfi import *

CMS2isoDepMuonWithCharged   = isoDepositReplace('pfAllMuons', 'pfAllChargedHadrons')
CMS2isoDepMuonWithNeutral   = isoDepositReplace('pfAllMuons', 'pfAllNeutralHadrons')
CMS2isoDepMuonWithPhotons   = isoDepositReplace('pfAllMuons', 'pfAllPhotons')

CMS2pfMuonIsoDepositsSequence = cms.Sequence(CMS2isoDepMuonWithCharged
                                             + CMS2isoDepMuonWithNeutral
                                             + CMS2isoDepMuonWithPhotons)

#
CMS2isoValMuonWithCharged = isoValMuonWithCharged.clone()
CMS2isoValMuonWithCharged.deposits[0].src = cms.InputTag("CMS2isoDepMuonWithCharged")
CMS2isoValMuonWithCharged.deposits[0].deltaR = cms.double(0.3)
CMS2isoValMuonWithNeutral = isoValMuonWithNeutral.clone()
CMS2isoValMuonWithNeutral.deposits[0].src = cms.InputTag("CMS2isoDepMuonWithNeutral")
CMS2isoValMuonWithNeutral.deposits[0].deltaR = cms.double(0.3)
CMS2isoValMuonWithPhotons = isoValMuonWithPhotons.clone()
CMS2isoValMuonWithPhotons.deposits[0].src = cms.InputTag("CMS2isoDepMuonWithPhotons")
CMS2isoValMuonWithPhotons.deposits[0].deltaR = cms.double(0.3)

#same, cone 0.4
CMS2isoValMuonWithCharged04 = isoValMuonWithCharged.clone()
CMS2isoValMuonWithCharged04.deposits[0].src = cms.InputTag("CMS2isoDepMuonWithCharged")
CMS2isoValMuonWithCharged04.deposits[0].deltaR = cms.double(0.4)
CMS2isoValMuonWithNeutral04 = isoValMuonWithNeutral.clone()
CMS2isoValMuonWithNeutral04.deposits[0].src = cms.InputTag("CMS2isoDepMuonWithNeutral")
CMS2isoValMuonWithNeutral04.deposits[0].deltaR = cms.double(0.4)
CMS2isoValMuonWithPhotons04 = isoValMuonWithPhotons.clone()
CMS2isoValMuonWithPhotons04.deposits[0].src = cms.InputTag("CMS2isoDepMuonWithPhotons")
CMS2isoValMuonWithPhotons04.deposits[0].deltaR = cms.double(0.4)

CMS2pfMuonIsolationFromDepositsSequence = cms.Sequence(CMS2isoValMuonWithCharged
                                                       + CMS2isoValMuonWithNeutral
                                                       + CMS2isoValMuonWithPhotons
                                                       + CMS2isoValMuonWithCharged04
                                                       + CMS2isoValMuonWithNeutral04
                                                       + CMS2isoValMuonWithPhotons04
                                                       )

CMS2pfMuonIsolationSequence = cms.Sequence(CMS2pfMuonIsoDepositsSequence
                                           + CMS2pfMuonIsolationFromDepositsSequence)
