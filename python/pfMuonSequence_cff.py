import FWCore.ParameterSet.Config as cms

from PhysicsTools.PFCandProducer.Isolation.tools_cfi import *
from PhysicsTools.PFCandProducer.Isolation.pfMuonIsolationFromDeposits_cff import *
from PhysicsTools.PFCandProducer.Isolation.isoValMuonWithCharged_cfi import *
from PhysicsTools.PFCandProducer.Isolation.isoValMuonWithNeutral_cfi import *
from PhysicsTools.PFCandProducer.Isolation.isoValMuonWithPhotons_cfi import *
from CMS2.NtupleMaker.pfMuonMaker_cfi import *

CMS2isoDepMuonWithCharged   = isoDepositReplace('pfAllMuons', 'pfAllChargedHadrons')
CMS2isoDepMuonWithNeutral   = isoDepositReplace('pfAllMuons', 'pfAllNeutralHadrons')
CMS2isoDepMuonWithPhotons   = isoDepositReplace('pfAllMuons', 'pfAllPhotons')
#CMS2isoDepMuonWithElectrons = isoDepositReplace('pfAllMuons', 'CMS2pfAllElectrons')
#CMS2isoDepMuonWithMuons     = isoDepositReplace('pfAllMuons', 'pfAllMuons')

CMS2pfMuonIsoDepositsSequence = cms.Sequence(CMS2isoDepMuonWithCharged
                                             + CMS2isoDepMuonWithNeutral
                                             + CMS2isoDepMuonWithPhotons)
#                                             + CMS2isoDepMuonWithElectrons
#                                             + CMS2isoDepMuonWithMuons)

CMS2isoValMuonWithCharged = isoValMuonWithCharged.clone()
CMS2isoValMuonWithCharged.deposits[0].src = cms.InputTag("CMS2isoDepMuonWithCharged")
CMS2isoValMuonWithCharged.deposits[0].deltaR = cms.double(0.3)
CMS2isoValMuonWithNeutral = isoValMuonWithNeutral.clone()
CMS2isoValMuonWithNeutral.deposits[0].src = cms.InputTag("CMS2isoDepMuonWithNeutral")
CMS2isoValMuonWithNeutral.deposits[0].deltaR = cms.double(0.3)
CMS2isoValMuonWithPhotons = isoValMuonWithPhotons.clone()
CMS2isoValMuonWithPhotons.deposits[0].src = cms.InputTag("CMS2isoDepMuonWithPhotons")
CMS2isoValMuonWithPhotons.deposits[0].deltaR = cms.double(0.3)
#CMS2isoValMuonWithElectrons = isoValMuonWithCharged.clone()
#CMS2isoValMuonWithElectrons.deposits[0].src = cms.InputTag("CMS2isoDepMuonWithElectrons")
#CMS2isoValMuonWithMuons = isoValMuonWithCharged.clone()
#CMS2isoValMuonWithMuons.deposits[0].src = cms.InputTag("CMS2isoDepMuonWithMuons")

CMS2pfMuonIsolationFromDepositsSequence = cms.Sequence(CMS2isoValMuonWithCharged
                                                       + CMS2isoValMuonWithNeutral
                                                       + CMS2isoValMuonWithPhotons)
#                                                       + CMS2isoValMuonWithElectrons
#                                                       + CMS2isoValMuonWithMuons)

CMS2pfMuonIsolationSequence = cms.Sequence(CMS2pfMuonIsoDepositsSequence
                                           + CMS2pfMuonIsolationFromDepositsSequence)
