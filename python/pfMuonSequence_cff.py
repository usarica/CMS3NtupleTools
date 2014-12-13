import FWCore.ParameterSet.Config as cms

from RecoMuon.MuonIsolation.muonPFIsolationValues_cff import *

from CommonTools.ParticleFlow.Isolation.tools_cfi import *
#from CommonTools.ParticleFlow.Isolation.pfMuonIsolationFromDeposits_cff import *
#from CommonTools.ParticleFlow.Isolation.isoValMuonWithCharged_cfi import *
#from CommonTools.ParticleFlow.Isolation.isoValMuonWithNeutral_cfi import *
#from CommonTools.ParticleFlow.Isolation.isoValMuonWithPhotons_cfi import *
from CMS3.NtupleMaker.pfMuonMaker_cfi import *

CMS2isoDepMuonWithCharged   = isoDepositReplace('pfAllMuons', 'pfAllChargedHadrons')
CMS2isoDepMuonWithNeutral   = isoDepositReplace('pfAllMuons', 'pfAllNeutralHadrons')
CMS2isoDepMuonWithPhotons   = isoDepositReplace('pfAllMuons', 'pfAllPhotons')

CMS2pfMuonIsoDepositsSequence = cms.Sequence(
  CMS2isoDepMuonWithCharged + 
  CMS2isoDepMuonWithNeutral + 
  CMS2isoDepMuonWithPhotons
)


# cone 0.3
CMS2isoValMuonWithCharged03                 = muPFIsoValueCharged03.clone()
CMS2isoValMuonWithCharged03.deposits[0].src = cms.InputTag("CMS2isoDepMuonWithCharged")

CMS2isoValMuonWithNeutral03                 = muPFIsoValueNeutral03.clone()
CMS2isoValMuonWithNeutral03.deposits[0].src = cms.InputTag("CMS2isoDepMuonWithNeutral")

CMS2isoValMuonWithPhotons03                 = muPFIsoValueGamma03.clone()
CMS2isoValMuonWithPhotons03.deposits[0].src = cms.InputTag("CMS2isoDepMuonWithPhotons")


## cone 0.4
CMS2isoValMuonWithCharged04                 = muPFIsoValueCharged04.clone()
CMS2isoValMuonWithCharged04.deposits[0].src = cms.InputTag("CMS2isoDepMuonWithCharged")

CMS2isoValMuonWithNeutral04                 = muPFIsoValueNeutral04.clone()
CMS2isoValMuonWithNeutral04.deposits[0].src = cms.InputTag("CMS2isoDepMuonWithNeutral")

CMS2isoValMuonWithPhotons04                 = muPFIsoValueGamma04.clone()
CMS2isoValMuonWithPhotons04.deposits[0].src = cms.InputTag("CMS2isoDepMuonWithPhotons")


CMS2pfMuonIsolationFromDepositsSequence = cms.Sequence(
  CMS2isoValMuonWithCharged03 +
  CMS2isoValMuonWithNeutral03 +
  CMS2isoValMuonWithPhotons03 +
  CMS2isoValMuonWithCharged04 +
  CMS2isoValMuonWithNeutral04 +
  CMS2isoValMuonWithPhotons04
)

CMS2pfMuonIsolationSequence = cms.Sequence(
  CMS2pfMuonIsoDepositsSequence +
  CMS2pfMuonIsolationFromDepositsSequence
)
