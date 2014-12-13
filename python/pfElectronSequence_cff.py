import FWCore.ParameterSet.Config as cms

from CommonTools.ParticleFlow.Isolation.tools_cfi import *
#from CommonTools.ParticleFlow.Isolation.pfElectronIsolationFromDeposits_cff import *
#from CommonTools.ParticleFlow.Isolation.isoValElectronWithCharged_cfi import *
#from CommonTools.ParticleFlow.Isolation.isoValElectronWithNeutral_cfi import *
#from CommonTools.ParticleFlow.Isolation.isoValElectronWithPhotons_cfi import *
from CMS3.NtupleMaker.pfElectronMaker_cfi import *

CMS2isoDepElectronWithCharged   = isoDepositReplace('CMS2pfAllElectrons', 'pfAllChargedHadrons')
CMS2isoDepElectronWithNeutral   = isoDepositReplace('CMS2pfAllElectrons', 'pfAllNeutralHadrons')
CMS2isoDepElectronWithPhotons   = isoDepositReplace('CMS2pfAllElectrons', 'pfAllPhotons')

CMS2pfElectronIsoDepositsSequence = cms.Sequence( 
  CMS2isoDepElectronWithCharged +
  CMS2isoDepElectronWithNeutral + 
  CMS2isoDepElectronWithPhotons
)

isoValElectronWithCharged = cms.EDProducer(
  "CandIsolatorFromDeposits",
  deposits = cms.VPSet(
    cms.PSet(
      src = cms.InputTag("isoDepElectronWithCharged"),
      deltaR = cms.double(0.4),
      weight = cms.string('1'),
      vetos = cms.vstring(),
      skipDefaultVeto = cms.bool(True),
      mode = cms.string('sum')
    )
  )
)

isoValElectronWithNeutral = cms.EDProducer(
  "CandIsolatorFromDeposits",
  deposits = cms.VPSet(
    cms.PSet(
      src = cms.InputTag("isoDepElectronWithNeutral"),
      deltaR = cms.double(0.4),
      weight = cms.string('1'), # 0.3333,
      vetos = cms.vstring('Threshold(0.5)'),
      skipDefaultVeto = cms.bool(True),
      mode = cms.string('sum')
    )
  )
)

isoValElectronWithPhotons = cms.EDProducer(
  "CandIsolatorFromDeposits",
  deposits = cms.VPSet(
    cms.PSet(
      src = cms.InputTag("isoDepElectronWithPhotons"),
      deltaR = cms.double(0.4),
      weight = cms.string('1'),
      vetos = cms.vstring('Threshold(0.5)'),
      skipDefaultVeto = cms.bool(True),
      mode = cms.string('sum')
    )
  )
)


#
CMS2isoValElectronWithCharged03                    = isoValElectronWithCharged.clone()
CMS2isoValElectronWithCharged03.deposits[0].src    = cms.InputTag("CMS2isoDepElectronWithCharged")
CMS2isoValElectronWithCharged03.deposits[0].deltaR = cms.double(0.3)
CMS2isoValElectronWithNeutral03                    = isoValElectronWithNeutral.clone()
CMS2isoValElectronWithNeutral03.deposits[0].src    = cms.InputTag("CMS2isoDepElectronWithNeutral")
CMS2isoValElectronWithNeutral03.deposits[0].deltaR = cms.double(0.3)
CMS2isoValElectronWithPhotons03                    = isoValElectronWithPhotons.clone()
CMS2isoValElectronWithPhotons03.deposits[0].src    = cms.InputTag("CMS2isoDepElectronWithPhotons")
CMS2isoValElectronWithPhotons03.deposits[0].deltaR = cms.double(0.3)

#
CMS2isoValElectronWithCharged04                    = isoValElectronWithCharged.clone()
CMS2isoValElectronWithCharged04.deposits[0].src    = cms.InputTag("CMS2isoDepElectronWithCharged")
CMS2isoValElectronWithCharged04.deposits[0].deltaR = cms.double(0.4)
CMS2isoValElectronWithNeutral04                    = isoValElectronWithNeutral.clone()
CMS2isoValElectronWithNeutral04.deposits[0].src    = cms.InputTag("CMS2isoDepElectronWithNeutral")
CMS2isoValElectronWithNeutral04.deposits[0].deltaR = cms.double(0.4)
CMS2isoValElectronWithPhotons04                    = isoValElectronWithPhotons.clone()
CMS2isoValElectronWithPhotons04.deposits[0].src    = cms.InputTag("CMS2isoDepElectronWithPhotons")
CMS2isoValElectronWithPhotons04.deposits[0].deltaR = cms.double(0.4)

CMS2pfElectronIsolationFromDepositsSequence = cms.Sequence(
  CMS2isoValElectronWithCharged03 +
  CMS2isoValElectronWithNeutral03 +
  CMS2isoValElectronWithPhotons03 +
  CMS2isoValElectronWithCharged04 +
  CMS2isoValElectronWithNeutral04 +
  CMS2isoValElectronWithPhotons04
)

CMS2pfElectronIsolationSequence = cms.Sequence(
  CMS2pfElectronIsoDepositsSequence +
  CMS2pfElectronIsolationFromDepositsSequence
)
