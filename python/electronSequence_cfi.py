import FWCore.ParameterSet.Config as cms

from RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequence_cff import *
from CMS2.NtupleMaker.electronDuplicateRemover_cfi import *
from CMS2.NtupleMaker.electronMaker_cfi import *

import FWCore.ParameterSet.Config as cms

from RecoEgamma.EgammaIsolationAlgos.eleTrackExtractorBlocks_cff import *

eleIsoDepositTkCMS2 = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("uniqueElectrons"),
    trackType = cms.string('candidate'),
    MultipleDepositsFlag = cms.bool(False),
    ExtractorPSet = cms.PSet(EleIsoTrackExtractorBlock)
)

from RecoEgamma.EgammaIsolationAlgos.eleEcalExtractorBlocks_cff import *

eleIsoDepositEcalFromHitsCMS2 = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("uniqueElectrons"),
    trackType = cms.string('candidate'),
    MultipleDepositsFlag = cms.bool(False),
    ExtractorPSet = cms.PSet(EleIsoEcalFromHitsExtractorBlock)
)


##in 21X, there doesn't seem to be an extractor from towers
from RecoEgamma.EgammaIsolationAlgos.eleHcalExtractorBlocks_cff import *
eleIsoDepositHcalFromTowersCMS2 = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("uniqueElectrons"),
    trackType = cms.string('candidate'),
    MultipleDepositsFlag = cms.bool(False),
    ExtractorPSet = cms.PSet(EleIsoHcalFromTowersExtractorBlock)
)

eleIsoDepositsCMS2 = cms.Sequence(eleIsoDepositTkCMS2+eleIsoDepositEcalFromHitsCMS2+eleIsoDepositHcalFromTowersCMS2)


#ValueMap Producers
eleIsoFromDepsTkCMS2 = cms.EDFilter("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        mode = cms.string('sum'),
        src = cms.InputTag("eleIsoDepositTkCMS2"),
        weight = cms.string('1'),
        deltaR = cms.double(0.3),
        vetos = cms.vstring('0.015','Threshold(1.0)'),
        skipDefaultVeto = cms.bool(True)
    ))
)
eleIsoFromDepsEcalFromHitsCMS2 = cms.EDFilter("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        mode = cms.string('sum'),
        src = cms.InputTag("eleIsoDepositEcalFromHitsCMS2"),
        weight = cms.string('1'),
        deltaR = cms.double(0.4),
        vetos = cms.vstring('EcalBarrel:0.045', 
                            'EcalBarrel:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)',
                            'EcalBarrel:ThresholdFromTransverse(0.08)',
                            'EcalEndcaps:ThresholdFromTransverse(0.3)',
                            'EcalEndcaps:0.070', 
                            'EcalEndcaps:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)'), 
        skipDefaultVeto = cms.bool(True)
    ))
)
eleIsoFromDepsHcalFromTowersCMS2 = cms.EDFilter("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("eleIsoDepositHcalFromTowersCMS2"),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('0.0'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)



eleIsoFromDepositsCMS2 = cms.Sequence(eleIsoFromDepsTkCMS2*eleIsoFromDepsEcalFromHitsCMS2*eleIsoFromDepsHcalFromTowersCMS2)
                                       
egammaIsolationSequenceCMS2 = cms.Sequence(eleIsoDepositsCMS2*eleIsoFromDepositsCMS2)

electronSequence = cms.Sequence(uniqueElectrons*egammaIsolationSequenceCMS2*electronMaker)
