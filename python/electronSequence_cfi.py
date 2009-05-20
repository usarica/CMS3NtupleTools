import FWCore.ParameterSet.Config as cms

from RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequence_cff import *
from CMS2.NtupleMaker.electronDuplicateRemover_cfi import *
from CMS2.NtupleMaker.electronMaker_cfi import *
from RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi import *

import FWCore.ParameterSet.Config as cms


###
### Isolation
###

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

###
### Electron ID
###

import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
eidRobustLooseCMS2 = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
eidRobustLooseCMS2.src = "uniqueElectrons"

import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
eidRobustTightCMS2 = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
eidRobustTightCMS2.src = "uniqueElectrons"
eidRobustTightCMS2.robustEleIDCuts.barrel = [0.015, 0.0092, 0.020, 0.0025]
eidRobustTightCMS2.robustEleIDCuts.endcap = [0.018, 0.025, 0.020, 0.0040]

import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
eidRobustHighEnergyCMS2 = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
eidRobustHighEnergyCMS2.src = "uniqueElectrons"
eidRobustHighEnergyCMS2.robustEleIDCuts.barrel = [0.050, 0.011, 0.090, 0.005]
eidRobustHighEnergyCMS2.robustEleIDCuts.endcap = [0.100, 0.0275, 0.090, 0.007]

import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
eidLooseCMS2 = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
eidLooseCMS2.src = "uniqueElectrons"
eidLooseCMS2.electronQuality = 'loose'

import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
eidTightCMS2 = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
eidTightCMS2.src = "uniqueElectrons"
eidTightCMS2.electronQuality = 'tight'

###
### Sequences
###

### Isolation sequences
eleIsoDepositsCMS2 = cms.Sequence(eleIsoDepositTkCMS2+eleIsoDepositEcalFromHitsCMS2+eleIsoDepositHcalFromTowersCMS2)
eleIsoFromDepositsCMS2 = cms.Sequence(eleIsoFromDepsTkCMS2*eleIsoFromDepsEcalFromHitsCMS2*eleIsoFromDepsHcalFromTowersCMS2)
egammaIsolationSequenceCMS2 = cms.Sequence(eleIsoDepositsCMS2*eleIsoFromDepositsCMS2)

### Electron ID sequences
egammaElectronIDCMS2 = cms.Sequence(eidRobustLooseCMS2+eidRobustTightCMS2+eidRobustHighEnergyCMS2+eidLooseCMS2+eidTightCMS2)

### Master sequence
electronSequence = cms.Sequence(uniqueElectrons*egammaIsolationSequenceCMS2*egammaElectronIDCMS2*electronMaker)

