import FWCore.ParameterSet.Config as cms
from RecoEgamma.EgammaIsolationAlgos.egammaIsolationSequence_cff import *

###
### Isolation
###

from RecoEgamma.EgammaIsolationAlgos.gamTrackExtractorBlocks_cff import *
gamIsoDepositTkCMS2 = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("photons"),
    trackType = cms.string('candidate'),
    MultipleDepositsFlag = cms.bool(False),
    ExtractorPSet = cms.PSet(GamIsoTrackExtractorBlock)
)

from RecoEgamma.EgammaIsolationAlgos.gamEcalExtractorBlocks_cff import *

gamIsoDepositEcalFromHitsCMS2 = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("photons"),
    trackType = cms.string('candidate'),
    MultipleDepositsFlag = cms.bool(False),
    ExtractorPSet = cms.PSet(GamIsoEcalFromHitsExtractorBlock)
)

## in the 22X gamIsoDeposits_cff, the HCAL uses an iso from hits
## extractor; modified by hand below to use towers instead
##
from RecoEgamma.EgammaIsolationAlgos.gamHcalExtractorBlocks_cff import *
gamIsoDepositHcalFromTowersCMS2 = cms.EDProducer("CandIsoDepositProducer",
    src = cms.InputTag("photons"),
    trackType = cms.string('candidate'),
    MultipleDepositsFlag = cms.bool(False),
    ExtractorPSet = cms.PSet(GamIsoHcalFromTowersExtractorBlock)
)

#ValueMap Producers
gamIsoFromDepsTkCMS2 = cms.EDFilter("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        mode = cms.string('sum'),
        src = cms.InputTag("gamIsoDepositTkCMS2"),
        weight = cms.string('1'),
        deltaR = cms.double(0.3),
        vetos = cms.vstring('0.015','Threshold(1.0)'),
        skipDefaultVeto = cms.bool(True)
    ))
)
gamIsoFromDepsEcalFromHitsCMS2 = cms.EDFilter("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        mode = cms.string('sum'),
        src = cms.InputTag("gamIsoDepositEcalFromHitsCMS2"),
        weight = cms.string('1'),
        deltaR = cms.double(0.4),
        vetos = cms.vstring('EcalBarrel:0.045', 
                            'EcalBarrel:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)', 
                            'EcalBarrel:Threshold(0.080)', 
                            'EcalEndcaps:0.070', 
                            'EcalEndcaps:RectangularEtaPhiVeto(-0.02,0.02,-0.5,0.5)', 
                            'EcalEndcaps:Threshold(0.30)'),
        skipDefaultVeto = cms.bool(True)
    ))
)
gamIsoFromDepsHcalFromTowersCMS2 = cms.EDFilter("CandIsolatorFromDeposits",
    deposits = cms.VPSet(cms.PSet(
        src = cms.InputTag("gamIsoDepositHcalFromTowersCMS2"),
        deltaR = cms.double(0.4),
        weight = cms.string('1'),
        vetos = cms.vstring('0.00'),
        skipDefaultVeto = cms.bool(True),
        mode = cms.string('sum')
    ))
)

###
### Sequences
###

### Isolation sequences
gamIsoDepositsCMS2 = cms.Sequence(gamIsoDepositTkCMS2+gamIsoDepositEcalFromHitsCMS2+gamIsoDepositHcalFromTowersCMS2)
gamIsoFromDepositsCMS2 = cms.Sequence(gamIsoFromDepsTkCMS2*gamIsoFromDepsEcalFromHitsCMS2*gamIsoFromDepsHcalFromTowersCMS2)

### Master sequence
gammaSequence = cms.Sequence(gamIsoDepositsCMS2 + gamIsoFromDepositsCMS2)
