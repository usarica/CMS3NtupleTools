import FWCore.ParameterSet.Config as cms

# New SeverityLevelComputer that forces RecHits with UserDefinedBit0 set to be excluded from new rechit collection
import JetMETAnalysis.HcalReflagging.RemoveAddSevLevel as cms2RemoveAddSevLevel
from RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi import *
from JetMETAnalysis.HcalReflagging.hbherechitreflaggerJETMET_cfi import *


hcalRecAlgos= cms2RemoveAddSevLevel.AddFlag(hcalRecAlgos,"UserDefinedBit0",10)

# HBHE RecHit reflagger
cms2hbherecoReflagged = hbherechitreflaggerJETMET.clone()
cms2hbherecoReflagged.debug=0


# Use the reflagged HF RecHits to make the CaloTowers
from RecoJets.JetProducers. CaloTowerSchemeB_cfi import *
from RecoMET.METProducers.CaloMET_cfi import met
from CMS2.NtupleMaker.deltaMETMaker_cfi import *

cms2towerMakerHCAL = towerMaker.clone()
cms2towerMakerHCAL.hbheInput = cms.InputTag("cms2hbherecoReflagged")
cms2towerMakerHCAL.UseEcalTiming = cms.bool(False)
cms2metHCAL = met.clone()
cms2metHCAL.src = cms.InputTag("cms2towerMakerHCAL")
deltaMETMakerHCAL = deltaMETMaker.clone()
deltaMETMakerHCAL.aliasPrefix = cms.untracked.string("evthcal")
deltaMETMakerHCAL.metInputTag_ = cms.InputTag("cms2metHCAL")


cms2HCALcleaningSequence = cms.Sequence(cms2hbherecoReflagged * cms2towerMakerHCAL * cms2metHCAL * deltaMETMakerHCAL)
