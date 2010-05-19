import FWCore.ParameterSet.Config as cms

isMC = False
#useHBHEcleaning = True
#useHBHEfiltering = False

HFPMTcleaningversion = 4   # version 1 = default (loose), version 2 = (medium), version 3 = (tight)
                            # version 4 = version 2 above+timing cut, version 5 = version 3 above+timing cut

#if useHBHEfiltering == True:
#    from CommonTools.RecoAlgos.HBHENoiseFilter_cfi import *
#    cms2HBHENoiseFilter = HBHENoiseFilter.clone()

# New SeverityLevelComputer that forces RecHits with UserDefinedBit0 set to be excluded from new rechit collection
import JetMETAnalysis.HcalReflagging.RemoveAddSevLevel as cms2RemoveAddSevLevel
from RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi import *
#cms2hcalRecAlgos = hcalRecAlgos.clone()
#cms2hcalRecAlgos = cms2RemoveAddSevLevel.RemoveFlag(cms2hcalRecAlgos,"HFLongShort")
hcalRecAlgos = cms2RemoveAddSevLevel.RemoveFlag(hcalRecAlgos,"HFLongShort")

# UserDefinedBit0 is used by both the HF and HBHE reflaggers
#cms2hcalRecAlgos = cms2RemoveAddSevLevel.AddFlag(cms2hcalRecAlgos,"UserDefinedBit0",10)
hcalRecAlgos = cms2RemoveAddSevLevel.AddFlag(hcalRecAlgos,"UserDefinedBit0",10)

# HF RecHit reflagger
from JetMETAnalysis.HcalReflagging.HFrechitreflaggerJETMET_cff import *
if HFPMTcleaningversion == 1:
    cms2hfrecoReflagged = HFrechitreflaggerJETMETv1.clone()
elif HFPMTcleaningversion == 2:
    cms2hfrecoReflagged = HFrechitreflaggerJETMETv2.clone()
elif HFPMTcleaningversion == 3:
    cms2hfrecoReflagged = HFrechitreflaggerJETMETv3.clone()
elif HFPMTcleaningversion == 4:
    if (isMC == False):
        cms2hfrecoReflagged = HFrechitreflaggerJETMETv4.clone()
    else:
        cms2hfrecoReflagged = HFrechitreflaggerJETMETv2.clone()
elif HFPMTcleaningversion == 5:
    if (isMC == False):
        cms2hfrecoReflagged = cms2HFrechitreflaggerJETMETv5.clone()
    else:
        cms2hfrecoReflagged = cms2HFrechitreflaggerJETMETv3.clone()

# HBHE RecHit reflagger
from JetMETAnalysis.HcalReflagging.hbherechitreflaggerJETMET_cfi import *
cms2hbherecoReflagged = hbherechitreflaggerJETMET.clone()
cms2hbherecoReflagged.debug = cms.untracked.int32(0)

# Use the reflagged HF RecHits to make the CaloTowers
from RecoJets.JetProducers. CaloTowerSchemeB_cfi import *
from RecoMET.METProducers.CaloMET_cfi import met
from CMS2.NtupleMaker.deltaMETMaker_cfi import *

cms2towerMakerHF = towerMaker.clone()
cms2towerMakerHF.hfInput = cms.InputTag("cms2hfrecoReflagged")
cms2metHF = met.clone()
cms2metHF.src = cms.InputTag("cms2towerMakerHF")
deltaMETMakerHF = deltaMETMaker.clone()
deltaMETMakerHF.aliasPrefix = cms.untracked.string("evthf")
deltaMETMakerHF.metInputTag_ = cms.InputTag("cms2metHF")

cms2towerMakerHCAL = towerMaker.clone()
cms2towerMakerHCAL.hbheInput = cms.InputTag("cms2hbherecoReflagged")
cms2metHCAL = met.clone()
cms2metHCAL.src = cms.InputTag("cms2towerMakerHCAL")
deltaMETMakerHCAL = deltaMETMaker.clone()
deltaMETMakerHCAL.aliasPrefix = cms.untracked.string("evthcal")
deltaMETMakerHCAL.metInputTag_ = cms.InputTag("cms2metHCAL")

# do the same thing for ECAL timing
cms2towerMakerECAL = towerMaker.clone()
cms2towerMakerECAL.UseEcalTiming = cms.bool(True)
cms2metECAL = met.clone()
cms2metECAL.src = cms.InputTag("cms2towerMakerECAL")
deltaMETMakerECAL = deltaMETMaker.clone()
deltaMETMakerECAL.aliasPrefix = cms.untracked.string("evtecal")
deltaMETMakerECAL.metInputTag_ = cms.InputTag("cms2metECAL")

cms2HFcleaningSequence   = cms.Sequence(cms2hfrecoReflagged * cms2towerMakerHF * cms2metHF * deltaMETMakerHF)
cms2HCALcleaningSequence = cms.Sequence(cms2hbherecoReflagged * cms2towerMakerHCAL * cms2metHCAL * deltaMETMakerHCAL)
cms2ECALcleaningSequence = cms.Sequence(cms2towerMakerECAL * cms2metECAL * deltaMETMakerECAL)
cms2CleaningSequence     = cms.Sequence(cms2HFcleaningSequence * cms2HCALcleaningSequence * cms2ECALcleaningSequence)
