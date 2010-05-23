import FWCore.ParameterSet.Config as cms

isMC = False

HFPMTcleaningversion = 2   # version 1 = default (loose), version 2 = (medium), version 3 = (tight)
                            # version 4 = version 2 above+timing cut, version 5 = version 3 above+timing cut

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
    cms2hfrecoReflaggedV2 = HFrechitreflaggerJETMETv1.clone()
elif HFPMTcleaningversion == 2:
    cms2hfrecoReflaggedV2 = HFrechitreflaggerJETMETv2.clone()
elif HFPMTcleaningversion == 3:
    cms2hfrecoReflaggedV2 = HFrechitreflaggerJETMETv3.clone()
elif HFPMTcleaningversion == 4:
    if (isMC == False):
        cms2hfrecoReflaggedV2 = HFrechitreflaggerJETMETv4.clone()
    else:
        cms2hfrecoReflaggedV2 = HFrechitreflaggerJETMETv2.clone()
elif HFPMTcleaningversion == 5:
    if (isMC == False):
        cms2hfrecoReflaggedV2 = cms2HFrechitreflaggerJETMETv5.clone()
    else:
        cms2hfrecoReflaggedV2 = cms2HFrechitreflaggerJETMETv3.clone()

# Use the reflagged HF RecHits to make the CaloTowers
from RecoJets.JetProducers. CaloTowerSchemeB_cfi import *
from RecoMET.METProducers.CaloMET_cfi import met
from CMS2.NtupleMaker.deltaMETMaker_cfi import *

cms2towerMakerECAL = towerMaker.clone()
cms2towerMakerECAL.hfInput = cms.InputTag("cms2hfrecoReflaggedV2")
cms2towerMakerECAL.UseEcalTiming = cms.bool(True)
cms2metECAL = met.clone()
cms2metECAL.src = cms.InputTag("cms2towerMakerECAL")
deltaMETMakerECAL = deltaMETMaker.clone()
deltaMETMakerECAL.aliasPrefix = cms.untracked.string("evtecal")
deltaMETMakerECAL.metInputTag_ = cms.InputTag("cms2metECAL")

cms2ECALcleaningSequence = cms.Sequence(cms2hfrecoReflaggedV2 * cms2towerMakerECAL * cms2metECAL * deltaMETMakerECAL)

