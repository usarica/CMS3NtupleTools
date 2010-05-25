import FWCore.ParameterSet.Config as cms

#isMC = False

#HFPMTcleaningversion = 2   # version 1 = default (loose), version 2 = (medium), version 3 = (tight)
                            # version 4 = version 2 above+timing cut, version 5 = version 3 above+timing cut

# New SeverityLevelComputer that forces RecHits with UserDefinedBit0 set to be excluded from new rechit collection
import JetMETAnalysis.HcalReflagging.RemoveAddSevLevel as cms2RemoveAddSevLevel
from RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi import *

#hcalRecAlgos = cms2RemoveAddSevLevel.RemoveFlag(hcalRecAlgos,"HFLongShort")

# UserDefinedBit0 is used by both the HF and HBHE reflaggers

hcalRecAlgos = cms2RemoveAddSevLevel.AddFlag(hcalRecAlgos,"UserDefinedBit0",10)

# HF RecHit reflagger
#from JetMETAnalysis.HcalReflagging.HFrechitreflaggerJETMET_cff import *
#if HFPMTcleaningversion == 1:
#    cms2hfrecoReflaggedv2 = HFrechitreflaggerJETMETv1.clone()
#elif HFPMTcleaningversion == 2:
#    cms2hfrecoReflaggedv2 = HFrechitreflaggerJETMETv2.clone()
#elif HFPMTcleaningversion == 3:
#    cms2hfrecoReflaggedv2 = HFrechitreflaggerJETMETv3.clone()
#elif HFPMTcleaningversion == 4:
#    if (isMC == False):
#        cms2hfrecoReflaggedv2 = HFrechitreflaggerJETMETv4.clone()
#    else:
#        cms2hfrecoReflaggedv2 = HFrechitreflaggerJETMETv2.clone()
#elif HFPMTcleaningversion == 5:
#    if (isMC == False):
#        cms2hfrecoReflaggedv2 = cms2HFrechitreflaggerJETMETv5.clone()
#    else:
#        cms2hfrecoReflaggedv2 = cms2HFrechitreflaggerJETMETv3.clone()

# HBHE RecHit reflagger
from JetMETAnalysis.HcalReflagging.hbherechitreflaggerJETMET_cfi import *
cms2hbherecoReflagged = hbherechitreflaggerJETMET.clone()
cms2hbherecoReflagged.debug = cms.untracked.int32(0)

# Use the reflagged HF RecHits to make the CaloTowers
from RecoJets.JetProducers. CaloTowerSchemeB_cfi import *
from RecoMET.METProducers.CaloMET_cfi import met
from CMS2.NtupleMaker.deltaMETMaker_cfi import *

cms2towerMakerHCAL = towerMaker.clone()
#cms2towerMakerHCAL.hfInput = cms.InputTag("cms2hfrecoReflaggedv2")
cms2towerMakerHCAL.hbheInput = cms.InputTag("cms2hbherecoReflagged")
cms2towerMakerHCAL.UseEcalTiming = cms.bool(False)
cms2metHCAL = met.clone()
cms2metHCAL.src = cms.InputTag("cms2towerMakerHCAL")
deltaMETMakerHCAL = deltaMETMaker.clone()
deltaMETMakerHCAL.aliasPrefix = cms.untracked.string("evthcal")
deltaMETMakerHCAL.metInputTag_ = cms.InputTag("cms2metHCAL")

#cms2HCALcleaningSequence = cms.Sequence(cms2hfrecoReflaggedv2 * cms2hbherecoReflagged * cms2towerMakerHCAL * cms2metHCAL * deltaMETMakerHCAL)
cms2HCALcleaningSequence = cms.Sequence(cms2hbherecoReflagged * cms2towerMakerHCAL * cms2metHCAL * deltaMETMakerHCAL)
