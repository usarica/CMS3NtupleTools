import FWCore.ParameterSet.Config as cms

isMC = False
version = 10


## HF RecHit reflagger -- specify type of HF cleaning to use
#from JetMETAnalysis.HcalReflagging.HFrechitreflaggerJETMET_cff import *
#if version==1:
#    cms2hfrecoReflagged = HFrechitreflaggerJETMETv1.clone()
#elif version==2:
#    cms2hfrecoReflagged = HFrechitreflaggerJETMETv2.clone()
#elif version==3:
#    cms2hfrecoReflagged = HFrechitreflaggerJETMETv3.clone()
#elif version==4: 
#    if (isMC==False):
#        cms2hfrecoReflagged = HFrechitreflaggerJETMETv4.clone()
#    else:
#        cms2hfrecoReflagged = HFrechitreflaggerJETMETv2.clone()  
#elif version==5:
#    if (isMC==False):
#        cms2hfrecoReflagged = HFrechitreflaggerJETMETv5.clone()
#    else:
#        cms2hfrecoReflagged = HFrechitreflaggerJETMETv3.clone()  
#
## CURRENT RECOMMENDATION
#elif version==10:
#    cms2hfrecoReflagged = HFrechitreflaggerJETMETv10.clone()
#    if (isMC==False):  # V10 cleaning uses results of prior flags when setting new flags; this is the current recommendation as of 21 July 2010 
#        cms2hfrecoReflagged.PETstat.flagsToSkip =string.atoi('10',2)
#        cms2hfrecoReflagged.S8S1stat.flagsToSkip=string.atoi('10010',2)
#        cms2hfrecoReflagged.S9S1stat.flagsToSkip=string.atoi('11010',2)
#        # Flag ordering
#        cms2hfrecoReflagged.FlagsToSet=(4,3,0)  # set flag 4 (HFPET -- also sets HFLongShort), then flag 3 (HFS8S1 -- also sets HFLongShort), then flag 0 (HFLongShort -- set directly via S9S1)



#import JetMETAnalysis.HcalReflagging.RemoveAddSevLevel as cms2RemoveAddSevLevel
from RecoLocalCalo.HcalRecAlgos.hcalRecAlgoESProd_cfi import *
#hcalRecAlgos=cms2RemoveAddSevLevel.AddFlag(hcalRecAlgos,"HFLongShort",11)

#if (isMC==False):  # Don't use HFDigiTime on MC !
#    hcalRecAlgos=cms2RemoveAddSevLevel.AddFlag(hcalRecAlgos,"HFDigiTime",11)



# Use the reflagged HF RecHits to make the CaloTowers
from RecoJets.JetProducers. CaloTowerSchemeB_cfi import *
from RecoMET.METProducers.CaloMET_cfi import met
from CMS2.NtupleMaker.deltaMETMaker_cfi import *



cms2towerMakerHF = towerMaker.clone()
cms2towerMakerHF.hfInput = cms.InputTag("cms2hfrecoReflagged")
cms2towerMakerHF.UseEcalTiming = cms.bool(False)
cms2metHF = met.clone()
cms2metHF.src = cms.InputTag("cms2towerMakerHF")
deltaMETMakerHF = deltaMETMaker.clone()
deltaMETMakerHF.aliasPrefix = cms.untracked.string("evthf")
deltaMETMakerHF.metInputTag_ = cms.InputTag("cms2metHF")

#cms2HFcleaningSequence   = cms.Sequence(cms2hfrecoReflagged * cms2towerMakerHF * cms2metHF * deltaMETMakerHF)
cms2HFcleaningSequence   = cms.Sequence(cms2towerMakerHF * cms2metHF * deltaMETMakerHF)

