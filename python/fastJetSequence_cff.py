import FWCore.ParameterSet.Config as cms
from CMS3.NtupleMaker.fastJetMaker_cfi import *

# miniAODrhoSequence = cms.Sequence ( fixedGridRhoAllMaker * fixedGridRhoFastJetAllMakerMETTools * fixedGridRhoFastJetAllMaker * fixedGridRhoFastJetAllCaloMaker * fixedGridRhoFastJetCentralCaloMaker * fixedGridRhoFastJetCentralChargedPileUpMaker * fixedGridRhoFastJetCentralNeutralMaker )
miniAODrhoSequence = cms.Sequence ( fixedGridRhoAllMaker * fixedGridRhoFastJetAllMaker * fixedGridRhoFastJetAllCaloMaker * fixedGridRhoFastJetCentralCaloMaker * fixedGridRhoFastJetCentralChargedPileUpMaker * fixedGridRhoFastJetCentralNeutralMaker )
