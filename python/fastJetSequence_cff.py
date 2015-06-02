import FWCore.ParameterSet.Config as cms

from RecoJets.Configuration.RecoPFJets_cff import *
from RecoJets.Configuration.RecoJets_cff import *
from CMS3.NtupleMaker.fastJetMaker_cfi import *

miniAODrhoSequence = cms.Sequence ( fixedGridRhoAllMaker * fixedGridRhoFastJetAllMaker * fixedGridRhoFastJetAllCaloMaker * fixedGridRhoFastJetCentralCaloMaker * fixedGridRhoFastJetCentralChargedPileUpMaker * fixedGridRhoFastJetCentralNeutralMaker )
