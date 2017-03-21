import FWCore.ParameterSet.Config as cms
from CMS3.NtupleMaker.fastJetMaker_cfi import *

miniAODrhoSequence = cms.Sequence ( fixedGridRhoFastJetCentralCaloMaker )
