import FWCore.ParameterSet.Config as cms

from CMS3.NtupleMaker.caloTowerMaker_cfi import *

cms2CaloTowerSequence = cms.Sequence( caloTowerMaker  )

