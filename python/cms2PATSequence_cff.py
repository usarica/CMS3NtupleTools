#Contains the CMS2 PAT Makers.
import FWCore.ParameterSet.Config as cms

from CMS2.NtupleMaker.patElectronMaker_cfi import *
from CMS2.NtupleMaker.patJetMaker_cfi import *
from CMS2.NtupleMaker.patMETMaker_cfi import *
from CMS2.NtupleMaker.patMuonMaker_cfi import *

cms2PATSequence     = cms.Sequence(patMuonMaker * patElectronMaker * patJetMaker * patMETMaker)
