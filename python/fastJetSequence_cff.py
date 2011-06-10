import FWCore.ParameterSet.Config as cms

from RecoJets.Configuration.RecoPFJets_cff import *
from CMS2.NtupleMaker.fastJetMaker_cfi import *

kt6PFJetsForRhoComputation = kt6PFJets.clone(doRhoFastjet = True)

fastJetSequence = cms.Sequence( kt6PFJets * kt6PFJetsForRhoComputation * fastJetMaker )
