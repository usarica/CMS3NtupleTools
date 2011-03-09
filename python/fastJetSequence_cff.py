import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.kt4PFJets_cfi import *
from CMS2.NtupleMaker.fastJetMaker_cfi import *

mykt6PFJets = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
mykt6PFJets.Rho_EtaMax = cms.double(2.5)
mykt6PFJets.Ghost_EtaMax = cms.double(2.5)
fastJetSequence = cms.Sequence( mykt6PFJets * fastJetMaker )
