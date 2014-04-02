import FWCore.ParameterSet.Config as cms

# sequence comes from:
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePFTauID

# tau-tagging general configuration
from RecoTauTag.Configuration.HPSPFTaus_cff import *
from RecoTauTag.Configuration.RecoPFTauTag_cff import *

CMS2PFtautagging = cms.Sequence(recoTauClassicHPSSequence
                                )
#CMS2PFtautagging.remove(kt6PFJetsForRhoComputationVoronoi)
