import FWCore.ParameterSet.Config as cms

######################
# Import MET Filters #
######################q

from RecoMET.METFilters.EcalDeadCellBoundaryEnergyFilter_cfi   import *
from RecoMET.METFilters.EcalDeadCellDeltaRFilter_cfi           import *
from RecoMET.METFilters.EcalDeadCellTriggerPrimitiveFilter_cfi import *
from RecoMET.METFilters.greedyMuonPFCandidateFilter_cfi        import *
from RecoMET.METFilters.hcalLaserEventFilter_cfi               import *
from RecoMET.METFilters.inconsistentMuonPFCandidateFilter_cfi  import *
from RecoMET.METFilters.jetIDFailureFilter_cfi                 import *
from RecoMET.METFilters.multiEventFilter_cfi                   import *
from RecoMET.METFilters.trackingFailureFilter_cfi              import *
from RecoMET.METFilters.eeBadScFilter_cfi                      import *
from RecoMET.METFilters.ecalLaserCorrFilter_cfi                import *

##################
# Clone Defaults #
##################

cms2EcalDeadCellBoundaryEnergyFilter   = EcalDeadCellBoundaryEnergyFilter.clone()
cms2EcalDeadCellDeltaRFilter           = simpleDRfilter.clone()
cms2EcalDeadCellTriggerPrimitiveFilter = EcalDeadCellTriggerPrimitiveFilter.clone()
cms2greedyMuonPFCandidateFilter        = greedyMuonPFCandidateFilter.clone()
cms2hcalLaserEventFilter               = hcalLaserEventFilter.clone()
cms2inconsistentMuonPFCandidateFilter  = inconsistentMuonPFCandidateFilter.clone()
cms2jetIDFailureFilter                 = jetIDFailure.clone()
cms2multiEventFailureFilter            = multiEventFilter.clone()
cms2trackingFailureFilter              = trackingFailureFilter.clone()
cms2eeBadScFilter                      = eeBadScFilter.clone()
cms2ecalLaserCorrFilter                = ecalLaserCorrFilter.clone()

####################
# Set Tagging Mode #
####################

cms2EcalDeadCellBoundaryEnergyFilter   .taggingMode = cms.bool(True)
cms2EcalDeadCellDeltaRFilter           .taggingMode = cms.bool(True)
cms2EcalDeadCellTriggerPrimitiveFilter .taggingMode = cms.bool(True)
cms2greedyMuonPFCandidateFilter        .taggingMode = cms.bool(True)
cms2hcalLaserEventFilter               .taggingMode = cms.bool(True)
cms2inconsistentMuonPFCandidateFilter  .taggingMode = cms.bool(True)
cms2jetIDFailureFilter                 .taggingMode = cms.bool(True)
cms2multiEventFailureFilter            .taggingMode = cms.bool(True)
cms2trackingFailureFilter              .taggingMode = cms.bool(True)
cms2eeBadScFilter                      .taggingMode = cms.bool(True)
cms2ecalLaserCorrFilter                .taggingMode = cms.bool(True)

cms2hcalLaserEventFilter.vetoByRunEventNumber = cms.untracked.bool(False)
cms2hcalLaserEventFilter.vetoByHBHEOccupancy  = cms.untracked.bool(True)

########################################
# Vertices for Tracking Failure Filter #
########################################

cms2trackingFailureFilter.VertexSource = cms.InputTag("offlinePrimaryVertices")


##################################
# Jets for jet ID Failure Filter #
##################################

#

#
cms2MetFilterSequence = cms.Sequence( 
  cms2EcalDeadCellBoundaryEnergyFilter * 
  cms2EcalDeadCellDeltaRFilter * 
  cms2EcalDeadCellTriggerPrimitiveFilter *
  cms2greedyMuonPFCandidateFilter *
  cms2hcalLaserEventFilter *      
  cms2inconsistentMuonPFCandidateFilter *
  #cms2jetIDFailureFilter *
  cms2multiEventFailureFilter *
  cms2trackingFailureFilter *
  cms2eeBadScFilter *
  cms2ecalLaserCorrFilter
)


