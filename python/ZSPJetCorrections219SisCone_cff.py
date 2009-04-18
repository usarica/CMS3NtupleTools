import FWCore.ParameterSet.Config as cms

# Jet corrections.
#
# Define the correction services for each algorithm.
#
ZSPJetCorrectorSisCone = cms.ESSource("ZSPJetCorrectionService",
    tagName = cms.string('ZSP_CMSSW219_Iterative_Cone_05'),
    label = cms.string('ZSPJetCorrectorSisCone')
#    tagName = cms.string('ZSP_CMSSW152_Iterative_Cone_05'),
)

#   
#   Define the producers of corrected jet collections for each algorithm.
#
ZSPJetCorJetSisCone = cms.EDProducer("CaloJetCorrectionProducer",
    src = cms.InputTag("sisCone5CaloJets"),
    correctors = cms.vstring('ZSPJetCorrectorSisCone'),
    alias = cms.untracked.string('ZSPJetCorSisCone')
)

#
#  Define a sequence to make all corrected jet collections at once.
#
ZSPJetCorrectionsSisCone = cms.Sequence(ZSPJetCorJetSisCone)

