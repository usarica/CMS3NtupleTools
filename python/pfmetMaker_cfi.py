
import FWCore.ParameterSet.Config as cms
import CMS3.NtupleMaker.configProcessName as configProcessName

pfmetMaker = cms.EDProducer("PFMETMaker",
                            aliasPrefix = cms.untracked.string("evt"),
                            pfMetInputTag_ = cms.InputTag("slimmedMETs","",configProcessName.name),
                            onlySaveTwoVector   = cms.bool(False),
                            doUncertainties   = cms.bool(True)
)

pfmetMakerEGClean = cms.EDProducer("PFMETMaker",
                            aliasPrefix = cms.untracked.string("evt_egclean"),
                            pfMetInputTag_ = cms.InputTag("slimmedMETsEGClean","",configProcessName.name),
                            onlySaveTwoVector   = cms.bool(False),
                            doUncertainties   = cms.bool(True)
)

pfmetMakerMuEGClean = cms.EDProducer("PFMETMaker",
                            aliasPrefix = cms.untracked.string("evt_muegclean"),
                            pfMetInputTag_ = cms.InputTag("slimmedMETsMuEGClean","",configProcessName.name),
                            onlySaveTwoVector   = cms.bool(False),
                            doUncertainties   = cms.bool(True)
)

pfmetMakerUncorr = cms.EDProducer("PFMETMaker",
                            aliasPrefix = cms.untracked.string("evt_uncorr"),
                            pfMetInputTag_ = cms.InputTag("slimmedMETsUncorrected","",configProcessName.name),
                            onlySaveTwoVector   = cms.bool(False),
                            doUncertainties   = cms.bool(True)
)

pfmetNoHFMaker = cms.EDProducer("PFMETMaker",
                            aliasPrefix = cms.untracked.string("evt_NoHF"),
                            pfMetInputTag_ = cms.InputTag("slimmedMETsNoHF","",configProcessName.name),
                                onlySaveTwoVector   = cms.bool(False),
                            doUncertainties   = cms.bool(False) # I don't know how uncertainties are treated in the NoHF case
)

pfmetpuppiMaker = cms.EDProducer("PFMETMaker",
                            aliasPrefix = cms.untracked.string("evt_puppi"),
                            pfMetInputTag_ = cms.InputTag("slimmedMETsPuppi","",configProcessName.name),
                            onlySaveTwoVector   = cms.bool(False),
                            doUncertainties   = cms.bool(True)
)

# this module gets the jets used by the MET tool box to recalculate T1 MET with OTF corrections
# T1pfmetMaker = cms.EDProducer("PFMETMaker",
#                               aliasPrefix = cms.untracked.string("evt_METToolbox"),
#                               pfMetInputTag_ = cms.InputTag("slimmedMETs","","CMS3"),
#                               onlySaveTwoVector   = cms.bool(True)
# )

# T1pfmetNoHFMaker = cms.EDProducer("PFMETMaker",
#                                   aliasPrefix = cms.untracked.string("evt_METToolboxNoHF"),
#                                   pfMetInputTag_ = cms.InputTag("slimmedMETsNoHF","","CMS3"),
#                                   onlySaveTwoVector   = cms.bool(True)
# )


