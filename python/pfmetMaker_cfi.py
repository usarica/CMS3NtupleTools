import FWCore.ParameterSet.Config as cms
import CMS3.NtupleMaker.configProcessName as configProcessName

pfmetMaker = cms.EDProducer("PFMETMaker",
                            aliasPrefix = cms.untracked.string("evt"),
                            pfMetInputTag_ = cms.InputTag("slimmedMETs","",configProcessName.name),
                            isData              = cms.bool(False),
                            onlySaveTwoVector   = cms.bool(False)
)

# this module gets the jets used by the MET tool box to recalculate T1 MET with OTF corrections
T1pfmetMaker = cms.EDProducer("PFMETMaker",
                              aliasPrefix = cms.untracked.string("evt_METToolbox"),
                              pfMetInputTag_ = cms.InputTag("slimmedMETs","","CMS3"),
                              isData              = cms.bool(False),
                              onlySaveTwoVector   = cms.bool(True)
)


