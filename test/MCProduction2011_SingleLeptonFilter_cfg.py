from MCProduction2011_cfg import *

# Single Lepton Filter
process.EventSelectionSingleFilt = cms.PSet(
  SelectEvents = cms.untracked.PSet(
    SelectEvents = cms.vstring('pWithRecoLepton', 'pWithGenLepton')
  )
)
process.out = cms.OutputModule(
        "PoolOutputModule",
        process.EventSelectionSingleFilt,
        dropMetaData = cms.untracked.string("NONE"),
        fileName = cms.untracked.string('ntuple.root')
)
process.load('CMS2.NtupleMaker.aSkimFilter_cfi')
process.load('CMS2.NtupleMaker.monolepGenFilter_cfi')
process.pWithRecoLepton = cms.Path(process.cms2WithEverything * process.aSkimFilter   )
process.pWithGenLepton  = cms.Path(process.cms2WithEverything * process.monolepGenFilter  )
