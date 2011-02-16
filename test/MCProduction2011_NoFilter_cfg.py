from MCProduction2011_cfg import *

process.out = cms.OutputModule(
        "PoolOutputModule",
        dropMetaData = cms.untracked.string("NONE"),
        fileName = cms.untracked.string('ntuple.root')
)
