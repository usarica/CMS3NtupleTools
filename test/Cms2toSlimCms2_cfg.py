import FWCore.ParameterSet.Config as cms

process = cms.Process("CMS2")

process.source = cms.Source("PoolSource",
    skipEvents = cms.untracked.uint32(0),
    fileNames  = cms.untracked.vstring(
        "file:/home/users/fgolf/devel/CMSSW_5_3_2_patch4/src/CMS2/NtupleMaker/test/full_cms2.root"
    ),
)

process.out = cms.OutputModule(
        "PoolOutputModule",
        fileName     = cms.untracked.string('test.root')
)
process.outpath      = cms.EndPath(process.out)

process.out.outputCommands = cms.untracked.vstring("keep *")

from CMS2.NtupleMaker.SlimCms2_cff import slimcms2
process.out.outputCommands.extend(slimcms2)
