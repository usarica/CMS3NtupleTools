from CMS2.NtupleMaker.RecoConfiguration2012_cfg import *

# Global Tag
## globaltag used for rereco of 2012D in 700
process.GlobalTag.globaltag = "FT_R_70_V1::All"

#Input
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:/nfs-7/userdata/olivito/miniaod/patTuple_mini_dimu.root')
)

# Output
process.out = cms.OutputModule(
  "PoolOutputModule",
  fileName     = cms.untracked.string('ntuple.root'),
  dropMetaData = cms.untracked.string("NONE")
)
process.outpath      = cms.EndPath(process.out)

# Branches
process.out.outputCommands = cms.untracked.vstring( 'drop *' )
process.out.outputCommands.extend(cms.untracked.vstring('keep *_*Maker*_*_CMS2*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop *_cms2towerMaker*_*_CMS2*'))
process.out.outputCommands.extend(cms.untracked.vstring('drop CaloTowers*_*_*_CMS2*'))

process.pfJetMETcorr.jetCorrLabel = cms.string("ak5PFL1FastL2L3Residual")

#
# process.cms2WithEverything = cms.Sequence( process.ak5PFJets * process.kt6PFJets * process.cms2CoreSequence * process.cms2PFNoTauSequence )
# process.cms2WithEverything.remove(process.jptMaker)
# process.cms2WithEverything.remove(process.hypTrilepMaker)
# process.cms2WithEverything.remove(process.hypQuadlepMaker)
# process.p                  = cms.Path( process.cms2WithEverything )

process.p                  = cms.Path( 
    #process.unpackedTracksAndVertices *
    process.beamSpotMaker *
    process.vertexMaker *
    process.pfCandidateMaker*
    process.eventMaker*
    process.electronMaker*
    process.muonMaker*
    process.pfJetMaker*
    process.pfmetMaker*
    process.hltMakerSequence*
    process.pftauMaker*
    process.photonMaker*
    process.muToTrigAssMaker*  # requires muonMaker
    process.elToTrigAssMaker*  # requires electronMaker
    process.puSummaryInfoMaker
    )

#
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.eventMaker.isData                        = cms.bool(True)
process.luminosityMaker.isData                   = process.eventMaker.isData

#process.maxEvents.input = 10
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck")
#process.Timing =cms.Service("Timing")        

#Slim CMS2
from CMS2.NtupleMaker.SlimCms2_cff import slimcms2
process.out.outputCommands.extend(slimcms2)
