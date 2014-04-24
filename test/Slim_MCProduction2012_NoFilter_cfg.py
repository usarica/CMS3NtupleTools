from CMS2.NtupleMaker.RecoConfiguration2012_cfg import *

# Global Tag
process.GlobalTag.globaltag = "START70_V6::All"

#Input
process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('file:/nfs-3/userdata/gzevi/patTuple_TTbar.root') # default file, old miniAOD version
#     fileNames = cms.untracked.vstring('file:patTuple_mini_eleClusFix.root') # just 50 events, but with prescales and fixed elecluster
     fileNames = cms.untracked.vstring('file:/nfs-3/userdata/gzevi/patTuple_mini_eleClusFix.root') # just 50 events, but with prescales and fixed elecluster
#    fileNames = cms.untracked.vstring('file:/nfs-3/userdata/gzevi/patTuple_mini_withL1.root') # just 35 events, but with prescales and L1GlobalTriggerReadoutRecord 
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

#
#process.cms2WithEverything = cms.Sequence( process.ak5PFJets * process.kt6PFJets * process.cms2CoreSequence * process.cms2PFNoTauSequence * process.cms2GENSequence )
#process.cms2WithEverything.remove(process.jptMaker)
#process.cms2WithEverything.remove(process.hypTrilepMaker)
#process.cms2WithEverything.remove(process.hypQuadlepMaker)

#process.cms2WithEverything = cms.Sequence( process.cms2PFNoTauSequence )
#process.p                  = cms.Path( process.cms2WithEverything )

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
    process.photonMaker 
    )


#
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.eventMaker.isData                        = cms.bool(False)
process.luminosityMaker.isData                   = process.eventMaker.isData

#Slim CMS2
from CMS2.NtupleMaker.SlimCms2_cff import slimcms2
process.out.outputCommands.extend(slimcms2)
