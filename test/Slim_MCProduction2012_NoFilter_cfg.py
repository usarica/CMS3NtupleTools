from CMS2.NtupleMaker.RecoConfiguration2012_cfg import *

# Global Tag
process.GlobalTag.globaltag = "START70_V6::All"

#Input
process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('file:/nfs-3/userdata/gzevi/DYJetsToLL_M-50_13TeV-madgraph-pythia8_PAT.root')
    fileNames = cms.untracked.vstring('file:/nfs-3/userdata/gzevi/TT_Tune4C_13TeV-pythia8-tauola_PAT.root')
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
    process.pftauMaker*
    process.photonMaker*
    process.genMaker*
    process.genJetMaker*
    process.muToTrigAssMaker*  # requires muonMaker
    process.elToTrigAssMaker*  # requires electronMaker
    process.candToGenAssMaker* # requires electronMaker, muonMaker, pfJetMaker, photonMaker
    process.pdfinfoMaker*
    process.puSummaryInfoMaker*
    process.recoConversionMaker*
    process.metFilterMaker*
    process.hcalNoiseSummaryMaker

    # Optional (filters)
    #process.dilepGenFilter*    # requires genMaker    

    # Optional (requires LHEEventProduct, need a SUSY sample to test)
    #process.sParmMakerSequence
    )


#
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.eventMaker.isData                        = cms.bool(False)
process.luminosityMaker.isData                   = process.eventMaker.isData

#Slim CMS2
from CMS2.NtupleMaker.SlimCms2_cff import slimcms2
process.out.outputCommands.extend(slimcms2)
