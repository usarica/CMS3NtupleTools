from CMS3.NtupleMaker.RecoConfiguration2012_cfg import *

# Global Tag
#process.GlobalTag.globaltag = "START70_V6::All"
#process.GlobalTag.globaltag = "MCRUN2_72_V1A::All"
process.GlobalTag.globaltag = "PHYS14_25_V1::All"
#process.GlobalTag.globaltag = "PLS170_V7AN1::All"

#Input
process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('file:/nfs-3/userdata/gzevi/DYJetsToLL_M-50_13TeV-madgraph-pythia8_PAT.root')
    fileNames = cms.untracked.vstring('file:/home/users/namin/stop/cms3/CMSSW_7_2_0/src/SMS-T1tttt_PU20bx25_tsg_PHYS14_25_V1.root')
)

# Output
process.out = cms.OutputModule(
        "PoolOutputModule",
        fileName     = cms.untracked.string('ntuple.root'),
        dropMetaData = cms.untracked.string("NONE")
)
process.outpath      = cms.EndPath(process.out)
process.maxEvents                     = cms.untracked.PSet( input = cms.untracked.int32(12000) )

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
    process.egmGsfElectronIDSequence *     
    process.beamSpotMaker *
    process.vertexMaker *
    process.secondaryVertexMaker *
    process.pfCandidateMaker*
    process.eventMaker*
    process.electronMaker*
    process.muonMaker*
    process.pfJetMaker*
    process.subJetMaker *
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
    process.hcalNoiseSummaryMaker*
    process.miniAODrhoSequence*
    process.hypDilepMaker

    # Optional (filters)
    #process.dilepGenFilter*    # requires genMaker    

    # Optional (requires LHEEventProduct, need a SUSY sample to test)
    #process.sParmMakerSequence
    )


#
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.eventMaker.isData                        = cms.bool(False)
process.luminosityMaker.isData                   = process.eventMaker.isData

#process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

#Slim CMS2
from CMS3.NtupleMaker.SlimCms2_cff import slimcms2
process.out.outputCommands.extend(slimcms2)
