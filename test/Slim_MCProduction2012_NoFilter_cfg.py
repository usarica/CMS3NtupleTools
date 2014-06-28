from CMS2.NtupleMaker.RecoConfiguration2012_cfg import *

# Global Tag
process.GlobalTag.globaltag = "START70_V6::All"
#process.GlobalTag.globaltag = "PLS170_V7AN1::All"

#Input
process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('file:/nfs-3/userdata/gzevi/DYJetsToLL_M-50_13TeV-madgraph-pythia8_PAT.root')
 fileNames = cms.untracked.vstring(
#       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/003E832C-8AFC-E311-B7AA-002590596490.root',
#       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/00907CE8-8CFC-E311-8AEB-0025905B8562.root',
#       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/00FC23A0-8AFC-E311-B5F1-0025905A6138.root',
#       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/029E0636-8CFC-E311-A0C8-0025905A610A.root',
#       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/042E80B7-8BFC-E311-9752-0025905A612A.root',
#       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/06591229-8DFC-E311-A896-002618943901.root',
#       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/065C10DA-8AFC-E311-98B4-003048D15DE0.root',
#       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/066046AA-8AFC-E311-9DA9-00261894396B.root',
#       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/0800059A-8AFC-E311-B8AC-0025905A48D0.root',
#       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/081F40D0-8BFC-E311-9145-002618FDA262.root',
#       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/081F8E04-8CFC-E311-95B0-00261894398C.root',
#       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/084B5BDF-8AFC-E311-8CF2-00261894398A.root',
       '/store/mc/Spring14miniaod/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_POSTLS170_V5-v1/00000/085E17FC-8AFC-E311-B079-0025905A6136.root') 

)

# Output
process.out = cms.OutputModule(
        "PoolOutputModule",
#        fileName     = cms.untracked.string('TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_PU20bx25_POSTLS170_V5-v1.root'),
        fileName     = cms.untracked.string('ntuple.root'),
        dropMetaData = cms.untracked.string("NONE")
)
process.outpath      = cms.EndPath(process.out)
process.maxEvents                     = cms.untracked.PSet( input = cms.untracked.int32(50) )

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
    process.hcalNoiseSummaryMaker*
    process.miniAODrhoSequence

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
