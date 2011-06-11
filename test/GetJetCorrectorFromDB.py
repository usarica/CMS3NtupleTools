import FWCore.ParameterSet.Config as cms
process = cms.Process("JEC")
#
# README
#
# This special Voodoo tool is "designed" to extract parameters of jet energy correction
# from a database. To run it you need:
# * input file (any that you can run on, don't ask why)
# * GlobalTag - it's critical. Looks like only one version of JEC is stored in
#   the database, so you need to get it right
# * "globalTag" parameter for JetCorrectorDBReader - it's nothing but prefix for output files
# 


# Voodoo
process.source = cms.Source("PoolSource",
      fileNames = cms.untracked.vstring('/store/mc/Spring11/GluGluToHToWWTo2L2Nu_M-160_7TeV-powheg-pythia6/AODSIM/PU_S1_START311_V1G1-v1/0013/068A6ABB-F350-E011-9DCF-00A0D1EEDDA8.root')
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

# Global Tag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# process.GlobalTag.globaltag = "START41_V0::All"
process.GlobalTag.globaltag = "START311_V2::All"

# jet correction extractor
process.readAK5PF = cms.EDAnalyzer('JetCorrectorDBReader',
                                   payloadName    = cms.untracked.string('AK5PF'),
                                   printScreen    = cms.untracked.bool(False),
                                   createTextFile = cms.untracked.bool(True),
                                   globalTag      = cms.untracked.string('Spring10') # voodoo parameter
                                   )
process.p = cms.Path( process.readAK5PF )
