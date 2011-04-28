from CMS2.NtupleMaker.RecoConfiguration2011_cfg import *

# Global Tag
process.GlobalTag.globaltag = "GR_R_311_V2::All"

# Load Filters
process.load("CMS2.NtupleMaker.sdFilter_cfi")
process.filter = cms.Path(process.sdFilter)
if "DoubleElectron" in str(process.source.fileNames):
    process.sdFilter.filterName_=cms.string("doubleElectron")
elif "DoubleMu" in str(process.source.fileNames):
    process.sdFilter.filterName_=cms.string("doubleMu")
elif "MuEG" in  str(process.source.fileNames):
    process.sdFilter.filterName_=cms.string("MuEG")
elif "SingleMu" in str(process.source.fileNames):
    process.sdFilter.filterName_=cms.string("SingleMu")
elif "Photon" in str(process.source.fileNames):
    process.sdFilter.filterName_=cms.string("Photon")
elif "ElectronHad" in str(process.source.fileNames):
    process.sdFilter.filterName_=cms.string("ElectronHad")
    process.hypDilepMaker.TightLepton_PtCut  = cms.double(5.0)
    process.hypDilepMaker.LooseLepton_PtCut  = cms.double(5.0)
    process.hypTrilepMaker.TightLepton_PtCut = cms.double(5.0)
    process.hypTrilepMaker.LooseLepton_PtCut = cms.double(5.0)
elif "MuHad" in str(process.source.fileNames):
    process.sdFilter.filterName_=cms.string("MuHad")
    process.hypDilepMaker.TightLepton_PtCut  = cms.double(5.0)
    process.hypDilepMaker.LooseLepton_PtCut  = cms.double(5.0)
    process.hypTrilepMaker.TightLepton_PtCut = cms.double(5.0)
    process.hypTrilepMaker.LooseLepton_PtCut = cms.double(5.0)
else:
    print 'filterName missing!'
# Output
process.out = cms.OutputModule(
        "PoolOutputModule",
        SelectEvents = cms.untracked.PSet( SelectEvents = cms.vstring('filter') ),
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
process.cms2WithEverything = cms.Sequence( process.sdFilter * process.ak5PFJets * process.kt6PFJets * process.cms2CoreSequence * process.cms2PFNoTauSequence )
process.p                  = cms.Path( process.cms2WithEverything )

#
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.eventMaker.isData                        = cms.bool(True)
