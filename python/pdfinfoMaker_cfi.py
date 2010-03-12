import FWCore.ParameterSet.Config as cms

pdfinfoMaker = cms.EDFilter( "PDFInfoMaker",
	aliasPrefix = cms.untracked.string("pdfinfo"),
)


