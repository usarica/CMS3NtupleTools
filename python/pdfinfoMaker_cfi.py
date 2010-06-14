import FWCore.ParameterSet.Config as cms

pdfinfoMaker = cms.EDProducer( "PDFInfoMaker",
	aliasPrefix = cms.untracked.string("pdfinfo"),
)


