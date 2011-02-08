import FWCore.ParameterSet.Config as cms

pdfinfoMaker = cms.EDProducer( "PDFInfoMaker",
        genEventInfoInputTag = cms.string("generator"),                               
	aliasPrefix = cms.untracked.string("pdfinfo")
)


