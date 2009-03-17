# The following comments couldn't be translated into the new config version:

# -*-sh-*-

import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("filterByEvent")

from Configuration.EventContent.EventContent_cff import *

process.configurationMetadata = cms.untracked.PSet(
        version = cms.untracked.string('$Revision: 1.1 $'),
        annotation = cms.untracked.string('CMS2'),
        name = cms.untracked.string('CMS2 test configuration')
)

##source 
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(os.environ.get("INFILE")),
                            eventsToProcess = cms.untracked.VEventID(cms.EventID(int(os.environ.get("RUN")), 
                                                                                 int(os.environ.get("EVENT"))))
)

## output
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string(os.environ.get("OUTFILE"))
)

#process.out.outputCommands = cms.untracked.vstring( 'drop *' )
#process.out.outputCommands.extend(AODSIMEventContent.outputCommands)

process.outpath = cms.EndPath(process.out)
