# The following comments couldn't be translated into the new config version:

# -*-sh-*-

import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("filterByLumi")

from Configuration.EventContent.EventContent_cff import *

process.configurationMetadata = cms.untracked.PSet(
        version = cms.untracked.string('$Revision: 1.1 $'),
        annotation = cms.untracked.string('CMS2'),
        name = cms.untracked.string('CMS2 test configuration')
)

##source 
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(os.environ.get("INFILE")),
                            # NB: the next line only makes sense if
                            # you've buggered PoolSource into
                            # inverting the meaning of "lumisToSkip",
                            # see http://omega.physics.ucsb.edu/twiki/bin/view/CMS/EventListToEdmFile
                            lumisToSkip = cms.untracked.VLuminosityBlockID(cms.LuminosityBlockID(int(os.environ.get("RUN")),
                                                                                                 int(os.environ.get("LUMI"))))
)

## output
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string(os.environ.get("OUTFILE"))
)

process.out.outputCommands = cms.untracked.vstring( 'drop *' )
process.out.outputCommands.extend(AODSIMEventContent.outputCommands)

process.outpath = cms.EndPath(process.out)
