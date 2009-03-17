# The following comments couldn't be translated into the new config version:

# -*-sh-*-

import os
import FWCore.ParameterSet.Config as cms

process = cms.Process("concat")

from Configuration.EventContent.EventContent_cff import *

process.configurationMetadata = cms.untracked.PSet(
        version = cms.untracked.string('$Revision: 1.1 $'),
        annotation = cms.untracked.string('CMS2'),
        name = cms.untracked.string('CMS2 test configuration')
)

##source 
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
 "file:RUN-2-EVENT-265-LUMI-1970-0659BAB4-F1E7-DD11-ABDE-0015170ACA8C.root", "file:RUN-2-EVENT-264-LUMI-1587-4C9AB13E-EEE7-DD11-9614-00151715BB94.root", "file:RUN-2-EVENT-235-LUMI-760-264126A9-B9E4-DD11-BABB-00E08178C0A5.root", "file:RUN-2-EVENT-119-LUMI-72-1CC76818-26E1-DD11-9C52-003048635DE8.root", "file:RUN-2-EVENT-32-LUMI-1176-B632118B-73E1-DD11-9FA4-0015170AC484.root", "file:RUN-2-EVENT-266-LUMI-2549-F64D33C5-0BE7-DD11-A8E5-0013D3DE26BD.root", "file:RUN-2-EVENT-6-LUMI-1140-685C59AC-E2E1-DD11-9A2D-00161725E4BB.root", "file:RUN-2-EVENT-381-LUMI-1914-06F23705-EFE7-DD11-9147-00E0812D8D70.root", "file:RUN-2-EVENT-381-LUMI-1914-A485DCEB-FBE6-DD11-9A9D-0013D3DE2655.root", "file:RUN-2-EVENT-269-LUMI-1374-6083B6EE-2CE8-DD11-86CB-00151715BB94.root", "file:RUN-2-EVENT-269-LUMI-1374-685C59AC-E2E1-DD11-9A2D-00161725E4BB.root", "file:RUN-2-EVENT-381-LUMI-1914-06F23705-EFE7-DD11-9147-00E0812D8D70.root", "file:RUN-2-EVENT-381-LUMI-1914-A485DCEB-FBE6-DD11-9A9D-0013D3DE2655.root", "file:RUN-2-EVENT-326-LUMI-1980-4C9AB13E-EEE7-DD11-9614-00151715BB94.root", "file:RUN-2-EVENT-326-LUMI-1980-D8253F67-EEE7-DD11-AE43-00E08178C111.root", "file:RUN-2-EVENT-350-LUMI-1808-F44DE2B7-B9E4-DD11-8CEC-00E0817917C3.root", "file:RUN-2-EVENT-55-LUMI-1075-78D3CEE3-EEE7-DD11-A305-00161725E4BF.root", "file:RUN-2-EVENT-55-LUMI-1075-E20DE4A0-E7E6-DD11-9182-0013D3DE269F.root", "file:RUN-2-EVENT-221-LUMI-97-18A223D5-2AE3-DD11-A4DA-00E08178C0C7.root", "file:RUN-2-EVENT-107-LUMI-2622-8A0164A3-E7E6-DD11-A68F-00161725E500.root",
        )
)

## output
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('zjetsCands.root')
)

#process.out.outputCommands = cms.untracked.vstring( 'drop *' )
#process.out.outputCommands.extend(AODSIMEventContent.outputCommands)

process.outpath = cms.EndPath(process.out)
