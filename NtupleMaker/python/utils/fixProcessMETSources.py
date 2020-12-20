import FWCore.ParameterSet.Config as cms

def fixProcessMETSources(process, year, metrecipe):
   if year == 2017 and metrecipe:
      if hasattr(process, "patPFMetModifiedMET") and hasattr(process, "pfMetModifiedMET"):
         _process_to_fix = getattr(process, "patPFMetModifiedMET")
         _process_to_fix.srcPFCands = process.pfMetModifiedMET.src
         print "fixProcessMETSources: Fixed process.patPFMetModifiedMET.srcPFCands"
      if hasattr(process, "patPFMetModifiedPuppiMET") and hasattr(process, "pfMetModifiedPuppiMET"):
         _process_to_fix = getattr(process, "patPFMetModifiedPuppiMET")
         _process_to_fix.srcPFCands = process.pfMetModifiedPuppiMET.src
         print "fixProcessMETSources: Fixed process.patPFMetModifiedPuppiMET.srcPFCands"
