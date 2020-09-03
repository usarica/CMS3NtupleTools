import FWCore.ParameterSet.Config as cms

def fixProcessPuppiSources(process, puppicollname, puppialgoname):
   if hasattr(process, "pfCandidateJetsWithEEnoiseModifiedPuppiMET"):
      _process_to_fix = getattr(process, "pfCandidateJetsWithEEnoiseModifiedPuppiMET")
      _process_to_fix.jetsrc = cms.InputTag(puppicollname)
      print "fixProcessPuppiSources: Fixed process.pfCandidateJetsWithEEnoiseModifiedPuppiMET"
   if hasattr(process, "pfcandidateClusteredModifiedPuppiMET"):
      _process_to_fix = getattr(process, "pfcandidateClusteredModifiedPuppiMET")
      _process_to_fix.src = cms.VInputTag(
         cms.InputTag("slimmedElectrons"), cms.InputTag("slimmedMuons"), cms.InputTag("slimmedTaus"), cms.InputTag("slimmedPhotons"), cms.InputTag(puppicollname)
      )
      print "fixProcessPuppiSources: Fixed process.pfCandidateJetsWithEEnoiseModifiedPuppiMET"
   if hasattr(process, "patCaloMet"):
      _process_to_fix = getattr(process, "patCaloMet")
      _process_to_fix.metSource = cms.InputTag("metrawCaloModifiedMET")
      print "fixProcessPuppiSources: Fixed process.pfCandidateJetsWithEEnoiseModifiedPuppiMET"
   if hasattr(process, "patSmearedJetsModifiedPuppiMET"):
      _process_to_fix = getattr(process, "patSmearedJetsModifiedPuppiMET")
      _process_to_fix.algo = cms.string(puppialgoname)
      _process_to_fix.algopt = cms.string(puppialgoname+"_pt")
      print "fixProcessPuppiSources: Fixed process.patSmearedJetsModifiedPuppiMET"
