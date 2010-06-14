import FWCore.ParameterSet.Config as cms

<<<<<<< pfElectronMaker_cfi.py
pfElectronMaker = cms.EDProducer("PFElectronMaker",
                               pfCandidatesTag = cms.InputTag("particleFlow","electrons")
=======
pfElectronMaker = cms.EDFilter("PFElectronMaker",
                               pfCandidatesTag     = cms.InputTag("particleFlow","electrons"),
                               isoc_vm_tag         = cms.InputTag("CMS2isoValElectronWithCharged"),                                  
                               ison_vm_tag         = cms.InputTag("CMS2isoValElectronWithNeutral"),       
                               isop_vm_tag         = cms.InputTag("CMS2isoValElectronWithPhotons"),       
                               pfAllElectrons_tag  = cms.InputTag("CMS2pfAllElectrons")
>>>>>>> 1.2
)
