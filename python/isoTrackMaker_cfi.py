import FWCore.ParameterSet.Config as cms
import CMS3.NtupleMaker.configProcessName as configProcessName

isoTrackMaker = cms.EDProducer("IsoTrackMaker",
                                  pfCandidatesTag     = cms.InputTag("packedPFCandidates","",configProcessName.name),
                                  isotrack_dz_cut = cms.double(0.1),
                                  isolation_dz_cut = cms.double(0.1),
                                  pflep_pt_cut = cms.double(5.0),
                                  pfhad_pt_cut = cms.double(10.0),
                                  coneR = cms.double(0.3)
)
