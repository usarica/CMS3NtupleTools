import FWCore.ParameterSet.Config as cms
import CMS3.NtupleMaker.configProcessName as configProcessName

fixedGridRhoFastJetCentralCaloMaker = cms.EDProducer("EnergyDensityMaker",
                                             input = cms.InputTag("fixedGridRhoFastjetCentralCalo","", configProcessName.name2),
                                             alias = cms.untracked.string("evt_fixgridfastjet_centralcalo_rho"))
