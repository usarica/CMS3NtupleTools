import FWCore.ParameterSet.Config as cms

fixedGridRhoFastJetCentralCaloMaker = cms.EDProducer("EnergyDensityMaker",
                                             input = cms.InputTag("fixedGridRhoFastjetCentralCalo"),
                                             alias = cms.untracked.string("evt_fixgridfastjet_centralcalo_rho"))

fixedGridRhoAllMaker = cms.EDProducer("EnergyDensityMaker",
                                      input = cms.InputTag("fixedGridRhoAll"),
                                      alias = cms.untracked.string("evt_fixgrid_all_rho"))

fixedGridRhoFastJetAllMaker = cms.EDProducer("EnergyDensityMaker",
                                             input = cms.InputTag("fixedGridRhoFastjetAll"),
                                             alias = cms.untracked.string("evt_fixgridfastjet_all_rho"))

fixedGridRhoFastJetAllCaloMaker = cms.EDProducer("EnergyDensityMaker",
                                             input = cms.InputTag("fixedGridRhoFastjetAllCalo"),
                                             alias = cms.untracked.string("evt_fixgridfastjet_allcalo_rho"))

fixedGridRhoFastJetCentralChargedPileUpMaker = cms.EDProducer("EnergyDensityMaker",
                                             input = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp"),
                                             alias = cms.untracked.string("evt_fixgridfastjet_centralchargedpileup_rho"))

fixedGridRhoFastJetCentralNeutralMaker = cms.EDProducer("EnergyDensityMaker",
                                             input = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
                                             alias = cms.untracked.string("evt_fixgridfastjet_centralneutral_rho"))

fixedGridRhoFastJetAllCentralMaker = cms.EDProducer("EnergyDensityMaker",
                                             input = cms.InputTag("fixedGridRhoFastjetCentral"),
                                             alias = cms.untracked.string("evt_fixgridfastjet_central_rho"))
