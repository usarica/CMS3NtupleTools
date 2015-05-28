import FWCore.ParameterSet.Config as cms
#from RecoJets.Configuration.RecoPFJets_cff import *
#from RecoJets.Configuration.RecoJets_cff import *
from CMS3.NtupleMaker.fastJetMaker_cfi import *
#from RecoTauTag.Configuration.HPSPFTaus_cff import kt6PFJetsForRhoComputationVoronoi

#kt6PFJetsDeterministicJEC = kt4PFJets.clone(
#    rParam = 0.6,
#    doAreaFastjet = True,
#    doRhoFastjet = True,
#    voronoiRfact = 0.9,
#    Rho_EtaMax = 5.0,
#    Ghost_EtaMax = 5.0
#)

#kt6PFJetsDeterministicIso = kt4PFJets.clone(
#    rParam = 0.6,
#    doAreaFastjet = True,
#    doRhoFastjet = True,
#    voronoiRfact = 0.9,
#    Rho_EtaMax = 2.5,
#    Ghost_EtaMax = 2.5
#)

#kt6PFJetsForRhoComputationVoronoi = kt6PFJets.clone(doRhoFastjet = True, voronoiRfact = 0.9)
#kt6PFJetsForRhoComputationActiveArea = kt6PFJets.clone(doRhoFastjet = True, voronoiRfact = -0.9)
#kt6PFJetsForRhoComputationRandom     = kt6PFJets.clone(doRhoFastjet = True, voronoiRfact = -1)

#kt6PFJetsForEGIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
#kt6PFJetsForEGIsolation.Rho_EtaMax = cms.double(2.5)

#kt6CaloJetsForMuHLT = kt6CaloJets.clone(doRhoFastjet = True, doAreaFastjet = False, voronoiRfact = 0.9, Ghost_EtaMax = 5.0, Rho_EtaMax = 2.5, doAreaDiskApprox = True, doPVCorrection = False)

#wwFastJetSequence = cms.Sequence( kt6PFJetsForRhoComputationVoronoi * kt6PFJetsForRhoComputationActiveArea * kt6PFJetsForRhoComputationRandom *
#                                  wwRhoDefaultMaker * wwRhoVoronoiMaker * wwRhoActiveAreaMaker * wwRhoRandomMaker )

#additionalFastJetSequence = cms.Sequence( fixedGridRhoAllMaker * fixedGridRhoFastJetAllMaker *
#                                          kt6CaloJetsRhoMaker * kt6CaloJetsCentralRhoMaker *
#                                          kt6PFJetsCentralChargedPileUpRhoMaker * kt6PFJetsCentralNeutralRhoMaker *
#                                          kt6PFJetsCentralNeutralTightRhoMaker * kt6PFJetsForEGIsolation * kt6PFJetsForEGIsolationRhoMaker *
#                                          kt6CaloJetsForMuHLT * kt6CaloJetsForMuHLTRhoMaker)

miniAODrhoSequence = cms.Sequence ( fixedGridRhoAllMaker * fixedGridRhoFastJetAllMaker * fixedGridRhoFastJetAllCaloMaker * fixedGridRhoFastJetCentralCaloMaker * fixedGridRhoFastJetCentralChargedPileUpMaker * fixedGridRhoFastJetCentralNeutralMaker )
#fastJetSequence = cms.Sequence( kt6PFJetsDeterministicJEC * kt6PFJetsDeterministicIso * fastJetMaker * wwFastJetSequence * additionalFastJetSequence)
#fastJetSequence = cms.Sequence( kt6PFJetsDeterministicJEC * kt6PFJetsDeterministicIso * wwFastJetSequence * additionalFastJetSequence)
