import FWCore.ParameterSet.Config as cms

process = cms.Process("CMS3")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(),
    noEventSort = cms.untracked.bool(True)
)
process.AODEventContent = cms.PSet(
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
    outputCommands = cms.untracked.vstring( ('drop *', 
        'keep ClusterSummary_clusterSummaryProducer_*_*', 
        'keep *_dt4DSegments_*_*', 
        'keep *_dt4DCosmicSegments_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*', 
        'keep *_castorreco_*_*', 
        'keep *_reducedHcalRecHits_*_*', 
        'keep HcalUnpackerReport_castorDigis_*_*', 
        'keep HcalUnpackerReport_hcalDigiAlCaMB_*_*', 
        'keep HcalUnpackerReport_hcalDigis_*_*', 
        'keep *_selectDigi_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEB_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEE_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsES_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep recoCaloClusters_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_hybridSuperClusters_uncleanOnlyHybridSuperClusters_*', 
        'keep recoCaloClusters_multi5x5SuperClusters_multi5x5EndcapBasicClusters_*', 
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_*', 
        'keep recoSuperClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoCaloClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoSuperClusters_particleFlowSuperClusterOOTECAL_*_*', 
        'keep recoCaloClusters_particleFlowSuperClusterOOTECAL_*_*', 
        'keep recoTracks_GsfGlobalElectronTest_*_*', 
        'keep recoGsfTracks_electronGsfTracks_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_conversionStepTracks_*_*', 
        'keep recoTracks_beamhaloTracks_*_*', 
        'keep recoTracks_ctfPixelLess_*_*', 
        'keep *_dedxHarmonic2_*_*', 
        'keep *_dedxPixelHarmonic2_*_*', 
        'keep *_dedxHitInfo_*_*', 
        'keep *_trackExtrapolator_*_*', 
        'keep *_generalTracks_MVAValues_*', 
        'keep *_generalTracks_MVAVals_*', 
        'keep *_ak4CaloJets_*_*', 
        'keep *_ak4PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHSSoftDrop_*_*', 
        'keep *_cmsTopTagPFJetsCHS_*_*', 
        'keep *_ak4PFJets_*_*', 
        'keep *_JetPlusTrackZSPCorJetAntiKt4_*_*', 
        'keep *_ak4TrackJets_*_*', 
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', 
        'keep *_caloTowers_*_*', 
        'keep *_CastorTowerReco_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertex_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertexPF_*_*', 
        'keep *_ak4JetTracksAssociatorExplicit_*_*', 
        'keep *_ak4JetExtender_*_*', 
        'keep *_ak4JetID_*_*', 
        'keep *_ak5CastorJets_*_*', 
        'keep *_ak5CastorJetID_*_*', 
        'keep *_ak7CastorJets_*_*', 
        'keep *_ak7CastorJetID_*_*', 
        'keep *_fixedGridRhoAll_*_*', 
        'keep *_fixedGridRhoFastjetAll_*_*', 
        'keep *_fixedGridRhoFastjetAllTmp_*_*', 
        'keep *_fixedGridRhoFastjetCentral_*_*', 
        'keep *_fixedGridRhoFastjetAllCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralChargedPileUp_*_*', 
        'keep *_fixedGridRhoFastjetCentralNeutral_*_*', 
        'drop doubles_*Jets_rhos_*', 
        'drop doubles_*Jets_sigmas_*', 
        'keep *_ak8PFJetsCHSSoftDropMass_*_*', 
        'keep recoCaloMETs_caloMet_*_*', 
        'keep recoCaloMETs_caloMetBE_*_*', 
        'keep recoCaloMETs_caloMetBEFO_*_*', 
        'keep recoCaloMETs_caloMetM_*_*', 
        'keep recoPFMETs_pfMet_*_*', 
        'keep recoPFMETs_pfChMet_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_*_*', 
        'drop recoHcalNoiseRBXs_*_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*', 
        'keep recoGlobalHaloData_GlobalHaloData_*_*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep *_muons_*_*', 
        'keep *_particleFlow_muons_*', 
        'drop *_muons_muons1stStep2muonsMap_*', 
        'drop recoIsoDepositedmValueMap_muons_*_*', 
        'drop doubleedmValueMap_muons_muPFIso*_*', 
        'keep recoTracks_standAloneMuons_*_*', 
        'keep recoTrackExtras_standAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_standAloneMuons_*_*', 
        'keep recoTracks_globalMuons_*_*', 
        'keep recoTrackExtras_globalMuons_*_*', 
        'keep recoTracks_tevMuons_*_*', 
        'keep recoTrackExtras_tevMuons_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_displacedTracks_*_*', 
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_*_*', 
        'keep recoTracks_displacedGlobalMuons_*_*', 
        'keep recoTrackExtras_displacedGlobalMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedGlobalMuons_*_*', 
        'keep recoTracks_cosmicMuons_*_*', 
        'keep recoMuons_muonsFromCosmics_*_*', 
        'keep recoTracks_cosmicMuons1Leg_*_*', 
        'keep recoMuons_muonsFromCosmics1Leg_*_*', 
        'keep recoTracks_refittedStandAloneMuons_*_*', 
        'keep recoTrackExtras_refittedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_refittedStandAloneMuons_*_*', 
        'keep recoTracks_displacedStandAloneMuons__*', 
        'keep recoTrackExtras_displacedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedStandAloneMuons_*_*', 
        'keep *_softPFElectronBJetTags_*_*', 
        'keep *_softPFMuonBJetTags_*_*', 
        'keep *_pfTrackCountingHighEffBJetTags_*_*', 
        'keep *_pfJetProbabilityBJetTags_*_*', 
        'keep *_pfJetBProbabilityBJetTags_*_*', 
        'keep *_pfSimpleSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfCombinedSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfGhostTrackBJetTags_*_*', 
        'keep *_pfCombinedMVAV2BJetTags_*_*', 
        'keep *_inclusiveCandidateSecondaryVertices_*_*', 
        'keep *_inclusiveCandidateSecondaryVerticesCvsL_*_*', 
        'keep *_pfCombinedCvsLJetTags_*_*', 
        'keep *_pfCombinedCvsBJetTags_*_*', 
        'keep *_pfChargeBJetTags_*_*', 
        'keep *_pfDeepCSVJetTags_*_*', 
        'keep *_pfDeepCMVAJetTags_*_*', 
        'keep recoRecoTauPiZeros_hpsPFTauProducer_pizeros_*', 
        'keep recoPFTaus_hpsPFTauProducer_*_*', 
        'keep *_hpsPFTauDiscriminationByDeadECALElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFinding_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingNewDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingOldDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseChargedIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByTightMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSum_*_*', 
        'keep *_hpsPFTauPUcorrPtSum_*_*', 
        'keep *_hpsPFTauChargedIsoPtSum_*_*', 
        'keep *_hpsPFTauTransverseImpactParameters_*_*', 
        'keep *_hpsPFTauFootprintCorrection_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeight_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6rawElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6LooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6MediumElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6TightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VTightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauChargedIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeightdR03_*_*', 
        'keep *_hpsPFTauFootprintCorrectiondR03_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalConedR03_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep  *_offlinePrimaryVertices__*', 
        'keep *_offlinePrimaryVerticesWithBS_*_*', 
        'keep *_offlinePrimaryVerticesFromCosmicTracks_*_*', 
        'keep *_nuclearInteractionMaker_*_*', 
        'keep *_generalV0Candidates_*_*', 
        'keep *_inclusiveSecondaryVertices_*_*', 
        'keep recoGsfElectronCores_gsfElectronCores_*_*', 
        'keep recoGsfElectronCores_gedGsfElectronCores_*_*', 
        'keep recoGsfElectrons_gsfElectrons_*_*', 
        'keep recoGsfElectrons_gedGsfElectrons_*_*', 
        'keep recoGsfElectronCores_uncleanedOnlyGsfElectronCores_*_*', 
        'keep recoGsfElectrons_uncleanedOnlyGsfElectrons_*_*', 
        'keep floatedmValueMap_eidRobustLoose_*_*', 
        'keep floatedmValueMap_eidRobustTight_*_*', 
        'keep floatedmValueMap_eidRobustHighEnergy_*_*', 
        'keep floatedmValueMap_eidLoose_*_*', 
        'keep floatedmValueMap_eidTight_*_*', 
        'keep *_egmGedGsfElectronPFIsolation_*_*', 
        'keep *_photonEcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonEcalPFClusterIsolationProducer_*_*', 
        'keep *_electronEcalPFClusterIsolationProducer_*_*', 
        'keep *_photonHcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonHcalPFClusterIsolationProducer_*_*', 
        'keep *_electronHcalPFClusterIsolationProducer_*_*', 
        'drop *_egmGsfElectronIDs_*_*', 
        'drop *_egmPhotonIDs_*_*', 
        'keep recoPhotonCores_gedPhotonCore_*_*', 
        'keep recoPhotons_gedPhotons_*_*', 
        'keep *_particleBasedIsolation_*_*', 
        'drop *_gedPhotons_valMapPFEgammaCandToPhoton_*', 
        'keep recoPhotonCores_photonCore_*_*', 
        'keep recoPhotons_photons_*_*', 
        'keep recoPhotonCores_ootPhotonCore_*_*', 
        'keep recoPhotons_ootPhotons_*_*', 
        'keep recoConversions_conversions_*_*', 
        'keep recoConversions_mustacheConversions_*_*', 
        'drop *_conversions_uncleanedConversions_*', 
        'keep recoConversions_allConversions_*_*', 
        'keep recoConversions_allConversionsOldEG_*_*', 
        'keep recoTracks_ckfOutInTracksFromConversions_*_*', 
        'keep recoTracks_ckfInOutTracksFromConversions_*_*', 
        'keep recoConversions_uncleanedOnlyAllConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep *_PhotonIDProdGED_*_*', 
        'keep *_hfRecoEcalCandidate_*_*', 
        'keep recoRecoEcalCandidates_hfRecoEcalCandidate_*_*', 
        'keep *_hfEMClusters_*_*', 
        'keep *_gedGsfElectronCores_*_*', 
        'keep *_gedGsfElectrons_*_*', 
        'drop CaloTowersSorted_towerMakerPF_*_*', 
        'drop *_pfElectronTranslator_*_*', 
        'keep recoPFRecHits_particleFlowClusterECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHCAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterPS_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHBHE_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitPS_Cleaned_*', 
        'keep recoCaloClusters_particleFlowEGamma_*_*', 
        'keep recoSuperClusters_particleFlowEGamma_*_*', 
        'keep recoCaloClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoSuperClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoConversions_particleFlowEGamma_*_*', 
        'keep recoPFCandidates_particleFlow_*_*', 
        'keep recoPFCandidates_particleFlowTmp_*_*', 
        'drop recoPFCandidates_particleFlowTmp__*', 
        'keep *_particleFlow_electrons_*', 
        'keep *_particleFlow_photons_*', 
        'keep *_particleFlow_muons_*', 
        'keep recoCaloClusters_pfElectronTranslator_*_*', 
        'keep recoPreshowerClusters_pfElectronTranslator_*_*', 
        'keep recoSuperClusters_pfElectronTranslator_*_*', 
        'keep recoCaloClusters_pfPhotonTranslator_*_*', 
        'keep recoPreshowerClusters_pfPhotonTranslator_*_*', 
        'keep recoSuperClusters_pfPhotonTranslator_*_*', 
        'keep recoPhotons_pfPhotonTranslator_*_*', 
        'keep recoPhotonCores_pfPhotonTranslator_*_*', 
        'keep recoConversions_pfPhotonTranslator_*_*', 
        'keep *_particleFlowPtrs_*_*', 
        'keep *_particleFlowTmpPtrs_*_*', 
        'keep *_chargedHadronPFTrackIsolation_*_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep LumiSummary_lumiProducer_*_*', 
        'drop *_hlt*_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep L1AcceptBunchCrossings_scalersRawToDigi_*_*', 
        'keep L1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep Level1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep BeamSpotOnlines_scalersRawToDigi_*_*', 
        'keep DcsStatuss_scalersRawToDigi_*_*', 
        'keep *_tcdsDigis_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'keep *_pfIsolatedElectronsEI_*_*', 
        'keep *_pfIsolatedMuonsEI_*_*', 
        'keep recoPFJets_pfJetsEI_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTagsEI_*_*', 
        'keep recoPFTaus_pfTausEI_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByDecayModeFinding_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByIsolation_*_*', 
        'keep recoPFMETs_pfMetEI_*_*' ) )
)

process.AODSIMEventContent = cms.PSet(
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    eventAutoFlushCompressedSize = cms.untracked.int32(31457280),
    outputCommands = cms.untracked.vstring( ('drop *', 
        'drop *', 
        'keep ClusterSummary_clusterSummaryProducer_*_*', 
        'keep *_dt4DSegments_*_*', 
        'keep *_dt4DCosmicSegments_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*', 
        'keep *_castorreco_*_*', 
        'keep *_reducedHcalRecHits_*_*', 
        'keep HcalUnpackerReport_castorDigis_*_*', 
        'keep HcalUnpackerReport_hcalDigiAlCaMB_*_*', 
        'keep HcalUnpackerReport_hcalDigis_*_*', 
        'keep *_selectDigi_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEB_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEE_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsES_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep recoCaloClusters_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_hybridSuperClusters_uncleanOnlyHybridSuperClusters_*', 
        'keep recoCaloClusters_multi5x5SuperClusters_multi5x5EndcapBasicClusters_*', 
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_*', 
        'keep recoSuperClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoCaloClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoSuperClusters_particleFlowSuperClusterOOTECAL_*_*', 
        'keep recoCaloClusters_particleFlowSuperClusterOOTECAL_*_*', 
        'keep recoTracks_GsfGlobalElectronTest_*_*', 
        'keep recoGsfTracks_electronGsfTracks_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_conversionStepTracks_*_*', 
        'keep recoTracks_beamhaloTracks_*_*', 
        'keep recoTracks_ctfPixelLess_*_*', 
        'keep *_dedxHarmonic2_*_*', 
        'keep *_dedxPixelHarmonic2_*_*', 
        'keep *_dedxHitInfo_*_*', 
        'keep *_trackExtrapolator_*_*', 
        'keep *_generalTracks_MVAValues_*', 
        'keep *_generalTracks_MVAVals_*', 
        'keep *_ak4CaloJets_*_*', 
        'keep *_ak4PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHSSoftDrop_*_*', 
        'keep *_cmsTopTagPFJetsCHS_*_*', 
        'keep *_ak4PFJets_*_*', 
        'keep *_JetPlusTrackZSPCorJetAntiKt4_*_*', 
        'keep *_ak4TrackJets_*_*', 
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', 
        'keep *_caloTowers_*_*', 
        'keep *_CastorTowerReco_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertex_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertexPF_*_*', 
        'keep *_ak4JetTracksAssociatorExplicit_*_*', 
        'keep *_ak4JetExtender_*_*', 
        'keep *_ak4JetID_*_*', 
        'keep *_ak5CastorJets_*_*', 
        'keep *_ak5CastorJetID_*_*', 
        'keep *_ak7CastorJets_*_*', 
        'keep *_ak7CastorJetID_*_*', 
        'keep *_fixedGridRhoAll_*_*', 
        'keep *_fixedGridRhoFastjetAll_*_*', 
        'keep *_fixedGridRhoFastjetAllTmp_*_*', 
        'keep *_fixedGridRhoFastjetCentral_*_*', 
        'keep *_fixedGridRhoFastjetAllCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralChargedPileUp_*_*', 
        'keep *_fixedGridRhoFastjetCentralNeutral_*_*', 
        'drop doubles_*Jets_rhos_*', 
        'drop doubles_*Jets_sigmas_*', 
        'keep *_ak8PFJetsCHSSoftDropMass_*_*', 
        'keep recoCaloMETs_caloMet_*_*', 
        'keep recoCaloMETs_caloMetBE_*_*', 
        'keep recoCaloMETs_caloMetBEFO_*_*', 
        'keep recoCaloMETs_caloMetM_*_*', 
        'keep recoPFMETs_pfMet_*_*', 
        'keep recoPFMETs_pfChMet_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_*_*', 
        'drop recoHcalNoiseRBXs_*_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*', 
        'keep recoGlobalHaloData_GlobalHaloData_*_*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep *_muons_*_*', 
        'keep *_particleFlow_muons_*', 
        'drop *_muons_muons1stStep2muonsMap_*', 
        'drop recoIsoDepositedmValueMap_muons_*_*', 
        'drop doubleedmValueMap_muons_muPFIso*_*', 
        'keep recoTracks_standAloneMuons_*_*', 
        'keep recoTrackExtras_standAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_standAloneMuons_*_*', 
        'keep recoTracks_globalMuons_*_*', 
        'keep recoTrackExtras_globalMuons_*_*', 
        'keep recoTracks_tevMuons_*_*', 
        'keep recoTrackExtras_tevMuons_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_displacedTracks_*_*', 
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_*_*', 
        'keep recoTracks_displacedGlobalMuons_*_*', 
        'keep recoTrackExtras_displacedGlobalMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedGlobalMuons_*_*', 
        'keep recoTracks_cosmicMuons_*_*', 
        'keep recoMuons_muonsFromCosmics_*_*', 
        'keep recoTracks_cosmicMuons1Leg_*_*', 
        'keep recoMuons_muonsFromCosmics1Leg_*_*', 
        'keep recoTracks_refittedStandAloneMuons_*_*', 
        'keep recoTrackExtras_refittedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_refittedStandAloneMuons_*_*', 
        'keep recoTracks_displacedStandAloneMuons__*', 
        'keep recoTrackExtras_displacedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedStandAloneMuons_*_*', 
        'keep *_softPFElectronBJetTags_*_*', 
        'keep *_softPFMuonBJetTags_*_*', 
        'keep *_pfTrackCountingHighEffBJetTags_*_*', 
        'keep *_pfJetProbabilityBJetTags_*_*', 
        'keep *_pfJetBProbabilityBJetTags_*_*', 
        'keep *_pfSimpleSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfCombinedSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfGhostTrackBJetTags_*_*', 
        'keep *_pfCombinedMVAV2BJetTags_*_*', 
        'keep *_inclusiveCandidateSecondaryVertices_*_*', 
        'keep *_inclusiveCandidateSecondaryVerticesCvsL_*_*', 
        'keep *_pfCombinedCvsLJetTags_*_*', 
        'keep *_pfCombinedCvsBJetTags_*_*', 
        'keep *_pfChargeBJetTags_*_*', 
        'keep *_pfDeepCSVJetTags_*_*', 
        'keep *_pfDeepCMVAJetTags_*_*', 
        'keep recoRecoTauPiZeros_hpsPFTauProducer_pizeros_*', 
        'keep recoPFTaus_hpsPFTauProducer_*_*', 
        'keep *_hpsPFTauDiscriminationByDeadECALElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFinding_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingNewDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingOldDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseChargedIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByTightMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSum_*_*', 
        'keep *_hpsPFTauPUcorrPtSum_*_*', 
        'keep *_hpsPFTauChargedIsoPtSum_*_*', 
        'keep *_hpsPFTauTransverseImpactParameters_*_*', 
        'keep *_hpsPFTauFootprintCorrection_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeight_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6rawElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6LooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6MediumElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6TightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VTightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauChargedIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeightdR03_*_*', 
        'keep *_hpsPFTauFootprintCorrectiondR03_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalConedR03_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep  *_offlinePrimaryVertices__*', 
        'keep *_offlinePrimaryVerticesWithBS_*_*', 
        'keep *_offlinePrimaryVerticesFromCosmicTracks_*_*', 
        'keep *_nuclearInteractionMaker_*_*', 
        'keep *_generalV0Candidates_*_*', 
        'keep *_inclusiveSecondaryVertices_*_*', 
        'keep recoGsfElectronCores_gsfElectronCores_*_*', 
        'keep recoGsfElectronCores_gedGsfElectronCores_*_*', 
        'keep recoGsfElectrons_gsfElectrons_*_*', 
        'keep recoGsfElectrons_gedGsfElectrons_*_*', 
        'keep recoGsfElectronCores_uncleanedOnlyGsfElectronCores_*_*', 
        'keep recoGsfElectrons_uncleanedOnlyGsfElectrons_*_*', 
        'keep floatedmValueMap_eidRobustLoose_*_*', 
        'keep floatedmValueMap_eidRobustTight_*_*', 
        'keep floatedmValueMap_eidRobustHighEnergy_*_*', 
        'keep floatedmValueMap_eidLoose_*_*', 
        'keep floatedmValueMap_eidTight_*_*', 
        'keep *_egmGedGsfElectronPFIsolation_*_*', 
        'keep *_photonEcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonEcalPFClusterIsolationProducer_*_*', 
        'keep *_electronEcalPFClusterIsolationProducer_*_*', 
        'keep *_photonHcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonHcalPFClusterIsolationProducer_*_*', 
        'keep *_electronHcalPFClusterIsolationProducer_*_*', 
        'drop *_egmGsfElectronIDs_*_*', 
        'drop *_egmPhotonIDs_*_*', 
        'keep recoPhotonCores_gedPhotonCore_*_*', 
        'keep recoPhotons_gedPhotons_*_*', 
        'keep *_particleBasedIsolation_*_*', 
        'drop *_gedPhotons_valMapPFEgammaCandToPhoton_*', 
        'keep recoPhotonCores_photonCore_*_*', 
        'keep recoPhotons_photons_*_*', 
        'keep recoPhotonCores_ootPhotonCore_*_*', 
        'keep recoPhotons_ootPhotons_*_*', 
        'keep recoConversions_conversions_*_*', 
        'keep recoConversions_mustacheConversions_*_*', 
        'drop *_conversions_uncleanedConversions_*', 
        'keep recoConversions_allConversions_*_*', 
        'keep recoConversions_allConversionsOldEG_*_*', 
        'keep recoTracks_ckfOutInTracksFromConversions_*_*', 
        'keep recoTracks_ckfInOutTracksFromConversions_*_*', 
        'keep recoConversions_uncleanedOnlyAllConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep *_PhotonIDProdGED_*_*', 
        'keep *_hfRecoEcalCandidate_*_*', 
        'keep recoRecoEcalCandidates_hfRecoEcalCandidate_*_*', 
        'keep *_hfEMClusters_*_*', 
        'keep *_gedGsfElectronCores_*_*', 
        'keep *_gedGsfElectrons_*_*', 
        'drop CaloTowersSorted_towerMakerPF_*_*', 
        'drop *_pfElectronTranslator_*_*', 
        'keep recoPFRecHits_particleFlowClusterECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHCAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterPS_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHBHE_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitPS_Cleaned_*', 
        'keep recoCaloClusters_particleFlowEGamma_*_*', 
        'keep recoSuperClusters_particleFlowEGamma_*_*', 
        'keep recoCaloClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoSuperClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoConversions_particleFlowEGamma_*_*', 
        'keep recoPFCandidates_particleFlow_*_*', 
        'keep recoPFCandidates_particleFlowTmp_*_*', 
        'drop recoPFCandidates_particleFlowTmp__*', 
        'keep *_particleFlow_electrons_*', 
        'keep *_particleFlow_photons_*', 
        'keep *_particleFlow_muons_*', 
        'keep recoCaloClusters_pfElectronTranslator_*_*', 
        'keep recoPreshowerClusters_pfElectronTranslator_*_*', 
        'keep recoSuperClusters_pfElectronTranslator_*_*', 
        'keep recoCaloClusters_pfPhotonTranslator_*_*', 
        'keep recoPreshowerClusters_pfPhotonTranslator_*_*', 
        'keep recoSuperClusters_pfPhotonTranslator_*_*', 
        'keep recoPhotons_pfPhotonTranslator_*_*', 
        'keep recoPhotonCores_pfPhotonTranslator_*_*', 
        'keep recoConversions_pfPhotonTranslator_*_*', 
        'keep *_particleFlowPtrs_*_*', 
        'keep *_particleFlowTmpPtrs_*_*', 
        'keep *_chargedHadronPFTrackIsolation_*_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep LumiSummary_lumiProducer_*_*', 
        'drop *_hlt*_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep L1AcceptBunchCrossings_scalersRawToDigi_*_*', 
        'keep L1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep Level1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep BeamSpotOnlines_scalersRawToDigi_*_*', 
        'keep DcsStatuss_scalersRawToDigi_*_*', 
        'keep *_tcdsDigis_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'keep *_pfIsolatedElectronsEI_*_*', 
        'keep *_pfIsolatedMuonsEI_*_*', 
        'keep recoPFJets_pfJetsEI_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTagsEI_*_*', 
        'keep recoPFTaus_pfTausEI_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByDecayModeFinding_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByIsolation_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep *_muonSimClassifier_*_*', 
        'keep *_ak4GenJets_*_*', 
        'keep *_ak8GenJets_*_*', 
        'keep *_ak4GenJetsNoNu_*_*', 
        'keep *_ak8GenJetsNoNu_*_*', 
        'keep *_genParticle_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*' ) )
)

process.BeamSpotAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_offlineBeamSpot_*_*')
)

process.BeamSpotFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_offlineBeamSpot_*_*')
)

process.BeamSpotRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_offlineBeamSpot_*_*')
)

process.CommonEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_logErrorHarvester_*_*')
)

process.CondDB = cms.PSet(
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS')
)

process.DATAMIXEREventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *', 
        'keep CSCDetIdCSCALCTDigiMuonDigiCollection_muonCSCDigis_MuonCSCALCTDigi_*', 
        'keep CSCDetIdCSCCLCTDigiMuonDigiCollection_muonCSCDigis_MuonCSCCLCTDigi_*', 
        'keep CSCDetIdCSCComparatorDigiMuonDigiCollection_muonCSCDigis_MuonCSCComparatorDigi_*', 
        'keep CSCDetIdCSCCorrelatedLCTDigiMuonDigiCollection_csctfDigis_*_*', 
        'keep CSCDetIdCSCCorrelatedLCTDigiMuonDigiCollection_muonCSCDigis_MuonCSCCorrelatedLCTDigi_*', 
        'keep CSCDetIdCSCRPCDigiMuonDigiCollection_muonCSCDigis_MuonCSCRPCDigi_*', 
        'keep CSCDetIdCSCStripDigiMuonDigiCollection_muonCSCDigis_MuonCSCStripDigi_*', 
        'keep CSCDetIdCSCWireDigiMuonDigiCollection_muonCSCDigis_MuonCSCWireDigi_*', 
        'keep DTLayerIdDTDigiMuonDigiCollection_muonDTDigis_*_*', 
        'keep PixelDigiedmDetSetVector_siPixelDigis_*_*', 
        'keep SiStripDigiedmDetSetVector_siStripDigis_*_*', 
        'keep RPCDetIdRPCDigiMuonDigiCollection_muonRPCDigis_*_*', 
        'keep HBHEDataFramesSorted_hcalDigis_*_*', 
        'keep HFDataFramesSorted_hcalDigis_*_*', 
        'keep HODataFramesSorted_hcalDigis_*_*', 
        'keep QIE10DataFrameHcalDataFrameContainer_hcalDigis_*_*', 
        'keep QIE11DataFrameHcalDataFrameContainer_hcalDigis_*_*', 
        'keep ZDCDataFramesSorted_hcalDigis_*_*', 
        'keep CastorDataFramesSorted_castorDigis_*_*', 
        'keep EBDigiCollection_ecalDigis_*_*', 
        'keep EEDigiCollection_ecalDigis_*_*', 
        'keep ESDigiCollection_ecalPreshowerDigis_*_*'),
    splitLevel = cms.untracked.int32(0)
)

process.DQMEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_MEtoEDMConverter_*_*'),
    splitLevel = cms.untracked.int32(0)
)

process.DigiToRawFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*')
)

process.EITopPAGEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_pfIsolatedElectronsEI_*_*', 
        'keep *_pfIsolatedMuonsEI_*_*', 
        'keep recoPFJets_pfJetsEI_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTagsEI_*_*', 
        'keep recoPFTaus_pfTausEI_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByDecayModeFinding_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByIsolation_*_*', 
        'keep recoPFMETs_pfMetEI_*_*')
)

process.EvtScalersAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep L1AcceptBunchCrossings_scalersRawToDigi_*_*', 
        'keep L1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep Level1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep BeamSpotOnlines_scalersRawToDigi_*_*', 
        'keep DcsStatuss_scalersRawToDigi_*_*')
)

process.EvtScalersRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep L1AcceptBunchCrossings_scalersRawToDigi_*_*', 
        'keep L1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep Level1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep BeamSpotOnlines_scalersRawToDigi_*_*', 
        'keep DcsStatuss_scalersRawToDigi_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*')
)

process.FASTPUEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_famosSimHits_*_*', 
        'keep *_MuonSimHits_*_*', 
        'drop *_famosSimHits_VertexTypes_*', 
        'keep *_generalTracksBeforeMixing_*_*', 
        'drop *_generalTracksBeforeMixing_MVAValues_*', 
        'drop *_generalTracksBeforeMixing_QualityMasks_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*')
)

process.FEVTDEBUGEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring( ('drop *', 
        'drop *', 
        'drop *', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep *_g4SimHits_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep CSCDetIdCSCComparatorDigiMuonDigiCollection_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*', 
        'keep EBSrFlagsSorted_simEcalDigis_*_*', 
        'keep EESrFlagsSorted_simEcalDigis_*_*', 
        'keep CrossingFramePlaybackInfoNew_*_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int6stdbitsetstdpairs_*_AffectedAPVList_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'keep *_mix_MergedTrackTruth_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep recoGenJets_*_*_*', 
        'keep *_genParticle_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep *_MEtoEDMConverter_*_*', 
        'keep *_randomEngineStateProducer_*_*', 
        'keep DetIdedmEDCollection_siStripDigis_*_*', 
        'keep DetIdedmEDCollection_siPixelDigis_*_*', 
        'keep PixelFEDChanneledmNewDetSetVector_siPixelDigis_*_*', 
        'keep *_siPixelClusters_*_*', 
        'keep *_siStripClusters_*_*', 
        'keep ClusterSummary_clusterSummaryProducer_*_*', 
        'keep *_dt1DRecHits_*_*', 
        'keep *_dt4DSegments_*_*', 
        'keep *_dt1DCosmicRecHits_*_*', 
        'keep *_dt4DCosmicSegments_*_*', 
        'keep *_csc2DRecHits_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*', 
        'keep *_hbhereco_*_*', 
        'keep *_hbheprereco_*_*', 
        'keep *_hfprereco_*_*', 
        'keep *_hfreco_*_*', 
        'keep *_horeco_*_*', 
        'keep HBHERecHitsSorted_hbherecoMB_*_*', 
        'keep HORecHitsSorted_horecoMB_*_*', 
        'keep HFRecHitsSorted_hfrecoMB_*_*', 
        'keep ZDCDataFramesSorted_hcalDigis_*_*', 
        'keep ZDCDataFramesSorted_castorDigis_*_*', 
        'keep QIE10DataFrameHcalDataFrameContainer_hcalDigis_ZDC_*', 
        'keep ZDCRecHitsSorted_zdcreco_*_*', 
        'keep *_reducedHcalRecHits_*_*', 
        'keep *_castorreco_*_*', 
        'keep HcalUnpackerReport_castorDigis_*_*', 
        'keep HcalUnpackerReport_hcalDigiAlCaMB_*_*', 
        'keep HcalUnpackerReport_hcalDigis_*_*', 
        'keep *_ecalPreshowerRecHit_*_*', 
        'keep *_ecalRecHit_*_*', 
        'keep *_ecalCompactTrigPrim_*_*', 
        'keep *_ecalTPSkim_*_*', 
        'keep EBSrFlagsSorted_ecalDigis__*', 
        'keep EESrFlagsSorted_ecalDigis__*', 
        'keep *_selectDigi_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEE_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEB_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsES_*_*', 
        'keep *_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep *_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_*', 
        'keep *_particleFlowSuperClusterECAL_*_*', 
        'keep *_particleFlowSuperClusterOOTECAL_*_*', 
        'drop recoClusterShapes_*_*_*', 
        'drop recoBasicClustersToOnerecoClusterShapesAssociation_*_*_*', 
        'drop recoBasicClusters_multi5x5BasicClusters_multi5x5BarrelBasicClusters_*', 
        'drop recoSuperClusters_multi5x5SuperClusters_multi5x5BarrelSuperClusters_*', 
        'keep *_CkfElectronCandidates_*_*', 
        'keep *_GsfGlobalElectronTest_*_*', 
        'keep *_electronMergedSeeds_*_*', 
        'keep recoGsfTracks_electronGsfTracks_*_*', 
        'keep recoGsfTrackExtras_electronGsfTracks_*_*', 
        'keep recoTrackExtras_electronGsfTracks_*_*', 
        'keep TrackingRecHitsOwned_electronGsfTracks_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTrackExtras_generalTracks_*_*', 
        'keep TrackingRecHitsOwned_generalTracks_*_*', 
        'keep *_generalTracks_MVAValues_*', 
        'keep *_generalTracks_MVAVals_*', 
        'keep TrackingRecHitsOwned_extraFromSeeds_*_*', 
        'keep uints_extraFromSeeds_*_*', 
        'keep recoTracks_beamhaloTracks_*_*', 
        'keep recoTrackExtras_beamhaloTracks_*_*', 
        'keep TrackingRecHitsOwned_beamhaloTracks_*_*', 
        'keep recoTracks_conversionStepTracks_*_*', 
        'keep recoTrackExtras_conversionStepTracks_*_*', 
        'keep TrackingRecHitsOwned_conversionStepTracks_*_*', 
        'keep *_ctfPixelLess_*_*', 
        'keep *_dedxTruncated40_*_*', 
        'keep *_dedxHitInfo_*_*', 
        'keep *_dedxHarmonic2_*_*', 
        'keep *_dedxPixelHarmonic2_*_*', 
        'keep *_trackExtrapolator_*_*', 
        'keep recoTracks_cosmicDCTracks_*_*', 
        'keep recoTrackExtras_cosmicDCTracks_*_*', 
        'keep TrackingRecHitsOwned_cosmicDCTracks_*_*', 
        'keep *_ak4CaloJets_*_*', 
        'keep *_ak4PFJets_*_*', 
        'keep *_ak4PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHSSoftDrop_*_*', 
        'keep *_cmsTopTagPFJetsCHS_*_*', 
        'keep *_JetPlusTrackZSPCorJetAntiKt4_*_*', 
        'keep *_ak4TrackJets_*_*', 
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', 
        'keep *_caloTowers_*_*', 
        'keep *_towerMaker_*_*', 
        'keep *_CastorTowerReco_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertex_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertexPF_*_*', 
        'keep *_ak4JetTracksAssociatorAtCaloFace_*_*', 
        'keep *_ak4JetTracksAssociatorExplicit_*_*', 
        'keep *_ak4JetExtender_*_*', 
        'keep *_ak4JetID_*_*', 
        'keep *_ak5CastorJets_*_*', 
        'keep *_ak5CastorJetID_*_*', 
        'keep *_ak7CastorJets_*_*', 
        'keep *_ak7CastorJetID_*_*', 
        'keep *_fixedGridRhoAll_*_*', 
        'keep *_fixedGridRhoFastjetAll_*_*', 
        'keep *_fixedGridRhoFastjetAllTmp_*_*', 
        'keep *_fixedGridRhoFastjetAllCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentral_*_*', 
        'keep *_fixedGridRhoFastjetCentralCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralChargedPileUp_*_*', 
        'keep *_fixedGridRhoFastjetCentralNeutral_*_*', 
        'keep *_ak8PFJetsCHSSoftDropMass_*_*', 
        'keep recoCaloMETs_caloMet_*_*', 
        'keep recoCaloMETs_caloMetBE_*_*', 
        'keep recoCaloMETs_caloMetBEFO_*_*', 
        'keep recoCaloMETs_caloMetM_*_*', 
        'keep recoPFMETs_pfMet_*_*', 
        'keep recoPFMETs_pfChMet_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_*_*', 
        'keep recoHcalNoiseRBXs_hcalnoise_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoEcalHaloData_EcalHaloData_*_*', 
        'keep recoGlobalHaloData_GlobalHaloData_*_*', 
        'keep recoHcalHaloData_HcalHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep *_MuonSeed_*_*', 
        'keep *_ancientMuonSeed_*_*', 
        'keep *_displacedMuonSeeds_*_*', 
        'keep TrackingRecHitsOwned_globalMuons_*_*', 
        'keep TrackingRecHitsOwned_tevMuons_*_*', 
        'keep *_CosmicMuonSeed_*_*', 
        'keep recoTrackExtras_cosmicMuons_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons_*_*', 
        'keep recoTrackExtras_cosmicMuons1Leg_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons1Leg_*_*', 
        'keep recoTracks_cosmicsVetoTracks_*_*', 
        'keep *_muons_*_*', 
        'keep *_particleFlow_muons_*', 
        'drop *_muons_muons1stStep2muonsMap_*', 
        'drop recoIsoDepositedmValueMap_muons_*_*', 
        'drop doubleedmValueMap_muons_muPFIso*_*', 
        'keep recoTracks_standAloneMuons_*_*', 
        'keep recoTrackExtras_standAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_standAloneMuons_*_*', 
        'keep recoTracks_globalMuons_*_*', 
        'keep recoTrackExtras_globalMuons_*_*', 
        'keep recoTracks_tevMuons_*_*', 
        'keep recoTrackExtras_tevMuons_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_displacedTracks_*_*', 
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_*_*', 
        'keep recoTracks_displacedGlobalMuons_*_*', 
        'keep recoTrackExtras_displacedGlobalMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedGlobalMuons_*_*', 
        'keep recoTracks_cosmicMuons_*_*', 
        'keep recoMuons_muonsFromCosmics_*_*', 
        'keep recoTracks_cosmicMuons1Leg_*_*', 
        'keep recoMuons_muonsFromCosmics1Leg_*_*', 
        'keep recoTracks_refittedStandAloneMuons_*_*', 
        'keep recoTrackExtras_refittedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_refittedStandAloneMuons_*_*', 
        'keep recoTracks_displacedStandAloneMuons__*', 
        'keep recoTrackExtras_displacedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedStandAloneMuons_*_*', 
        'keep *_muIsoDepositTk_*_*', 
        'keep *_muIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muIsoDepositJets_*_*', 
        'keep *_muGlobalIsoDepositCtfTk_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muGlobalIsoDepositJets_*_*', 
        'keep *_softPFMuonsTagInfos_*_*', 
        'keep *_softPFElectronsTagInfos_*_*', 
        'keep *_softPFElectronBJetTags_*_*', 
        'keep *_softPFMuonBJetTags_*_*', 
        'keep *_pfImpactParameterTagInfos_*_*', 
        'keep *_pfTrackCountingHighEffBJetTags_*_*', 
        'keep *_pfJetProbabilityBJetTags_*_*', 
        'keep *_pfJetBProbabilityBJetTags_*_*', 
        'keep *_pfSecondaryVertexTagInfos_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderTagInfos_*_*', 
        'keep *_pfSimpleSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfCombinedSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfGhostTrackVertexTagInfos_*_*', 
        'keep *_pfGhostTrackBJetTags_*_*', 
        'keep *_pfCombinedMVAV2BJetTags_*_*', 
        'keep *_inclusiveCandidateSecondaryVertices_*_*', 
        'keep *_inclusiveCandidateSecondaryVerticesCvsL_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderCvsLTagInfos_*_*', 
        'keep *_pfCombinedCvsLJetTags_*_*', 
        'keep *_pfCombinedCvsBJetTags_*_*', 
        'keep *_pfChargeBJetTags_*_*', 
        'keep *_pfDeepCSVJetTags_*_*', 
        'keep *_pfDeepCMVAJetTags_*_*', 
        'keep recoRecoTauPiZeros_hpsPFTauProducer_pizeros_*', 
        'keep recoPFTaus_hpsPFTauProducer_*_*', 
        'keep *_hpsPFTauDiscriminationByDeadECALElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFinding_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingNewDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingOldDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseChargedIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByTightMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSum_*_*', 
        'keep *_hpsPFTauPUcorrPtSum_*_*', 
        'keep *_hpsPFTauChargedIsoPtSum_*_*', 
        'keep *_hpsPFTauTransverseImpactParameters_*_*', 
        'keep *_hpsPFTauFootprintCorrection_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeight_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6rawElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6LooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6MediumElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6TightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VTightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauChargedIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeightdR03_*_*', 
        'keep *_hpsPFTauFootprintCorrectiondR03_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalConedR03_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep  *_offlinePrimaryVertices__*', 
        'keep *_offlinePrimaryVerticesWithBS_*_*', 
        'keep *_offlinePrimaryVerticesFromCosmicTracks_*_*', 
        'keep *_nuclearInteractionMaker_*_*', 
        'keep *_generalV0Candidates_*_*', 
        'keep *_inclusiveSecondaryVertices_*_*', 
        'keep recoGsfElectronCores_gsfElectronCores_*_*', 
        'keep recoGsfElectronCores_gedGsfElectronCores_*_*', 
        'keep recoGsfElectrons_gsfElectrons_*_*', 
        'keep recoGsfElectrons_gedGsfElectrons_*_*', 
        'keep recoGsfElectronCores_uncleanedOnlyGsfElectronCores_*_*', 
        'keep recoGsfElectrons_uncleanedOnlyGsfElectrons_*_*', 
        'keep floatedmValueMap_eidRobustLoose_*_*', 
        'keep floatedmValueMap_eidRobustTight_*_*', 
        'keep floatedmValueMap_eidRobustHighEnergy_*_*', 
        'keep floatedmValueMap_eidLoose_*_*', 
        'keep floatedmValueMap_eidTight_*_*', 
        'keep *_egmGedGsfElectronPFIsolation_*_*', 
        'keep *_photonEcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonEcalPFClusterIsolationProducer_*_*', 
        'keep *_electronEcalPFClusterIsolationProducer_*_*', 
        'keep *_photonHcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonHcalPFClusterIsolationProducer_*_*', 
        'keep *_electronHcalPFClusterIsolationProducer_*_*', 
        'drop *_egmGsfElectronIDs_*_*', 
        'drop *_egmPhotonIDs_*_*', 
        'keep *_gedPhotonCore_*_*', 
        'keep *_gedPhotons_*_*', 
        'keep *_particleBasedIsolation_*_*', 
        'keep recoPhotons_mustachePhotons_*_*', 
        'keep recoPhotonCores_mustachePhotonCore_*_*', 
        'keep recoPhotons_ootPhotons_*_*', 
        'keep recoPhotonCores_ootPhotonCore_*_*', 
        'keep recoPhotons_photons_*_*', 
        'keep recoPhotonCores_photonCore_*_*', 
        'keep recoConversions_conversions_*_*', 
        'keep recoConversions_mustacheConversions_*_*', 
        'drop *_conversions_uncleanedConversions_*', 
        'drop *_gedPhotonsTmp_valMapPFEgammaCandToPhoton_*', 
        'keep recoConversions_allConversions_*_*', 
        'keep recoConversions_allConversionsOldEG_*_*', 
        'keep recoTracks_ckfOutInTracksFromConversions_*_*', 
        'keep recoTracks_ckfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfInOutTracksFromConversions_*_*', 
        'keep recoConversions_uncleanedOnlyAllConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep *_PhotonIDProdGED_*_*', 
        'keep recoRecoEcalCandidates_hfRecoEcalCandidate_*_*', 
        'keep *_hfEMClusters_*_*', 
        'keep *_gedGsfElectronCores_*_*', 
        'keep *_gedGsfElectrons_*_*', 
        'keep *_pixelTracks_*_*', 
        'keep *_pixelVertices_*_*', 
        'drop CaloTowersSorted_towerMakerPF_*_*', 
        'keep recoPFRecHits_particleFlowClusterECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHCAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterPS_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHBHE_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitPS_Cleaned_*', 
        'keep recoPFClusters_particleFlowClusterECAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHCAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHO_*_*', 
        'keep recoPFClusters_particleFlowClusterHF_*_*', 
        'keep recoPFClusters_particleFlowClusterPS_*_*', 
        'keep recoPFBlocks_particleFlowBlock_*_*', 
        'keep recoPFCandidates_particleFlowEGamma_*_*', 
        'keep recoCaloClusters_particleFlowEGamma_*_*', 
        'keep recoSuperClusters_particleFlowEGamma_*_*', 
        'keep recoConversions_particleFlowEGamma_*_*', 
        'keep recoPFCandidates_particleFlow_*_*', 
        'keep recoPFCandidates_particleFlowTmp_electrons_*', 
        'keep recoPFCandidates_particleFlowTmp_*_*', 
        'drop recoPFCandidates_particleFlowTmp__*', 
        'keep recoPFDisplacedVertexs_particleFlowDisplacedVertex_*_*', 
        'keep *_pfElectronTranslator_*_*', 
        'keep *_pfPhotonTranslator_*_*', 
        'keep *_particleFlow_electrons_*', 
        'keep *_particleFlow_photons_*', 
        'keep *_particleFlow_muons_*', 
        'keep *_trackerDrivenElectronSeeds_preid_*', 
        'keep *_particleFlowPtrs_*_*', 
        'keep *_particleFlowTmpPtrs_*_*', 
        'keep *_chargedHadronPFTrackIsolation_*_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep L1MuGMTReadoutCollection_gtDigis_*_*', 
        'keep L1GctEmCand*_gctDigis_*_*', 
        'keep L1GctJetCand*_gctDigis_*_*', 
        'keep L1GctEtHad*_gctDigis_*_*', 
        'keep L1GctEtMiss*_gctDigis_*_*', 
        'keep L1GctEtTotal*_gctDigis_*_*', 
        'keep L1GctHtMiss*_gctDigis_*_*', 
        'keep L1GctJetCounts*_gctDigis_*_*', 
        'keep L1GctHFRingEtSums*_gctDigis_*_*', 
        'keep L1GctHFBitCounts*_gctDigis_*_*', 
        'keep LumiDetails_lumiProducer_*_*', 
        'keep LumiSummary_lumiProducer_*_*', 
        'drop *_hlt*_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep *_ak4GenJets_*_*', 
        'keep *_ak8GenJets_*_*', 
        'keep *_ak4GenJetsNoNu_*_*', 
        'keep *_ak8GenJetsNoNu_*_*', 
        'keep *_genParticle_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep SimTracks_g4SimHits_*_*', 
        'keep SimVertexs_g4SimHits_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*', 
        'keep *_muonSimClassifier_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'keep L1AcceptBunchCrossings_scalersRawToDigi_*_*', 
        'keep L1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep Level1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep BeamSpotOnlines_scalersRawToDigi_*_*', 
        'keep DcsStatuss_scalersRawToDigi_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'keep *_pfIsolatedElectronsEI_*_*', 
        'keep *_pfIsolatedMuonsEI_*_*', 
        'keep recoPFJets_pfJetsEI_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTagsEI_*_*', 
        'keep recoPFTaus_pfTausEI_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByDecayModeFinding_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByIsolation_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep *_tcdsDigis_*_*', 
        'keep *_simCscTriggerPrimitiveDigis_*_*', 
        'keep *_simDtTriggerPrimitiveDigis_*_*', 
        'keep *_simRpcTriggerDigis_*_*', 
        'keep *_simRctDigis_*_*', 
        'keep *_simCsctfDigis_*_*', 
        'keep *_simCsctfTrackDigis_*_*', 
        'keep *_simDttfDigis_*_*', 
        'keep *_simGctDigis_*_*', 
        'keep *_simCaloStage1Digis_*_*', 
        'keep *_simCaloStage1FinalDigis_*_*', 
        'keep *_simCaloStage2Layer1Digis_*_*', 
        'keep *_simCaloStage2Digis_*_*', 
        'keep *_simGmtDigis_*_*', 
        'keep *_simBmtfDigis_*_*', 
        'keep *_simOmtfDigis_*_*', 
        'keep *_simEmtfDigis_*_*', 
        'keep *_simGmtStage2Digis_*_*', 
        'keep *_simGtDigis_*_*', 
        'keep *_simGtStage2Digis_*_*', 
        'keep *_cscTriggerPrimitiveDigis_*_*', 
        'keep *_dtTriggerPrimitiveDigis_*_*', 
        'keep *_rpcTriggerDigis_*_*', 
        'keep *_rctDigis_*_*', 
        'keep *_csctfDigis_*_*', 
        'keep *_csctfTrackDigis_*_*', 
        'keep *_dttfDigis_*_*', 
        'keep *_gctDigis_*_*', 
        'keep *_gmtDigis_*_*', 
        'keep *_gtDigis_*_*', 
        'keep *_gtEvmDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep LumiDetails_lumiProducer_*_*', 
        'keep LumiSummary_lumiProducer_*_*', 
        'drop *_trackingtruthprod_*_*', 
        'drop *_electrontruth_*_*', 
        'keep *_mix_MergedTrackTruth_*', 
        'keep CrossingFramePlaybackInfoNew_*_*_*', 
        'keep *_simSiPixelDigis_*_*', 
        'keep *_simSiStripDigis_*_*', 
        'drop *_mix_simSiPixelDigis*_*', 
        'drop *_mix_simSiStripDigis*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep *_trackingParticleRecoTrackAsssociation_*_*', 
        'keep *_assoc2secStepTk_*_*', 
        'keep *_assoc2thStepTk_*_*', 
        'keep *_assoc2GsfTracks_*_*', 
        'keep *_assocOutInConversionTracks_*_*', 
        'keep *_assocInOutConversionTracks_*_*', 
        'keep *_simMuonCSCDigis_*_*', 
        'keep *_simMuonDTDigis_*_*', 
        'keep *_simMuonRPCDigis_*_*', 
        'keep *_simEcalDigis_*_*', 
        'keep *_simEcalPreshowerDigis_*_*', 
        'keep *_simEcalTriggerPrimitiveDigis_*_*', 
        'keep *_simEcalEBTriggerPrimitiveDigis_*_*', 
        'keep *_simHcalDigis_*_*', 
        'keep ZDCDataFramesSorted_simHcalUnsuppressedDigis_*_*', 
        'drop ZDCDataFramesSorted_mix_simHcalUnsuppressedDigis*_*', 
        'keep *_simHcalTriggerPrimitiveDigis_*_*', 
        'keep *_mix_HcalSamples_*', 
        'keep *_mixData_HcalSamples_*', 
        'keep *_mix_HcalHits_*', 
        'keep *_mixData_HcalHits_*', 
        'keep *_DMHcalTriggerPrimitiveDigis_*_*' ) ),
    splitLevel = cms.untracked.int32(0)
)

process.FEVTDEBUGHLTEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring( ('drop *', 
        'drop *', 
        'drop *', 
        'drop *', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep *_g4SimHits_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep CSCDetIdCSCComparatorDigiMuonDigiCollection_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*', 
        'keep EBSrFlagsSorted_simEcalDigis_*_*', 
        'keep EESrFlagsSorted_simEcalDigis_*_*', 
        'keep CrossingFramePlaybackInfoNew_*_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int6stdbitsetstdpairs_*_AffectedAPVList_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'keep *_mix_MergedTrackTruth_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep recoGenJets_*_*_*', 
        'keep *_genParticle_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep *_MEtoEDMConverter_*_*', 
        'keep *_randomEngineStateProducer_*_*', 
        'keep DetIdedmEDCollection_siStripDigis_*_*', 
        'keep DetIdedmEDCollection_siPixelDigis_*_*', 
        'keep PixelFEDChanneledmNewDetSetVector_siPixelDigis_*_*', 
        'keep *_siPixelClusters_*_*', 
        'keep *_siStripClusters_*_*', 
        'keep ClusterSummary_clusterSummaryProducer_*_*', 
        'keep *_dt1DRecHits_*_*', 
        'keep *_dt4DSegments_*_*', 
        'keep *_dt1DCosmicRecHits_*_*', 
        'keep *_dt4DCosmicSegments_*_*', 
        'keep *_csc2DRecHits_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*', 
        'keep *_hbhereco_*_*', 
        'keep *_hbheprereco_*_*', 
        'keep *_hfprereco_*_*', 
        'keep *_hfreco_*_*', 
        'keep *_horeco_*_*', 
        'keep HBHERecHitsSorted_hbherecoMB_*_*', 
        'keep HORecHitsSorted_horecoMB_*_*', 
        'keep HFRecHitsSorted_hfrecoMB_*_*', 
        'keep ZDCDataFramesSorted_hcalDigis_*_*', 
        'keep ZDCDataFramesSorted_castorDigis_*_*', 
        'keep QIE10DataFrameHcalDataFrameContainer_hcalDigis_ZDC_*', 
        'keep ZDCRecHitsSorted_zdcreco_*_*', 
        'keep *_reducedHcalRecHits_*_*', 
        'keep *_castorreco_*_*', 
        'keep HcalUnpackerReport_castorDigis_*_*', 
        'keep HcalUnpackerReport_hcalDigiAlCaMB_*_*', 
        'keep HcalUnpackerReport_hcalDigis_*_*', 
        'keep *_ecalPreshowerRecHit_*_*', 
        'keep *_ecalRecHit_*_*', 
        'keep *_ecalCompactTrigPrim_*_*', 
        'keep *_ecalTPSkim_*_*', 
        'keep EBSrFlagsSorted_ecalDigis__*', 
        'keep EESrFlagsSorted_ecalDigis__*', 
        'keep *_selectDigi_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEE_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEB_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsES_*_*', 
        'keep *_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep *_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_*', 
        'keep *_particleFlowSuperClusterECAL_*_*', 
        'keep *_particleFlowSuperClusterOOTECAL_*_*', 
        'drop recoClusterShapes_*_*_*', 
        'drop recoBasicClustersToOnerecoClusterShapesAssociation_*_*_*', 
        'drop recoBasicClusters_multi5x5BasicClusters_multi5x5BarrelBasicClusters_*', 
        'drop recoSuperClusters_multi5x5SuperClusters_multi5x5BarrelSuperClusters_*', 
        'keep *_CkfElectronCandidates_*_*', 
        'keep *_GsfGlobalElectronTest_*_*', 
        'keep *_electronMergedSeeds_*_*', 
        'keep recoGsfTracks_electronGsfTracks_*_*', 
        'keep recoGsfTrackExtras_electronGsfTracks_*_*', 
        'keep recoTrackExtras_electronGsfTracks_*_*', 
        'keep TrackingRecHitsOwned_electronGsfTracks_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTrackExtras_generalTracks_*_*', 
        'keep TrackingRecHitsOwned_generalTracks_*_*', 
        'keep *_generalTracks_MVAValues_*', 
        'keep *_generalTracks_MVAVals_*', 
        'keep TrackingRecHitsOwned_extraFromSeeds_*_*', 
        'keep uints_extraFromSeeds_*_*', 
        'keep recoTracks_beamhaloTracks_*_*', 
        'keep recoTrackExtras_beamhaloTracks_*_*', 
        'keep TrackingRecHitsOwned_beamhaloTracks_*_*', 
        'keep recoTracks_conversionStepTracks_*_*', 
        'keep recoTrackExtras_conversionStepTracks_*_*', 
        'keep TrackingRecHitsOwned_conversionStepTracks_*_*', 
        'keep *_ctfPixelLess_*_*', 
        'keep *_dedxTruncated40_*_*', 
        'keep *_dedxHitInfo_*_*', 
        'keep *_dedxHarmonic2_*_*', 
        'keep *_dedxPixelHarmonic2_*_*', 
        'keep *_trackExtrapolator_*_*', 
        'keep recoTracks_cosmicDCTracks_*_*', 
        'keep recoTrackExtras_cosmicDCTracks_*_*', 
        'keep TrackingRecHitsOwned_cosmicDCTracks_*_*', 
        'keep *_ak4CaloJets_*_*', 
        'keep *_ak4PFJets_*_*', 
        'keep *_ak4PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHSSoftDrop_*_*', 
        'keep *_cmsTopTagPFJetsCHS_*_*', 
        'keep *_JetPlusTrackZSPCorJetAntiKt4_*_*', 
        'keep *_ak4TrackJets_*_*', 
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', 
        'keep *_caloTowers_*_*', 
        'keep *_towerMaker_*_*', 
        'keep *_CastorTowerReco_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertex_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertexPF_*_*', 
        'keep *_ak4JetTracksAssociatorAtCaloFace_*_*', 
        'keep *_ak4JetTracksAssociatorExplicit_*_*', 
        'keep *_ak4JetExtender_*_*', 
        'keep *_ak4JetID_*_*', 
        'keep *_ak5CastorJets_*_*', 
        'keep *_ak5CastorJetID_*_*', 
        'keep *_ak7CastorJets_*_*', 
        'keep *_ak7CastorJetID_*_*', 
        'keep *_fixedGridRhoAll_*_*', 
        'keep *_fixedGridRhoFastjetAll_*_*', 
        'keep *_fixedGridRhoFastjetAllTmp_*_*', 
        'keep *_fixedGridRhoFastjetAllCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentral_*_*', 
        'keep *_fixedGridRhoFastjetCentralCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralChargedPileUp_*_*', 
        'keep *_fixedGridRhoFastjetCentralNeutral_*_*', 
        'keep *_ak8PFJetsCHSSoftDropMass_*_*', 
        'keep recoCaloMETs_caloMet_*_*', 
        'keep recoCaloMETs_caloMetBE_*_*', 
        'keep recoCaloMETs_caloMetBEFO_*_*', 
        'keep recoCaloMETs_caloMetM_*_*', 
        'keep recoPFMETs_pfMet_*_*', 
        'keep recoPFMETs_pfChMet_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_*_*', 
        'keep recoHcalNoiseRBXs_hcalnoise_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoEcalHaloData_EcalHaloData_*_*', 
        'keep recoGlobalHaloData_GlobalHaloData_*_*', 
        'keep recoHcalHaloData_HcalHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep *_MuonSeed_*_*', 
        'keep *_ancientMuonSeed_*_*', 
        'keep *_displacedMuonSeeds_*_*', 
        'keep TrackingRecHitsOwned_globalMuons_*_*', 
        'keep TrackingRecHitsOwned_tevMuons_*_*', 
        'keep *_CosmicMuonSeed_*_*', 
        'keep recoTrackExtras_cosmicMuons_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons_*_*', 
        'keep recoTrackExtras_cosmicMuons1Leg_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons1Leg_*_*', 
        'keep recoTracks_cosmicsVetoTracks_*_*', 
        'keep *_muons_*_*', 
        'keep *_particleFlow_muons_*', 
        'drop *_muons_muons1stStep2muonsMap_*', 
        'drop recoIsoDepositedmValueMap_muons_*_*', 
        'drop doubleedmValueMap_muons_muPFIso*_*', 
        'keep recoTracks_standAloneMuons_*_*', 
        'keep recoTrackExtras_standAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_standAloneMuons_*_*', 
        'keep recoTracks_globalMuons_*_*', 
        'keep recoTrackExtras_globalMuons_*_*', 
        'keep recoTracks_tevMuons_*_*', 
        'keep recoTrackExtras_tevMuons_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_displacedTracks_*_*', 
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_*_*', 
        'keep recoTracks_displacedGlobalMuons_*_*', 
        'keep recoTrackExtras_displacedGlobalMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedGlobalMuons_*_*', 
        'keep recoTracks_cosmicMuons_*_*', 
        'keep recoMuons_muonsFromCosmics_*_*', 
        'keep recoTracks_cosmicMuons1Leg_*_*', 
        'keep recoMuons_muonsFromCosmics1Leg_*_*', 
        'keep recoTracks_refittedStandAloneMuons_*_*', 
        'keep recoTrackExtras_refittedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_refittedStandAloneMuons_*_*', 
        'keep recoTracks_displacedStandAloneMuons__*', 
        'keep recoTrackExtras_displacedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedStandAloneMuons_*_*', 
        'keep *_muIsoDepositTk_*_*', 
        'keep *_muIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muIsoDepositJets_*_*', 
        'keep *_muGlobalIsoDepositCtfTk_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muGlobalIsoDepositJets_*_*', 
        'keep *_softPFMuonsTagInfos_*_*', 
        'keep *_softPFElectronsTagInfos_*_*', 
        'keep *_softPFElectronBJetTags_*_*', 
        'keep *_softPFMuonBJetTags_*_*', 
        'keep *_pfImpactParameterTagInfos_*_*', 
        'keep *_pfTrackCountingHighEffBJetTags_*_*', 
        'keep *_pfJetProbabilityBJetTags_*_*', 
        'keep *_pfJetBProbabilityBJetTags_*_*', 
        'keep *_pfSecondaryVertexTagInfos_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderTagInfos_*_*', 
        'keep *_pfSimpleSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfCombinedSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfGhostTrackVertexTagInfos_*_*', 
        'keep *_pfGhostTrackBJetTags_*_*', 
        'keep *_pfCombinedMVAV2BJetTags_*_*', 
        'keep *_inclusiveCandidateSecondaryVertices_*_*', 
        'keep *_inclusiveCandidateSecondaryVerticesCvsL_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderCvsLTagInfos_*_*', 
        'keep *_pfCombinedCvsLJetTags_*_*', 
        'keep *_pfCombinedCvsBJetTags_*_*', 
        'keep *_pfChargeBJetTags_*_*', 
        'keep *_pfDeepCSVJetTags_*_*', 
        'keep *_pfDeepCMVAJetTags_*_*', 
        'keep recoRecoTauPiZeros_hpsPFTauProducer_pizeros_*', 
        'keep recoPFTaus_hpsPFTauProducer_*_*', 
        'keep *_hpsPFTauDiscriminationByDeadECALElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFinding_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingNewDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingOldDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseChargedIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByTightMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSum_*_*', 
        'keep *_hpsPFTauPUcorrPtSum_*_*', 
        'keep *_hpsPFTauChargedIsoPtSum_*_*', 
        'keep *_hpsPFTauTransverseImpactParameters_*_*', 
        'keep *_hpsPFTauFootprintCorrection_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeight_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6rawElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6LooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6MediumElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6TightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VTightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauChargedIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeightdR03_*_*', 
        'keep *_hpsPFTauFootprintCorrectiondR03_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalConedR03_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep  *_offlinePrimaryVertices__*', 
        'keep *_offlinePrimaryVerticesWithBS_*_*', 
        'keep *_offlinePrimaryVerticesFromCosmicTracks_*_*', 
        'keep *_nuclearInteractionMaker_*_*', 
        'keep *_generalV0Candidates_*_*', 
        'keep *_inclusiveSecondaryVertices_*_*', 
        'keep recoGsfElectronCores_gsfElectronCores_*_*', 
        'keep recoGsfElectronCores_gedGsfElectronCores_*_*', 
        'keep recoGsfElectrons_gsfElectrons_*_*', 
        'keep recoGsfElectrons_gedGsfElectrons_*_*', 
        'keep recoGsfElectronCores_uncleanedOnlyGsfElectronCores_*_*', 
        'keep recoGsfElectrons_uncleanedOnlyGsfElectrons_*_*', 
        'keep floatedmValueMap_eidRobustLoose_*_*', 
        'keep floatedmValueMap_eidRobustTight_*_*', 
        'keep floatedmValueMap_eidRobustHighEnergy_*_*', 
        'keep floatedmValueMap_eidLoose_*_*', 
        'keep floatedmValueMap_eidTight_*_*', 
        'keep *_egmGedGsfElectronPFIsolation_*_*', 
        'keep *_photonEcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonEcalPFClusterIsolationProducer_*_*', 
        'keep *_electronEcalPFClusterIsolationProducer_*_*', 
        'keep *_photonHcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonHcalPFClusterIsolationProducer_*_*', 
        'keep *_electronHcalPFClusterIsolationProducer_*_*', 
        'drop *_egmGsfElectronIDs_*_*', 
        'drop *_egmPhotonIDs_*_*', 
        'keep *_gedPhotonCore_*_*', 
        'keep *_gedPhotons_*_*', 
        'keep *_particleBasedIsolation_*_*', 
        'keep recoPhotons_mustachePhotons_*_*', 
        'keep recoPhotonCores_mustachePhotonCore_*_*', 
        'keep recoPhotons_ootPhotons_*_*', 
        'keep recoPhotonCores_ootPhotonCore_*_*', 
        'keep recoPhotons_photons_*_*', 
        'keep recoPhotonCores_photonCore_*_*', 
        'keep recoConversions_conversions_*_*', 
        'keep recoConversions_mustacheConversions_*_*', 
        'drop *_conversions_uncleanedConversions_*', 
        'drop *_gedPhotonsTmp_valMapPFEgammaCandToPhoton_*', 
        'keep recoConversions_allConversions_*_*', 
        'keep recoConversions_allConversionsOldEG_*_*', 
        'keep recoTracks_ckfOutInTracksFromConversions_*_*', 
        'keep recoTracks_ckfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfInOutTracksFromConversions_*_*', 
        'keep recoConversions_uncleanedOnlyAllConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep *_PhotonIDProdGED_*_*', 
        'keep recoRecoEcalCandidates_hfRecoEcalCandidate_*_*', 
        'keep *_hfEMClusters_*_*', 
        'keep *_gedGsfElectronCores_*_*', 
        'keep *_gedGsfElectrons_*_*', 
        'keep *_pixelTracks_*_*', 
        'keep *_pixelVertices_*_*', 
        'drop CaloTowersSorted_towerMakerPF_*_*', 
        'keep recoPFRecHits_particleFlowClusterECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHCAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterPS_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHBHE_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitPS_Cleaned_*', 
        'keep recoPFClusters_particleFlowClusterECAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHCAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHO_*_*', 
        'keep recoPFClusters_particleFlowClusterHF_*_*', 
        'keep recoPFClusters_particleFlowClusterPS_*_*', 
        'keep recoPFBlocks_particleFlowBlock_*_*', 
        'keep recoPFCandidates_particleFlowEGamma_*_*', 
        'keep recoCaloClusters_particleFlowEGamma_*_*', 
        'keep recoSuperClusters_particleFlowEGamma_*_*', 
        'keep recoConversions_particleFlowEGamma_*_*', 
        'keep recoPFCandidates_particleFlow_*_*', 
        'keep recoPFCandidates_particleFlowTmp_electrons_*', 
        'keep recoPFCandidates_particleFlowTmp_*_*', 
        'drop recoPFCandidates_particleFlowTmp__*', 
        'keep recoPFDisplacedVertexs_particleFlowDisplacedVertex_*_*', 
        'keep *_pfElectronTranslator_*_*', 
        'keep *_pfPhotonTranslator_*_*', 
        'keep *_particleFlow_electrons_*', 
        'keep *_particleFlow_photons_*', 
        'keep *_particleFlow_muons_*', 
        'keep *_trackerDrivenElectronSeeds_preid_*', 
        'keep *_particleFlowPtrs_*_*', 
        'keep *_particleFlowTmpPtrs_*_*', 
        'keep *_chargedHadronPFTrackIsolation_*_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep L1MuGMTReadoutCollection_gtDigis_*_*', 
        'keep L1GctEmCand*_gctDigis_*_*', 
        'keep L1GctJetCand*_gctDigis_*_*', 
        'keep L1GctEtHad*_gctDigis_*_*', 
        'keep L1GctEtMiss*_gctDigis_*_*', 
        'keep L1GctEtTotal*_gctDigis_*_*', 
        'keep L1GctHtMiss*_gctDigis_*_*', 
        'keep L1GctJetCounts*_gctDigis_*_*', 
        'keep L1GctHFRingEtSums*_gctDigis_*_*', 
        'keep L1GctHFBitCounts*_gctDigis_*_*', 
        'keep LumiDetails_lumiProducer_*_*', 
        'keep LumiSummary_lumiProducer_*_*', 
        'drop *_hlt*_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep *_ak4GenJets_*_*', 
        'keep *_ak8GenJets_*_*', 
        'keep *_ak4GenJetsNoNu_*_*', 
        'keep *_ak8GenJetsNoNu_*_*', 
        'keep *_genParticle_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep SimTracks_g4SimHits_*_*', 
        'keep SimVertexs_g4SimHits_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*', 
        'keep *_muonSimClassifier_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'keep L1AcceptBunchCrossings_scalersRawToDigi_*_*', 
        'keep L1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep Level1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep BeamSpotOnlines_scalersRawToDigi_*_*', 
        'keep DcsStatuss_scalersRawToDigi_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'keep *_pfIsolatedElectronsEI_*_*', 
        'keep *_pfIsolatedMuonsEI_*_*', 
        'keep recoPFJets_pfJetsEI_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTagsEI_*_*', 
        'keep recoPFTaus_pfTausEI_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByDecayModeFinding_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByIsolation_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep *_tcdsDigis_*_*', 
        'keep *_simCscTriggerPrimitiveDigis_*_*', 
        'keep *_simDtTriggerPrimitiveDigis_*_*', 
        'keep *_simRpcTriggerDigis_*_*', 
        'keep *_simRctDigis_*_*', 
        'keep *_simCsctfDigis_*_*', 
        'keep *_simCsctfTrackDigis_*_*', 
        'keep *_simDttfDigis_*_*', 
        'keep *_simGctDigis_*_*', 
        'keep *_simCaloStage1Digis_*_*', 
        'keep *_simCaloStage1FinalDigis_*_*', 
        'keep *_simCaloStage2Layer1Digis_*_*', 
        'keep *_simCaloStage2Digis_*_*', 
        'keep *_simGmtDigis_*_*', 
        'keep *_simBmtfDigis_*_*', 
        'keep *_simOmtfDigis_*_*', 
        'keep *_simEmtfDigis_*_*', 
        'keep *_simGmtStage2Digis_*_*', 
        'keep *_simGtDigis_*_*', 
        'keep *_simGtStage2Digis_*_*', 
        'keep *_cscTriggerPrimitiveDigis_*_*', 
        'keep *_dtTriggerPrimitiveDigis_*_*', 
        'keep *_rpcTriggerDigis_*_*', 
        'keep *_rctDigis_*_*', 
        'keep *_csctfDigis_*_*', 
        'keep *_csctfTrackDigis_*_*', 
        'keep *_dttfDigis_*_*', 
        'keep *_gctDigis_*_*', 
        'keep *_gmtDigis_*_*', 
        'keep *_gtDigis_*_*', 
        'keep *_gtEvmDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep LumiDetails_lumiProducer_*_*', 
        'keep LumiSummary_lumiProducer_*_*', 
        'drop *_trackingtruthprod_*_*', 
        'drop *_electrontruth_*_*', 
        'keep *_mix_MergedTrackTruth_*', 
        'keep CrossingFramePlaybackInfoNew_*_*_*', 
        'keep *_simSiPixelDigis_*_*', 
        'keep *_simSiStripDigis_*_*', 
        'drop *_mix_simSiPixelDigis*_*', 
        'drop *_mix_simSiStripDigis*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep *_trackingParticleRecoTrackAsssociation_*_*', 
        'keep *_assoc2secStepTk_*_*', 
        'keep *_assoc2thStepTk_*_*', 
        'keep *_assoc2GsfTracks_*_*', 
        'keep *_assocOutInConversionTracks_*_*', 
        'keep *_assocInOutConversionTracks_*_*', 
        'keep *_simMuonCSCDigis_*_*', 
        'keep *_simMuonDTDigis_*_*', 
        'keep *_simMuonRPCDigis_*_*', 
        'keep *_simEcalDigis_*_*', 
        'keep *_simEcalPreshowerDigis_*_*', 
        'keep *_simEcalTriggerPrimitiveDigis_*_*', 
        'keep *_simEcalEBTriggerPrimitiveDigis_*_*', 
        'keep *_simHcalDigis_*_*', 
        'keep ZDCDataFramesSorted_simHcalUnsuppressedDigis_*_*', 
        'drop ZDCDataFramesSorted_mix_simHcalUnsuppressedDigis*_*', 
        'keep *_simHcalTriggerPrimitiveDigis_*_*', 
        'keep *_mix_HcalSamples_*', 
        'keep *_mixData_HcalSamples_*', 
        'keep *_mix_HcalHits_*', 
        'keep *_mixData_HcalHits_*', 
        'keep *_DMHcalTriggerPrimitiveDigis_*_*', 
        'drop *_hlt*_*_*', 
        'keep *_hltAK4CaloJetsCorrectedIDPassed_*_*', 
        'keep *_hltAK4CaloJetsIDPassed_*_*', 
        'keep *_hltAK4CaloJets_*_*', 
        'keep *_hltAK4PFJetsCorrected_*_*', 
        'keep *_hltAK4PFJetsForTaus_*_*', 
        'keep *_hltAK4PFJets_*_*', 
        'keep *_hltAlCaEtaEBRechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaEtaEBRechitsToDigis_*_*', 
        'keep *_hltAlCaEtaEERechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaEtaEERechitsToDigis_*_*', 
        'keep *_hltAlCaEtaRecHitsFilterEBonly_*_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegionalLowPU_etaEcalRecHitsES_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegional_etaEcalRecHitsES_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonly_*_*', 
        'keep *_hltAlCaEtaRecHitsFilter_*_*', 
        'keep *_hltAlCaPhiSymStream_*_*', 
        'keep *_hltAlCaPi0EBRechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaPi0EBRechitsToDigis_*_*', 
        'keep *_hltAlCaPi0EERechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaPi0EERechitsToDigis_*_*', 
        'keep *_hltAlCaPi0RecHitsFilterEBonly_*_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegionalLowPU_pi0EcalRecHitsES_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegional_pi0EcalRecHitsES_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonly_*_*', 
        'keep *_hltAlCaPi0RecHitsFilter_*_*', 
        'keep *_hltBLifetimeL25AssociatorbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL25BJetTagsbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL25TagInfosbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3AssociatorbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3BJetTagsbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3TagInfosbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeRegionalCtfWithMaterialTracksbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBSoftMuonDiJet110Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet110Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet20Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet20Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet40Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet40Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet70Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet70Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonMu5L3_*_*', 
        'keep *_hltCSVJetTagSingleTopEle27_*_*', 
        'keep *_hltCSVJetTagSingleTopIsoMu24_*_*', 
        'keep *_hltCaloJetCorrectedRegional_*_*', 
        'keep *_hltCaloJetCorrected_*_*', 
        'keep *_hltCaloJetL1FastJetCorrected_*_*', 
        'keep *_hltCaloStage2Digis_*_*', 
        'keep *_hltCleanedCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCleanedHiCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCombinedSecondaryVertexBJetTagsCalo_*_*', 
        'keep *_hltCombinedSecondaryVertexBJetTagsPF_*_*', 
        'keep *_hltConvPFTausTightIsoTrackFindingIsolation_*_*', 
        'keep *_hltConvPFTausTightIsoTrackFinding_*_*', 
        'keep *_hltConvPFTausTightIsoTrackPt5Isolation_*_*', 
        'keep *_hltConvPFTausTightIsoTrackPt5_*_*', 
        'keep *_hltConvPFTausTightIso_*_*', 
        'keep *_hltConvPFTausTrackFindingLooseIsolation_*_*', 
        'keep *_hltConvPFTausTrackFinding_*_*', 
        'keep *_hltConvPFTaus_*_*', 
        'keep *_hltCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCorrectedIslandEndcapSuperClustersHI_*_*', 
        'keep *_hltCsc2DRecHits_*_*', 
        'keep *_hltCscSegments_*_*', 
        'keep *_hltDisplacedhltIter4PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolationProng4L1HLTMatched_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolationProng4_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolation_*_*', 
        'keep *_hltDoublePFTau25TrackPt5_*_*', 
        'keep *_hltDoublePFTau25_*_*', 
        'keep *_hltDoublePFTauTightIso45Track5_*_*', 
        'keep *_hltDoublePFTauTightIso45Track_*_*', 
        'keep *_hltDt4DSegments_*_*', 
        'keep *_hltEcalPhiSymFilter_*_*', 
        'keep *_hltEcalRecHitAll_*_*', 
        'keep *_hltEcalRecHit_*_*', 
        'keep *_hltEgammaCandidates_*_*', 
        'keep *_hltEgammaGsfElectrons_*_*', 
        'keep *_hltEgammaGsfTracks_*_*', 
        'keep *_hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1SingleEG18orL1SingleEG20_*_*', 
        'keep *_hltEle20CaloIdVTTrkIdTDphiFilter_*_*', 
        'keep *_hltEle27WP85PixelMatchElectronsL1Seeded_*_*', 
        'keep *_hltElectronsVertex_*_*', 
        'keep *_hltFEDSelectorLumiPixels_*_*', 
        'keep *_hltFastPVPixelTracksMerger_*_*', 
        'keep *_hltFastPVPixelTracksRecover_*_*', 
        'keep *_hltFastPVPixelTracks_*_*', 
        'keep *_hltFastPVPixelVertices3D_*_*', 
        'keep *_hltFastPVPixelVertices_*_*', 
        'keep *_hltFastPixelBLifetimeL3TagInfos_*_*', 
        'keep *_hltFastPrimaryVertex_*_*', 
        'keep *_hltFilterDoubleIsoPFTau45Trk5LeadTrack5IsolationL1HLTMatched_*_*', 
        'keep *_hltFilterL2EtCutDoublePFIsoTau45Trk5_*_*', 
        'keep *_hltFilterL2EtCutSingleIsoPFTau35Trk20MET70_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20LeadTrackPt20_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20MET60LeadTrack20IsolationL1HLTMatched_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20MET70LeadTrack20IsolationL1HLTMatched_*_*', 
        'keep *_hltGmtStage2Digis_*_*', 
        'keep *_hltGoodOnlinePVs_*_*', 
        'keep *_hltGtStage2Digis_*_*', 
        'keep *_hltHICaloJetCorrected_*_*', 
        'keep *_hltHICaloJetIDPassed_*_*', 
        'keep *_hltHIGoodLooseTracks_*_*', 
        'keep *_hltHIPixel3PrimTracks_*_*', 
        'keep *_hltHISelectedVertex_*_*', 
        'keep *_hltHISiPixelClusters_*_*', 
        'keep *_hltHITIPTCorrectorHB_*_*', 
        'keep *_hltHITIPTCorrectorHE_*_*', 
        'keep *_hltHbhereco_*_*', 
        'keep *_hltHfreco_*_*', 
        'keep *_hltHiCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltHiCorrectedIslandEndcapSuperClustersHI_*_*', 
        'keep *_hltHiIslandSuperClustersHI_*_*', 
        'keep *_hltHoreco_*_*', 
        'keep *_hltIsolPixelTrackProdHB_*_*', 
        'keep *_hltIsolPixelTrackProdHE_*_*', 
        'keep *_hltIter0ElectronsTrackSelectionHighPurity_*_*', 
        'keep *_hltIter0HighPtTkMuPixelTracks_*_*', 
        'keep *_hltIter0HighPtTkMuTrackSelectionHighPurity_*_*', 
        'keep *_hltIter0PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter1Merged_*_*', 
        'keep *_hltIter1PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter2HighPtTkMuMerged_*_*', 
        'keep *_hltIter2HighPtTkMuTrackSelectionHighPurity_*_*', 
        'keep *_hltIter2MergedForBTag_*_*', 
        'keep *_hltIter2MergedForElectrons_*_*', 
        'keep *_hltIter2Merged_*_*', 
        'keep *_hltIter2PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter3Merged_*_*', 
        'keep *_hltIter4Merged_*_*', 
        'keep *_hltIterativeCone5PileupSubtractionCaloJets_*_*', 
        'keep *_hltL1HLTSingleIsoPFTau35Trk20Met60JetsMatch_*_*', 
        'keep *_hltL1IsoElectronTrackIsol_*_*', 
        'keep *_hltL1NonIsoElectronTrackIsol_*_*', 
        'keep *_hltL1SeededRecoEcalCandidate_*_*', 
        'keep *_hltL1sDoubleTauJet44Eta2p17orDoubleJet64Central_*_*', 
        'keep *_hltL1sDoubleTauJet44erorDoubleJetC64_*_*', 
        'keep *_hltL1sL1EG18er_*_*', 
        'keep *_hltL1sL1ETM36ORETM40_*_*', 
        'keep *_hltL1sL1Jet52ETM30_*_*', 
        'keep *_hltL1sL1SingleEG12_*_*', 
        'keep *_hltL1sL1SingleEG15_*_*', 
        'keep *_hltL1sL1SingleEG18orL1SingleEG20_*_*', 
        'keep *_hltL1sL1SingleMu10_*_*', 
        'keep *_hltL1sL1SingleMu14Eta2p1_*_*', 
        'keep *_hltL1sMu16Eta2p1_*_*', 
        'keep *_hltL2MuonCandidatesNoVtx_*_*', 
        'keep *_hltL2MuonCandidates_*_*', 
        'keep *_hltL2MuonSeeds_*_*', 
        'keep *_hltL2Muons_*_*', 
        'keep *_hltL2TauJets_*_*', 
        'keep *_hltL3MuonCandidates_*_*', 
        'keep *_hltL3MuonsIOHit_*_*', 
        'keep *_hltL3MuonsLinksCombination_*_*', 
        'keep *_hltL3MuonsOIHit_*_*', 
        'keep *_hltL3MuonsOIState_*_*', 
        'keep *_hltL3Muons_*_*', 
        'keep *_hltL3NoFiltersNoVtxMuonCandidates_*_*', 
        'keep *_hltL3NoFiltersNoVtxMuons_*_*', 
        'keep *_hltL3SecondaryVertexTagInfos_*_*', 
        'keep *_hltL3TkFromL2OICombination_*_*', 
        'keep *_hltL3TkTracksFromL2IOHit_*_*', 
        'keep *_hltL3TkTracksFromL2OIHit_*_*', 
        'keep *_hltL3TkTracksFromL2OIState_*_*', 
        'keep *_hltL3TkTracksFromL2_*_*', 
        'keep *_hltL3TrackCandidateFromL2IOHit_*_*', 
        'keep *_hltL3TrackCandidateFromL2OIHit_*_*', 
        'keep *_hltL3TrackCandidateFromL2OIState_*_*', 
        'keep *_hltL3TrajSeedIOHit_*_*', 
        'keep *_hltL3TrajSeedOIHit_*_*', 
        'keep *_hltL3TrajSeedOIState_*_*', 
        'keep *_hltL3TrajectorySeed_*_*', 
        'keep *_hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f18QL3crIsoRhoFiltered0p15_*_*', 
        'keep *_hltLeadingCentralJets30SingleTopEle27_*_*', 
        'keep *_hltLeadingCentralJets30SingleTopIsoMu24_*_*', 
        'keep *_hltMergedTracks_*_*', 
        'keep *_hltMet_*_*', 
        'keep *_hltMu8Ele17CaloIdTCaloIsoVLPixelMatchFilter_*_*', 
        'keep *_hltMuTrackJpsiCtfTrackCands_*_*', 
        'keep *_hltMuTrackJpsiPixelTrackCands_*_*', 
        'keep *_hltMuonCSCDigis_*_*', 
        'keep *_hltMuonCSCDigis_MuonCSCStripDigi_*', 
        'keep *_hltMuonCSCDigis_MuonCSCWireDigi_*', 
        'keep *_hltMuonDTDigis_*_*', 
        'keep *_hltMuonRPCDigis_*_*', 
        'keep *_hltOnlineBeamSpot_*_*', 
        'keep *_hltOnlinePrimaryVertices_*_*', 
        'keep *_hltOverlapFilterEle20LooseIsoPFTau20OldVersion_*_*', 
        'keep *_hltOverlapFilterIsoEle20MediumIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15IsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15MediumIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15TightIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu18LooseIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu18PFTau25TrackPt5Prong4_*_*', 
        'keep *_hltPFJetForBtag_*_*', 
        'keep *_hltPFMuonMerging_*_*', 
        'keep *_hltPFTau15TrackLooseIso_*_*', 
        'keep *_hltPFTau15Track_*_*', 
        'keep *_hltPFTau15_*_*', 
        'keep *_hltPFTau20IsoMuVertex_*_*', 
        'keep *_hltPFTau20TrackLooseIso_*_*', 
        'keep *_hltPFTau20Track_*_*', 
        'keep *_hltPFTau20_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolationProng4IsoMuVertex_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolationProng4_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolation_*_*', 
        'keep *_hltPFTau25TrackPt5_*_*', 
        'keep *_hltPFTau25_*_*', 
        'keep *_hltPFTau35TrackPt20LooseIsoProng2_*_*', 
        'keep *_hltPFTau35TrackPt20LooseIso_*_*', 
        'keep *_hltPFTau35TrackPt20_*_*', 
        'keep *_hltPFTau35Track_*_*', 
        'keep *_hltPFTau35_*_*', 
        'keep *_hltPFTauEleVertex20_*_*', 
        'keep *_hltPFTauJetTracksAssociator_*_*', 
        'keep *_hltPFTauMediumIso20TrackMediumIso_*_*', 
        'keep *_hltPFTauMediumIso20Track_*_*', 
        'keep *_hltPFTauMediumIso20_*_*', 
        'keep *_hltPFTauMediumIso35Track_*_*', 
        'keep *_hltPFTauMediumIso35_*_*', 
        'keep *_hltPFTauTagInfo_*_*', 
        'keep *_hltPFTauTightIso20TrackTightIso_*_*', 
        'keep *_hltPFTauTightIso20Track_*_*', 
        'keep *_hltPFTauTightIso20_*_*', 
        'keep *_hltPFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltParticleFlowForTaus_*_*', 
        'keep *_hltParticleFlow_*_*', 
        'keep *_hltPixelMatch3HitElectronsActivity_*_*', 
        'keep *_hltPixelMatch3HitElectronsL1Seeded_*_*', 
        'keep *_hltPixelMatchCleanElectronsL1Seeded_*_*', 
        'keep *_hltPixelMatchElectronsActivity_*_*', 
        'keep *_hltPixelMatchElectronsL1Iso_*_*', 
        'keep *_hltPixelMatchElectronsL1NonIso_*_*', 
        'keep *_hltPixelMatchElectronsL1Seeded_*_*', 
        'keep *_hltPixelTracksElectrons_*_*', 
        'keep *_hltPixelTracksMerged_*_*', 
        'keep *_hltPixelTracks_*_*', 
        'keep *_hltPixelVertices3DbbPhi_*_*', 
        'keep *_hltPixelVertices_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidateSC4_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidateSC5_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidate_*_*', 
        'keep *_hltRpcRecHits_*_*', 
        'keep *_hltSelector4CentralJetsL1FastJet_*_*', 
        'keep *_hltSelector8CentralJetsL1FastJet_*_*', 
        'keep *_hltSelectorJets20L1FastJet_*_*', 
        'keep *_hltSiPixelCluster_*_*', 
        'keep *_hltSiPixelClusters_*_*', 
        'keep *_hltSiStripClusters_*_*', 
        'keep *_hltSiStripRawToClustersFacility_*_*', 
        'keep *_hltSingleMu15L3Filtered15_*_*', 
        'keep *_hltSingleMuIsoL1s14L3IsoFiltered15eta2p1_*_*', 
        'keep *_hltSingleMuIsoL3IsoFiltered15_*_*', 
        'keep *_hltTowerMakerForAll_*_*', 
        'keep *_hltTowerMakerForMuons_*_*', 
        'keep *_hltTriggerSummaryAOD_*_*', 
        'keep *_hltTriggerSummaryRAW_*_*', 
        'keep *_hltTrimmedPixelVertices_*_*', 
        'keep *_hltVerticesL3_*_*', 
        'keep *_hltVerticesPFFilter_*_*', 
        'keep *_hltVerticesPFSelector_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_rawDataRepacker_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_virginRawDataRepacker_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep L2MuonTrajectorySeeds_hltL2MuonSeeds_*_*', 
        'keep L3MuonTrajectorySeeds_hltHIL3TrajSeedOIHit_*_*', 
        'keep L3MuonTrajectorySeeds_hltHIL3TrajectorySeed_*_*', 
        'keep L3MuonTrajectorySeeds_hltL3TrajSeedOIState_*_*', 
        'keep LumiScalerss_hltScalersRawToDigi_*_*', 
        'keep SiPixelClusteredmNewDetSetVector_hltSiPixelClusters_*_*', 
        'keep TrackCandidates_hltHIL3TrackCandidateFromL2OIHit_*_*', 
        'keep TrackCandidates_hltHIL3TrackCandidateFromL2OIState_*_*', 
        'keep TrackingRecHitsOwned_hltL3Muons_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep recoCaloJets_*_*_*', 
        'keep recoCaloMETs_*_*_*', 
        'keep recoCaloMETs_hltMet_*_*', 
        'keep recoCompositeCandidates_*_*_*', 
        'keep recoElectrons_*_*_*', 
        'keep recoIsolatedPixelTrackCandidates_*_*_*', 
        'keep recoMETs_*_*_*', 
        'keep recoPFJets_*_*_*', 
        'keep recoPFTaus_*_*_*', 
        'keep recoRecoChargedCandidates_*_*_*', 
        'keep recoRecoChargedCandidates_hltHIL3MuonCandidates_*_*', 
        'keep recoRecoChargedCandidates_hltL2MuonCandidates_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsoHLTClusterShape_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonEcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalForHE_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsoHLTClusterShape_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonEcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalForHE_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalIsol_*_*', 
        'keep recoRecoEcalCandidates_*_*_*', 
        'keep recoRecoEcalCandidates_hltL1IsoRecoEcalCandidate_*_*', 
        'keep recoRecoEcalCandidates_hltL1NonIsoRecoEcalCandidate_*_*', 
        'keep recoTrackExtras_hltHIL3MuonsOIHit_*_*', 
        'keep recoTrackExtras_hltHIL3MuonsOIState_*_*', 
        'keep recoTracks_hltHIL3MuonsOIHit_*_*', 
        'keep recoTracks_hltHIL3MuonsOIState_*_*', 
        'keep recoTracks_hltHIL3Muons_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2OIHit_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2OIState_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2_*_*', 
        'keep triggerTriggerEventWithRefs_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep triggerTriggerFilterObjectWithRefs_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep *_*_MergedTrackTruth_*', 
        'keep *_*_StripDigiSimLink_*', 
        'keep *_*_PixelDigiSimLink_*', 
        'keep *_*_MuonCSCStripDigiSimLinks_*', 
        'keep *_*_MuonCSCWireDigiSimLinks_*', 
        'keep *_*_RPCDigiSimLink_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_*_*_*' ) ),
    splitLevel = cms.untracked.int32(0)
)

process.FEVTEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring( ('drop *', 
        'drop *', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep DetIdedmEDCollection_siStripDigis_*_*', 
        'keep DetIdedmEDCollection_siPixelDigis_*_*', 
        'keep PixelFEDChanneledmNewDetSetVector_siPixelDigis_*_*', 
        'keep *_siPixelClusters_*_*', 
        'keep *_siStripClusters_*_*', 
        'keep ClusterSummary_clusterSummaryProducer_*_*', 
        'keep *_dt1DRecHits_*_*', 
        'keep *_dt4DSegments_*_*', 
        'keep *_dt1DCosmicRecHits_*_*', 
        'keep *_dt4DCosmicSegments_*_*', 
        'keep *_csc2DRecHits_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*', 
        'keep *_hbhereco_*_*', 
        'keep *_hbheprereco_*_*', 
        'keep *_hfprereco_*_*', 
        'keep *_hfreco_*_*', 
        'keep *_horeco_*_*', 
        'keep HBHERecHitsSorted_hbherecoMB_*_*', 
        'keep HORecHitsSorted_horecoMB_*_*', 
        'keep HFRecHitsSorted_hfrecoMB_*_*', 
        'keep ZDCDataFramesSorted_hcalDigis_*_*', 
        'keep ZDCDataFramesSorted_castorDigis_*_*', 
        'keep QIE10DataFrameHcalDataFrameContainer_hcalDigis_ZDC_*', 
        'keep ZDCRecHitsSorted_zdcreco_*_*', 
        'keep *_reducedHcalRecHits_*_*', 
        'keep *_castorreco_*_*', 
        'keep HcalUnpackerReport_castorDigis_*_*', 
        'keep HcalUnpackerReport_hcalDigiAlCaMB_*_*', 
        'keep HcalUnpackerReport_hcalDigis_*_*', 
        'keep *_ecalPreshowerRecHit_*_*', 
        'keep *_ecalRecHit_*_*', 
        'keep *_ecalCompactTrigPrim_*_*', 
        'keep *_ecalTPSkim_*_*', 
        'keep EBSrFlagsSorted_ecalDigis__*', 
        'keep EESrFlagsSorted_ecalDigis__*', 
        'keep *_selectDigi_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEE_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEB_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsES_*_*', 
        'keep *_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep *_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_*', 
        'keep *_particleFlowSuperClusterECAL_*_*', 
        'keep *_particleFlowSuperClusterOOTECAL_*_*', 
        'drop recoClusterShapes_*_*_*', 
        'drop recoBasicClustersToOnerecoClusterShapesAssociation_*_*_*', 
        'drop recoBasicClusters_multi5x5BasicClusters_multi5x5BarrelBasicClusters_*', 
        'drop recoSuperClusters_multi5x5SuperClusters_multi5x5BarrelSuperClusters_*', 
        'keep *_CkfElectronCandidates_*_*', 
        'keep *_GsfGlobalElectronTest_*_*', 
        'keep *_electronMergedSeeds_*_*', 
        'keep recoGsfTracks_electronGsfTracks_*_*', 
        'keep recoGsfTrackExtras_electronGsfTracks_*_*', 
        'keep recoTrackExtras_electronGsfTracks_*_*', 
        'keep TrackingRecHitsOwned_electronGsfTracks_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTrackExtras_generalTracks_*_*', 
        'keep TrackingRecHitsOwned_generalTracks_*_*', 
        'keep *_generalTracks_MVAValues_*', 
        'keep *_generalTracks_MVAVals_*', 
        'keep TrackingRecHitsOwned_extraFromSeeds_*_*', 
        'keep uints_extraFromSeeds_*_*', 
        'keep recoTracks_beamhaloTracks_*_*', 
        'keep recoTrackExtras_beamhaloTracks_*_*', 
        'keep TrackingRecHitsOwned_beamhaloTracks_*_*', 
        'keep recoTracks_conversionStepTracks_*_*', 
        'keep recoTrackExtras_conversionStepTracks_*_*', 
        'keep TrackingRecHitsOwned_conversionStepTracks_*_*', 
        'keep *_ctfPixelLess_*_*', 
        'keep *_dedxTruncated40_*_*', 
        'keep *_dedxHitInfo_*_*', 
        'keep *_dedxHarmonic2_*_*', 
        'keep *_dedxPixelHarmonic2_*_*', 
        'keep *_trackExtrapolator_*_*', 
        'keep recoTracks_cosmicDCTracks_*_*', 
        'keep recoTrackExtras_cosmicDCTracks_*_*', 
        'keep TrackingRecHitsOwned_cosmicDCTracks_*_*', 
        'keep *_ak4CaloJets_*_*', 
        'keep *_ak4PFJets_*_*', 
        'keep *_ak4PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHSSoftDrop_*_*', 
        'keep *_cmsTopTagPFJetsCHS_*_*', 
        'keep *_JetPlusTrackZSPCorJetAntiKt4_*_*', 
        'keep *_ak4TrackJets_*_*', 
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', 
        'keep *_caloTowers_*_*', 
        'keep *_towerMaker_*_*', 
        'keep *_CastorTowerReco_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertex_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertexPF_*_*', 
        'keep *_ak4JetTracksAssociatorAtCaloFace_*_*', 
        'keep *_ak4JetTracksAssociatorExplicit_*_*', 
        'keep *_ak4JetExtender_*_*', 
        'keep *_ak4JetID_*_*', 
        'keep *_ak5CastorJets_*_*', 
        'keep *_ak5CastorJetID_*_*', 
        'keep *_ak7CastorJets_*_*', 
        'keep *_ak7CastorJetID_*_*', 
        'keep *_fixedGridRhoAll_*_*', 
        'keep *_fixedGridRhoFastjetAll_*_*', 
        'keep *_fixedGridRhoFastjetAllTmp_*_*', 
        'keep *_fixedGridRhoFastjetAllCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentral_*_*', 
        'keep *_fixedGridRhoFastjetCentralCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralChargedPileUp_*_*', 
        'keep *_fixedGridRhoFastjetCentralNeutral_*_*', 
        'keep *_ak8PFJetsCHSSoftDropMass_*_*', 
        'keep recoCaloMETs_caloMet_*_*', 
        'keep recoCaloMETs_caloMetBE_*_*', 
        'keep recoCaloMETs_caloMetBEFO_*_*', 
        'keep recoCaloMETs_caloMetM_*_*', 
        'keep recoPFMETs_pfMet_*_*', 
        'keep recoPFMETs_pfChMet_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_*_*', 
        'keep recoHcalNoiseRBXs_hcalnoise_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoEcalHaloData_EcalHaloData_*_*', 
        'keep recoGlobalHaloData_GlobalHaloData_*_*', 
        'keep recoHcalHaloData_HcalHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep *_MuonSeed_*_*', 
        'keep *_ancientMuonSeed_*_*', 
        'keep *_displacedMuonSeeds_*_*', 
        'keep TrackingRecHitsOwned_globalMuons_*_*', 
        'keep TrackingRecHitsOwned_tevMuons_*_*', 
        'keep *_CosmicMuonSeed_*_*', 
        'keep recoTrackExtras_cosmicMuons_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons_*_*', 
        'keep recoTrackExtras_cosmicMuons1Leg_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons1Leg_*_*', 
        'keep recoTracks_cosmicsVetoTracks_*_*', 
        'keep *_muons_*_*', 
        'keep *_particleFlow_muons_*', 
        'drop *_muons_muons1stStep2muonsMap_*', 
        'drop recoIsoDepositedmValueMap_muons_*_*', 
        'drop doubleedmValueMap_muons_muPFIso*_*', 
        'keep recoTracks_standAloneMuons_*_*', 
        'keep recoTrackExtras_standAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_standAloneMuons_*_*', 
        'keep recoTracks_globalMuons_*_*', 
        'keep recoTrackExtras_globalMuons_*_*', 
        'keep recoTracks_tevMuons_*_*', 
        'keep recoTrackExtras_tevMuons_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_displacedTracks_*_*', 
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_*_*', 
        'keep recoTracks_displacedGlobalMuons_*_*', 
        'keep recoTrackExtras_displacedGlobalMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedGlobalMuons_*_*', 
        'keep recoTracks_cosmicMuons_*_*', 
        'keep recoMuons_muonsFromCosmics_*_*', 
        'keep recoTracks_cosmicMuons1Leg_*_*', 
        'keep recoMuons_muonsFromCosmics1Leg_*_*', 
        'keep recoTracks_refittedStandAloneMuons_*_*', 
        'keep recoTrackExtras_refittedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_refittedStandAloneMuons_*_*', 
        'keep recoTracks_displacedStandAloneMuons__*', 
        'keep recoTrackExtras_displacedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedStandAloneMuons_*_*', 
        'keep *_muIsoDepositTk_*_*', 
        'keep *_muIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muIsoDepositJets_*_*', 
        'keep *_muGlobalIsoDepositCtfTk_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muGlobalIsoDepositJets_*_*', 
        'keep *_softPFMuonsTagInfos_*_*', 
        'keep *_softPFElectronsTagInfos_*_*', 
        'keep *_softPFElectronBJetTags_*_*', 
        'keep *_softPFMuonBJetTags_*_*', 
        'keep *_pfImpactParameterTagInfos_*_*', 
        'keep *_pfTrackCountingHighEffBJetTags_*_*', 
        'keep *_pfJetProbabilityBJetTags_*_*', 
        'keep *_pfJetBProbabilityBJetTags_*_*', 
        'keep *_pfSecondaryVertexTagInfos_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderTagInfos_*_*', 
        'keep *_pfSimpleSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfCombinedSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfGhostTrackVertexTagInfos_*_*', 
        'keep *_pfGhostTrackBJetTags_*_*', 
        'keep *_pfCombinedMVAV2BJetTags_*_*', 
        'keep *_inclusiveCandidateSecondaryVertices_*_*', 
        'keep *_inclusiveCandidateSecondaryVerticesCvsL_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderCvsLTagInfos_*_*', 
        'keep *_pfCombinedCvsLJetTags_*_*', 
        'keep *_pfCombinedCvsBJetTags_*_*', 
        'keep *_pfChargeBJetTags_*_*', 
        'keep *_pfDeepCSVJetTags_*_*', 
        'keep *_pfDeepCMVAJetTags_*_*', 
        'keep recoRecoTauPiZeros_hpsPFTauProducer_pizeros_*', 
        'keep recoPFTaus_hpsPFTauProducer_*_*', 
        'keep *_hpsPFTauDiscriminationByDeadECALElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFinding_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingNewDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingOldDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseChargedIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByTightMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSum_*_*', 
        'keep *_hpsPFTauPUcorrPtSum_*_*', 
        'keep *_hpsPFTauChargedIsoPtSum_*_*', 
        'keep *_hpsPFTauTransverseImpactParameters_*_*', 
        'keep *_hpsPFTauFootprintCorrection_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeight_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6rawElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6LooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6MediumElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6TightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VTightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauChargedIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeightdR03_*_*', 
        'keep *_hpsPFTauFootprintCorrectiondR03_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalConedR03_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep  *_offlinePrimaryVertices__*', 
        'keep *_offlinePrimaryVerticesWithBS_*_*', 
        'keep *_offlinePrimaryVerticesFromCosmicTracks_*_*', 
        'keep *_nuclearInteractionMaker_*_*', 
        'keep *_generalV0Candidates_*_*', 
        'keep *_inclusiveSecondaryVertices_*_*', 
        'keep recoGsfElectronCores_gsfElectronCores_*_*', 
        'keep recoGsfElectronCores_gedGsfElectronCores_*_*', 
        'keep recoGsfElectrons_gsfElectrons_*_*', 
        'keep recoGsfElectrons_gedGsfElectrons_*_*', 
        'keep recoGsfElectronCores_uncleanedOnlyGsfElectronCores_*_*', 
        'keep recoGsfElectrons_uncleanedOnlyGsfElectrons_*_*', 
        'keep floatedmValueMap_eidRobustLoose_*_*', 
        'keep floatedmValueMap_eidRobustTight_*_*', 
        'keep floatedmValueMap_eidRobustHighEnergy_*_*', 
        'keep floatedmValueMap_eidLoose_*_*', 
        'keep floatedmValueMap_eidTight_*_*', 
        'keep *_egmGedGsfElectronPFIsolation_*_*', 
        'keep *_photonEcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonEcalPFClusterIsolationProducer_*_*', 
        'keep *_electronEcalPFClusterIsolationProducer_*_*', 
        'keep *_photonHcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonHcalPFClusterIsolationProducer_*_*', 
        'keep *_electronHcalPFClusterIsolationProducer_*_*', 
        'drop *_egmGsfElectronIDs_*_*', 
        'drop *_egmPhotonIDs_*_*', 
        'keep *_gedPhotonCore_*_*', 
        'keep *_gedPhotons_*_*', 
        'keep *_particleBasedIsolation_*_*', 
        'keep recoPhotons_mustachePhotons_*_*', 
        'keep recoPhotonCores_mustachePhotonCore_*_*', 
        'keep recoPhotons_ootPhotons_*_*', 
        'keep recoPhotonCores_ootPhotonCore_*_*', 
        'keep recoPhotons_photons_*_*', 
        'keep recoPhotonCores_photonCore_*_*', 
        'keep recoConversions_conversions_*_*', 
        'keep recoConversions_mustacheConversions_*_*', 
        'drop *_conversions_uncleanedConversions_*', 
        'drop *_gedPhotonsTmp_valMapPFEgammaCandToPhoton_*', 
        'keep recoConversions_allConversions_*_*', 
        'keep recoConversions_allConversionsOldEG_*_*', 
        'keep recoTracks_ckfOutInTracksFromConversions_*_*', 
        'keep recoTracks_ckfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfInOutTracksFromConversions_*_*', 
        'keep recoConversions_uncleanedOnlyAllConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep *_PhotonIDProdGED_*_*', 
        'keep recoRecoEcalCandidates_hfRecoEcalCandidate_*_*', 
        'keep *_hfEMClusters_*_*', 
        'keep *_gedGsfElectronCores_*_*', 
        'keep *_gedGsfElectrons_*_*', 
        'keep *_pixelTracks_*_*', 
        'keep *_pixelVertices_*_*', 
        'drop CaloTowersSorted_towerMakerPF_*_*', 
        'keep recoPFRecHits_particleFlowClusterECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHCAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterPS_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHBHE_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitPS_Cleaned_*', 
        'keep recoPFClusters_particleFlowClusterECAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHCAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHO_*_*', 
        'keep recoPFClusters_particleFlowClusterHF_*_*', 
        'keep recoPFClusters_particleFlowClusterPS_*_*', 
        'keep recoPFBlocks_particleFlowBlock_*_*', 
        'keep recoPFCandidates_particleFlowEGamma_*_*', 
        'keep recoCaloClusters_particleFlowEGamma_*_*', 
        'keep recoSuperClusters_particleFlowEGamma_*_*', 
        'keep recoConversions_particleFlowEGamma_*_*', 
        'keep recoPFCandidates_particleFlow_*_*', 
        'keep recoPFCandidates_particleFlowTmp_electrons_*', 
        'keep recoPFCandidates_particleFlowTmp_*_*', 
        'drop recoPFCandidates_particleFlowTmp__*', 
        'keep recoPFDisplacedVertexs_particleFlowDisplacedVertex_*_*', 
        'keep *_pfElectronTranslator_*_*', 
        'keep *_pfPhotonTranslator_*_*', 
        'keep *_particleFlow_electrons_*', 
        'keep *_particleFlow_photons_*', 
        'keep *_particleFlow_muons_*', 
        'keep *_trackerDrivenElectronSeeds_preid_*', 
        'keep *_particleFlowPtrs_*_*', 
        'keep *_particleFlowTmpPtrs_*_*', 
        'keep *_chargedHadronPFTrackIsolation_*_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep L1MuGMTReadoutCollection_gtDigis_*_*', 
        'keep L1GctEmCand*_gctDigis_*_*', 
        'keep L1GctJetCand*_gctDigis_*_*', 
        'keep L1GctEtHad*_gctDigis_*_*', 
        'keep L1GctEtMiss*_gctDigis_*_*', 
        'keep L1GctEtTotal*_gctDigis_*_*', 
        'keep L1GctHtMiss*_gctDigis_*_*', 
        'keep L1GctJetCounts*_gctDigis_*_*', 
        'keep L1GctHFRingEtSums*_gctDigis_*_*', 
        'keep L1GctHFBitCounts*_gctDigis_*_*', 
        'keep LumiDetails_lumiProducer_*_*', 
        'keep LumiSummary_lumiProducer_*_*', 
        'drop *_hlt*_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep L1AcceptBunchCrossings_scalersRawToDigi_*_*', 
        'keep L1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep Level1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep BeamSpotOnlines_scalersRawToDigi_*_*', 
        'keep DcsStatuss_scalersRawToDigi_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep *_tcdsDigis_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'keep *_pfIsolatedElectronsEI_*_*', 
        'keep *_pfIsolatedMuonsEI_*_*', 
        'keep recoPFJets_pfJetsEI_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTagsEI_*_*', 
        'keep recoPFTaus_pfTausEI_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByDecayModeFinding_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByIsolation_*_*', 
        'keep recoPFMETs_pfMetEI_*_*' ) ),
    splitLevel = cms.untracked.int32(0)
)

process.FEVTHLTALLEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring( ('drop *', 
        'drop *', 
        'drop *', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep DetIdedmEDCollection_siStripDigis_*_*', 
        'keep DetIdedmEDCollection_siPixelDigis_*_*', 
        'keep PixelFEDChanneledmNewDetSetVector_siPixelDigis_*_*', 
        'keep *_siPixelClusters_*_*', 
        'keep *_siStripClusters_*_*', 
        'keep ClusterSummary_clusterSummaryProducer_*_*', 
        'keep *_dt1DRecHits_*_*', 
        'keep *_dt4DSegments_*_*', 
        'keep *_dt1DCosmicRecHits_*_*', 
        'keep *_dt4DCosmicSegments_*_*', 
        'keep *_csc2DRecHits_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*', 
        'keep *_hbhereco_*_*', 
        'keep *_hbheprereco_*_*', 
        'keep *_hfprereco_*_*', 
        'keep *_hfreco_*_*', 
        'keep *_horeco_*_*', 
        'keep HBHERecHitsSorted_hbherecoMB_*_*', 
        'keep HORecHitsSorted_horecoMB_*_*', 
        'keep HFRecHitsSorted_hfrecoMB_*_*', 
        'keep ZDCDataFramesSorted_hcalDigis_*_*', 
        'keep ZDCDataFramesSorted_castorDigis_*_*', 
        'keep QIE10DataFrameHcalDataFrameContainer_hcalDigis_ZDC_*', 
        'keep ZDCRecHitsSorted_zdcreco_*_*', 
        'keep *_reducedHcalRecHits_*_*', 
        'keep *_castorreco_*_*', 
        'keep HcalUnpackerReport_castorDigis_*_*', 
        'keep HcalUnpackerReport_hcalDigiAlCaMB_*_*', 
        'keep HcalUnpackerReport_hcalDigis_*_*', 
        'keep *_ecalPreshowerRecHit_*_*', 
        'keep *_ecalRecHit_*_*', 
        'keep *_ecalCompactTrigPrim_*_*', 
        'keep *_ecalTPSkim_*_*', 
        'keep EBSrFlagsSorted_ecalDigis__*', 
        'keep EESrFlagsSorted_ecalDigis__*', 
        'keep *_selectDigi_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEE_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEB_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsES_*_*', 
        'keep *_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep *_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_*', 
        'keep *_particleFlowSuperClusterECAL_*_*', 
        'keep *_particleFlowSuperClusterOOTECAL_*_*', 
        'drop recoClusterShapes_*_*_*', 
        'drop recoBasicClustersToOnerecoClusterShapesAssociation_*_*_*', 
        'drop recoBasicClusters_multi5x5BasicClusters_multi5x5BarrelBasicClusters_*', 
        'drop recoSuperClusters_multi5x5SuperClusters_multi5x5BarrelSuperClusters_*', 
        'keep *_CkfElectronCandidates_*_*', 
        'keep *_GsfGlobalElectronTest_*_*', 
        'keep *_electronMergedSeeds_*_*', 
        'keep recoGsfTracks_electronGsfTracks_*_*', 
        'keep recoGsfTrackExtras_electronGsfTracks_*_*', 
        'keep recoTrackExtras_electronGsfTracks_*_*', 
        'keep TrackingRecHitsOwned_electronGsfTracks_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTrackExtras_generalTracks_*_*', 
        'keep TrackingRecHitsOwned_generalTracks_*_*', 
        'keep *_generalTracks_MVAValues_*', 
        'keep *_generalTracks_MVAVals_*', 
        'keep TrackingRecHitsOwned_extraFromSeeds_*_*', 
        'keep uints_extraFromSeeds_*_*', 
        'keep recoTracks_beamhaloTracks_*_*', 
        'keep recoTrackExtras_beamhaloTracks_*_*', 
        'keep TrackingRecHitsOwned_beamhaloTracks_*_*', 
        'keep recoTracks_conversionStepTracks_*_*', 
        'keep recoTrackExtras_conversionStepTracks_*_*', 
        'keep TrackingRecHitsOwned_conversionStepTracks_*_*', 
        'keep *_ctfPixelLess_*_*', 
        'keep *_dedxTruncated40_*_*', 
        'keep *_dedxHitInfo_*_*', 
        'keep *_dedxHarmonic2_*_*', 
        'keep *_dedxPixelHarmonic2_*_*', 
        'keep *_trackExtrapolator_*_*', 
        'keep recoTracks_cosmicDCTracks_*_*', 
        'keep recoTrackExtras_cosmicDCTracks_*_*', 
        'keep TrackingRecHitsOwned_cosmicDCTracks_*_*', 
        'keep *_ak4CaloJets_*_*', 
        'keep *_ak4PFJets_*_*', 
        'keep *_ak4PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHSSoftDrop_*_*', 
        'keep *_cmsTopTagPFJetsCHS_*_*', 
        'keep *_JetPlusTrackZSPCorJetAntiKt4_*_*', 
        'keep *_ak4TrackJets_*_*', 
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', 
        'keep *_caloTowers_*_*', 
        'keep *_towerMaker_*_*', 
        'keep *_CastorTowerReco_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertex_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertexPF_*_*', 
        'keep *_ak4JetTracksAssociatorAtCaloFace_*_*', 
        'keep *_ak4JetTracksAssociatorExplicit_*_*', 
        'keep *_ak4JetExtender_*_*', 
        'keep *_ak4JetID_*_*', 
        'keep *_ak5CastorJets_*_*', 
        'keep *_ak5CastorJetID_*_*', 
        'keep *_ak7CastorJets_*_*', 
        'keep *_ak7CastorJetID_*_*', 
        'keep *_fixedGridRhoAll_*_*', 
        'keep *_fixedGridRhoFastjetAll_*_*', 
        'keep *_fixedGridRhoFastjetAllTmp_*_*', 
        'keep *_fixedGridRhoFastjetAllCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentral_*_*', 
        'keep *_fixedGridRhoFastjetCentralCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralChargedPileUp_*_*', 
        'keep *_fixedGridRhoFastjetCentralNeutral_*_*', 
        'keep *_ak8PFJetsCHSSoftDropMass_*_*', 
        'keep recoCaloMETs_caloMet_*_*', 
        'keep recoCaloMETs_caloMetBE_*_*', 
        'keep recoCaloMETs_caloMetBEFO_*_*', 
        'keep recoCaloMETs_caloMetM_*_*', 
        'keep recoPFMETs_pfMet_*_*', 
        'keep recoPFMETs_pfChMet_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_*_*', 
        'keep recoHcalNoiseRBXs_hcalnoise_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoEcalHaloData_EcalHaloData_*_*', 
        'keep recoGlobalHaloData_GlobalHaloData_*_*', 
        'keep recoHcalHaloData_HcalHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep *_MuonSeed_*_*', 
        'keep *_ancientMuonSeed_*_*', 
        'keep *_displacedMuonSeeds_*_*', 
        'keep TrackingRecHitsOwned_globalMuons_*_*', 
        'keep TrackingRecHitsOwned_tevMuons_*_*', 
        'keep *_CosmicMuonSeed_*_*', 
        'keep recoTrackExtras_cosmicMuons_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons_*_*', 
        'keep recoTrackExtras_cosmicMuons1Leg_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons1Leg_*_*', 
        'keep recoTracks_cosmicsVetoTracks_*_*', 
        'keep *_muons_*_*', 
        'keep *_particleFlow_muons_*', 
        'drop *_muons_muons1stStep2muonsMap_*', 
        'drop recoIsoDepositedmValueMap_muons_*_*', 
        'drop doubleedmValueMap_muons_muPFIso*_*', 
        'keep recoTracks_standAloneMuons_*_*', 
        'keep recoTrackExtras_standAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_standAloneMuons_*_*', 
        'keep recoTracks_globalMuons_*_*', 
        'keep recoTrackExtras_globalMuons_*_*', 
        'keep recoTracks_tevMuons_*_*', 
        'keep recoTrackExtras_tevMuons_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_displacedTracks_*_*', 
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_*_*', 
        'keep recoTracks_displacedGlobalMuons_*_*', 
        'keep recoTrackExtras_displacedGlobalMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedGlobalMuons_*_*', 
        'keep recoTracks_cosmicMuons_*_*', 
        'keep recoMuons_muonsFromCosmics_*_*', 
        'keep recoTracks_cosmicMuons1Leg_*_*', 
        'keep recoMuons_muonsFromCosmics1Leg_*_*', 
        'keep recoTracks_refittedStandAloneMuons_*_*', 
        'keep recoTrackExtras_refittedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_refittedStandAloneMuons_*_*', 
        'keep recoTracks_displacedStandAloneMuons__*', 
        'keep recoTrackExtras_displacedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedStandAloneMuons_*_*', 
        'keep *_muIsoDepositTk_*_*', 
        'keep *_muIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muIsoDepositJets_*_*', 
        'keep *_muGlobalIsoDepositCtfTk_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muGlobalIsoDepositJets_*_*', 
        'keep *_softPFMuonsTagInfos_*_*', 
        'keep *_softPFElectronsTagInfos_*_*', 
        'keep *_softPFElectronBJetTags_*_*', 
        'keep *_softPFMuonBJetTags_*_*', 
        'keep *_pfImpactParameterTagInfos_*_*', 
        'keep *_pfTrackCountingHighEffBJetTags_*_*', 
        'keep *_pfJetProbabilityBJetTags_*_*', 
        'keep *_pfJetBProbabilityBJetTags_*_*', 
        'keep *_pfSecondaryVertexTagInfos_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderTagInfos_*_*', 
        'keep *_pfSimpleSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfCombinedSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfGhostTrackVertexTagInfos_*_*', 
        'keep *_pfGhostTrackBJetTags_*_*', 
        'keep *_pfCombinedMVAV2BJetTags_*_*', 
        'keep *_inclusiveCandidateSecondaryVertices_*_*', 
        'keep *_inclusiveCandidateSecondaryVerticesCvsL_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderCvsLTagInfos_*_*', 
        'keep *_pfCombinedCvsLJetTags_*_*', 
        'keep *_pfCombinedCvsBJetTags_*_*', 
        'keep *_pfChargeBJetTags_*_*', 
        'keep *_pfDeepCSVJetTags_*_*', 
        'keep *_pfDeepCMVAJetTags_*_*', 
        'keep recoRecoTauPiZeros_hpsPFTauProducer_pizeros_*', 
        'keep recoPFTaus_hpsPFTauProducer_*_*', 
        'keep *_hpsPFTauDiscriminationByDeadECALElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFinding_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingNewDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingOldDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseChargedIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByTightMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSum_*_*', 
        'keep *_hpsPFTauPUcorrPtSum_*_*', 
        'keep *_hpsPFTauChargedIsoPtSum_*_*', 
        'keep *_hpsPFTauTransverseImpactParameters_*_*', 
        'keep *_hpsPFTauFootprintCorrection_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeight_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6rawElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6LooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6MediumElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6TightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VTightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauChargedIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeightdR03_*_*', 
        'keep *_hpsPFTauFootprintCorrectiondR03_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalConedR03_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep  *_offlinePrimaryVertices__*', 
        'keep *_offlinePrimaryVerticesWithBS_*_*', 
        'keep *_offlinePrimaryVerticesFromCosmicTracks_*_*', 
        'keep *_nuclearInteractionMaker_*_*', 
        'keep *_generalV0Candidates_*_*', 
        'keep *_inclusiveSecondaryVertices_*_*', 
        'keep recoGsfElectronCores_gsfElectronCores_*_*', 
        'keep recoGsfElectronCores_gedGsfElectronCores_*_*', 
        'keep recoGsfElectrons_gsfElectrons_*_*', 
        'keep recoGsfElectrons_gedGsfElectrons_*_*', 
        'keep recoGsfElectronCores_uncleanedOnlyGsfElectronCores_*_*', 
        'keep recoGsfElectrons_uncleanedOnlyGsfElectrons_*_*', 
        'keep floatedmValueMap_eidRobustLoose_*_*', 
        'keep floatedmValueMap_eidRobustTight_*_*', 
        'keep floatedmValueMap_eidRobustHighEnergy_*_*', 
        'keep floatedmValueMap_eidLoose_*_*', 
        'keep floatedmValueMap_eidTight_*_*', 
        'keep *_egmGedGsfElectronPFIsolation_*_*', 
        'keep *_photonEcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonEcalPFClusterIsolationProducer_*_*', 
        'keep *_electronEcalPFClusterIsolationProducer_*_*', 
        'keep *_photonHcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonHcalPFClusterIsolationProducer_*_*', 
        'keep *_electronHcalPFClusterIsolationProducer_*_*', 
        'drop *_egmGsfElectronIDs_*_*', 
        'drop *_egmPhotonIDs_*_*', 
        'keep *_gedPhotonCore_*_*', 
        'keep *_gedPhotons_*_*', 
        'keep *_particleBasedIsolation_*_*', 
        'keep recoPhotons_mustachePhotons_*_*', 
        'keep recoPhotonCores_mustachePhotonCore_*_*', 
        'keep recoPhotons_ootPhotons_*_*', 
        'keep recoPhotonCores_ootPhotonCore_*_*', 
        'keep recoPhotons_photons_*_*', 
        'keep recoPhotonCores_photonCore_*_*', 
        'keep recoConversions_conversions_*_*', 
        'keep recoConversions_mustacheConversions_*_*', 
        'drop *_conversions_uncleanedConversions_*', 
        'drop *_gedPhotonsTmp_valMapPFEgammaCandToPhoton_*', 
        'keep recoConversions_allConversions_*_*', 
        'keep recoConversions_allConversionsOldEG_*_*', 
        'keep recoTracks_ckfOutInTracksFromConversions_*_*', 
        'keep recoTracks_ckfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfInOutTracksFromConversions_*_*', 
        'keep recoConversions_uncleanedOnlyAllConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep *_PhotonIDProdGED_*_*', 
        'keep recoRecoEcalCandidates_hfRecoEcalCandidate_*_*', 
        'keep *_hfEMClusters_*_*', 
        'keep *_gedGsfElectronCores_*_*', 
        'keep *_gedGsfElectrons_*_*', 
        'keep *_pixelTracks_*_*', 
        'keep *_pixelVertices_*_*', 
        'drop CaloTowersSorted_towerMakerPF_*_*', 
        'keep recoPFRecHits_particleFlowClusterECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHCAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterPS_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHBHE_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitPS_Cleaned_*', 
        'keep recoPFClusters_particleFlowClusterECAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHCAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHO_*_*', 
        'keep recoPFClusters_particleFlowClusterHF_*_*', 
        'keep recoPFClusters_particleFlowClusterPS_*_*', 
        'keep recoPFBlocks_particleFlowBlock_*_*', 
        'keep recoPFCandidates_particleFlowEGamma_*_*', 
        'keep recoCaloClusters_particleFlowEGamma_*_*', 
        'keep recoSuperClusters_particleFlowEGamma_*_*', 
        'keep recoConversions_particleFlowEGamma_*_*', 
        'keep recoPFCandidates_particleFlow_*_*', 
        'keep recoPFCandidates_particleFlowTmp_electrons_*', 
        'keep recoPFCandidates_particleFlowTmp_*_*', 
        'drop recoPFCandidates_particleFlowTmp__*', 
        'keep recoPFDisplacedVertexs_particleFlowDisplacedVertex_*_*', 
        'keep *_pfElectronTranslator_*_*', 
        'keep *_pfPhotonTranslator_*_*', 
        'keep *_particleFlow_electrons_*', 
        'keep *_particleFlow_photons_*', 
        'keep *_particleFlow_muons_*', 
        'keep *_trackerDrivenElectronSeeds_preid_*', 
        'keep *_particleFlowPtrs_*_*', 
        'keep *_particleFlowTmpPtrs_*_*', 
        'keep *_chargedHadronPFTrackIsolation_*_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep L1MuGMTReadoutCollection_gtDigis_*_*', 
        'keep L1GctEmCand*_gctDigis_*_*', 
        'keep L1GctJetCand*_gctDigis_*_*', 
        'keep L1GctEtHad*_gctDigis_*_*', 
        'keep L1GctEtMiss*_gctDigis_*_*', 
        'keep L1GctEtTotal*_gctDigis_*_*', 
        'keep L1GctHtMiss*_gctDigis_*_*', 
        'keep L1GctJetCounts*_gctDigis_*_*', 
        'keep L1GctHFRingEtSums*_gctDigis_*_*', 
        'keep L1GctHFBitCounts*_gctDigis_*_*', 
        'keep LumiDetails_lumiProducer_*_*', 
        'keep LumiSummary_lumiProducer_*_*', 
        'drop *_hlt*_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep L1AcceptBunchCrossings_scalersRawToDigi_*_*', 
        'keep L1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep Level1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep BeamSpotOnlines_scalersRawToDigi_*_*', 
        'keep DcsStatuss_scalersRawToDigi_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep *_tcdsDigis_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'keep *_pfIsolatedElectronsEI_*_*', 
        'keep *_pfIsolatedMuonsEI_*_*', 
        'keep recoPFJets_pfJetsEI_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTagsEI_*_*', 
        'keep recoPFTaus_pfTausEI_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByDecayModeFinding_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByIsolation_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep *_*_*_HLT' ) ),
    splitLevel = cms.untracked.int32(0)
)

process.FEVTSIMEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring( ('drop *', 
        'drop *', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep *_g4SimHits_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep CSCDetIdCSCComparatorDigiMuonDigiCollection_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*', 
        'keep EBSrFlagsSorted_simEcalDigis_*_*', 
        'keep EESrFlagsSorted_simEcalDigis_*_*', 
        'keep CrossingFramePlaybackInfoNew_*_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int6stdbitsetstdpairs_*_AffectedAPVList_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'keep *_mix_MergedTrackTruth_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep recoGenJets_*_*_*', 
        'keep *_genParticle_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep *_MEtoEDMConverter_*_*', 
        'keep *_randomEngineStateProducer_*_*', 
        'keep DetIdedmEDCollection_siStripDigis_*_*', 
        'keep DetIdedmEDCollection_siPixelDigis_*_*', 
        'keep PixelFEDChanneledmNewDetSetVector_siPixelDigis_*_*', 
        'keep *_siPixelClusters_*_*', 
        'keep *_siStripClusters_*_*', 
        'keep ClusterSummary_clusterSummaryProducer_*_*', 
        'keep *_dt1DRecHits_*_*', 
        'keep *_dt4DSegments_*_*', 
        'keep *_dt1DCosmicRecHits_*_*', 
        'keep *_dt4DCosmicSegments_*_*', 
        'keep *_csc2DRecHits_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*', 
        'keep *_hbhereco_*_*', 
        'keep *_hbheprereco_*_*', 
        'keep *_hfprereco_*_*', 
        'keep *_hfreco_*_*', 
        'keep *_horeco_*_*', 
        'keep HBHERecHitsSorted_hbherecoMB_*_*', 
        'keep HORecHitsSorted_horecoMB_*_*', 
        'keep HFRecHitsSorted_hfrecoMB_*_*', 
        'keep ZDCDataFramesSorted_hcalDigis_*_*', 
        'keep ZDCDataFramesSorted_castorDigis_*_*', 
        'keep QIE10DataFrameHcalDataFrameContainer_hcalDigis_ZDC_*', 
        'keep ZDCRecHitsSorted_zdcreco_*_*', 
        'keep *_reducedHcalRecHits_*_*', 
        'keep *_castorreco_*_*', 
        'keep HcalUnpackerReport_castorDigis_*_*', 
        'keep HcalUnpackerReport_hcalDigiAlCaMB_*_*', 
        'keep HcalUnpackerReport_hcalDigis_*_*', 
        'keep *_ecalPreshowerRecHit_*_*', 
        'keep *_ecalRecHit_*_*', 
        'keep *_ecalCompactTrigPrim_*_*', 
        'keep *_ecalTPSkim_*_*', 
        'keep EBSrFlagsSorted_ecalDigis__*', 
        'keep EESrFlagsSorted_ecalDigis__*', 
        'keep *_selectDigi_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEE_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEB_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsES_*_*', 
        'keep *_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep *_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_*', 
        'keep *_particleFlowSuperClusterECAL_*_*', 
        'keep *_particleFlowSuperClusterOOTECAL_*_*', 
        'drop recoClusterShapes_*_*_*', 
        'drop recoBasicClustersToOnerecoClusterShapesAssociation_*_*_*', 
        'drop recoBasicClusters_multi5x5BasicClusters_multi5x5BarrelBasicClusters_*', 
        'drop recoSuperClusters_multi5x5SuperClusters_multi5x5BarrelSuperClusters_*', 
        'keep *_CkfElectronCandidates_*_*', 
        'keep *_GsfGlobalElectronTest_*_*', 
        'keep *_electronMergedSeeds_*_*', 
        'keep recoGsfTracks_electronGsfTracks_*_*', 
        'keep recoGsfTrackExtras_electronGsfTracks_*_*', 
        'keep recoTrackExtras_electronGsfTracks_*_*', 
        'keep TrackingRecHitsOwned_electronGsfTracks_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTrackExtras_generalTracks_*_*', 
        'keep TrackingRecHitsOwned_generalTracks_*_*', 
        'keep *_generalTracks_MVAValues_*', 
        'keep *_generalTracks_MVAVals_*', 
        'keep TrackingRecHitsOwned_extraFromSeeds_*_*', 
        'keep uints_extraFromSeeds_*_*', 
        'keep recoTracks_beamhaloTracks_*_*', 
        'keep recoTrackExtras_beamhaloTracks_*_*', 
        'keep TrackingRecHitsOwned_beamhaloTracks_*_*', 
        'keep recoTracks_conversionStepTracks_*_*', 
        'keep recoTrackExtras_conversionStepTracks_*_*', 
        'keep TrackingRecHitsOwned_conversionStepTracks_*_*', 
        'keep *_ctfPixelLess_*_*', 
        'keep *_dedxTruncated40_*_*', 
        'keep *_dedxHitInfo_*_*', 
        'keep *_dedxHarmonic2_*_*', 
        'keep *_dedxPixelHarmonic2_*_*', 
        'keep *_trackExtrapolator_*_*', 
        'keep recoTracks_cosmicDCTracks_*_*', 
        'keep recoTrackExtras_cosmicDCTracks_*_*', 
        'keep TrackingRecHitsOwned_cosmicDCTracks_*_*', 
        'keep *_ak4CaloJets_*_*', 
        'keep *_ak4PFJets_*_*', 
        'keep *_ak4PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHSSoftDrop_*_*', 
        'keep *_cmsTopTagPFJetsCHS_*_*', 
        'keep *_JetPlusTrackZSPCorJetAntiKt4_*_*', 
        'keep *_ak4TrackJets_*_*', 
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', 
        'keep *_caloTowers_*_*', 
        'keep *_towerMaker_*_*', 
        'keep *_CastorTowerReco_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertex_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertexPF_*_*', 
        'keep *_ak4JetTracksAssociatorAtCaloFace_*_*', 
        'keep *_ak4JetTracksAssociatorExplicit_*_*', 
        'keep *_ak4JetExtender_*_*', 
        'keep *_ak4JetID_*_*', 
        'keep *_ak5CastorJets_*_*', 
        'keep *_ak5CastorJetID_*_*', 
        'keep *_ak7CastorJets_*_*', 
        'keep *_ak7CastorJetID_*_*', 
        'keep *_fixedGridRhoAll_*_*', 
        'keep *_fixedGridRhoFastjetAll_*_*', 
        'keep *_fixedGridRhoFastjetAllTmp_*_*', 
        'keep *_fixedGridRhoFastjetAllCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentral_*_*', 
        'keep *_fixedGridRhoFastjetCentralCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralChargedPileUp_*_*', 
        'keep *_fixedGridRhoFastjetCentralNeutral_*_*', 
        'keep *_ak8PFJetsCHSSoftDropMass_*_*', 
        'keep recoCaloMETs_caloMet_*_*', 
        'keep recoCaloMETs_caloMetBE_*_*', 
        'keep recoCaloMETs_caloMetBEFO_*_*', 
        'keep recoCaloMETs_caloMetM_*_*', 
        'keep recoPFMETs_pfMet_*_*', 
        'keep recoPFMETs_pfChMet_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_*_*', 
        'keep recoHcalNoiseRBXs_hcalnoise_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoEcalHaloData_EcalHaloData_*_*', 
        'keep recoGlobalHaloData_GlobalHaloData_*_*', 
        'keep recoHcalHaloData_HcalHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep *_MuonSeed_*_*', 
        'keep *_ancientMuonSeed_*_*', 
        'keep *_displacedMuonSeeds_*_*', 
        'keep TrackingRecHitsOwned_globalMuons_*_*', 
        'keep TrackingRecHitsOwned_tevMuons_*_*', 
        'keep *_CosmicMuonSeed_*_*', 
        'keep recoTrackExtras_cosmicMuons_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons_*_*', 
        'keep recoTrackExtras_cosmicMuons1Leg_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons1Leg_*_*', 
        'keep recoTracks_cosmicsVetoTracks_*_*', 
        'keep *_muons_*_*', 
        'keep *_particleFlow_muons_*', 
        'drop *_muons_muons1stStep2muonsMap_*', 
        'drop recoIsoDepositedmValueMap_muons_*_*', 
        'drop doubleedmValueMap_muons_muPFIso*_*', 
        'keep recoTracks_standAloneMuons_*_*', 
        'keep recoTrackExtras_standAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_standAloneMuons_*_*', 
        'keep recoTracks_globalMuons_*_*', 
        'keep recoTrackExtras_globalMuons_*_*', 
        'keep recoTracks_tevMuons_*_*', 
        'keep recoTrackExtras_tevMuons_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_displacedTracks_*_*', 
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_*_*', 
        'keep recoTracks_displacedGlobalMuons_*_*', 
        'keep recoTrackExtras_displacedGlobalMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedGlobalMuons_*_*', 
        'keep recoTracks_cosmicMuons_*_*', 
        'keep recoMuons_muonsFromCosmics_*_*', 
        'keep recoTracks_cosmicMuons1Leg_*_*', 
        'keep recoMuons_muonsFromCosmics1Leg_*_*', 
        'keep recoTracks_refittedStandAloneMuons_*_*', 
        'keep recoTrackExtras_refittedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_refittedStandAloneMuons_*_*', 
        'keep recoTracks_displacedStandAloneMuons__*', 
        'keep recoTrackExtras_displacedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedStandAloneMuons_*_*', 
        'keep *_muIsoDepositTk_*_*', 
        'keep *_muIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muIsoDepositJets_*_*', 
        'keep *_muGlobalIsoDepositCtfTk_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muGlobalIsoDepositJets_*_*', 
        'keep *_softPFMuonsTagInfos_*_*', 
        'keep *_softPFElectronsTagInfos_*_*', 
        'keep *_softPFElectronBJetTags_*_*', 
        'keep *_softPFMuonBJetTags_*_*', 
        'keep *_pfImpactParameterTagInfos_*_*', 
        'keep *_pfTrackCountingHighEffBJetTags_*_*', 
        'keep *_pfJetProbabilityBJetTags_*_*', 
        'keep *_pfJetBProbabilityBJetTags_*_*', 
        'keep *_pfSecondaryVertexTagInfos_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderTagInfos_*_*', 
        'keep *_pfSimpleSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfCombinedSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfGhostTrackVertexTagInfos_*_*', 
        'keep *_pfGhostTrackBJetTags_*_*', 
        'keep *_pfCombinedMVAV2BJetTags_*_*', 
        'keep *_inclusiveCandidateSecondaryVertices_*_*', 
        'keep *_inclusiveCandidateSecondaryVerticesCvsL_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderCvsLTagInfos_*_*', 
        'keep *_pfCombinedCvsLJetTags_*_*', 
        'keep *_pfCombinedCvsBJetTags_*_*', 
        'keep *_pfChargeBJetTags_*_*', 
        'keep *_pfDeepCSVJetTags_*_*', 
        'keep *_pfDeepCMVAJetTags_*_*', 
        'keep recoRecoTauPiZeros_hpsPFTauProducer_pizeros_*', 
        'keep recoPFTaus_hpsPFTauProducer_*_*', 
        'keep *_hpsPFTauDiscriminationByDeadECALElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFinding_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingNewDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingOldDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseChargedIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByTightMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSum_*_*', 
        'keep *_hpsPFTauPUcorrPtSum_*_*', 
        'keep *_hpsPFTauChargedIsoPtSum_*_*', 
        'keep *_hpsPFTauTransverseImpactParameters_*_*', 
        'keep *_hpsPFTauFootprintCorrection_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeight_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6rawElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6LooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6MediumElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6TightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VTightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauChargedIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeightdR03_*_*', 
        'keep *_hpsPFTauFootprintCorrectiondR03_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalConedR03_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep  *_offlinePrimaryVertices__*', 
        'keep *_offlinePrimaryVerticesWithBS_*_*', 
        'keep *_offlinePrimaryVerticesFromCosmicTracks_*_*', 
        'keep *_nuclearInteractionMaker_*_*', 
        'keep *_generalV0Candidates_*_*', 
        'keep *_inclusiveSecondaryVertices_*_*', 
        'keep recoGsfElectronCores_gsfElectronCores_*_*', 
        'keep recoGsfElectronCores_gedGsfElectronCores_*_*', 
        'keep recoGsfElectrons_gsfElectrons_*_*', 
        'keep recoGsfElectrons_gedGsfElectrons_*_*', 
        'keep recoGsfElectronCores_uncleanedOnlyGsfElectronCores_*_*', 
        'keep recoGsfElectrons_uncleanedOnlyGsfElectrons_*_*', 
        'keep floatedmValueMap_eidRobustLoose_*_*', 
        'keep floatedmValueMap_eidRobustTight_*_*', 
        'keep floatedmValueMap_eidRobustHighEnergy_*_*', 
        'keep floatedmValueMap_eidLoose_*_*', 
        'keep floatedmValueMap_eidTight_*_*', 
        'keep *_egmGedGsfElectronPFIsolation_*_*', 
        'keep *_photonEcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonEcalPFClusterIsolationProducer_*_*', 
        'keep *_electronEcalPFClusterIsolationProducer_*_*', 
        'keep *_photonHcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonHcalPFClusterIsolationProducer_*_*', 
        'keep *_electronHcalPFClusterIsolationProducer_*_*', 
        'drop *_egmGsfElectronIDs_*_*', 
        'drop *_egmPhotonIDs_*_*', 
        'keep *_gedPhotonCore_*_*', 
        'keep *_gedPhotons_*_*', 
        'keep *_particleBasedIsolation_*_*', 
        'keep recoPhotons_mustachePhotons_*_*', 
        'keep recoPhotonCores_mustachePhotonCore_*_*', 
        'keep recoPhotons_ootPhotons_*_*', 
        'keep recoPhotonCores_ootPhotonCore_*_*', 
        'keep recoPhotons_photons_*_*', 
        'keep recoPhotonCores_photonCore_*_*', 
        'keep recoConversions_conversions_*_*', 
        'keep recoConversions_mustacheConversions_*_*', 
        'drop *_conversions_uncleanedConversions_*', 
        'drop *_gedPhotonsTmp_valMapPFEgammaCandToPhoton_*', 
        'keep recoConversions_allConversions_*_*', 
        'keep recoConversions_allConversionsOldEG_*_*', 
        'keep recoTracks_ckfOutInTracksFromConversions_*_*', 
        'keep recoTracks_ckfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfInOutTracksFromConversions_*_*', 
        'keep recoConversions_uncleanedOnlyAllConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep *_PhotonIDProdGED_*_*', 
        'keep recoRecoEcalCandidates_hfRecoEcalCandidate_*_*', 
        'keep *_hfEMClusters_*_*', 
        'keep *_gedGsfElectronCores_*_*', 
        'keep *_gedGsfElectrons_*_*', 
        'keep *_pixelTracks_*_*', 
        'keep *_pixelVertices_*_*', 
        'drop CaloTowersSorted_towerMakerPF_*_*', 
        'keep recoPFRecHits_particleFlowClusterECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHCAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterPS_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHBHE_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitPS_Cleaned_*', 
        'keep recoPFClusters_particleFlowClusterECAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHCAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHO_*_*', 
        'keep recoPFClusters_particleFlowClusterHF_*_*', 
        'keep recoPFClusters_particleFlowClusterPS_*_*', 
        'keep recoPFBlocks_particleFlowBlock_*_*', 
        'keep recoPFCandidates_particleFlowEGamma_*_*', 
        'keep recoCaloClusters_particleFlowEGamma_*_*', 
        'keep recoSuperClusters_particleFlowEGamma_*_*', 
        'keep recoConversions_particleFlowEGamma_*_*', 
        'keep recoPFCandidates_particleFlow_*_*', 
        'keep recoPFCandidates_particleFlowTmp_electrons_*', 
        'keep recoPFCandidates_particleFlowTmp_*_*', 
        'drop recoPFCandidates_particleFlowTmp__*', 
        'keep recoPFDisplacedVertexs_particleFlowDisplacedVertex_*_*', 
        'keep *_pfElectronTranslator_*_*', 
        'keep *_pfPhotonTranslator_*_*', 
        'keep *_particleFlow_electrons_*', 
        'keep *_particleFlow_photons_*', 
        'keep *_particleFlow_muons_*', 
        'keep *_trackerDrivenElectronSeeds_preid_*', 
        'keep *_particleFlowPtrs_*_*', 
        'keep *_particleFlowTmpPtrs_*_*', 
        'keep *_chargedHadronPFTrackIsolation_*_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep L1MuGMTReadoutCollection_gtDigis_*_*', 
        'keep L1GctEmCand*_gctDigis_*_*', 
        'keep L1GctJetCand*_gctDigis_*_*', 
        'keep L1GctEtHad*_gctDigis_*_*', 
        'keep L1GctEtMiss*_gctDigis_*_*', 
        'keep L1GctEtTotal*_gctDigis_*_*', 
        'keep L1GctHtMiss*_gctDigis_*_*', 
        'keep L1GctJetCounts*_gctDigis_*_*', 
        'keep L1GctHFRingEtSums*_gctDigis_*_*', 
        'keep L1GctHFBitCounts*_gctDigis_*_*', 
        'keep LumiDetails_lumiProducer_*_*', 
        'keep LumiSummary_lumiProducer_*_*', 
        'drop *_hlt*_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep *_ak4GenJets_*_*', 
        'keep *_ak8GenJets_*_*', 
        'keep *_ak4GenJetsNoNu_*_*', 
        'keep *_ak8GenJetsNoNu_*_*', 
        'keep *_genParticle_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep SimTracks_g4SimHits_*_*', 
        'keep SimVertexs_g4SimHits_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*', 
        'keep *_muonSimClassifier_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'keep L1AcceptBunchCrossings_scalersRawToDigi_*_*', 
        'keep L1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep Level1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep BeamSpotOnlines_scalersRawToDigi_*_*', 
        'keep DcsStatuss_scalersRawToDigi_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'keep *_pfIsolatedElectronsEI_*_*', 
        'keep *_pfIsolatedMuonsEI_*_*', 
        'keep recoPFJets_pfJetsEI_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTagsEI_*_*', 
        'keep recoPFTaus_pfTausEI_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByDecayModeFinding_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByIsolation_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep *_tcdsDigis_*_*' ) ),
    splitLevel = cms.untracked.int32(0)
)

process.GENRAWEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *', 
        'drop *', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep SimTracks_g4SimHits_*_*', 
        'keep SimVertexs_g4SimHits_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*', 
        'keep *_muonSimClassifier_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep recoGenJets_*_*_*', 
        'keep *_genParticle_*_*', 
        'keep *_MEtoEDMConverter_*_*', 
        'keep *_randomEngineStateProducer_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep *_logErrorHarvester_*_*'),
    splitLevel = cms.untracked.int32(0)
)

process.GeneratorInterfaceAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*')
)

process.GeneratorInterfaceLHE = cms.PSet(
    outputCommands = cms.untracked.vstring('keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep *_externalLHEProducer_LHEScriptOutput_*')
)

process.GeneratorInterfaceRAW = cms.PSet(
    outputCommands = cms.untracked.vstring('keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*')
)

process.GeneratorInterfaceRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*')
)

process.HFRecalParameterBlock = cms.PSet(
    HFdepthOneParameterA = cms.vdouble(0.004123, 0.00602, 0.008201, 0.010489, 0.013379, 
        0.016997, 0.021464, 0.027371, 0.034195, 0.044807, 
        0.058939, 0.125497),
    HFdepthOneParameterB = cms.vdouble(-4e-06, -2e-06, 0.0, 4e-06, 1.5e-05, 
        2.6e-05, 6.3e-05, 8.4e-05, 0.00016, 0.000107, 
        0.000425, 0.000209),
    HFdepthTwoParameterA = cms.vdouble(0.002861, 0.004168, 0.0064, 0.008388, 0.011601, 
        0.014425, 0.018633, 0.023232, 0.028274, 0.035447, 
        0.051579, 0.086593),
    HFdepthTwoParameterB = cms.vdouble(-2e-06, -0.0, -7e-06, -6e-06, -2e-06, 
        1e-06, 1.9e-05, 3.1e-05, 6.7e-05, 1.2e-05, 
        0.000157, -3e-06)
)

process.HLTDEBUGEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring( ('drop *', 
        'keep *_logErrorHarvester_*_*', 
        'drop *_hlt*_*_*', 
        'keep *_hltAK4CaloJetsCorrectedIDPassed_*_*', 
        'keep *_hltAK4CaloJetsIDPassed_*_*', 
        'keep *_hltAK4CaloJets_*_*', 
        'keep *_hltAK4PFJetsCorrected_*_*', 
        'keep *_hltAK4PFJetsForTaus_*_*', 
        'keep *_hltAK4PFJets_*_*', 
        'keep *_hltAlCaEtaEBRechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaEtaEBRechitsToDigis_*_*', 
        'keep *_hltAlCaEtaEERechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaEtaEERechitsToDigis_*_*', 
        'keep *_hltAlCaEtaRecHitsFilterEBonly_*_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegionalLowPU_etaEcalRecHitsES_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegional_etaEcalRecHitsES_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonly_*_*', 
        'keep *_hltAlCaEtaRecHitsFilter_*_*', 
        'keep *_hltAlCaPhiSymStream_*_*', 
        'keep *_hltAlCaPi0EBRechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaPi0EBRechitsToDigis_*_*', 
        'keep *_hltAlCaPi0EERechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaPi0EERechitsToDigis_*_*', 
        'keep *_hltAlCaPi0RecHitsFilterEBonly_*_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegionalLowPU_pi0EcalRecHitsES_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegional_pi0EcalRecHitsES_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonly_*_*', 
        'keep *_hltAlCaPi0RecHitsFilter_*_*', 
        'keep *_hltBLifetimeL25AssociatorbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL25BJetTagsbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL25TagInfosbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3AssociatorbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3BJetTagsbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3TagInfosbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeRegionalCtfWithMaterialTracksbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBSoftMuonDiJet110Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet110Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet20Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet20Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet40Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet40Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet70Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet70Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonMu5L3_*_*', 
        'keep *_hltCSVJetTagSingleTopEle27_*_*', 
        'keep *_hltCSVJetTagSingleTopIsoMu24_*_*', 
        'keep *_hltCaloJetCorrectedRegional_*_*', 
        'keep *_hltCaloJetCorrected_*_*', 
        'keep *_hltCaloJetL1FastJetCorrected_*_*', 
        'keep *_hltCaloStage2Digis_*_*', 
        'keep *_hltCleanedCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCleanedHiCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCombinedSecondaryVertexBJetTagsCalo_*_*', 
        'keep *_hltCombinedSecondaryVertexBJetTagsPF_*_*', 
        'keep *_hltConvPFTausTightIsoTrackFindingIsolation_*_*', 
        'keep *_hltConvPFTausTightIsoTrackFinding_*_*', 
        'keep *_hltConvPFTausTightIsoTrackPt5Isolation_*_*', 
        'keep *_hltConvPFTausTightIsoTrackPt5_*_*', 
        'keep *_hltConvPFTausTightIso_*_*', 
        'keep *_hltConvPFTausTrackFindingLooseIsolation_*_*', 
        'keep *_hltConvPFTausTrackFinding_*_*', 
        'keep *_hltConvPFTaus_*_*', 
        'keep *_hltCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCorrectedIslandEndcapSuperClustersHI_*_*', 
        'keep *_hltCsc2DRecHits_*_*', 
        'keep *_hltCscSegments_*_*', 
        'keep *_hltDisplacedhltIter4PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolationProng4L1HLTMatched_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolationProng4_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolation_*_*', 
        'keep *_hltDoublePFTau25TrackPt5_*_*', 
        'keep *_hltDoublePFTau25_*_*', 
        'keep *_hltDoublePFTauTightIso45Track5_*_*', 
        'keep *_hltDoublePFTauTightIso45Track_*_*', 
        'keep *_hltDt4DSegments_*_*', 
        'keep *_hltEcalPhiSymFilter_*_*', 
        'keep *_hltEcalRecHitAll_*_*', 
        'keep *_hltEcalRecHit_*_*', 
        'keep *_hltEgammaCandidates_*_*', 
        'keep *_hltEgammaGsfElectrons_*_*', 
        'keep *_hltEgammaGsfTracks_*_*', 
        'keep *_hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1SingleEG18orL1SingleEG20_*_*', 
        'keep *_hltEle20CaloIdVTTrkIdTDphiFilter_*_*', 
        'keep *_hltEle27WP85PixelMatchElectronsL1Seeded_*_*', 
        'keep *_hltElectronsVertex_*_*', 
        'keep *_hltFEDSelectorLumiPixels_*_*', 
        'keep *_hltFastPVPixelTracksMerger_*_*', 
        'keep *_hltFastPVPixelTracksRecover_*_*', 
        'keep *_hltFastPVPixelTracks_*_*', 
        'keep *_hltFastPVPixelVertices3D_*_*', 
        'keep *_hltFastPVPixelVertices_*_*', 
        'keep *_hltFastPixelBLifetimeL3TagInfos_*_*', 
        'keep *_hltFastPrimaryVertex_*_*', 
        'keep *_hltFilterDoubleIsoPFTau45Trk5LeadTrack5IsolationL1HLTMatched_*_*', 
        'keep *_hltFilterL2EtCutDoublePFIsoTau45Trk5_*_*', 
        'keep *_hltFilterL2EtCutSingleIsoPFTau35Trk20MET70_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20LeadTrackPt20_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20MET60LeadTrack20IsolationL1HLTMatched_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20MET70LeadTrack20IsolationL1HLTMatched_*_*', 
        'keep *_hltGmtStage2Digis_*_*', 
        'keep *_hltGoodOnlinePVs_*_*', 
        'keep *_hltGtStage2Digis_*_*', 
        'keep *_hltHICaloJetCorrected_*_*', 
        'keep *_hltHICaloJetIDPassed_*_*', 
        'keep *_hltHIGoodLooseTracks_*_*', 
        'keep *_hltHIPixel3PrimTracks_*_*', 
        'keep *_hltHISelectedVertex_*_*', 
        'keep *_hltHISiPixelClusters_*_*', 
        'keep *_hltHITIPTCorrectorHB_*_*', 
        'keep *_hltHITIPTCorrectorHE_*_*', 
        'keep *_hltHbhereco_*_*', 
        'keep *_hltHfreco_*_*', 
        'keep *_hltHiCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltHiCorrectedIslandEndcapSuperClustersHI_*_*', 
        'keep *_hltHiIslandSuperClustersHI_*_*', 
        'keep *_hltHoreco_*_*', 
        'keep *_hltIsolPixelTrackProdHB_*_*', 
        'keep *_hltIsolPixelTrackProdHE_*_*', 
        'keep *_hltIter0ElectronsTrackSelectionHighPurity_*_*', 
        'keep *_hltIter0HighPtTkMuPixelTracks_*_*', 
        'keep *_hltIter0HighPtTkMuTrackSelectionHighPurity_*_*', 
        'keep *_hltIter0PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter1Merged_*_*', 
        'keep *_hltIter1PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter2HighPtTkMuMerged_*_*', 
        'keep *_hltIter2HighPtTkMuTrackSelectionHighPurity_*_*', 
        'keep *_hltIter2MergedForBTag_*_*', 
        'keep *_hltIter2MergedForElectrons_*_*', 
        'keep *_hltIter2Merged_*_*', 
        'keep *_hltIter2PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter3Merged_*_*', 
        'keep *_hltIter4Merged_*_*', 
        'keep *_hltIterativeCone5PileupSubtractionCaloJets_*_*', 
        'keep *_hltL1HLTSingleIsoPFTau35Trk20Met60JetsMatch_*_*', 
        'keep *_hltL1IsoElectronTrackIsol_*_*', 
        'keep *_hltL1NonIsoElectronTrackIsol_*_*', 
        'keep *_hltL1SeededRecoEcalCandidate_*_*', 
        'keep *_hltL1sDoubleTauJet44Eta2p17orDoubleJet64Central_*_*', 
        'keep *_hltL1sDoubleTauJet44erorDoubleJetC64_*_*', 
        'keep *_hltL1sL1EG18er_*_*', 
        'keep *_hltL1sL1ETM36ORETM40_*_*', 
        'keep *_hltL1sL1Jet52ETM30_*_*', 
        'keep *_hltL1sL1SingleEG12_*_*', 
        'keep *_hltL1sL1SingleEG15_*_*', 
        'keep *_hltL1sL1SingleEG18orL1SingleEG20_*_*', 
        'keep *_hltL1sL1SingleMu10_*_*', 
        'keep *_hltL1sL1SingleMu14Eta2p1_*_*', 
        'keep *_hltL1sMu16Eta2p1_*_*', 
        'keep *_hltL2MuonCandidatesNoVtx_*_*', 
        'keep *_hltL2MuonCandidates_*_*', 
        'keep *_hltL2MuonSeeds_*_*', 
        'keep *_hltL2Muons_*_*', 
        'keep *_hltL2TauJets_*_*', 
        'keep *_hltL3MuonCandidates_*_*', 
        'keep *_hltL3MuonsIOHit_*_*', 
        'keep *_hltL3MuonsLinksCombination_*_*', 
        'keep *_hltL3MuonsOIHit_*_*', 
        'keep *_hltL3MuonsOIState_*_*', 
        'keep *_hltL3Muons_*_*', 
        'keep *_hltL3NoFiltersNoVtxMuonCandidates_*_*', 
        'keep *_hltL3NoFiltersNoVtxMuons_*_*', 
        'keep *_hltL3SecondaryVertexTagInfos_*_*', 
        'keep *_hltL3TkFromL2OICombination_*_*', 
        'keep *_hltL3TkTracksFromL2IOHit_*_*', 
        'keep *_hltL3TkTracksFromL2OIHit_*_*', 
        'keep *_hltL3TkTracksFromL2OIState_*_*', 
        'keep *_hltL3TkTracksFromL2_*_*', 
        'keep *_hltL3TrackCandidateFromL2IOHit_*_*', 
        'keep *_hltL3TrackCandidateFromL2OIHit_*_*', 
        'keep *_hltL3TrackCandidateFromL2OIState_*_*', 
        'keep *_hltL3TrajSeedIOHit_*_*', 
        'keep *_hltL3TrajSeedOIHit_*_*', 
        'keep *_hltL3TrajSeedOIState_*_*', 
        'keep *_hltL3TrajectorySeed_*_*', 
        'keep *_hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f18QL3crIsoRhoFiltered0p15_*_*', 
        'keep *_hltLeadingCentralJets30SingleTopEle27_*_*', 
        'keep *_hltLeadingCentralJets30SingleTopIsoMu24_*_*', 
        'keep *_hltMergedTracks_*_*', 
        'keep *_hltMet_*_*', 
        'keep *_hltMu8Ele17CaloIdTCaloIsoVLPixelMatchFilter_*_*', 
        'keep *_hltMuTrackJpsiCtfTrackCands_*_*', 
        'keep *_hltMuTrackJpsiPixelTrackCands_*_*', 
        'keep *_hltMuonCSCDigis_*_*', 
        'keep *_hltMuonCSCDigis_MuonCSCStripDigi_*', 
        'keep *_hltMuonCSCDigis_MuonCSCWireDigi_*', 
        'keep *_hltMuonDTDigis_*_*', 
        'keep *_hltMuonRPCDigis_*_*', 
        'keep *_hltOnlineBeamSpot_*_*', 
        'keep *_hltOnlinePrimaryVertices_*_*', 
        'keep *_hltOverlapFilterEle20LooseIsoPFTau20OldVersion_*_*', 
        'keep *_hltOverlapFilterIsoEle20MediumIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15IsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15MediumIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15TightIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu18LooseIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu18PFTau25TrackPt5Prong4_*_*', 
        'keep *_hltPFJetForBtag_*_*', 
        'keep *_hltPFMuonMerging_*_*', 
        'keep *_hltPFTau15TrackLooseIso_*_*', 
        'keep *_hltPFTau15Track_*_*', 
        'keep *_hltPFTau15_*_*', 
        'keep *_hltPFTau20IsoMuVertex_*_*', 
        'keep *_hltPFTau20TrackLooseIso_*_*', 
        'keep *_hltPFTau20Track_*_*', 
        'keep *_hltPFTau20_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolationProng4IsoMuVertex_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolationProng4_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolation_*_*', 
        'keep *_hltPFTau25TrackPt5_*_*', 
        'keep *_hltPFTau25_*_*', 
        'keep *_hltPFTau35TrackPt20LooseIsoProng2_*_*', 
        'keep *_hltPFTau35TrackPt20LooseIso_*_*', 
        'keep *_hltPFTau35TrackPt20_*_*', 
        'keep *_hltPFTau35Track_*_*', 
        'keep *_hltPFTau35_*_*', 
        'keep *_hltPFTauEleVertex20_*_*', 
        'keep *_hltPFTauJetTracksAssociator_*_*', 
        'keep *_hltPFTauMediumIso20TrackMediumIso_*_*', 
        'keep *_hltPFTauMediumIso20Track_*_*', 
        'keep *_hltPFTauMediumIso20_*_*', 
        'keep *_hltPFTauMediumIso35Track_*_*', 
        'keep *_hltPFTauMediumIso35_*_*', 
        'keep *_hltPFTauTagInfo_*_*', 
        'keep *_hltPFTauTightIso20TrackTightIso_*_*', 
        'keep *_hltPFTauTightIso20Track_*_*', 
        'keep *_hltPFTauTightIso20_*_*', 
        'keep *_hltPFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltParticleFlowForTaus_*_*', 
        'keep *_hltParticleFlow_*_*', 
        'keep *_hltPixelMatch3HitElectronsActivity_*_*', 
        'keep *_hltPixelMatch3HitElectronsL1Seeded_*_*', 
        'keep *_hltPixelMatchCleanElectronsL1Seeded_*_*', 
        'keep *_hltPixelMatchElectronsActivity_*_*', 
        'keep *_hltPixelMatchElectronsL1Iso_*_*', 
        'keep *_hltPixelMatchElectronsL1NonIso_*_*', 
        'keep *_hltPixelMatchElectronsL1Seeded_*_*', 
        'keep *_hltPixelTracksElectrons_*_*', 
        'keep *_hltPixelTracksMerged_*_*', 
        'keep *_hltPixelTracks_*_*', 
        'keep *_hltPixelVertices3DbbPhi_*_*', 
        'keep *_hltPixelVertices_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidateSC4_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidateSC5_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidate_*_*', 
        'keep *_hltRpcRecHits_*_*', 
        'keep *_hltSelector4CentralJetsL1FastJet_*_*', 
        'keep *_hltSelector8CentralJetsL1FastJet_*_*', 
        'keep *_hltSelectorJets20L1FastJet_*_*', 
        'keep *_hltSiPixelCluster_*_*', 
        'keep *_hltSiPixelClusters_*_*', 
        'keep *_hltSiStripClusters_*_*', 
        'keep *_hltSiStripRawToClustersFacility_*_*', 
        'keep *_hltSingleMu15L3Filtered15_*_*', 
        'keep *_hltSingleMuIsoL1s14L3IsoFiltered15eta2p1_*_*', 
        'keep *_hltSingleMuIsoL3IsoFiltered15_*_*', 
        'keep *_hltTowerMakerForAll_*_*', 
        'keep *_hltTowerMakerForMuons_*_*', 
        'keep *_hltTriggerSummaryAOD_*_*', 
        'keep *_hltTriggerSummaryRAW_*_*', 
        'keep *_hltTrimmedPixelVertices_*_*', 
        'keep *_hltVerticesL3_*_*', 
        'keep *_hltVerticesPFFilter_*_*', 
        'keep *_hltVerticesPFSelector_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_rawDataRepacker_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_virginRawDataRepacker_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep L2MuonTrajectorySeeds_hltL2MuonSeeds_*_*', 
        'keep L3MuonTrajectorySeeds_hltHIL3TrajSeedOIHit_*_*', 
        'keep L3MuonTrajectorySeeds_hltHIL3TrajectorySeed_*_*', 
        'keep L3MuonTrajectorySeeds_hltL3TrajSeedOIState_*_*', 
        'keep LumiScalerss_hltScalersRawToDigi_*_*', 
        'keep SiPixelClusteredmNewDetSetVector_hltSiPixelClusters_*_*', 
        'keep TrackCandidates_hltHIL3TrackCandidateFromL2OIHit_*_*', 
        'keep TrackCandidates_hltHIL3TrackCandidateFromL2OIState_*_*', 
        'keep TrackingRecHitsOwned_hltL3Muons_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep recoCaloJets_*_*_*', 
        'keep recoCaloMETs_*_*_*', 
        'keep recoCaloMETs_hltMet_*_*', 
        'keep recoCompositeCandidates_*_*_*', 
        'keep recoElectrons_*_*_*', 
        'keep recoIsolatedPixelTrackCandidates_*_*_*', 
        'keep recoMETs_*_*_*', 
        'keep recoPFJets_*_*_*', 
        'keep recoPFTaus_*_*_*', 
        'keep recoRecoChargedCandidates_*_*_*', 
        'keep recoRecoChargedCandidates_hltHIL3MuonCandidates_*_*', 
        'keep recoRecoChargedCandidates_hltL2MuonCandidates_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsoHLTClusterShape_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonEcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalForHE_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsoHLTClusterShape_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonEcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalForHE_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalIsol_*_*', 
        'keep recoRecoEcalCandidates_*_*_*', 
        'keep recoRecoEcalCandidates_hltL1IsoRecoEcalCandidate_*_*', 
        'keep recoRecoEcalCandidates_hltL1NonIsoRecoEcalCandidate_*_*', 
        'keep recoTrackExtras_hltHIL3MuonsOIHit_*_*', 
        'keep recoTrackExtras_hltHIL3MuonsOIState_*_*', 
        'keep recoTracks_hltHIL3MuonsOIHit_*_*', 
        'keep recoTracks_hltHIL3MuonsOIState_*_*', 
        'keep recoTracks_hltHIL3Muons_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2OIHit_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2OIState_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2_*_*', 
        'keep triggerTriggerEventWithRefs_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep triggerTriggerFilterObjectWithRefs_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*' ) ),
    splitLevel = cms.untracked.int32(0)
)

process.HLTDebugFEVT = cms.PSet(
    outputCommands = cms.vstring( ('drop *_hlt*_*_*', 
        'keep *_hltAK4CaloJetsCorrectedIDPassed_*_*', 
        'keep *_hltAK4CaloJetsIDPassed_*_*', 
        'keep *_hltAK4CaloJets_*_*', 
        'keep *_hltAK4PFJetsCorrected_*_*', 
        'keep *_hltAK4PFJetsForTaus_*_*', 
        'keep *_hltAK4PFJets_*_*', 
        'keep *_hltAlCaEtaEBRechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaEtaEBRechitsToDigis_*_*', 
        'keep *_hltAlCaEtaEERechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaEtaEERechitsToDigis_*_*', 
        'keep *_hltAlCaEtaRecHitsFilterEBonly_*_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegionalLowPU_etaEcalRecHitsES_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegional_etaEcalRecHitsES_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonly_*_*', 
        'keep *_hltAlCaEtaRecHitsFilter_*_*', 
        'keep *_hltAlCaPhiSymStream_*_*', 
        'keep *_hltAlCaPi0EBRechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaPi0EBRechitsToDigis_*_*', 
        'keep *_hltAlCaPi0EERechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaPi0EERechitsToDigis_*_*', 
        'keep *_hltAlCaPi0RecHitsFilterEBonly_*_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegionalLowPU_pi0EcalRecHitsES_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegional_pi0EcalRecHitsES_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonly_*_*', 
        'keep *_hltAlCaPi0RecHitsFilter_*_*', 
        'keep *_hltBLifetimeL25AssociatorbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL25BJetTagsbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL25TagInfosbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3AssociatorbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3BJetTagsbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3TagInfosbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeRegionalCtfWithMaterialTracksbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBSoftMuonDiJet110Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet110Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet20Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet20Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet40Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet40Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet70Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet70Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonMu5L3_*_*', 
        'keep *_hltCSVJetTagSingleTopEle27_*_*', 
        'keep *_hltCSVJetTagSingleTopIsoMu24_*_*', 
        'keep *_hltCaloJetCorrectedRegional_*_*', 
        'keep *_hltCaloJetCorrected_*_*', 
        'keep *_hltCaloJetL1FastJetCorrected_*_*', 
        'keep *_hltCaloStage2Digis_*_*', 
        'keep *_hltCleanedCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCleanedHiCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCombinedSecondaryVertexBJetTagsCalo_*_*', 
        'keep *_hltCombinedSecondaryVertexBJetTagsPF_*_*', 
        'keep *_hltConvPFTausTightIsoTrackFindingIsolation_*_*', 
        'keep *_hltConvPFTausTightIsoTrackFinding_*_*', 
        'keep *_hltConvPFTausTightIsoTrackPt5Isolation_*_*', 
        'keep *_hltConvPFTausTightIsoTrackPt5_*_*', 
        'keep *_hltConvPFTausTightIso_*_*', 
        'keep *_hltConvPFTausTrackFindingLooseIsolation_*_*', 
        'keep *_hltConvPFTausTrackFinding_*_*', 
        'keep *_hltConvPFTaus_*_*', 
        'keep *_hltCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCorrectedIslandEndcapSuperClustersHI_*_*', 
        'keep *_hltCsc2DRecHits_*_*', 
        'keep *_hltCscSegments_*_*', 
        'keep *_hltDisplacedhltIter4PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolationProng4L1HLTMatched_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolationProng4_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolation_*_*', 
        'keep *_hltDoublePFTau25TrackPt5_*_*', 
        'keep *_hltDoublePFTau25_*_*', 
        'keep *_hltDoublePFTauTightIso45Track5_*_*', 
        'keep *_hltDoublePFTauTightIso45Track_*_*', 
        'keep *_hltDt4DSegments_*_*', 
        'keep *_hltEcalPhiSymFilter_*_*', 
        'keep *_hltEcalRecHitAll_*_*', 
        'keep *_hltEcalRecHit_*_*', 
        'keep *_hltEgammaCandidates_*_*', 
        'keep *_hltEgammaGsfElectrons_*_*', 
        'keep *_hltEgammaGsfTracks_*_*', 
        'keep *_hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1SingleEG18orL1SingleEG20_*_*', 
        'keep *_hltEle20CaloIdVTTrkIdTDphiFilter_*_*', 
        'keep *_hltEle27WP85PixelMatchElectronsL1Seeded_*_*', 
        'keep *_hltElectronsVertex_*_*', 
        'keep *_hltFEDSelectorLumiPixels_*_*', 
        'keep *_hltFastPVPixelTracksMerger_*_*', 
        'keep *_hltFastPVPixelTracksRecover_*_*', 
        'keep *_hltFastPVPixelTracks_*_*', 
        'keep *_hltFastPVPixelVertices3D_*_*', 
        'keep *_hltFastPVPixelVertices_*_*', 
        'keep *_hltFastPixelBLifetimeL3TagInfos_*_*', 
        'keep *_hltFastPrimaryVertex_*_*', 
        'keep *_hltFilterDoubleIsoPFTau45Trk5LeadTrack5IsolationL1HLTMatched_*_*', 
        'keep *_hltFilterL2EtCutDoublePFIsoTau45Trk5_*_*', 
        'keep *_hltFilterL2EtCutSingleIsoPFTau35Trk20MET70_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20LeadTrackPt20_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20MET60LeadTrack20IsolationL1HLTMatched_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20MET70LeadTrack20IsolationL1HLTMatched_*_*', 
        'keep *_hltGmtStage2Digis_*_*', 
        'keep *_hltGoodOnlinePVs_*_*', 
        'keep *_hltGtStage2Digis_*_*', 
        'keep *_hltHICaloJetCorrected_*_*', 
        'keep *_hltHICaloJetIDPassed_*_*', 
        'keep *_hltHIGoodLooseTracks_*_*', 
        'keep *_hltHIPixel3PrimTracks_*_*', 
        'keep *_hltHISelectedVertex_*_*', 
        'keep *_hltHISiPixelClusters_*_*', 
        'keep *_hltHITIPTCorrectorHB_*_*', 
        'keep *_hltHITIPTCorrectorHE_*_*', 
        'keep *_hltHbhereco_*_*', 
        'keep *_hltHfreco_*_*', 
        'keep *_hltHiCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltHiCorrectedIslandEndcapSuperClustersHI_*_*', 
        'keep *_hltHiIslandSuperClustersHI_*_*', 
        'keep *_hltHoreco_*_*', 
        'keep *_hltIsolPixelTrackProdHB_*_*', 
        'keep *_hltIsolPixelTrackProdHE_*_*', 
        'keep *_hltIter0ElectronsTrackSelectionHighPurity_*_*', 
        'keep *_hltIter0HighPtTkMuPixelTracks_*_*', 
        'keep *_hltIter0HighPtTkMuTrackSelectionHighPurity_*_*', 
        'keep *_hltIter0PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter1Merged_*_*', 
        'keep *_hltIter1PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter2HighPtTkMuMerged_*_*', 
        'keep *_hltIter2HighPtTkMuTrackSelectionHighPurity_*_*', 
        'keep *_hltIter2MergedForBTag_*_*', 
        'keep *_hltIter2MergedForElectrons_*_*', 
        'keep *_hltIter2Merged_*_*', 
        'keep *_hltIter2PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter3Merged_*_*', 
        'keep *_hltIter4Merged_*_*', 
        'keep *_hltIterativeCone5PileupSubtractionCaloJets_*_*', 
        'keep *_hltL1HLTSingleIsoPFTau35Trk20Met60JetsMatch_*_*', 
        'keep *_hltL1IsoElectronTrackIsol_*_*', 
        'keep *_hltL1NonIsoElectronTrackIsol_*_*', 
        'keep *_hltL1SeededRecoEcalCandidate_*_*', 
        'keep *_hltL1sDoubleTauJet44Eta2p17orDoubleJet64Central_*_*', 
        'keep *_hltL1sDoubleTauJet44erorDoubleJetC64_*_*', 
        'keep *_hltL1sL1EG18er_*_*', 
        'keep *_hltL1sL1ETM36ORETM40_*_*', 
        'keep *_hltL1sL1Jet52ETM30_*_*', 
        'keep *_hltL1sL1SingleEG12_*_*', 
        'keep *_hltL1sL1SingleEG15_*_*', 
        'keep *_hltL1sL1SingleEG18orL1SingleEG20_*_*', 
        'keep *_hltL1sL1SingleMu10_*_*', 
        'keep *_hltL1sL1SingleMu14Eta2p1_*_*', 
        'keep *_hltL1sMu16Eta2p1_*_*', 
        'keep *_hltL2MuonCandidatesNoVtx_*_*', 
        'keep *_hltL2MuonCandidates_*_*', 
        'keep *_hltL2MuonSeeds_*_*', 
        'keep *_hltL2Muons_*_*', 
        'keep *_hltL2TauJets_*_*', 
        'keep *_hltL3MuonCandidates_*_*', 
        'keep *_hltL3MuonsIOHit_*_*', 
        'keep *_hltL3MuonsLinksCombination_*_*', 
        'keep *_hltL3MuonsOIHit_*_*', 
        'keep *_hltL3MuonsOIState_*_*', 
        'keep *_hltL3Muons_*_*', 
        'keep *_hltL3NoFiltersNoVtxMuonCandidates_*_*', 
        'keep *_hltL3NoFiltersNoVtxMuons_*_*', 
        'keep *_hltL3SecondaryVertexTagInfos_*_*', 
        'keep *_hltL3TkFromL2OICombination_*_*', 
        'keep *_hltL3TkTracksFromL2IOHit_*_*', 
        'keep *_hltL3TkTracksFromL2OIHit_*_*', 
        'keep *_hltL3TkTracksFromL2OIState_*_*', 
        'keep *_hltL3TkTracksFromL2_*_*', 
        'keep *_hltL3TrackCandidateFromL2IOHit_*_*', 
        'keep *_hltL3TrackCandidateFromL2OIHit_*_*', 
        'keep *_hltL3TrackCandidateFromL2OIState_*_*', 
        'keep *_hltL3TrajSeedIOHit_*_*', 
        'keep *_hltL3TrajSeedOIHit_*_*', 
        'keep *_hltL3TrajSeedOIState_*_*', 
        'keep *_hltL3TrajectorySeed_*_*', 
        'keep *_hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f18QL3crIsoRhoFiltered0p15_*_*', 
        'keep *_hltLeadingCentralJets30SingleTopEle27_*_*', 
        'keep *_hltLeadingCentralJets30SingleTopIsoMu24_*_*', 
        'keep *_hltMergedTracks_*_*', 
        'keep *_hltMet_*_*', 
        'keep *_hltMu8Ele17CaloIdTCaloIsoVLPixelMatchFilter_*_*', 
        'keep *_hltMuTrackJpsiCtfTrackCands_*_*', 
        'keep *_hltMuTrackJpsiPixelTrackCands_*_*', 
        'keep *_hltMuonCSCDigis_*_*', 
        'keep *_hltMuonCSCDigis_MuonCSCStripDigi_*', 
        'keep *_hltMuonCSCDigis_MuonCSCWireDigi_*', 
        'keep *_hltMuonDTDigis_*_*', 
        'keep *_hltMuonRPCDigis_*_*', 
        'keep *_hltOnlineBeamSpot_*_*', 
        'keep *_hltOnlinePrimaryVertices_*_*', 
        'keep *_hltOverlapFilterEle20LooseIsoPFTau20OldVersion_*_*', 
        'keep *_hltOverlapFilterIsoEle20MediumIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15IsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15MediumIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15TightIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu18LooseIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu18PFTau25TrackPt5Prong4_*_*', 
        'keep *_hltPFJetForBtag_*_*', 
        'keep *_hltPFMuonMerging_*_*', 
        'keep *_hltPFTau15TrackLooseIso_*_*', 
        'keep *_hltPFTau15Track_*_*', 
        'keep *_hltPFTau15_*_*', 
        'keep *_hltPFTau20IsoMuVertex_*_*', 
        'keep *_hltPFTau20TrackLooseIso_*_*', 
        'keep *_hltPFTau20Track_*_*', 
        'keep *_hltPFTau20_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolationProng4IsoMuVertex_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolationProng4_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolation_*_*', 
        'keep *_hltPFTau25TrackPt5_*_*', 
        'keep *_hltPFTau25_*_*', 
        'keep *_hltPFTau35TrackPt20LooseIsoProng2_*_*', 
        'keep *_hltPFTau35TrackPt20LooseIso_*_*', 
        'keep *_hltPFTau35TrackPt20_*_*', 
        'keep *_hltPFTau35Track_*_*', 
        'keep *_hltPFTau35_*_*', 
        'keep *_hltPFTauEleVertex20_*_*', 
        'keep *_hltPFTauJetTracksAssociator_*_*', 
        'keep *_hltPFTauMediumIso20TrackMediumIso_*_*', 
        'keep *_hltPFTauMediumIso20Track_*_*', 
        'keep *_hltPFTauMediumIso20_*_*', 
        'keep *_hltPFTauMediumIso35Track_*_*', 
        'keep *_hltPFTauMediumIso35_*_*', 
        'keep *_hltPFTauTagInfo_*_*', 
        'keep *_hltPFTauTightIso20TrackTightIso_*_*', 
        'keep *_hltPFTauTightIso20Track_*_*', 
        'keep *_hltPFTauTightIso20_*_*', 
        'keep *_hltPFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltParticleFlowForTaus_*_*', 
        'keep *_hltParticleFlow_*_*', 
        'keep *_hltPixelMatch3HitElectronsActivity_*_*', 
        'keep *_hltPixelMatch3HitElectronsL1Seeded_*_*', 
        'keep *_hltPixelMatchCleanElectronsL1Seeded_*_*', 
        'keep *_hltPixelMatchElectronsActivity_*_*', 
        'keep *_hltPixelMatchElectronsL1Iso_*_*', 
        'keep *_hltPixelMatchElectronsL1NonIso_*_*', 
        'keep *_hltPixelMatchElectronsL1Seeded_*_*', 
        'keep *_hltPixelTracksElectrons_*_*', 
        'keep *_hltPixelTracksMerged_*_*', 
        'keep *_hltPixelTracks_*_*', 
        'keep *_hltPixelVertices3DbbPhi_*_*', 
        'keep *_hltPixelVertices_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidateSC4_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidateSC5_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidate_*_*', 
        'keep *_hltRpcRecHits_*_*', 
        'keep *_hltSelector4CentralJetsL1FastJet_*_*', 
        'keep *_hltSelector8CentralJetsL1FastJet_*_*', 
        'keep *_hltSelectorJets20L1FastJet_*_*', 
        'keep *_hltSiPixelCluster_*_*', 
        'keep *_hltSiPixelClusters_*_*', 
        'keep *_hltSiStripClusters_*_*', 
        'keep *_hltSiStripRawToClustersFacility_*_*', 
        'keep *_hltSingleMu15L3Filtered15_*_*', 
        'keep *_hltSingleMuIsoL1s14L3IsoFiltered15eta2p1_*_*', 
        'keep *_hltSingleMuIsoL3IsoFiltered15_*_*', 
        'keep *_hltTowerMakerForAll_*_*', 
        'keep *_hltTowerMakerForMuons_*_*', 
        'keep *_hltTriggerSummaryAOD_*_*', 
        'keep *_hltTriggerSummaryRAW_*_*', 
        'keep *_hltTrimmedPixelVertices_*_*', 
        'keep *_hltVerticesL3_*_*', 
        'keep *_hltVerticesPFFilter_*_*', 
        'keep *_hltVerticesPFSelector_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_rawDataRepacker_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_virginRawDataRepacker_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep L2MuonTrajectorySeeds_hltL2MuonSeeds_*_*', 
        'keep L3MuonTrajectorySeeds_hltHIL3TrajSeedOIHit_*_*', 
        'keep L3MuonTrajectorySeeds_hltHIL3TrajectorySeed_*_*', 
        'keep L3MuonTrajectorySeeds_hltL3TrajSeedOIState_*_*', 
        'keep LumiScalerss_hltScalersRawToDigi_*_*', 
        'keep SiPixelClusteredmNewDetSetVector_hltSiPixelClusters_*_*', 
        'keep TrackCandidates_hltHIL3TrackCandidateFromL2OIHit_*_*', 
        'keep TrackCandidates_hltHIL3TrackCandidateFromL2OIState_*_*', 
        'keep TrackingRecHitsOwned_hltL3Muons_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep recoCaloJets_*_*_*', 
        'keep recoCaloMETs_*_*_*', 
        'keep recoCaloMETs_hltMet_*_*', 
        'keep recoCompositeCandidates_*_*_*', 
        'keep recoElectrons_*_*_*', 
        'keep recoIsolatedPixelTrackCandidates_*_*_*', 
        'keep recoMETs_*_*_*', 
        'keep recoPFJets_*_*_*', 
        'keep recoPFTaus_*_*_*', 
        'keep recoRecoChargedCandidates_*_*_*', 
        'keep recoRecoChargedCandidates_hltHIL3MuonCandidates_*_*', 
        'keep recoRecoChargedCandidates_hltL2MuonCandidates_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsoHLTClusterShape_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonEcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalForHE_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsoHLTClusterShape_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonEcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalForHE_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalIsol_*_*', 
        'keep recoRecoEcalCandidates_*_*_*', 
        'keep recoRecoEcalCandidates_hltL1IsoRecoEcalCandidate_*_*', 
        'keep recoRecoEcalCandidates_hltL1NonIsoRecoEcalCandidate_*_*', 
        'keep recoTrackExtras_hltHIL3MuonsOIHit_*_*', 
        'keep recoTrackExtras_hltHIL3MuonsOIState_*_*', 
        'keep recoTracks_hltHIL3MuonsOIHit_*_*', 
        'keep recoTracks_hltHIL3MuonsOIState_*_*', 
        'keep recoTracks_hltHIL3Muons_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2OIHit_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2OIState_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2_*_*', 
        'keep triggerTriggerEventWithRefs_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep triggerTriggerFilterObjectWithRefs_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*' ) )
)

process.HLTDebugRAW = cms.PSet(
    outputCommands = cms.vstring( ('drop *_hlt*_*_*', 
        'keep *_hltAK4CaloJetsCorrectedIDPassed_*_*', 
        'keep *_hltAK4CaloJetsIDPassed_*_*', 
        'keep *_hltAK4CaloJets_*_*', 
        'keep *_hltAK4PFJetsCorrected_*_*', 
        'keep *_hltAK4PFJetsForTaus_*_*', 
        'keep *_hltAK4PFJets_*_*', 
        'keep *_hltAlCaEtaEBRechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaEtaEBRechitsToDigis_*_*', 
        'keep *_hltAlCaEtaEERechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaEtaEERechitsToDigis_*_*', 
        'keep *_hltAlCaEtaRecHitsFilterEBonly_*_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegionalLowPU_etaEcalRecHitsES_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegional_etaEcalRecHitsES_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonly_*_*', 
        'keep *_hltAlCaEtaRecHitsFilter_*_*', 
        'keep *_hltAlCaPhiSymStream_*_*', 
        'keep *_hltAlCaPi0EBRechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaPi0EBRechitsToDigis_*_*', 
        'keep *_hltAlCaPi0EERechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaPi0EERechitsToDigis_*_*', 
        'keep *_hltAlCaPi0RecHitsFilterEBonly_*_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegionalLowPU_pi0EcalRecHitsES_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegional_pi0EcalRecHitsES_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonly_*_*', 
        'keep *_hltAlCaPi0RecHitsFilter_*_*', 
        'keep *_hltBLifetimeL25AssociatorbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL25BJetTagsbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL25TagInfosbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3AssociatorbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3BJetTagsbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3TagInfosbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeRegionalCtfWithMaterialTracksbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBSoftMuonDiJet110Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet110Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet20Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet20Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet40Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet40Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet70Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet70Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonMu5L3_*_*', 
        'keep *_hltCSVJetTagSingleTopEle27_*_*', 
        'keep *_hltCSVJetTagSingleTopIsoMu24_*_*', 
        'keep *_hltCaloJetCorrectedRegional_*_*', 
        'keep *_hltCaloJetCorrected_*_*', 
        'keep *_hltCaloJetL1FastJetCorrected_*_*', 
        'keep *_hltCaloStage2Digis_*_*', 
        'keep *_hltCleanedCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCleanedHiCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCombinedSecondaryVertexBJetTagsCalo_*_*', 
        'keep *_hltCombinedSecondaryVertexBJetTagsPF_*_*', 
        'keep *_hltConvPFTausTightIsoTrackFindingIsolation_*_*', 
        'keep *_hltConvPFTausTightIsoTrackFinding_*_*', 
        'keep *_hltConvPFTausTightIsoTrackPt5Isolation_*_*', 
        'keep *_hltConvPFTausTightIsoTrackPt5_*_*', 
        'keep *_hltConvPFTausTightIso_*_*', 
        'keep *_hltConvPFTausTrackFindingLooseIsolation_*_*', 
        'keep *_hltConvPFTausTrackFinding_*_*', 
        'keep *_hltConvPFTaus_*_*', 
        'keep *_hltCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCorrectedIslandEndcapSuperClustersHI_*_*', 
        'keep *_hltCsc2DRecHits_*_*', 
        'keep *_hltCscSegments_*_*', 
        'keep *_hltDisplacedhltIter4PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolationProng4L1HLTMatched_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolationProng4_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolation_*_*', 
        'keep *_hltDoublePFTau25TrackPt5_*_*', 
        'keep *_hltDoublePFTau25_*_*', 
        'keep *_hltDoublePFTauTightIso45Track5_*_*', 
        'keep *_hltDoublePFTauTightIso45Track_*_*', 
        'keep *_hltDt4DSegments_*_*', 
        'keep *_hltEcalPhiSymFilter_*_*', 
        'keep *_hltEcalRecHitAll_*_*', 
        'keep *_hltEcalRecHit_*_*', 
        'keep *_hltEgammaCandidates_*_*', 
        'keep *_hltEgammaGsfElectrons_*_*', 
        'keep *_hltEgammaGsfTracks_*_*', 
        'keep *_hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1SingleEG18orL1SingleEG20_*_*', 
        'keep *_hltEle20CaloIdVTTrkIdTDphiFilter_*_*', 
        'keep *_hltEle27WP85PixelMatchElectronsL1Seeded_*_*', 
        'keep *_hltElectronsVertex_*_*', 
        'keep *_hltFEDSelectorLumiPixels_*_*', 
        'keep *_hltFastPVPixelTracksMerger_*_*', 
        'keep *_hltFastPVPixelTracksRecover_*_*', 
        'keep *_hltFastPVPixelTracks_*_*', 
        'keep *_hltFastPVPixelVertices3D_*_*', 
        'keep *_hltFastPVPixelVertices_*_*', 
        'keep *_hltFastPixelBLifetimeL3TagInfos_*_*', 
        'keep *_hltFastPrimaryVertex_*_*', 
        'keep *_hltFilterDoubleIsoPFTau45Trk5LeadTrack5IsolationL1HLTMatched_*_*', 
        'keep *_hltFilterL2EtCutDoublePFIsoTau45Trk5_*_*', 
        'keep *_hltFilterL2EtCutSingleIsoPFTau35Trk20MET70_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20LeadTrackPt20_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20MET60LeadTrack20IsolationL1HLTMatched_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20MET70LeadTrack20IsolationL1HLTMatched_*_*', 
        'keep *_hltGmtStage2Digis_*_*', 
        'keep *_hltGoodOnlinePVs_*_*', 
        'keep *_hltGtStage2Digis_*_*', 
        'keep *_hltHICaloJetCorrected_*_*', 
        'keep *_hltHICaloJetIDPassed_*_*', 
        'keep *_hltHIGoodLooseTracks_*_*', 
        'keep *_hltHIPixel3PrimTracks_*_*', 
        'keep *_hltHISelectedVertex_*_*', 
        'keep *_hltHISiPixelClusters_*_*', 
        'keep *_hltHITIPTCorrectorHB_*_*', 
        'keep *_hltHITIPTCorrectorHE_*_*', 
        'keep *_hltHbhereco_*_*', 
        'keep *_hltHfreco_*_*', 
        'keep *_hltHiCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltHiCorrectedIslandEndcapSuperClustersHI_*_*', 
        'keep *_hltHiIslandSuperClustersHI_*_*', 
        'keep *_hltHoreco_*_*', 
        'keep *_hltIsolPixelTrackProdHB_*_*', 
        'keep *_hltIsolPixelTrackProdHE_*_*', 
        'keep *_hltIter0ElectronsTrackSelectionHighPurity_*_*', 
        'keep *_hltIter0HighPtTkMuPixelTracks_*_*', 
        'keep *_hltIter0HighPtTkMuTrackSelectionHighPurity_*_*', 
        'keep *_hltIter0PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter1Merged_*_*', 
        'keep *_hltIter1PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter2HighPtTkMuMerged_*_*', 
        'keep *_hltIter2HighPtTkMuTrackSelectionHighPurity_*_*', 
        'keep *_hltIter2MergedForBTag_*_*', 
        'keep *_hltIter2MergedForElectrons_*_*', 
        'keep *_hltIter2Merged_*_*', 
        'keep *_hltIter2PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter3Merged_*_*', 
        'keep *_hltIter4Merged_*_*', 
        'keep *_hltIterativeCone5PileupSubtractionCaloJets_*_*', 
        'keep *_hltL1HLTSingleIsoPFTau35Trk20Met60JetsMatch_*_*', 
        'keep *_hltL1IsoElectronTrackIsol_*_*', 
        'keep *_hltL1NonIsoElectronTrackIsol_*_*', 
        'keep *_hltL1SeededRecoEcalCandidate_*_*', 
        'keep *_hltL1sDoubleTauJet44Eta2p17orDoubleJet64Central_*_*', 
        'keep *_hltL1sDoubleTauJet44erorDoubleJetC64_*_*', 
        'keep *_hltL1sL1EG18er_*_*', 
        'keep *_hltL1sL1ETM36ORETM40_*_*', 
        'keep *_hltL1sL1Jet52ETM30_*_*', 
        'keep *_hltL1sL1SingleEG12_*_*', 
        'keep *_hltL1sL1SingleEG15_*_*', 
        'keep *_hltL1sL1SingleEG18orL1SingleEG20_*_*', 
        'keep *_hltL1sL1SingleMu10_*_*', 
        'keep *_hltL1sL1SingleMu14Eta2p1_*_*', 
        'keep *_hltL1sMu16Eta2p1_*_*', 
        'keep *_hltL2MuonCandidatesNoVtx_*_*', 
        'keep *_hltL2MuonCandidates_*_*', 
        'keep *_hltL2MuonSeeds_*_*', 
        'keep *_hltL2Muons_*_*', 
        'keep *_hltL2TauJets_*_*', 
        'keep *_hltL3MuonCandidates_*_*', 
        'keep *_hltL3MuonsIOHit_*_*', 
        'keep *_hltL3MuonsLinksCombination_*_*', 
        'keep *_hltL3MuonsOIHit_*_*', 
        'keep *_hltL3MuonsOIState_*_*', 
        'keep *_hltL3Muons_*_*', 
        'keep *_hltL3NoFiltersNoVtxMuonCandidates_*_*', 
        'keep *_hltL3NoFiltersNoVtxMuons_*_*', 
        'keep *_hltL3SecondaryVertexTagInfos_*_*', 
        'keep *_hltL3TkFromL2OICombination_*_*', 
        'keep *_hltL3TkTracksFromL2IOHit_*_*', 
        'keep *_hltL3TkTracksFromL2OIHit_*_*', 
        'keep *_hltL3TkTracksFromL2OIState_*_*', 
        'keep *_hltL3TkTracksFromL2_*_*', 
        'keep *_hltL3TrackCandidateFromL2IOHit_*_*', 
        'keep *_hltL3TrackCandidateFromL2OIHit_*_*', 
        'keep *_hltL3TrackCandidateFromL2OIState_*_*', 
        'keep *_hltL3TrajSeedIOHit_*_*', 
        'keep *_hltL3TrajSeedOIHit_*_*', 
        'keep *_hltL3TrajSeedOIState_*_*', 
        'keep *_hltL3TrajectorySeed_*_*', 
        'keep *_hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f18QL3crIsoRhoFiltered0p15_*_*', 
        'keep *_hltLeadingCentralJets30SingleTopEle27_*_*', 
        'keep *_hltLeadingCentralJets30SingleTopIsoMu24_*_*', 
        'keep *_hltMergedTracks_*_*', 
        'keep *_hltMet_*_*', 
        'keep *_hltMu8Ele17CaloIdTCaloIsoVLPixelMatchFilter_*_*', 
        'keep *_hltMuTrackJpsiCtfTrackCands_*_*', 
        'keep *_hltMuTrackJpsiPixelTrackCands_*_*', 
        'keep *_hltMuonCSCDigis_*_*', 
        'keep *_hltMuonCSCDigis_MuonCSCStripDigi_*', 
        'keep *_hltMuonCSCDigis_MuonCSCWireDigi_*', 
        'keep *_hltMuonDTDigis_*_*', 
        'keep *_hltMuonRPCDigis_*_*', 
        'keep *_hltOnlineBeamSpot_*_*', 
        'keep *_hltOnlinePrimaryVertices_*_*', 
        'keep *_hltOverlapFilterEle20LooseIsoPFTau20OldVersion_*_*', 
        'keep *_hltOverlapFilterIsoEle20MediumIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15IsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15MediumIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15TightIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu18LooseIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu18PFTau25TrackPt5Prong4_*_*', 
        'keep *_hltPFJetForBtag_*_*', 
        'keep *_hltPFMuonMerging_*_*', 
        'keep *_hltPFTau15TrackLooseIso_*_*', 
        'keep *_hltPFTau15Track_*_*', 
        'keep *_hltPFTau15_*_*', 
        'keep *_hltPFTau20IsoMuVertex_*_*', 
        'keep *_hltPFTau20TrackLooseIso_*_*', 
        'keep *_hltPFTau20Track_*_*', 
        'keep *_hltPFTau20_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolationProng4IsoMuVertex_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolationProng4_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolation_*_*', 
        'keep *_hltPFTau25TrackPt5_*_*', 
        'keep *_hltPFTau25_*_*', 
        'keep *_hltPFTau35TrackPt20LooseIsoProng2_*_*', 
        'keep *_hltPFTau35TrackPt20LooseIso_*_*', 
        'keep *_hltPFTau35TrackPt20_*_*', 
        'keep *_hltPFTau35Track_*_*', 
        'keep *_hltPFTau35_*_*', 
        'keep *_hltPFTauEleVertex20_*_*', 
        'keep *_hltPFTauJetTracksAssociator_*_*', 
        'keep *_hltPFTauMediumIso20TrackMediumIso_*_*', 
        'keep *_hltPFTauMediumIso20Track_*_*', 
        'keep *_hltPFTauMediumIso20_*_*', 
        'keep *_hltPFTauMediumIso35Track_*_*', 
        'keep *_hltPFTauMediumIso35_*_*', 
        'keep *_hltPFTauTagInfo_*_*', 
        'keep *_hltPFTauTightIso20TrackTightIso_*_*', 
        'keep *_hltPFTauTightIso20Track_*_*', 
        'keep *_hltPFTauTightIso20_*_*', 
        'keep *_hltPFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltParticleFlowForTaus_*_*', 
        'keep *_hltParticleFlow_*_*', 
        'keep *_hltPixelMatch3HitElectronsActivity_*_*', 
        'keep *_hltPixelMatch3HitElectronsL1Seeded_*_*', 
        'keep *_hltPixelMatchCleanElectronsL1Seeded_*_*', 
        'keep *_hltPixelMatchElectronsActivity_*_*', 
        'keep *_hltPixelMatchElectronsL1Iso_*_*', 
        'keep *_hltPixelMatchElectronsL1NonIso_*_*', 
        'keep *_hltPixelMatchElectronsL1Seeded_*_*', 
        'keep *_hltPixelTracksElectrons_*_*', 
        'keep *_hltPixelTracksMerged_*_*', 
        'keep *_hltPixelTracks_*_*', 
        'keep *_hltPixelVertices3DbbPhi_*_*', 
        'keep *_hltPixelVertices_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidateSC4_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidateSC5_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidate_*_*', 
        'keep *_hltRpcRecHits_*_*', 
        'keep *_hltSelector4CentralJetsL1FastJet_*_*', 
        'keep *_hltSelector8CentralJetsL1FastJet_*_*', 
        'keep *_hltSelectorJets20L1FastJet_*_*', 
        'keep *_hltSiPixelCluster_*_*', 
        'keep *_hltSiPixelClusters_*_*', 
        'keep *_hltSiStripClusters_*_*', 
        'keep *_hltSiStripRawToClustersFacility_*_*', 
        'keep *_hltSingleMu15L3Filtered15_*_*', 
        'keep *_hltSingleMuIsoL1s14L3IsoFiltered15eta2p1_*_*', 
        'keep *_hltSingleMuIsoL3IsoFiltered15_*_*', 
        'keep *_hltTowerMakerForAll_*_*', 
        'keep *_hltTowerMakerForMuons_*_*', 
        'keep *_hltTriggerSummaryAOD_*_*', 
        'keep *_hltTriggerSummaryRAW_*_*', 
        'keep *_hltTrimmedPixelVertices_*_*', 
        'keep *_hltVerticesL3_*_*', 
        'keep *_hltVerticesPFFilter_*_*', 
        'keep *_hltVerticesPFSelector_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_rawDataRepacker_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_virginRawDataRepacker_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep L2MuonTrajectorySeeds_hltL2MuonSeeds_*_*', 
        'keep L3MuonTrajectorySeeds_hltHIL3TrajSeedOIHit_*_*', 
        'keep L3MuonTrajectorySeeds_hltHIL3TrajectorySeed_*_*', 
        'keep L3MuonTrajectorySeeds_hltL3TrajSeedOIState_*_*', 
        'keep LumiScalerss_hltScalersRawToDigi_*_*', 
        'keep SiPixelClusteredmNewDetSetVector_hltSiPixelClusters_*_*', 
        'keep TrackCandidates_hltHIL3TrackCandidateFromL2OIHit_*_*', 
        'keep TrackCandidates_hltHIL3TrackCandidateFromL2OIState_*_*', 
        'keep TrackingRecHitsOwned_hltL3Muons_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep recoCaloJets_*_*_*', 
        'keep recoCaloMETs_*_*_*', 
        'keep recoCaloMETs_hltMet_*_*', 
        'keep recoCompositeCandidates_*_*_*', 
        'keep recoElectrons_*_*_*', 
        'keep recoIsolatedPixelTrackCandidates_*_*_*', 
        'keep recoMETs_*_*_*', 
        'keep recoPFJets_*_*_*', 
        'keep recoPFTaus_*_*_*', 
        'keep recoRecoChargedCandidates_*_*_*', 
        'keep recoRecoChargedCandidates_hltHIL3MuonCandidates_*_*', 
        'keep recoRecoChargedCandidates_hltL2MuonCandidates_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsoHLTClusterShape_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonEcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalForHE_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsoHLTClusterShape_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonEcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalForHE_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalIsol_*_*', 
        'keep recoRecoEcalCandidates_*_*_*', 
        'keep recoRecoEcalCandidates_hltL1IsoRecoEcalCandidate_*_*', 
        'keep recoRecoEcalCandidates_hltL1NonIsoRecoEcalCandidate_*_*', 
        'keep recoTrackExtras_hltHIL3MuonsOIHit_*_*', 
        'keep recoTrackExtras_hltHIL3MuonsOIState_*_*', 
        'keep recoTracks_hltHIL3MuonsOIHit_*_*', 
        'keep recoTracks_hltHIL3MuonsOIState_*_*', 
        'keep recoTracks_hltHIL3Muons_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2OIHit_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2OIState_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2_*_*', 
        'keep triggerTriggerEventWithRefs_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep triggerTriggerFilterObjectWithRefs_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*' ) )
)

process.HLTScouting = cms.PSet(
    outputCommands = cms.vstring('keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*')
)

process.HLTriggerAOD = cms.PSet(
    outputCommands = cms.vstring('drop *_hlt*_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*')
)

process.HLTriggerRAW = cms.PSet(
    outputCommands = cms.vstring('drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*')
)

process.HLTriggerRECO = cms.PSet(
    outputCommands = cms.vstring('drop *_hlt*_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*')
)

process.IOMCRAW = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_randomEngineStateProducer_*_*')
)

process.L1TriggerAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep LumiSummary_lumiProducer_*_*')
)

process.L1TriggerFEVTDEBUG = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_simCscTriggerPrimitiveDigis_*_*', 
        'keep *_simDtTriggerPrimitiveDigis_*_*', 
        'keep *_simRpcTriggerDigis_*_*', 
        'keep *_simRctDigis_*_*', 
        'keep *_simCsctfDigis_*_*', 
        'keep *_simCsctfTrackDigis_*_*', 
        'keep *_simDttfDigis_*_*', 
        'keep *_simGctDigis_*_*', 
        'keep *_simCaloStage1Digis_*_*', 
        'keep *_simCaloStage1FinalDigis_*_*', 
        'keep *_simCaloStage2Layer1Digis_*_*', 
        'keep *_simCaloStage2Digis_*_*', 
        'keep *_simGmtDigis_*_*', 
        'keep *_simBmtfDigis_*_*', 
        'keep *_simOmtfDigis_*_*', 
        'keep *_simEmtfDigis_*_*', 
        'keep *_simGmtStage2Digis_*_*', 
        'keep *_simGtDigis_*_*', 
        'keep *_simGtStage2Digis_*_*', 
        'keep *_cscTriggerPrimitiveDigis_*_*', 
        'keep *_dtTriggerPrimitiveDigis_*_*', 
        'keep *_rpcTriggerDigis_*_*', 
        'keep *_rctDigis_*_*', 
        'keep *_csctfDigis_*_*', 
        'keep *_csctfTrackDigis_*_*', 
        'keep *_dttfDigis_*_*', 
        'keep *_gctDigis_*_*', 
        'keep *_gmtDigis_*_*', 
        'keep *_gtDigis_*_*', 
        'keep *_gtEvmDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep LumiDetails_lumiProducer_*_*', 
        'keep LumiSummary_lumiProducer_*_*')
)

process.L1TriggerRAW = cms.PSet(
    outputCommands = cms.untracked.vstring('keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*')
)

process.L1TriggerRAWDEBUG = cms.PSet(
    outputCommands = cms.untracked.vstring('keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*')
)

process.L1TriggerRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep L1MuGMTReadoutCollection_gtDigis_*_*', 
        'keep L1GctEmCand*_gctDigis_*_*', 
        'keep L1GctJetCand*_gctDigis_*_*', 
        'keep L1GctEtHad*_gctDigis_*_*', 
        'keep L1GctEtMiss*_gctDigis_*_*', 
        'keep L1GctEtTotal*_gctDigis_*_*', 
        'keep L1GctHtMiss*_gctDigis_*_*', 
        'keep L1GctJetCounts*_gctDigis_*_*', 
        'keep L1GctHFRingEtSums*_gctDigis_*_*', 
        'keep L1GctHFBitCounts*_gctDigis_*_*', 
        'keep LumiDetails_lumiProducer_*_*', 
        'keep LumiSummary_lumiProducer_*_*')
)

process.LHEEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep *_externalLHEProducer_LHEScriptOutput_*'),
    splitLevel = cms.untracked.int32(0)
)

process.METSignificance_params = cms.PSet(
    EB_EtResPar = cms.vdouble(0.2, 0.03, 0.005),
    EB_PhiResPar = cms.vdouble(0.00502),
    EE_EtResPar = cms.vdouble(0.2, 0.03, 0.005),
    EE_PhiResPar = cms.vdouble(0.02511),
    HB_EtResPar = cms.vdouble(0.0, 1.22, 0.05),
    HB_PhiResPar = cms.vdouble(0.02511),
    HE_EtResPar = cms.vdouble(0.0, 1.3, 0.05),
    HE_PhiResPar = cms.vdouble(0.02511),
    HF_EtResPar = cms.vdouble(0.0, 1.82, 0.09),
    HF_PhiResPar = cms.vdouble(0.05022),
    HO_EtResPar = cms.vdouble(0.0, 1.3, 0.005),
    HO_PhiResPar = cms.vdouble(0.02511),
    PF_EtResType1 = cms.vdouble(0.05, 0, 0),
    PF_EtResType2 = cms.vdouble(0.05, 0, 0),
    PF_EtResType3 = cms.vdouble(0.05, 0, 0),
    PF_EtResType4 = cms.vdouble(0.042, 0.1, 0.0),
    PF_EtResType5 = cms.vdouble(0.41, 0.52, 0.25),
    PF_EtResType6 = cms.vdouble(0.0, 1.22, 0.05),
    PF_EtResType7 = cms.vdouble(0.0, 1.22, 0.05),
    PF_PhiResType1 = cms.vdouble(0.002),
    PF_PhiResType2 = cms.vdouble(0.002),
    PF_PhiResType3 = cms.vdouble(0.002),
    PF_PhiResType4 = cms.vdouble(0.0028, 0.0, 0.0022),
    PF_PhiResType5 = cms.vdouble(0.1, 0.1, 0.13),
    PF_PhiResType6 = cms.vdouble(0.02511),
    PF_PhiResType7 = cms.vdouble(0.02511),
    jdphi0 = cms.vdouble(0.034, 0.034, 0.034, 0.034, 0.032, 
        0.031, 0.028, 0.027, 0.027, 0.027),
    jdphi1 = cms.vdouble(0.034, 0.035, 0.035, 0.035, 0.035, 
        0.034, 0.031, 0.03, 0.029, 0.027),
    jdphi2 = cms.vdouble(0.04, 0.04, 0.04, 0.04, 0.04, 
        0.038, 0.036, 0.035, 0.034, 0.033),
    jdphi3 = cms.vdouble(0.042, 0.043, 0.044, 0.043, 0.041, 
        0.039, 0.039, 0.036, 0.034, 0.031),
    jdphi4 = cms.vdouble(0.042, 0.042, 0.043, 0.042, 0.038, 
        0.036, 0.036, 0.033, 0.031, 0.031),
    jdphi5 = cms.vdouble(0.069, 0.069, 0.064, 0.058, 0.053, 
        0.049, 0.049, 0.043, 0.039, 0.04),
    jdphi6 = cms.vdouble(0.084, 0.08, 0.072, 0.065, 0.066, 
        0.06, 0.051, 0.049, 0.045, 0.045),
    jdphi7 = cms.vdouble(0.077, 0.072, 0.059, 0.05, 0.045, 
        0.042, 0.039, 0.039, 0.037, 0.031),
    jdphi8 = cms.vdouble(0.059, 0.057, 0.051, 0.044, 0.038, 
        0.035, 0.037, 0.032, 0.028, 0.028),
    jdphi9 = cms.vdouble(0.062, 0.059, 0.053, 0.047, 0.042, 
        0.045, 0.036, 0.032, 0.034, 0.044),
    jdpt0 = cms.vdouble(0.749, 0.829, 1.099, 1.355, 1.584, 
        1.807, 2.035, 2.217, 2.378, 2.591),
    jdpt1 = cms.vdouble(0.718, 0.813, 1.133, 1.384, 1.588, 
        1.841, 2.115, 2.379, 2.508, 2.772),
    jdpt2 = cms.vdouble(0.841, 0.937, 1.316, 1.605, 1.919, 
        2.295, 2.562, 2.722, 2.943, 3.293),
    jdpt3 = cms.vdouble(0.929, 1.04, 1.46, 1.74, 2.042, 
        2.289, 2.639, 2.837, 2.946, 2.971),
    jdpt4 = cms.vdouble(0.85, 0.961, 1.337, 1.593, 1.854, 
        2.005, 2.209, 2.533, 2.812, 3.047),
    jdpt5 = cms.vdouble(1.049, 1.149, 1.607, 1.869, 2.012, 
        2.219, 2.289, 2.412, 2.695, 2.865),
    jdpt6 = cms.vdouble(1.213, 1.298, 1.716, 2.015, 2.191, 
        2.612, 2.863, 2.879, 2.925, 2.902),
    jdpt7 = cms.vdouble(1.094, 1.139, 1.436, 1.672, 1.831, 
        2.05, 2.267, 2.549, 2.785, 2.86),
    jdpt8 = cms.vdouble(0.889, 0.939, 1.166, 1.365, 1.553, 
        1.805, 2.06, 2.22, 2.268, 2.247),
    jdpt9 = cms.vdouble(0.843, 0.885, 1.245, 1.665, 1.944, 
        1.981, 1.972, 2.875, 3.923, 7.51),
    ptresolthreshold = cms.double(10.0),
    resolutionsAlgo = cms.string('AK5PF'),
    resolutionsEra = cms.string('Spring10')
)

process.MEtoEDMConverterAOD = cms.PSet(
    outputCommands = cms.untracked.vstring()
)

process.MEtoEDMConverterFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_MEtoEDMConverter_*_*')
)

process.MEtoEDMConverterRECO = cms.PSet(
    outputCommands = cms.untracked.vstring()
)

process.MINIAODEventContent = cms.PSet(
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    outputCommands = cms.untracked.vstring('drop *', 
        'drop *', 
        'keep *_slimmedPhotons_*_*', 
        'keep *_slimmedOOTPhotons_*_*', 
        'keep *_slimmedElectrons_*_*', 
        'keep *_slimmedMuons_*_*', 
        'keep *_slimmedTaus_*_*', 
        'keep *_slimmedTausBoosted_*_*', 
        'keep *_slimmedCaloJets_*_*', 
        'keep *_slimmedJets_*_*', 
        'drop recoBaseTagInfosOwned_slimmedJets_*_*', 
        'keep *_slimmedJetsAK8_*_*', 
        'keep *_slimmedJetsPuppi_*_*', 
        'keep *_slimmedMETs_*_*', 
        'keep *_slimmedMETsNoHF_*_*', 
        'keep *_slimmedMETsPuppi_*_*', 
        'keep *_slimmedSecondaryVertices_*_*', 
        'keep *_slimmedLambdaVertices_*_*', 
        'keep *_slimmedKshortVertices_*_*', 
        'keep *_slimmedJetsAK8PFPuppiSoftDropPacked_SubJets_*', 
        'keep recoPhotonCores_reducedEgamma_*_*', 
        'keep recoGsfElectronCores_reducedEgamma_*_*', 
        'keep recoConversions_reducedEgamma_*_*', 
        'keep recoSuperClusters_reducedEgamma_*_*', 
        'keep recoCaloClusters_reducedEgamma_*_*', 
        'keep EcalRecHitsSorted_reducedEgamma_*_*', 
        'keep recoGsfTracks_reducedEgamma_*_*', 
        'drop *_*_caloTowers_*', 
        'drop *_*_pfCandidates_*', 
        'drop *_*_genJets_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep *_offlineSlimmedPrimaryVertices_*_*', 
        'keep patPackedCandidates_packedPFCandidates_*_*', 
        'keep *_isolatedTracks_*_*', 
        'keep *_oniaPhotonCandidates_*_*', 
        'keep *_bunchSpacingProducer_*_*', 
        'keep double_fixedGridRhoAll__*', 
        'keep double_fixedGridRhoFastjetAll__*', 
        'keep double_fixedGridRhoFastjetAllCalo__*', 
        'keep double_fixedGridRhoFastjetCentral_*_*', 
        'keep double_fixedGridRhoFastjetCentralCalo__*', 
        'keep double_fixedGridRhoFastjetCentralChargedPileUp__*', 
        'keep double_fixedGridRhoFastjetCentralNeutral__*', 
        'keep *_slimmedPatTrigger_*_*', 
        'keep patPackedTriggerPrescales_patTrigger__*', 
        'keep patPackedTriggerPrescales_patTrigger_l1max_*', 
        'keep patPackedTriggerPrescales_patTrigger_l1min_*', 
        'keep *_l1extraParticles_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_gtStage2Digis__*', 
        'keep *_gmtStage2Digis_Muon_*', 
        'keep *_caloStage2Digis_Jet_*', 
        'keep *_caloStage2Digis_Tau_*', 
        'keep *_caloStage2Digis_EGamma_*', 
        'keep *_caloStage2Digis_EtSum_*', 
        'keep *_TriggerResults_*_HLT', 
        'keep *_TriggerResults_*_*', 
        'keep patPackedCandidates_lostTracks_*_*', 
        'keep HcalNoiseSummary_hcalnoise__*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep CTPPSLocalTrackLites_ctppsLocalTrackLiteProducer_*_*')
)

process.MINIAODSIMEventContent = cms.PSet(
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    outputCommands = cms.untracked.vstring('drop *', 
        'drop *', 
        'keep *_slimmedPhotons_*_*', 
        'keep *_slimmedOOTPhotons_*_*', 
        'keep *_slimmedElectrons_*_*', 
        'keep *_slimmedMuons_*_*', 
        'keep *_slimmedTaus_*_*', 
        'keep *_slimmedTausBoosted_*_*', 
        'keep *_slimmedCaloJets_*_*', 
        'keep *_slimmedJets_*_*', 
        'drop recoBaseTagInfosOwned_slimmedJets_*_*', 
        'keep *_slimmedJetsAK8_*_*', 
        'keep *_slimmedJetsPuppi_*_*', 
        'keep *_slimmedMETs_*_*', 
        'keep *_slimmedMETsNoHF_*_*', 
        'keep *_slimmedMETsPuppi_*_*', 
        'keep *_slimmedSecondaryVertices_*_*', 
        'keep *_slimmedLambdaVertices_*_*', 
        'keep *_slimmedKshortVertices_*_*', 
        'keep *_slimmedJetsAK8PFPuppiSoftDropPacked_SubJets_*', 
        'keep recoPhotonCores_reducedEgamma_*_*', 
        'keep recoGsfElectronCores_reducedEgamma_*_*', 
        'keep recoConversions_reducedEgamma_*_*', 
        'keep recoSuperClusters_reducedEgamma_*_*', 
        'keep recoCaloClusters_reducedEgamma_*_*', 
        'keep EcalRecHitsSorted_reducedEgamma_*_*', 
        'keep recoGsfTracks_reducedEgamma_*_*', 
        'drop *_*_caloTowers_*', 
        'drop *_*_pfCandidates_*', 
        'drop *_*_genJets_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep *_offlineSlimmedPrimaryVertices_*_*', 
        'keep patPackedCandidates_packedPFCandidates_*_*', 
        'keep *_isolatedTracks_*_*', 
        'keep *_oniaPhotonCandidates_*_*', 
        'keep *_bunchSpacingProducer_*_*', 
        'keep double_fixedGridRhoAll__*', 
        'keep double_fixedGridRhoFastjetAll__*', 
        'keep double_fixedGridRhoFastjetAllCalo__*', 
        'keep double_fixedGridRhoFastjetCentral_*_*', 
        'keep double_fixedGridRhoFastjetCentralCalo__*', 
        'keep double_fixedGridRhoFastjetCentralChargedPileUp__*', 
        'keep double_fixedGridRhoFastjetCentralNeutral__*', 
        'keep *_slimmedPatTrigger_*_*', 
        'keep patPackedTriggerPrescales_patTrigger__*', 
        'keep patPackedTriggerPrescales_patTrigger_l1max_*', 
        'keep patPackedTriggerPrescales_patTrigger_l1min_*', 
        'keep *_l1extraParticles_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_gtStage2Digis__*', 
        'keep *_gmtStage2Digis_Muon_*', 
        'keep *_caloStage2Digis_Jet_*', 
        'keep *_caloStage2Digis_Tau_*', 
        'keep *_caloStage2Digis_EGamma_*', 
        'keep *_caloStage2Digis_EtSum_*', 
        'keep *_TriggerResults_*_HLT', 
        'keep *_TriggerResults_*_*', 
        'keep patPackedCandidates_lostTracks_*_*', 
        'keep HcalNoiseSummary_hcalnoise__*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep CTPPSLocalTrackLites_ctppsLocalTrackLiteProducer_*_*', 
        'keep patPackedGenParticles_packedGenParticles_*_*', 
        'keep recoGenParticles_prunedGenParticles_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_*_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep recoGenParticles_genPUProtons_*_*', 
        'keep *_slimmedGenJetsFlavourInfos_*_*', 
        'keep *_slimmedGenJets__*', 
        'keep *_slimmedGenJetsAK8__*', 
        'keep *_slimmedGenJetsAK8SoftDropSubJets__*', 
        'keep *_genMetTrue_*_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep GenRunInfoProduct_*_*_*', 
        'keep *_genParticles_xyz0_*', 
        'keep *_genParticles_t0_*', 
        'keep PileupSummaryInfos_slimmedAddPileupInfo_*_*', 
        'keep L1GtTriggerMenuLite_l1GtTriggerMenuLite__*')
)

process.MINIGENEventContent = cms.PSet(
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep patPackedGenParticles_packedGenParticles_*_*', 
        'keep recoGenParticles_prunedGenParticles_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_*_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep recoGenParticles_genPUProtons_*_*', 
        'keep *_slimmedGenJetsFlavourInfos_*_*', 
        'keep *_slimmedGenJets__*', 
        'keep *_slimmedGenJetsAK8__*', 
        'keep *_slimmedGenJetsAK8SoftDropSubJets__*', 
        'keep *_genMetTrue_*_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep GenRunInfoProduct_*_*_*', 
        'keep *_genParticles_xyz0_*', 
        'keep *_genParticles_t0_*')
)

process.MIXINGMODULEEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_cfWriter_*_*'),
    splitLevel = cms.untracked.int32(0)
)

process.MicroEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_slimmedPhotons_*_*', 
        'keep *_slimmedOOTPhotons_*_*', 
        'keep *_slimmedElectrons_*_*', 
        'keep *_slimmedMuons_*_*', 
        'keep *_slimmedTaus_*_*', 
        'keep *_slimmedTausBoosted_*_*', 
        'keep *_slimmedCaloJets_*_*', 
        'keep *_slimmedJets_*_*', 
        'drop recoBaseTagInfosOwned_slimmedJets_*_*', 
        'keep *_slimmedJetsAK8_*_*', 
        'keep *_slimmedJetsPuppi_*_*', 
        'keep *_slimmedMETs_*_*', 
        'keep *_slimmedMETsNoHF_*_*', 
        'keep *_slimmedMETsPuppi_*_*', 
        'keep *_slimmedSecondaryVertices_*_*', 
        'keep *_slimmedLambdaVertices_*_*', 
        'keep *_slimmedKshortVertices_*_*', 
        'keep *_slimmedJetsAK8PFPuppiSoftDropPacked_SubJets_*', 
        'keep recoPhotonCores_reducedEgamma_*_*', 
        'keep recoGsfElectronCores_reducedEgamma_*_*', 
        'keep recoConversions_reducedEgamma_*_*', 
        'keep recoSuperClusters_reducedEgamma_*_*', 
        'keep recoCaloClusters_reducedEgamma_*_*', 
        'keep EcalRecHitsSorted_reducedEgamma_*_*', 
        'keep recoGsfTracks_reducedEgamma_*_*', 
        'drop *_*_caloTowers_*', 
        'drop *_*_pfCandidates_*', 
        'drop *_*_genJets_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep *_offlineSlimmedPrimaryVertices_*_*', 
        'keep patPackedCandidates_packedPFCandidates_*_*', 
        'keep *_isolatedTracks_*_*', 
        'keep *_oniaPhotonCandidates_*_*', 
        'keep *_bunchSpacingProducer_*_*', 
        'keep double_fixedGridRhoAll__*', 
        'keep double_fixedGridRhoFastjetAll__*', 
        'keep double_fixedGridRhoFastjetAllCalo__*', 
        'keep double_fixedGridRhoFastjetCentral_*_*', 
        'keep double_fixedGridRhoFastjetCentralCalo__*', 
        'keep double_fixedGridRhoFastjetCentralChargedPileUp__*', 
        'keep double_fixedGridRhoFastjetCentralNeutral__*', 
        'keep *_slimmedPatTrigger_*_*', 
        'keep patPackedTriggerPrescales_patTrigger__*', 
        'keep patPackedTriggerPrescales_patTrigger_l1max_*', 
        'keep patPackedTriggerPrescales_patTrigger_l1min_*', 
        'keep *_l1extraParticles_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_gtStage2Digis__*', 
        'keep *_gmtStage2Digis_Muon_*', 
        'keep *_caloStage2Digis_Jet_*', 
        'keep *_caloStage2Digis_Tau_*', 
        'keep *_caloStage2Digis_EGamma_*', 
        'keep *_caloStage2Digis_EtSum_*', 
        'keep *_TriggerResults_*_HLT', 
        'keep *_TriggerResults_*_*', 
        'keep patPackedCandidates_lostTracks_*_*', 
        'keep HcalNoiseSummary_hcalnoise__*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep CTPPSLocalTrackLites_ctppsLocalTrackLiteProducer_*_*')
)

process.MicroEventContentGEN = cms.PSet(
    outputCommands = cms.untracked.vstring('keep patPackedGenParticles_packedGenParticles_*_*', 
        'keep recoGenParticles_prunedGenParticles_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_*_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep recoGenParticles_genPUProtons_*_*', 
        'keep *_slimmedGenJetsFlavourInfos_*_*', 
        'keep *_slimmedGenJets__*', 
        'keep *_slimmedGenJetsAK8__*', 
        'keep *_slimmedGenJetsAK8SoftDropSubJets__*', 
        'keep *_genMetTrue_*_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep GenRunInfoProduct_*_*_*', 
        'keep *_genParticles_xyz0_*', 
        'keep *_genParticles_t0_*')
)

process.MicroEventContentMC = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_slimmedPhotons_*_*', 
        'keep *_slimmedOOTPhotons_*_*', 
        'keep *_slimmedElectrons_*_*', 
        'keep *_slimmedMuons_*_*', 
        'keep *_slimmedTaus_*_*', 
        'keep *_slimmedTausBoosted_*_*', 
        'keep *_slimmedCaloJets_*_*', 
        'keep *_slimmedJets_*_*', 
        'drop recoBaseTagInfosOwned_slimmedJets_*_*', 
        'keep *_slimmedJetsAK8_*_*', 
        'keep *_slimmedJetsPuppi_*_*', 
        'keep *_slimmedMETs_*_*', 
        'keep *_slimmedMETsNoHF_*_*', 
        'keep *_slimmedMETsPuppi_*_*', 
        'keep *_slimmedSecondaryVertices_*_*', 
        'keep *_slimmedLambdaVertices_*_*', 
        'keep *_slimmedKshortVertices_*_*', 
        'keep *_slimmedJetsAK8PFPuppiSoftDropPacked_SubJets_*', 
        'keep recoPhotonCores_reducedEgamma_*_*', 
        'keep recoGsfElectronCores_reducedEgamma_*_*', 
        'keep recoConversions_reducedEgamma_*_*', 
        'keep recoSuperClusters_reducedEgamma_*_*', 
        'keep recoCaloClusters_reducedEgamma_*_*', 
        'keep EcalRecHitsSorted_reducedEgamma_*_*', 
        'keep recoGsfTracks_reducedEgamma_*_*', 
        'drop *_*_caloTowers_*', 
        'drop *_*_pfCandidates_*', 
        'drop *_*_genJets_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep *_offlineSlimmedPrimaryVertices_*_*', 
        'keep patPackedCandidates_packedPFCandidates_*_*', 
        'keep *_isolatedTracks_*_*', 
        'keep *_oniaPhotonCandidates_*_*', 
        'keep *_bunchSpacingProducer_*_*', 
        'keep double_fixedGridRhoAll__*', 
        'keep double_fixedGridRhoFastjetAll__*', 
        'keep double_fixedGridRhoFastjetAllCalo__*', 
        'keep double_fixedGridRhoFastjetCentral_*_*', 
        'keep double_fixedGridRhoFastjetCentralCalo__*', 
        'keep double_fixedGridRhoFastjetCentralChargedPileUp__*', 
        'keep double_fixedGridRhoFastjetCentralNeutral__*', 
        'keep *_slimmedPatTrigger_*_*', 
        'keep patPackedTriggerPrescales_patTrigger__*', 
        'keep patPackedTriggerPrescales_patTrigger_l1max_*', 
        'keep patPackedTriggerPrescales_patTrigger_l1min_*', 
        'keep *_l1extraParticles_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_gtStage2Digis__*', 
        'keep *_gmtStage2Digis_Muon_*', 
        'keep *_caloStage2Digis_Jet_*', 
        'keep *_caloStage2Digis_Tau_*', 
        'keep *_caloStage2Digis_EGamma_*', 
        'keep *_caloStage2Digis_EtSum_*', 
        'keep *_TriggerResults_*_HLT', 
        'keep *_TriggerResults_*_*', 
        'keep patPackedCandidates_lostTracks_*_*', 
        'keep HcalNoiseSummary_hcalnoise__*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep CTPPSLocalTrackLites_ctppsLocalTrackLiteProducer_*_*', 
        'keep patPackedGenParticles_packedGenParticles_*_*', 
        'keep recoGenParticles_prunedGenParticles_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_*_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep recoGenParticles_genPUProtons_*_*', 
        'keep *_slimmedGenJetsFlavourInfos_*_*', 
        'keep *_slimmedGenJets__*', 
        'keep *_slimmedGenJetsAK8__*', 
        'keep *_slimmedGenJetsAK8SoftDropSubJets__*', 
        'keep *_genMetTrue_*_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep GenRunInfoProduct_*_*_*', 
        'keep *_genParticles_xyz0_*', 
        'keep *_genParticles_t0_*', 
        'keep PileupSummaryInfos_slimmedAddPileupInfo_*_*', 
        'keep L1GtTriggerMenuLite_l1GtTriggerMenuLite__*')
)

process.NANOAODEventContent = cms.PSet(
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep nanoaodFlatTable_*Table_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep nanoaodMergeableCounterTable_*Table_*_*', 
        'keep nanoaodUniqueString_nanoMetadata_*_*')
)

process.NANOAODSIMEventContent = cms.PSet(
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep nanoaodFlatTable_*Table_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep nanoaodMergeableCounterTable_*Table_*_*', 
        'keep nanoaodUniqueString_nanoMetadata_*_*')
)

process.NanoAODEDMEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *', 
        'keep nanoaodFlatTable_*Table_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep nanoaodMergeableCounterTable_*Table_*_*', 
        'keep nanoaodUniqueString_nanoMetadata_*_*')
)

process.PREMIXEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *', 
        'drop *', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep SimTracks_g4SimHits_*_*', 
        'keep SimVertexs_g4SimHits_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep CrossingFramePlaybackInfoNew_*_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int6stdbitsetstdpairs_*_AffectedAPVList_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'keep *_mix_MergedTrackTruth_*', 
        'keep recoGenMETs_*_*_*', 
        'keep recoGenJets_*_*_*', 
        'keep *_genParticle_*_*', 
        'keep *_MEtoEDMConverter_*_*', 
        'keep *_randomEngineStateProducer_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'keep RPCDetIdRPCDigiMuonDigiCollection_simMuonRPCDigis_*_*', 
        'keep *_mix_MergedTrackTruth_*', 
        'keep StripDigiSimLinkedmDetSetVector_simSiStripDigis_*_*', 
        'keep PixelDigiSimLinkedmDetSetVector_simSiPixelDigis_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_*_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_*_*_*'),
    splitLevel = cms.untracked.int32(0)
)

process.PREMIXRAWEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *', 
        'drop *', 
        'drop *', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep *_g4SimHits_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep CSCDetIdCSCComparatorDigiMuonDigiCollection_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*', 
        'keep EBSrFlagsSorted_simEcalDigis_*_*', 
        'keep EESrFlagsSorted_simEcalDigis_*_*', 
        'keep CrossingFramePlaybackInfoNew_*_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int6stdbitsetstdpairs_*_AffectedAPVList_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'keep *_mix_MergedTrackTruth_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep recoGenJets_*_*_*', 
        'keep *_genParticle_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep *_MEtoEDMConverter_*_*', 
        'keep *_randomEngineStateProducer_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'keep CrossingFramePlaybackInfoNew_*_*_*', 
        'drop CrossingFramePlaybackInfoNew_mix_*_*', 
        'keep *_*_MergedTrackTruth_*', 
        'keep *_*_StripDigiSimLink_*', 
        'keep *_*_PixelDigiSimLink_*', 
        'keep *_*_MuonCSCStripDigiSimLinks_*', 
        'keep *_*_MuonCSCWireDigiSimLinks_*', 
        'keep *_*_RPCDigiSimLink_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_*_*_*'),
    splitLevel = cms.untracked.int32(0)
)

process.RAWAODSIMEventContent = cms.PSet(
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
    outputCommands = cms.untracked.vstring( ('drop *', 
        'drop *', 
        'drop *', 
        'keep ClusterSummary_clusterSummaryProducer_*_*', 
        'keep *_dt4DSegments_*_*', 
        'keep *_dt4DCosmicSegments_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*', 
        'keep *_castorreco_*_*', 
        'keep *_reducedHcalRecHits_*_*', 
        'keep HcalUnpackerReport_castorDigis_*_*', 
        'keep HcalUnpackerReport_hcalDigiAlCaMB_*_*', 
        'keep HcalUnpackerReport_hcalDigis_*_*', 
        'keep *_selectDigi_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEB_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEE_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsES_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep recoCaloClusters_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_hybridSuperClusters_uncleanOnlyHybridSuperClusters_*', 
        'keep recoCaloClusters_multi5x5SuperClusters_multi5x5EndcapBasicClusters_*', 
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_*', 
        'keep recoSuperClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoCaloClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoSuperClusters_particleFlowSuperClusterOOTECAL_*_*', 
        'keep recoCaloClusters_particleFlowSuperClusterOOTECAL_*_*', 
        'keep recoTracks_GsfGlobalElectronTest_*_*', 
        'keep recoGsfTracks_electronGsfTracks_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_conversionStepTracks_*_*', 
        'keep recoTracks_beamhaloTracks_*_*', 
        'keep recoTracks_ctfPixelLess_*_*', 
        'keep *_dedxHarmonic2_*_*', 
        'keep *_dedxPixelHarmonic2_*_*', 
        'keep *_dedxHitInfo_*_*', 
        'keep *_trackExtrapolator_*_*', 
        'keep *_generalTracks_MVAValues_*', 
        'keep *_generalTracks_MVAVals_*', 
        'keep *_ak4CaloJets_*_*', 
        'keep *_ak4PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHSSoftDrop_*_*', 
        'keep *_cmsTopTagPFJetsCHS_*_*', 
        'keep *_ak4PFJets_*_*', 
        'keep *_JetPlusTrackZSPCorJetAntiKt4_*_*', 
        'keep *_ak4TrackJets_*_*', 
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', 
        'keep *_caloTowers_*_*', 
        'keep *_CastorTowerReco_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertex_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertexPF_*_*', 
        'keep *_ak4JetTracksAssociatorExplicit_*_*', 
        'keep *_ak4JetExtender_*_*', 
        'keep *_ak4JetID_*_*', 
        'keep *_ak5CastorJets_*_*', 
        'keep *_ak5CastorJetID_*_*', 
        'keep *_ak7CastorJets_*_*', 
        'keep *_ak7CastorJetID_*_*', 
        'keep *_fixedGridRhoAll_*_*', 
        'keep *_fixedGridRhoFastjetAll_*_*', 
        'keep *_fixedGridRhoFastjetAllTmp_*_*', 
        'keep *_fixedGridRhoFastjetCentral_*_*', 
        'keep *_fixedGridRhoFastjetAllCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralChargedPileUp_*_*', 
        'keep *_fixedGridRhoFastjetCentralNeutral_*_*', 
        'drop doubles_*Jets_rhos_*', 
        'drop doubles_*Jets_sigmas_*', 
        'keep *_ak8PFJetsCHSSoftDropMass_*_*', 
        'keep recoCaloMETs_caloMet_*_*', 
        'keep recoCaloMETs_caloMetBE_*_*', 
        'keep recoCaloMETs_caloMetBEFO_*_*', 
        'keep recoCaloMETs_caloMetM_*_*', 
        'keep recoPFMETs_pfMet_*_*', 
        'keep recoPFMETs_pfChMet_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_*_*', 
        'drop recoHcalNoiseRBXs_*_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*', 
        'keep recoGlobalHaloData_GlobalHaloData_*_*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep *_muons_*_*', 
        'keep *_particleFlow_muons_*', 
        'drop *_muons_muons1stStep2muonsMap_*', 
        'drop recoIsoDepositedmValueMap_muons_*_*', 
        'drop doubleedmValueMap_muons_muPFIso*_*', 
        'keep recoTracks_standAloneMuons_*_*', 
        'keep recoTrackExtras_standAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_standAloneMuons_*_*', 
        'keep recoTracks_globalMuons_*_*', 
        'keep recoTrackExtras_globalMuons_*_*', 
        'keep recoTracks_tevMuons_*_*', 
        'keep recoTrackExtras_tevMuons_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_displacedTracks_*_*', 
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_*_*', 
        'keep recoTracks_displacedGlobalMuons_*_*', 
        'keep recoTrackExtras_displacedGlobalMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedGlobalMuons_*_*', 
        'keep recoTracks_cosmicMuons_*_*', 
        'keep recoMuons_muonsFromCosmics_*_*', 
        'keep recoTracks_cosmicMuons1Leg_*_*', 
        'keep recoMuons_muonsFromCosmics1Leg_*_*', 
        'keep recoTracks_refittedStandAloneMuons_*_*', 
        'keep recoTrackExtras_refittedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_refittedStandAloneMuons_*_*', 
        'keep recoTracks_displacedStandAloneMuons__*', 
        'keep recoTrackExtras_displacedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedStandAloneMuons_*_*', 
        'keep *_softPFElectronBJetTags_*_*', 
        'keep *_softPFMuonBJetTags_*_*', 
        'keep *_pfTrackCountingHighEffBJetTags_*_*', 
        'keep *_pfJetProbabilityBJetTags_*_*', 
        'keep *_pfJetBProbabilityBJetTags_*_*', 
        'keep *_pfSimpleSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfCombinedSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfGhostTrackBJetTags_*_*', 
        'keep *_pfCombinedMVAV2BJetTags_*_*', 
        'keep *_inclusiveCandidateSecondaryVertices_*_*', 
        'keep *_inclusiveCandidateSecondaryVerticesCvsL_*_*', 
        'keep *_pfCombinedCvsLJetTags_*_*', 
        'keep *_pfCombinedCvsBJetTags_*_*', 
        'keep *_pfChargeBJetTags_*_*', 
        'keep *_pfDeepCSVJetTags_*_*', 
        'keep *_pfDeepCMVAJetTags_*_*', 
        'keep recoRecoTauPiZeros_hpsPFTauProducer_pizeros_*', 
        'keep recoPFTaus_hpsPFTauProducer_*_*', 
        'keep *_hpsPFTauDiscriminationByDeadECALElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFinding_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingNewDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingOldDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseChargedIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByTightMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSum_*_*', 
        'keep *_hpsPFTauPUcorrPtSum_*_*', 
        'keep *_hpsPFTauChargedIsoPtSum_*_*', 
        'keep *_hpsPFTauTransverseImpactParameters_*_*', 
        'keep *_hpsPFTauFootprintCorrection_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeight_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6rawElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6LooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6MediumElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6TightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VTightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauChargedIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeightdR03_*_*', 
        'keep *_hpsPFTauFootprintCorrectiondR03_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalConedR03_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep  *_offlinePrimaryVertices__*', 
        'keep *_offlinePrimaryVerticesWithBS_*_*', 
        'keep *_offlinePrimaryVerticesFromCosmicTracks_*_*', 
        'keep *_nuclearInteractionMaker_*_*', 
        'keep *_generalV0Candidates_*_*', 
        'keep *_inclusiveSecondaryVertices_*_*', 
        'keep recoGsfElectronCores_gsfElectronCores_*_*', 
        'keep recoGsfElectronCores_gedGsfElectronCores_*_*', 
        'keep recoGsfElectrons_gsfElectrons_*_*', 
        'keep recoGsfElectrons_gedGsfElectrons_*_*', 
        'keep recoGsfElectronCores_uncleanedOnlyGsfElectronCores_*_*', 
        'keep recoGsfElectrons_uncleanedOnlyGsfElectrons_*_*', 
        'keep floatedmValueMap_eidRobustLoose_*_*', 
        'keep floatedmValueMap_eidRobustTight_*_*', 
        'keep floatedmValueMap_eidRobustHighEnergy_*_*', 
        'keep floatedmValueMap_eidLoose_*_*', 
        'keep floatedmValueMap_eidTight_*_*', 
        'keep *_egmGedGsfElectronPFIsolation_*_*', 
        'keep *_photonEcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonEcalPFClusterIsolationProducer_*_*', 
        'keep *_electronEcalPFClusterIsolationProducer_*_*', 
        'keep *_photonHcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonHcalPFClusterIsolationProducer_*_*', 
        'keep *_electronHcalPFClusterIsolationProducer_*_*', 
        'drop *_egmGsfElectronIDs_*_*', 
        'drop *_egmPhotonIDs_*_*', 
        'keep recoPhotonCores_gedPhotonCore_*_*', 
        'keep recoPhotons_gedPhotons_*_*', 
        'keep *_particleBasedIsolation_*_*', 
        'drop *_gedPhotons_valMapPFEgammaCandToPhoton_*', 
        'keep recoPhotonCores_photonCore_*_*', 
        'keep recoPhotons_photons_*_*', 
        'keep recoPhotonCores_ootPhotonCore_*_*', 
        'keep recoPhotons_ootPhotons_*_*', 
        'keep recoConversions_conversions_*_*', 
        'keep recoConversions_mustacheConversions_*_*', 
        'drop *_conversions_uncleanedConversions_*', 
        'keep recoConversions_allConversions_*_*', 
        'keep recoConversions_allConversionsOldEG_*_*', 
        'keep recoTracks_ckfOutInTracksFromConversions_*_*', 
        'keep recoTracks_ckfInOutTracksFromConversions_*_*', 
        'keep recoConversions_uncleanedOnlyAllConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep *_PhotonIDProdGED_*_*', 
        'keep *_hfRecoEcalCandidate_*_*', 
        'keep recoRecoEcalCandidates_hfRecoEcalCandidate_*_*', 
        'keep *_hfEMClusters_*_*', 
        'keep *_gedGsfElectronCores_*_*', 
        'keep *_gedGsfElectrons_*_*', 
        'drop CaloTowersSorted_towerMakerPF_*_*', 
        'drop *_pfElectronTranslator_*_*', 
        'keep recoPFRecHits_particleFlowClusterECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHCAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterPS_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHBHE_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitPS_Cleaned_*', 
        'keep recoCaloClusters_particleFlowEGamma_*_*', 
        'keep recoSuperClusters_particleFlowEGamma_*_*', 
        'keep recoCaloClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoSuperClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoConversions_particleFlowEGamma_*_*', 
        'keep recoPFCandidates_particleFlow_*_*', 
        'keep recoPFCandidates_particleFlowTmp_*_*', 
        'drop recoPFCandidates_particleFlowTmp__*', 
        'keep *_particleFlow_electrons_*', 
        'keep *_particleFlow_photons_*', 
        'keep *_particleFlow_muons_*', 
        'keep recoCaloClusters_pfElectronTranslator_*_*', 
        'keep recoPreshowerClusters_pfElectronTranslator_*_*', 
        'keep recoSuperClusters_pfElectronTranslator_*_*', 
        'keep recoCaloClusters_pfPhotonTranslator_*_*', 
        'keep recoPreshowerClusters_pfPhotonTranslator_*_*', 
        'keep recoSuperClusters_pfPhotonTranslator_*_*', 
        'keep recoPhotons_pfPhotonTranslator_*_*', 
        'keep recoPhotonCores_pfPhotonTranslator_*_*', 
        'keep recoConversions_pfPhotonTranslator_*_*', 
        'keep *_particleFlowPtrs_*_*', 
        'keep *_particleFlowTmpPtrs_*_*', 
        'keep *_chargedHadronPFTrackIsolation_*_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep LumiSummary_lumiProducer_*_*', 
        'drop *_hlt*_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep L1AcceptBunchCrossings_scalersRawToDigi_*_*', 
        'keep L1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep Level1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep BeamSpotOnlines_scalersRawToDigi_*_*', 
        'keep DcsStatuss_scalersRawToDigi_*_*', 
        'keep *_tcdsDigis_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'keep *_pfIsolatedElectronsEI_*_*', 
        'keep *_pfIsolatedMuonsEI_*_*', 
        'keep recoPFJets_pfJetsEI_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTagsEI_*_*', 
        'keep recoPFTaus_pfTausEI_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByDecayModeFinding_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByIsolation_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep *_muonSimClassifier_*_*', 
        'keep *_ak4GenJets_*_*', 
        'keep *_ak8GenJets_*_*', 
        'keep *_ak4GenJetsNoNu_*_*', 
        'keep *_ak8GenJetsNoNu_*_*', 
        'keep *_genParticle_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep SimVertexs_g4SimHits_*_*' ) )
)

process.RAWDEBUGEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *', 
        'drop *', 
        'drop *', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep *_g4SimHits_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep CSCDetIdCSCComparatorDigiMuonDigiCollection_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*', 
        'keep EBSrFlagsSorted_simEcalDigis_*_*', 
        'keep EESrFlagsSorted_simEcalDigis_*_*', 
        'keep CrossingFramePlaybackInfoNew_*_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int6stdbitsetstdpairs_*_AffectedAPVList_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'keep *_mix_MergedTrackTruth_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep recoGenJets_*_*_*', 
        'keep *_genParticle_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep *_MEtoEDMConverter_*_*', 
        'keep *_randomEngineStateProducer_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'keep PixelDigiSimLinkedmDetSetVector_simSiPixelDigis_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simSiStripDigis_*_*', 
        'drop *_mix_simSiPixelDigis*_*', 
        'drop *_mix_simSiStripDigis*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'drop *_trackingtruthprod_*_*', 
        'drop *_electrontruth_*_*', 
        'keep *_mix_MergedTrackTruth_*', 
        'keep CrossingFramePlaybackInfoNew_*_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*'),
    splitLevel = cms.untracked.int32(0)
)

process.RAWDEBUGHLTEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring( ('drop *', 
        'drop *', 
        'drop *', 
        'drop *', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep *_g4SimHits_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep CSCDetIdCSCComparatorDigiMuonDigiCollection_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*', 
        'keep EBSrFlagsSorted_simEcalDigis_*_*', 
        'keep EESrFlagsSorted_simEcalDigis_*_*', 
        'keep CrossingFramePlaybackInfoNew_*_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int6stdbitsetstdpairs_*_AffectedAPVList_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'keep *_mix_MergedTrackTruth_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep recoGenJets_*_*_*', 
        'keep *_genParticle_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep *_MEtoEDMConverter_*_*', 
        'keep *_randomEngineStateProducer_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'keep PixelDigiSimLinkedmDetSetVector_simSiPixelDigis_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simSiStripDigis_*_*', 
        'drop *_mix_simSiPixelDigis*_*', 
        'drop *_mix_simSiStripDigis*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'drop *_trackingtruthprod_*_*', 
        'drop *_electrontruth_*_*', 
        'keep *_mix_MergedTrackTruth_*', 
        'keep CrossingFramePlaybackInfoNew_*_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'drop *_hlt*_*_*', 
        'keep *_hltAK4CaloJetsCorrectedIDPassed_*_*', 
        'keep *_hltAK4CaloJetsIDPassed_*_*', 
        'keep *_hltAK4CaloJets_*_*', 
        'keep *_hltAK4PFJetsCorrected_*_*', 
        'keep *_hltAK4PFJetsForTaus_*_*', 
        'keep *_hltAK4PFJets_*_*', 
        'keep *_hltAlCaEtaEBRechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaEtaEBRechitsToDigis_*_*', 
        'keep *_hltAlCaEtaEERechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaEtaEERechitsToDigis_*_*', 
        'keep *_hltAlCaEtaRecHitsFilterEBonly_*_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegionalLowPU_etaEcalRecHitsES_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegional_etaEcalRecHitsES_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonly_*_*', 
        'keep *_hltAlCaEtaRecHitsFilter_*_*', 
        'keep *_hltAlCaPhiSymStream_*_*', 
        'keep *_hltAlCaPi0EBRechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaPi0EBRechitsToDigis_*_*', 
        'keep *_hltAlCaPi0EERechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaPi0EERechitsToDigis_*_*', 
        'keep *_hltAlCaPi0RecHitsFilterEBonly_*_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegionalLowPU_pi0EcalRecHitsES_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegional_pi0EcalRecHitsES_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonly_*_*', 
        'keep *_hltAlCaPi0RecHitsFilter_*_*', 
        'keep *_hltBLifetimeL25AssociatorbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL25BJetTagsbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL25TagInfosbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3AssociatorbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3BJetTagsbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3TagInfosbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeRegionalCtfWithMaterialTracksbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBSoftMuonDiJet110Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet110Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet20Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet20Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet40Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet40Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet70Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet70Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonMu5L3_*_*', 
        'keep *_hltCSVJetTagSingleTopEle27_*_*', 
        'keep *_hltCSVJetTagSingleTopIsoMu24_*_*', 
        'keep *_hltCaloJetCorrectedRegional_*_*', 
        'keep *_hltCaloJetCorrected_*_*', 
        'keep *_hltCaloJetL1FastJetCorrected_*_*', 
        'keep *_hltCaloStage2Digis_*_*', 
        'keep *_hltCleanedCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCleanedHiCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCombinedSecondaryVertexBJetTagsCalo_*_*', 
        'keep *_hltCombinedSecondaryVertexBJetTagsPF_*_*', 
        'keep *_hltConvPFTausTightIsoTrackFindingIsolation_*_*', 
        'keep *_hltConvPFTausTightIsoTrackFinding_*_*', 
        'keep *_hltConvPFTausTightIsoTrackPt5Isolation_*_*', 
        'keep *_hltConvPFTausTightIsoTrackPt5_*_*', 
        'keep *_hltConvPFTausTightIso_*_*', 
        'keep *_hltConvPFTausTrackFindingLooseIsolation_*_*', 
        'keep *_hltConvPFTausTrackFinding_*_*', 
        'keep *_hltConvPFTaus_*_*', 
        'keep *_hltCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCorrectedIslandEndcapSuperClustersHI_*_*', 
        'keep *_hltCsc2DRecHits_*_*', 
        'keep *_hltCscSegments_*_*', 
        'keep *_hltDisplacedhltIter4PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolationProng4L1HLTMatched_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolationProng4_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolation_*_*', 
        'keep *_hltDoublePFTau25TrackPt5_*_*', 
        'keep *_hltDoublePFTau25_*_*', 
        'keep *_hltDoublePFTauTightIso45Track5_*_*', 
        'keep *_hltDoublePFTauTightIso45Track_*_*', 
        'keep *_hltDt4DSegments_*_*', 
        'keep *_hltEcalPhiSymFilter_*_*', 
        'keep *_hltEcalRecHitAll_*_*', 
        'keep *_hltEcalRecHit_*_*', 
        'keep *_hltEgammaCandidates_*_*', 
        'keep *_hltEgammaGsfElectrons_*_*', 
        'keep *_hltEgammaGsfTracks_*_*', 
        'keep *_hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1SingleEG18orL1SingleEG20_*_*', 
        'keep *_hltEle20CaloIdVTTrkIdTDphiFilter_*_*', 
        'keep *_hltEle27WP85PixelMatchElectronsL1Seeded_*_*', 
        'keep *_hltElectronsVertex_*_*', 
        'keep *_hltFEDSelectorLumiPixels_*_*', 
        'keep *_hltFastPVPixelTracksMerger_*_*', 
        'keep *_hltFastPVPixelTracksRecover_*_*', 
        'keep *_hltFastPVPixelTracks_*_*', 
        'keep *_hltFastPVPixelVertices3D_*_*', 
        'keep *_hltFastPVPixelVertices_*_*', 
        'keep *_hltFastPixelBLifetimeL3TagInfos_*_*', 
        'keep *_hltFastPrimaryVertex_*_*', 
        'keep *_hltFilterDoubleIsoPFTau45Trk5LeadTrack5IsolationL1HLTMatched_*_*', 
        'keep *_hltFilterL2EtCutDoublePFIsoTau45Trk5_*_*', 
        'keep *_hltFilterL2EtCutSingleIsoPFTau35Trk20MET70_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20LeadTrackPt20_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20MET60LeadTrack20IsolationL1HLTMatched_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20MET70LeadTrack20IsolationL1HLTMatched_*_*', 
        'keep *_hltGmtStage2Digis_*_*', 
        'keep *_hltGoodOnlinePVs_*_*', 
        'keep *_hltGtStage2Digis_*_*', 
        'keep *_hltHICaloJetCorrected_*_*', 
        'keep *_hltHICaloJetIDPassed_*_*', 
        'keep *_hltHIGoodLooseTracks_*_*', 
        'keep *_hltHIPixel3PrimTracks_*_*', 
        'keep *_hltHISelectedVertex_*_*', 
        'keep *_hltHISiPixelClusters_*_*', 
        'keep *_hltHITIPTCorrectorHB_*_*', 
        'keep *_hltHITIPTCorrectorHE_*_*', 
        'keep *_hltHbhereco_*_*', 
        'keep *_hltHfreco_*_*', 
        'keep *_hltHiCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltHiCorrectedIslandEndcapSuperClustersHI_*_*', 
        'keep *_hltHiIslandSuperClustersHI_*_*', 
        'keep *_hltHoreco_*_*', 
        'keep *_hltIsolPixelTrackProdHB_*_*', 
        'keep *_hltIsolPixelTrackProdHE_*_*', 
        'keep *_hltIter0ElectronsTrackSelectionHighPurity_*_*', 
        'keep *_hltIter0HighPtTkMuPixelTracks_*_*', 
        'keep *_hltIter0HighPtTkMuTrackSelectionHighPurity_*_*', 
        'keep *_hltIter0PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter1Merged_*_*', 
        'keep *_hltIter1PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter2HighPtTkMuMerged_*_*', 
        'keep *_hltIter2HighPtTkMuTrackSelectionHighPurity_*_*', 
        'keep *_hltIter2MergedForBTag_*_*', 
        'keep *_hltIter2MergedForElectrons_*_*', 
        'keep *_hltIter2Merged_*_*', 
        'keep *_hltIter2PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter3Merged_*_*', 
        'keep *_hltIter4Merged_*_*', 
        'keep *_hltIterativeCone5PileupSubtractionCaloJets_*_*', 
        'keep *_hltL1HLTSingleIsoPFTau35Trk20Met60JetsMatch_*_*', 
        'keep *_hltL1IsoElectronTrackIsol_*_*', 
        'keep *_hltL1NonIsoElectronTrackIsol_*_*', 
        'keep *_hltL1SeededRecoEcalCandidate_*_*', 
        'keep *_hltL1sDoubleTauJet44Eta2p17orDoubleJet64Central_*_*', 
        'keep *_hltL1sDoubleTauJet44erorDoubleJetC64_*_*', 
        'keep *_hltL1sL1EG18er_*_*', 
        'keep *_hltL1sL1ETM36ORETM40_*_*', 
        'keep *_hltL1sL1Jet52ETM30_*_*', 
        'keep *_hltL1sL1SingleEG12_*_*', 
        'keep *_hltL1sL1SingleEG15_*_*', 
        'keep *_hltL1sL1SingleEG18orL1SingleEG20_*_*', 
        'keep *_hltL1sL1SingleMu10_*_*', 
        'keep *_hltL1sL1SingleMu14Eta2p1_*_*', 
        'keep *_hltL1sMu16Eta2p1_*_*', 
        'keep *_hltL2MuonCandidatesNoVtx_*_*', 
        'keep *_hltL2MuonCandidates_*_*', 
        'keep *_hltL2MuonSeeds_*_*', 
        'keep *_hltL2Muons_*_*', 
        'keep *_hltL2TauJets_*_*', 
        'keep *_hltL3MuonCandidates_*_*', 
        'keep *_hltL3MuonsIOHit_*_*', 
        'keep *_hltL3MuonsLinksCombination_*_*', 
        'keep *_hltL3MuonsOIHit_*_*', 
        'keep *_hltL3MuonsOIState_*_*', 
        'keep *_hltL3Muons_*_*', 
        'keep *_hltL3NoFiltersNoVtxMuonCandidates_*_*', 
        'keep *_hltL3NoFiltersNoVtxMuons_*_*', 
        'keep *_hltL3SecondaryVertexTagInfos_*_*', 
        'keep *_hltL3TkFromL2OICombination_*_*', 
        'keep *_hltL3TkTracksFromL2IOHit_*_*', 
        'keep *_hltL3TkTracksFromL2OIHit_*_*', 
        'keep *_hltL3TkTracksFromL2OIState_*_*', 
        'keep *_hltL3TkTracksFromL2_*_*', 
        'keep *_hltL3TrackCandidateFromL2IOHit_*_*', 
        'keep *_hltL3TrackCandidateFromL2OIHit_*_*', 
        'keep *_hltL3TrackCandidateFromL2OIState_*_*', 
        'keep *_hltL3TrajSeedIOHit_*_*', 
        'keep *_hltL3TrajSeedOIHit_*_*', 
        'keep *_hltL3TrajSeedOIState_*_*', 
        'keep *_hltL3TrajectorySeed_*_*', 
        'keep *_hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f18QL3crIsoRhoFiltered0p15_*_*', 
        'keep *_hltLeadingCentralJets30SingleTopEle27_*_*', 
        'keep *_hltLeadingCentralJets30SingleTopIsoMu24_*_*', 
        'keep *_hltMergedTracks_*_*', 
        'keep *_hltMet_*_*', 
        'keep *_hltMu8Ele17CaloIdTCaloIsoVLPixelMatchFilter_*_*', 
        'keep *_hltMuTrackJpsiCtfTrackCands_*_*', 
        'keep *_hltMuTrackJpsiPixelTrackCands_*_*', 
        'keep *_hltMuonCSCDigis_*_*', 
        'keep *_hltMuonCSCDigis_MuonCSCStripDigi_*', 
        'keep *_hltMuonCSCDigis_MuonCSCWireDigi_*', 
        'keep *_hltMuonDTDigis_*_*', 
        'keep *_hltMuonRPCDigis_*_*', 
        'keep *_hltOnlineBeamSpot_*_*', 
        'keep *_hltOnlinePrimaryVertices_*_*', 
        'keep *_hltOverlapFilterEle20LooseIsoPFTau20OldVersion_*_*', 
        'keep *_hltOverlapFilterIsoEle20MediumIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15IsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15MediumIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15TightIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu18LooseIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu18PFTau25TrackPt5Prong4_*_*', 
        'keep *_hltPFJetForBtag_*_*', 
        'keep *_hltPFMuonMerging_*_*', 
        'keep *_hltPFTau15TrackLooseIso_*_*', 
        'keep *_hltPFTau15Track_*_*', 
        'keep *_hltPFTau15_*_*', 
        'keep *_hltPFTau20IsoMuVertex_*_*', 
        'keep *_hltPFTau20TrackLooseIso_*_*', 
        'keep *_hltPFTau20Track_*_*', 
        'keep *_hltPFTau20_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolationProng4IsoMuVertex_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolationProng4_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolation_*_*', 
        'keep *_hltPFTau25TrackPt5_*_*', 
        'keep *_hltPFTau25_*_*', 
        'keep *_hltPFTau35TrackPt20LooseIsoProng2_*_*', 
        'keep *_hltPFTau35TrackPt20LooseIso_*_*', 
        'keep *_hltPFTau35TrackPt20_*_*', 
        'keep *_hltPFTau35Track_*_*', 
        'keep *_hltPFTau35_*_*', 
        'keep *_hltPFTauEleVertex20_*_*', 
        'keep *_hltPFTauJetTracksAssociator_*_*', 
        'keep *_hltPFTauMediumIso20TrackMediumIso_*_*', 
        'keep *_hltPFTauMediumIso20Track_*_*', 
        'keep *_hltPFTauMediumIso20_*_*', 
        'keep *_hltPFTauMediumIso35Track_*_*', 
        'keep *_hltPFTauMediumIso35_*_*', 
        'keep *_hltPFTauTagInfo_*_*', 
        'keep *_hltPFTauTightIso20TrackTightIso_*_*', 
        'keep *_hltPFTauTightIso20Track_*_*', 
        'keep *_hltPFTauTightIso20_*_*', 
        'keep *_hltPFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltParticleFlowForTaus_*_*', 
        'keep *_hltParticleFlow_*_*', 
        'keep *_hltPixelMatch3HitElectronsActivity_*_*', 
        'keep *_hltPixelMatch3HitElectronsL1Seeded_*_*', 
        'keep *_hltPixelMatchCleanElectronsL1Seeded_*_*', 
        'keep *_hltPixelMatchElectronsActivity_*_*', 
        'keep *_hltPixelMatchElectronsL1Iso_*_*', 
        'keep *_hltPixelMatchElectronsL1NonIso_*_*', 
        'keep *_hltPixelMatchElectronsL1Seeded_*_*', 
        'keep *_hltPixelTracksElectrons_*_*', 
        'keep *_hltPixelTracksMerged_*_*', 
        'keep *_hltPixelTracks_*_*', 
        'keep *_hltPixelVertices3DbbPhi_*_*', 
        'keep *_hltPixelVertices_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidateSC4_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidateSC5_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidate_*_*', 
        'keep *_hltRpcRecHits_*_*', 
        'keep *_hltSelector4CentralJetsL1FastJet_*_*', 
        'keep *_hltSelector8CentralJetsL1FastJet_*_*', 
        'keep *_hltSelectorJets20L1FastJet_*_*', 
        'keep *_hltSiPixelCluster_*_*', 
        'keep *_hltSiPixelClusters_*_*', 
        'keep *_hltSiStripClusters_*_*', 
        'keep *_hltSiStripRawToClustersFacility_*_*', 
        'keep *_hltSingleMu15L3Filtered15_*_*', 
        'keep *_hltSingleMuIsoL1s14L3IsoFiltered15eta2p1_*_*', 
        'keep *_hltSingleMuIsoL3IsoFiltered15_*_*', 
        'keep *_hltTowerMakerForAll_*_*', 
        'keep *_hltTowerMakerForMuons_*_*', 
        'keep *_hltTriggerSummaryAOD_*_*', 
        'keep *_hltTriggerSummaryRAW_*_*', 
        'keep *_hltTrimmedPixelVertices_*_*', 
        'keep *_hltVerticesL3_*_*', 
        'keep *_hltVerticesPFFilter_*_*', 
        'keep *_hltVerticesPFSelector_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_rawDataRepacker_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_virginRawDataRepacker_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep L2MuonTrajectorySeeds_hltL2MuonSeeds_*_*', 
        'keep L3MuonTrajectorySeeds_hltHIL3TrajSeedOIHit_*_*', 
        'keep L3MuonTrajectorySeeds_hltHIL3TrajectorySeed_*_*', 
        'keep L3MuonTrajectorySeeds_hltL3TrajSeedOIState_*_*', 
        'keep LumiScalerss_hltScalersRawToDigi_*_*', 
        'keep SiPixelClusteredmNewDetSetVector_hltSiPixelClusters_*_*', 
        'keep TrackCandidates_hltHIL3TrackCandidateFromL2OIHit_*_*', 
        'keep TrackCandidates_hltHIL3TrackCandidateFromL2OIState_*_*', 
        'keep TrackingRecHitsOwned_hltL3Muons_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep recoCaloJets_*_*_*', 
        'keep recoCaloMETs_*_*_*', 
        'keep recoCaloMETs_hltMet_*_*', 
        'keep recoCompositeCandidates_*_*_*', 
        'keep recoElectrons_*_*_*', 
        'keep recoIsolatedPixelTrackCandidates_*_*_*', 
        'keep recoMETs_*_*_*', 
        'keep recoPFJets_*_*_*', 
        'keep recoPFTaus_*_*_*', 
        'keep recoRecoChargedCandidates_*_*_*', 
        'keep recoRecoChargedCandidates_hltHIL3MuonCandidates_*_*', 
        'keep recoRecoChargedCandidates_hltL2MuonCandidates_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsoHLTClusterShape_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonEcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalForHE_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsoHLTClusterShape_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonEcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalForHE_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalIsol_*_*', 
        'keep recoRecoEcalCandidates_*_*_*', 
        'keep recoRecoEcalCandidates_hltL1IsoRecoEcalCandidate_*_*', 
        'keep recoRecoEcalCandidates_hltL1NonIsoRecoEcalCandidate_*_*', 
        'keep recoTrackExtras_hltHIL3MuonsOIHit_*_*', 
        'keep recoTrackExtras_hltHIL3MuonsOIState_*_*', 
        'keep recoTracks_hltHIL3MuonsOIHit_*_*', 
        'keep recoTracks_hltHIL3MuonsOIState_*_*', 
        'keep recoTracks_hltHIL3Muons_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2OIHit_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2OIState_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2_*_*', 
        'keep triggerTriggerEventWithRefs_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep triggerTriggerFilterObjectWithRefs_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*' ) ),
    splitLevel = cms.untracked.int32(0)
)

process.RAWEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*'),
    splitLevel = cms.untracked.int32(0)
)

process.RAWMINIAODEventContent = cms.PSet(
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
    outputCommands = cms.untracked.vstring('drop *', 
        'drop *', 
        'keep *_slimmedPhotons_*_*', 
        'keep *_slimmedOOTPhotons_*_*', 
        'keep *_slimmedElectrons_*_*', 
        'keep *_slimmedMuons_*_*', 
        'keep *_slimmedTaus_*_*', 
        'keep *_slimmedTausBoosted_*_*', 
        'keep *_slimmedCaloJets_*_*', 
        'keep *_slimmedJets_*_*', 
        'drop recoBaseTagInfosOwned_slimmedJets_*_*', 
        'keep *_slimmedJetsAK8_*_*', 
        'keep *_slimmedJetsPuppi_*_*', 
        'keep *_slimmedMETs_*_*', 
        'keep *_slimmedMETsNoHF_*_*', 
        'keep *_slimmedMETsPuppi_*_*', 
        'keep *_slimmedSecondaryVertices_*_*', 
        'keep *_slimmedLambdaVertices_*_*', 
        'keep *_slimmedKshortVertices_*_*', 
        'keep *_slimmedJetsAK8PFPuppiSoftDropPacked_SubJets_*', 
        'keep recoPhotonCores_reducedEgamma_*_*', 
        'keep recoGsfElectronCores_reducedEgamma_*_*', 
        'keep recoConversions_reducedEgamma_*_*', 
        'keep recoSuperClusters_reducedEgamma_*_*', 
        'keep recoCaloClusters_reducedEgamma_*_*', 
        'keep EcalRecHitsSorted_reducedEgamma_*_*', 
        'keep recoGsfTracks_reducedEgamma_*_*', 
        'drop *_*_caloTowers_*', 
        'drop *_*_pfCandidates_*', 
        'drop *_*_genJets_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep *_offlineSlimmedPrimaryVertices_*_*', 
        'keep patPackedCandidates_packedPFCandidates_*_*', 
        'keep *_isolatedTracks_*_*', 
        'keep *_oniaPhotonCandidates_*_*', 
        'keep *_bunchSpacingProducer_*_*', 
        'keep double_fixedGridRhoAll__*', 
        'keep double_fixedGridRhoFastjetAll__*', 
        'keep double_fixedGridRhoFastjetAllCalo__*', 
        'keep double_fixedGridRhoFastjetCentral_*_*', 
        'keep double_fixedGridRhoFastjetCentralCalo__*', 
        'keep double_fixedGridRhoFastjetCentralChargedPileUp__*', 
        'keep double_fixedGridRhoFastjetCentralNeutral__*', 
        'keep *_slimmedPatTrigger_*_*', 
        'keep patPackedTriggerPrescales_patTrigger__*', 
        'keep patPackedTriggerPrescales_patTrigger_l1max_*', 
        'keep patPackedTriggerPrescales_patTrigger_l1min_*', 
        'keep *_l1extraParticles_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_gtStage2Digis__*', 
        'keep *_gmtStage2Digis_Muon_*', 
        'keep *_caloStage2Digis_Jet_*', 
        'keep *_caloStage2Digis_Tau_*', 
        'keep *_caloStage2Digis_EGamma_*', 
        'keep *_caloStage2Digis_EtSum_*', 
        'keep *_TriggerResults_*_HLT', 
        'keep *_TriggerResults_*_*', 
        'keep patPackedCandidates_lostTracks_*_*', 
        'keep HcalNoiseSummary_hcalnoise__*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep CTPPSLocalTrackLites_ctppsLocalTrackLiteProducer_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*')
)

process.RAWMINIAODSIMEventContent = cms.PSet(
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
    outputCommands = cms.untracked.vstring('drop *', 
        'drop *', 
        'keep *_slimmedPhotons_*_*', 
        'keep *_slimmedOOTPhotons_*_*', 
        'keep *_slimmedElectrons_*_*', 
        'keep *_slimmedMuons_*_*', 
        'keep *_slimmedTaus_*_*', 
        'keep *_slimmedTausBoosted_*_*', 
        'keep *_slimmedCaloJets_*_*', 
        'keep *_slimmedJets_*_*', 
        'drop recoBaseTagInfosOwned_slimmedJets_*_*', 
        'keep *_slimmedJetsAK8_*_*', 
        'keep *_slimmedJetsPuppi_*_*', 
        'keep *_slimmedMETs_*_*', 
        'keep *_slimmedMETsNoHF_*_*', 
        'keep *_slimmedMETsPuppi_*_*', 
        'keep *_slimmedSecondaryVertices_*_*', 
        'keep *_slimmedLambdaVertices_*_*', 
        'keep *_slimmedKshortVertices_*_*', 
        'keep *_slimmedJetsAK8PFPuppiSoftDropPacked_SubJets_*', 
        'keep recoPhotonCores_reducedEgamma_*_*', 
        'keep recoGsfElectronCores_reducedEgamma_*_*', 
        'keep recoConversions_reducedEgamma_*_*', 
        'keep recoSuperClusters_reducedEgamma_*_*', 
        'keep recoCaloClusters_reducedEgamma_*_*', 
        'keep EcalRecHitsSorted_reducedEgamma_*_*', 
        'keep recoGsfTracks_reducedEgamma_*_*', 
        'drop *_*_caloTowers_*', 
        'drop *_*_pfCandidates_*', 
        'drop *_*_genJets_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep *_offlineSlimmedPrimaryVertices_*_*', 
        'keep patPackedCandidates_packedPFCandidates_*_*', 
        'keep *_isolatedTracks_*_*', 
        'keep *_oniaPhotonCandidates_*_*', 
        'keep *_bunchSpacingProducer_*_*', 
        'keep double_fixedGridRhoAll__*', 
        'keep double_fixedGridRhoFastjetAll__*', 
        'keep double_fixedGridRhoFastjetAllCalo__*', 
        'keep double_fixedGridRhoFastjetCentral_*_*', 
        'keep double_fixedGridRhoFastjetCentralCalo__*', 
        'keep double_fixedGridRhoFastjetCentralChargedPileUp__*', 
        'keep double_fixedGridRhoFastjetCentralNeutral__*', 
        'keep *_slimmedPatTrigger_*_*', 
        'keep patPackedTriggerPrescales_patTrigger__*', 
        'keep patPackedTriggerPrescales_patTrigger_l1max_*', 
        'keep patPackedTriggerPrescales_patTrigger_l1min_*', 
        'keep *_l1extraParticles_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_gtStage2Digis__*', 
        'keep *_gmtStage2Digis_Muon_*', 
        'keep *_caloStage2Digis_Jet_*', 
        'keep *_caloStage2Digis_Tau_*', 
        'keep *_caloStage2Digis_EGamma_*', 
        'keep *_caloStage2Digis_EtSum_*', 
        'keep *_TriggerResults_*_HLT', 
        'keep *_TriggerResults_*_*', 
        'keep patPackedCandidates_lostTracks_*_*', 
        'keep HcalNoiseSummary_hcalnoise__*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep CTPPSLocalTrackLites_ctppsLocalTrackLiteProducer_*_*', 
        'keep patPackedGenParticles_packedGenParticles_*_*', 
        'keep recoGenParticles_prunedGenParticles_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_*_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep recoGenParticles_genPUProtons_*_*', 
        'keep *_slimmedGenJetsFlavourInfos_*_*', 
        'keep *_slimmedGenJets__*', 
        'keep *_slimmedGenJetsAK8__*', 
        'keep *_slimmedGenJetsAK8SoftDropSubJets__*', 
        'keep *_genMetTrue_*_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep GenRunInfoProduct_*_*_*', 
        'keep *_genParticles_xyz0_*', 
        'keep *_genParticles_t0_*', 
        'keep PileupSummaryInfos_slimmedAddPileupInfo_*_*', 
        'keep L1GtTriggerMenuLite_l1GtTriggerMenuLite__*', 
        'keep SimVertexs_g4SimHits_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*')
)

process.RAWRECODEBUGHLTEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring( ('drop *', 
        'drop *', 
        'drop *', 
        'drop *', 
        'keep DetIdedmEDCollection_siStripDigis_*_*', 
        'keep DetIdedmEDCollection_siPixelDigis_*_*', 
        'keep PixelFEDChanneledmNewDetSetVector_siPixelDigis_*_*', 
        'keep *_siPixelClusters_*_*', 
        'keep *_siStripClusters_*_*', 
        'keep ClusterSummary_clusterSummaryProducer_*_*', 
        'keep *_dt1DRecHits_*_*', 
        'keep *_dt4DSegments_*_*', 
        'keep *_dt1DCosmicRecHits_*_*', 
        'keep *_dt4DCosmicSegments_*_*', 
        'keep *_csc2DRecHits_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*', 
        'keep *_hbhereco_*_*', 
        'keep *_hbheprereco_*_*', 
        'keep *_hfprereco_*_*', 
        'keep *_hfreco_*_*', 
        'keep *_horeco_*_*', 
        'keep HBHERecHitsSorted_hbherecoMB_*_*', 
        'keep HORecHitsSorted_horecoMB_*_*', 
        'keep HFRecHitsSorted_hfrecoMB_*_*', 
        'keep ZDCDataFramesSorted_hcalDigis_*_*', 
        'keep ZDCDataFramesSorted_castorDigis_*_*', 
        'keep QIE10DataFrameHcalDataFrameContainer_hcalDigis_ZDC_*', 
        'keep ZDCRecHitsSorted_zdcreco_*_*', 
        'keep *_reducedHcalRecHits_*_*', 
        'keep *_castorreco_*_*', 
        'keep HcalUnpackerReport_castorDigis_*_*', 
        'keep HcalUnpackerReport_hcalDigiAlCaMB_*_*', 
        'keep HcalUnpackerReport_hcalDigis_*_*', 
        'keep *_ecalPreshowerRecHit_*_*', 
        'keep *_ecalRecHit_*_*', 
        'keep *_ecalCompactTrigPrim_*_*', 
        'keep *_ecalTPSkim_*_*', 
        'keep EBSrFlagsSorted_ecalDigis__*', 
        'keep EESrFlagsSorted_ecalDigis__*', 
        'keep *_selectDigi_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEE_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEB_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsES_*_*', 
        'keep *_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep *_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_*', 
        'keep *_particleFlowSuperClusterECAL_*_*', 
        'keep *_particleFlowSuperClusterOOTECAL_*_*', 
        'drop recoClusterShapes_*_*_*', 
        'drop recoBasicClustersToOnerecoClusterShapesAssociation_*_*_*', 
        'drop recoBasicClusters_multi5x5BasicClusters_multi5x5BarrelBasicClusters_*', 
        'drop recoSuperClusters_multi5x5SuperClusters_multi5x5BarrelSuperClusters_*', 
        'keep *_CkfElectronCandidates_*_*', 
        'keep *_GsfGlobalElectronTest_*_*', 
        'keep *_electronMergedSeeds_*_*', 
        'keep recoGsfTracks_electronGsfTracks_*_*', 
        'keep recoGsfTrackExtras_electronGsfTracks_*_*', 
        'keep recoTrackExtras_electronGsfTracks_*_*', 
        'keep TrackingRecHitsOwned_electronGsfTracks_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTrackExtras_generalTracks_*_*', 
        'keep TrackingRecHitsOwned_generalTracks_*_*', 
        'keep *_generalTracks_MVAValues_*', 
        'keep *_generalTracks_MVAVals_*', 
        'keep TrackingRecHitsOwned_extraFromSeeds_*_*', 
        'keep uints_extraFromSeeds_*_*', 
        'keep recoTracks_beamhaloTracks_*_*', 
        'keep recoTrackExtras_beamhaloTracks_*_*', 
        'keep TrackingRecHitsOwned_beamhaloTracks_*_*', 
        'keep recoTracks_conversionStepTracks_*_*', 
        'keep recoTrackExtras_conversionStepTracks_*_*', 
        'keep TrackingRecHitsOwned_conversionStepTracks_*_*', 
        'keep *_ctfPixelLess_*_*', 
        'keep *_dedxTruncated40_*_*', 
        'keep *_dedxHitInfo_*_*', 
        'keep *_dedxHarmonic2_*_*', 
        'keep *_dedxPixelHarmonic2_*_*', 
        'keep *_trackExtrapolator_*_*', 
        'keep recoTracks_cosmicDCTracks_*_*', 
        'keep recoTrackExtras_cosmicDCTracks_*_*', 
        'keep TrackingRecHitsOwned_cosmicDCTracks_*_*', 
        'keep *_ak4CaloJets_*_*', 
        'keep *_ak4PFJets_*_*', 
        'keep *_ak4PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHSSoftDrop_*_*', 
        'keep *_cmsTopTagPFJetsCHS_*_*', 
        'keep *_JetPlusTrackZSPCorJetAntiKt4_*_*', 
        'keep *_ak4TrackJets_*_*', 
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', 
        'keep *_caloTowers_*_*', 
        'keep *_towerMaker_*_*', 
        'keep *_CastorTowerReco_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertex_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertexPF_*_*', 
        'keep *_ak4JetTracksAssociatorAtCaloFace_*_*', 
        'keep *_ak4JetTracksAssociatorExplicit_*_*', 
        'keep *_ak4JetExtender_*_*', 
        'keep *_ak4JetID_*_*', 
        'keep *_ak5CastorJets_*_*', 
        'keep *_ak5CastorJetID_*_*', 
        'keep *_ak7CastorJets_*_*', 
        'keep *_ak7CastorJetID_*_*', 
        'keep *_fixedGridRhoAll_*_*', 
        'keep *_fixedGridRhoFastjetAll_*_*', 
        'keep *_fixedGridRhoFastjetAllTmp_*_*', 
        'keep *_fixedGridRhoFastjetAllCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentral_*_*', 
        'keep *_fixedGridRhoFastjetCentralCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralChargedPileUp_*_*', 
        'keep *_fixedGridRhoFastjetCentralNeutral_*_*', 
        'keep *_ak8PFJetsCHSSoftDropMass_*_*', 
        'keep recoCaloMETs_caloMet_*_*', 
        'keep recoCaloMETs_caloMetBE_*_*', 
        'keep recoCaloMETs_caloMetBEFO_*_*', 
        'keep recoCaloMETs_caloMetM_*_*', 
        'keep recoPFMETs_pfMet_*_*', 
        'keep recoPFMETs_pfChMet_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_*_*', 
        'keep recoHcalNoiseRBXs_hcalnoise_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoEcalHaloData_EcalHaloData_*_*', 
        'keep recoGlobalHaloData_GlobalHaloData_*_*', 
        'keep recoHcalHaloData_HcalHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep *_MuonSeed_*_*', 
        'keep *_ancientMuonSeed_*_*', 
        'keep *_displacedMuonSeeds_*_*', 
        'keep TrackingRecHitsOwned_globalMuons_*_*', 
        'keep TrackingRecHitsOwned_tevMuons_*_*', 
        'keep *_CosmicMuonSeed_*_*', 
        'keep recoTrackExtras_cosmicMuons_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons_*_*', 
        'keep recoTrackExtras_cosmicMuons1Leg_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons1Leg_*_*', 
        'keep recoTracks_cosmicsVetoTracks_*_*', 
        'keep *_muons_*_*', 
        'keep *_particleFlow_muons_*', 
        'drop *_muons_muons1stStep2muonsMap_*', 
        'drop recoIsoDepositedmValueMap_muons_*_*', 
        'drop doubleedmValueMap_muons_muPFIso*_*', 
        'keep recoTracks_standAloneMuons_*_*', 
        'keep recoTrackExtras_standAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_standAloneMuons_*_*', 
        'keep recoTracks_globalMuons_*_*', 
        'keep recoTrackExtras_globalMuons_*_*', 
        'keep recoTracks_tevMuons_*_*', 
        'keep recoTrackExtras_tevMuons_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_displacedTracks_*_*', 
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_*_*', 
        'keep recoTracks_displacedGlobalMuons_*_*', 
        'keep recoTrackExtras_displacedGlobalMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedGlobalMuons_*_*', 
        'keep recoTracks_cosmicMuons_*_*', 
        'keep recoMuons_muonsFromCosmics_*_*', 
        'keep recoTracks_cosmicMuons1Leg_*_*', 
        'keep recoMuons_muonsFromCosmics1Leg_*_*', 
        'keep recoTracks_refittedStandAloneMuons_*_*', 
        'keep recoTrackExtras_refittedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_refittedStandAloneMuons_*_*', 
        'keep recoTracks_displacedStandAloneMuons__*', 
        'keep recoTrackExtras_displacedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedStandAloneMuons_*_*', 
        'keep *_muIsoDepositTk_*_*', 
        'keep *_muIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muIsoDepositJets_*_*', 
        'keep *_muGlobalIsoDepositCtfTk_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muGlobalIsoDepositJets_*_*', 
        'keep *_softPFMuonsTagInfos_*_*', 
        'keep *_softPFElectronsTagInfos_*_*', 
        'keep *_softPFElectronBJetTags_*_*', 
        'keep *_softPFMuonBJetTags_*_*', 
        'keep *_pfImpactParameterTagInfos_*_*', 
        'keep *_pfTrackCountingHighEffBJetTags_*_*', 
        'keep *_pfJetProbabilityBJetTags_*_*', 
        'keep *_pfJetBProbabilityBJetTags_*_*', 
        'keep *_pfSecondaryVertexTagInfos_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderTagInfos_*_*', 
        'keep *_pfSimpleSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfCombinedSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfGhostTrackVertexTagInfos_*_*', 
        'keep *_pfGhostTrackBJetTags_*_*', 
        'keep *_pfCombinedMVAV2BJetTags_*_*', 
        'keep *_inclusiveCandidateSecondaryVertices_*_*', 
        'keep *_inclusiveCandidateSecondaryVerticesCvsL_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderCvsLTagInfos_*_*', 
        'keep *_pfCombinedCvsLJetTags_*_*', 
        'keep *_pfCombinedCvsBJetTags_*_*', 
        'keep *_pfChargeBJetTags_*_*', 
        'keep *_pfDeepCSVJetTags_*_*', 
        'keep *_pfDeepCMVAJetTags_*_*', 
        'keep recoRecoTauPiZeros_hpsPFTauProducer_pizeros_*', 
        'keep recoPFTaus_hpsPFTauProducer_*_*', 
        'keep *_hpsPFTauDiscriminationByDeadECALElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFinding_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingNewDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingOldDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseChargedIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByTightMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSum_*_*', 
        'keep *_hpsPFTauPUcorrPtSum_*_*', 
        'keep *_hpsPFTauChargedIsoPtSum_*_*', 
        'keep *_hpsPFTauTransverseImpactParameters_*_*', 
        'keep *_hpsPFTauFootprintCorrection_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeight_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6rawElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6LooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6MediumElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6TightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VTightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauChargedIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeightdR03_*_*', 
        'keep *_hpsPFTauFootprintCorrectiondR03_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalConedR03_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep  *_offlinePrimaryVertices__*', 
        'keep *_offlinePrimaryVerticesWithBS_*_*', 
        'keep *_offlinePrimaryVerticesFromCosmicTracks_*_*', 
        'keep *_nuclearInteractionMaker_*_*', 
        'keep *_generalV0Candidates_*_*', 
        'keep *_inclusiveSecondaryVertices_*_*', 
        'keep recoGsfElectronCores_gsfElectronCores_*_*', 
        'keep recoGsfElectronCores_gedGsfElectronCores_*_*', 
        'keep recoGsfElectrons_gsfElectrons_*_*', 
        'keep recoGsfElectrons_gedGsfElectrons_*_*', 
        'keep recoGsfElectronCores_uncleanedOnlyGsfElectronCores_*_*', 
        'keep recoGsfElectrons_uncleanedOnlyGsfElectrons_*_*', 
        'keep floatedmValueMap_eidRobustLoose_*_*', 
        'keep floatedmValueMap_eidRobustTight_*_*', 
        'keep floatedmValueMap_eidRobustHighEnergy_*_*', 
        'keep floatedmValueMap_eidLoose_*_*', 
        'keep floatedmValueMap_eidTight_*_*', 
        'keep *_egmGedGsfElectronPFIsolation_*_*', 
        'keep *_photonEcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonEcalPFClusterIsolationProducer_*_*', 
        'keep *_electronEcalPFClusterIsolationProducer_*_*', 
        'keep *_photonHcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonHcalPFClusterIsolationProducer_*_*', 
        'keep *_electronHcalPFClusterIsolationProducer_*_*', 
        'drop *_egmGsfElectronIDs_*_*', 
        'drop *_egmPhotonIDs_*_*', 
        'keep *_gedPhotonCore_*_*', 
        'keep *_gedPhotons_*_*', 
        'keep *_particleBasedIsolation_*_*', 
        'keep recoPhotons_mustachePhotons_*_*', 
        'keep recoPhotonCores_mustachePhotonCore_*_*', 
        'keep recoPhotons_ootPhotons_*_*', 
        'keep recoPhotonCores_ootPhotonCore_*_*', 
        'keep recoPhotons_photons_*_*', 
        'keep recoPhotonCores_photonCore_*_*', 
        'keep recoConversions_conversions_*_*', 
        'keep recoConversions_mustacheConversions_*_*', 
        'drop *_conversions_uncleanedConversions_*', 
        'drop *_gedPhotonsTmp_valMapPFEgammaCandToPhoton_*', 
        'keep recoConversions_allConversions_*_*', 
        'keep recoConversions_allConversionsOldEG_*_*', 
        'keep recoTracks_ckfOutInTracksFromConversions_*_*', 
        'keep recoTracks_ckfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfInOutTracksFromConversions_*_*', 
        'keep recoConversions_uncleanedOnlyAllConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep *_PhotonIDProdGED_*_*', 
        'keep recoRecoEcalCandidates_hfRecoEcalCandidate_*_*', 
        'keep *_hfEMClusters_*_*', 
        'keep *_gedGsfElectronCores_*_*', 
        'keep *_gedGsfElectrons_*_*', 
        'keep *_pixelTracks_*_*', 
        'keep *_pixelVertices_*_*', 
        'drop CaloTowersSorted_towerMakerPF_*_*', 
        'keep recoPFRecHits_particleFlowClusterECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHCAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterPS_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHBHE_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitPS_Cleaned_*', 
        'keep recoPFClusters_particleFlowClusterECAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHCAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHO_*_*', 
        'keep recoPFClusters_particleFlowClusterHF_*_*', 
        'keep recoPFClusters_particleFlowClusterPS_*_*', 
        'keep recoPFBlocks_particleFlowBlock_*_*', 
        'keep recoPFCandidates_particleFlowEGamma_*_*', 
        'keep recoCaloClusters_particleFlowEGamma_*_*', 
        'keep recoSuperClusters_particleFlowEGamma_*_*', 
        'keep recoConversions_particleFlowEGamma_*_*', 
        'keep recoPFCandidates_particleFlow_*_*', 
        'keep recoPFCandidates_particleFlowTmp_electrons_*', 
        'keep recoPFCandidates_particleFlowTmp_*_*', 
        'drop recoPFCandidates_particleFlowTmp__*', 
        'keep recoPFDisplacedVertexs_particleFlowDisplacedVertex_*_*', 
        'keep *_pfElectronTranslator_*_*', 
        'keep *_pfPhotonTranslator_*_*', 
        'keep *_particleFlow_electrons_*', 
        'keep *_particleFlow_photons_*', 
        'keep *_particleFlow_muons_*', 
        'keep *_trackerDrivenElectronSeeds_preid_*', 
        'keep *_particleFlowPtrs_*_*', 
        'keep *_particleFlowTmpPtrs_*_*', 
        'keep *_chargedHadronPFTrackIsolation_*_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep L1MuGMTReadoutCollection_gtDigis_*_*', 
        'keep L1GctEmCand*_gctDigis_*_*', 
        'keep L1GctJetCand*_gctDigis_*_*', 
        'keep L1GctEtHad*_gctDigis_*_*', 
        'keep L1GctEtMiss*_gctDigis_*_*', 
        'keep L1GctEtTotal*_gctDigis_*_*', 
        'keep L1GctHtMiss*_gctDigis_*_*', 
        'keep L1GctJetCounts*_gctDigis_*_*', 
        'keep L1GctHFRingEtSums*_gctDigis_*_*', 
        'keep L1GctHFBitCounts*_gctDigis_*_*', 
        'keep LumiDetails_lumiProducer_*_*', 
        'keep LumiSummary_lumiProducer_*_*', 
        'drop *_hlt*_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep L1AcceptBunchCrossings_scalersRawToDigi_*_*', 
        'keep L1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep Level1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep BeamSpotOnlines_scalersRawToDigi_*_*', 
        'keep DcsStatuss_scalersRawToDigi_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep *_tcdsDigis_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'keep *_pfIsolatedElectronsEI_*_*', 
        'keep *_pfIsolatedMuonsEI_*_*', 
        'keep recoPFJets_pfJetsEI_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTagsEI_*_*', 
        'keep recoPFTaus_pfTausEI_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByDecayModeFinding_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByIsolation_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep *_ak4GenJets_*_*', 
        'keep *_ak8GenJets_*_*', 
        'keep *_ak4GenJetsNoNu_*_*', 
        'keep *_ak8GenJetsNoNu_*_*', 
        'keep *_genParticle_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep SimTracks_g4SimHits_*_*', 
        'keep SimVertexs_g4SimHits_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*', 
        'keep *_muonSimClassifier_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'drop *_hlt*_*_*', 
        'keep *_hltAK4CaloJetsCorrectedIDPassed_*_*', 
        'keep *_hltAK4CaloJetsIDPassed_*_*', 
        'keep *_hltAK4CaloJets_*_*', 
        'keep *_hltAK4PFJetsCorrected_*_*', 
        'keep *_hltAK4PFJetsForTaus_*_*', 
        'keep *_hltAK4PFJets_*_*', 
        'keep *_hltAlCaEtaEBRechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaEtaEBRechitsToDigis_*_*', 
        'keep *_hltAlCaEtaEERechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaEtaEERechitsToDigis_*_*', 
        'keep *_hltAlCaEtaRecHitsFilterEBonly_*_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegionalLowPU_etaEcalRecHitsES_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegional_etaEcalRecHitsES_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonly_*_*', 
        'keep *_hltAlCaEtaRecHitsFilter_*_*', 
        'keep *_hltAlCaPhiSymStream_*_*', 
        'keep *_hltAlCaPi0EBRechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaPi0EBRechitsToDigis_*_*', 
        'keep *_hltAlCaPi0EERechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaPi0EERechitsToDigis_*_*', 
        'keep *_hltAlCaPi0RecHitsFilterEBonly_*_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegionalLowPU_pi0EcalRecHitsES_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegional_pi0EcalRecHitsES_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonly_*_*', 
        'keep *_hltAlCaPi0RecHitsFilter_*_*', 
        'keep *_hltBLifetimeL25AssociatorbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL25BJetTagsbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL25TagInfosbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3AssociatorbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3BJetTagsbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3TagInfosbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeRegionalCtfWithMaterialTracksbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBSoftMuonDiJet110Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet110Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet20Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet20Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet40Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet40Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet70Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet70Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonMu5L3_*_*', 
        'keep *_hltCSVJetTagSingleTopEle27_*_*', 
        'keep *_hltCSVJetTagSingleTopIsoMu24_*_*', 
        'keep *_hltCaloJetCorrectedRegional_*_*', 
        'keep *_hltCaloJetCorrected_*_*', 
        'keep *_hltCaloJetL1FastJetCorrected_*_*', 
        'keep *_hltCaloStage2Digis_*_*', 
        'keep *_hltCleanedCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCleanedHiCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCombinedSecondaryVertexBJetTagsCalo_*_*', 
        'keep *_hltCombinedSecondaryVertexBJetTagsPF_*_*', 
        'keep *_hltConvPFTausTightIsoTrackFindingIsolation_*_*', 
        'keep *_hltConvPFTausTightIsoTrackFinding_*_*', 
        'keep *_hltConvPFTausTightIsoTrackPt5Isolation_*_*', 
        'keep *_hltConvPFTausTightIsoTrackPt5_*_*', 
        'keep *_hltConvPFTausTightIso_*_*', 
        'keep *_hltConvPFTausTrackFindingLooseIsolation_*_*', 
        'keep *_hltConvPFTausTrackFinding_*_*', 
        'keep *_hltConvPFTaus_*_*', 
        'keep *_hltCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCorrectedIslandEndcapSuperClustersHI_*_*', 
        'keep *_hltCsc2DRecHits_*_*', 
        'keep *_hltCscSegments_*_*', 
        'keep *_hltDisplacedhltIter4PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolationProng4L1HLTMatched_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolationProng4_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolation_*_*', 
        'keep *_hltDoublePFTau25TrackPt5_*_*', 
        'keep *_hltDoublePFTau25_*_*', 
        'keep *_hltDoublePFTauTightIso45Track5_*_*', 
        'keep *_hltDoublePFTauTightIso45Track_*_*', 
        'keep *_hltDt4DSegments_*_*', 
        'keep *_hltEcalPhiSymFilter_*_*', 
        'keep *_hltEcalRecHitAll_*_*', 
        'keep *_hltEcalRecHit_*_*', 
        'keep *_hltEgammaCandidates_*_*', 
        'keep *_hltEgammaGsfElectrons_*_*', 
        'keep *_hltEgammaGsfTracks_*_*', 
        'keep *_hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1SingleEG18orL1SingleEG20_*_*', 
        'keep *_hltEle20CaloIdVTTrkIdTDphiFilter_*_*', 
        'keep *_hltEle27WP85PixelMatchElectronsL1Seeded_*_*', 
        'keep *_hltElectronsVertex_*_*', 
        'keep *_hltFEDSelectorLumiPixels_*_*', 
        'keep *_hltFastPVPixelTracksMerger_*_*', 
        'keep *_hltFastPVPixelTracksRecover_*_*', 
        'keep *_hltFastPVPixelTracks_*_*', 
        'keep *_hltFastPVPixelVertices3D_*_*', 
        'keep *_hltFastPVPixelVertices_*_*', 
        'keep *_hltFastPixelBLifetimeL3TagInfos_*_*', 
        'keep *_hltFastPrimaryVertex_*_*', 
        'keep *_hltFilterDoubleIsoPFTau45Trk5LeadTrack5IsolationL1HLTMatched_*_*', 
        'keep *_hltFilterL2EtCutDoublePFIsoTau45Trk5_*_*', 
        'keep *_hltFilterL2EtCutSingleIsoPFTau35Trk20MET70_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20LeadTrackPt20_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20MET60LeadTrack20IsolationL1HLTMatched_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20MET70LeadTrack20IsolationL1HLTMatched_*_*', 
        'keep *_hltGmtStage2Digis_*_*', 
        'keep *_hltGoodOnlinePVs_*_*', 
        'keep *_hltGtStage2Digis_*_*', 
        'keep *_hltHICaloJetCorrected_*_*', 
        'keep *_hltHICaloJetIDPassed_*_*', 
        'keep *_hltHIGoodLooseTracks_*_*', 
        'keep *_hltHIPixel3PrimTracks_*_*', 
        'keep *_hltHISelectedVertex_*_*', 
        'keep *_hltHISiPixelClusters_*_*', 
        'keep *_hltHITIPTCorrectorHB_*_*', 
        'keep *_hltHITIPTCorrectorHE_*_*', 
        'keep *_hltHbhereco_*_*', 
        'keep *_hltHfreco_*_*', 
        'keep *_hltHiCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltHiCorrectedIslandEndcapSuperClustersHI_*_*', 
        'keep *_hltHiIslandSuperClustersHI_*_*', 
        'keep *_hltHoreco_*_*', 
        'keep *_hltIsolPixelTrackProdHB_*_*', 
        'keep *_hltIsolPixelTrackProdHE_*_*', 
        'keep *_hltIter0ElectronsTrackSelectionHighPurity_*_*', 
        'keep *_hltIter0HighPtTkMuPixelTracks_*_*', 
        'keep *_hltIter0HighPtTkMuTrackSelectionHighPurity_*_*', 
        'keep *_hltIter0PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter1Merged_*_*', 
        'keep *_hltIter1PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter2HighPtTkMuMerged_*_*', 
        'keep *_hltIter2HighPtTkMuTrackSelectionHighPurity_*_*', 
        'keep *_hltIter2MergedForBTag_*_*', 
        'keep *_hltIter2MergedForElectrons_*_*', 
        'keep *_hltIter2Merged_*_*', 
        'keep *_hltIter2PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter3Merged_*_*', 
        'keep *_hltIter4Merged_*_*', 
        'keep *_hltIterativeCone5PileupSubtractionCaloJets_*_*', 
        'keep *_hltL1HLTSingleIsoPFTau35Trk20Met60JetsMatch_*_*', 
        'keep *_hltL1IsoElectronTrackIsol_*_*', 
        'keep *_hltL1NonIsoElectronTrackIsol_*_*', 
        'keep *_hltL1SeededRecoEcalCandidate_*_*', 
        'keep *_hltL1sDoubleTauJet44Eta2p17orDoubleJet64Central_*_*', 
        'keep *_hltL1sDoubleTauJet44erorDoubleJetC64_*_*', 
        'keep *_hltL1sL1EG18er_*_*', 
        'keep *_hltL1sL1ETM36ORETM40_*_*', 
        'keep *_hltL1sL1Jet52ETM30_*_*', 
        'keep *_hltL1sL1SingleEG12_*_*', 
        'keep *_hltL1sL1SingleEG15_*_*', 
        'keep *_hltL1sL1SingleEG18orL1SingleEG20_*_*', 
        'keep *_hltL1sL1SingleMu10_*_*', 
        'keep *_hltL1sL1SingleMu14Eta2p1_*_*', 
        'keep *_hltL1sMu16Eta2p1_*_*', 
        'keep *_hltL2MuonCandidatesNoVtx_*_*', 
        'keep *_hltL2MuonCandidates_*_*', 
        'keep *_hltL2MuonSeeds_*_*', 
        'keep *_hltL2Muons_*_*', 
        'keep *_hltL2TauJets_*_*', 
        'keep *_hltL3MuonCandidates_*_*', 
        'keep *_hltL3MuonsIOHit_*_*', 
        'keep *_hltL3MuonsLinksCombination_*_*', 
        'keep *_hltL3MuonsOIHit_*_*', 
        'keep *_hltL3MuonsOIState_*_*', 
        'keep *_hltL3Muons_*_*', 
        'keep *_hltL3NoFiltersNoVtxMuonCandidates_*_*', 
        'keep *_hltL3NoFiltersNoVtxMuons_*_*', 
        'keep *_hltL3SecondaryVertexTagInfos_*_*', 
        'keep *_hltL3TkFromL2OICombination_*_*', 
        'keep *_hltL3TkTracksFromL2IOHit_*_*', 
        'keep *_hltL3TkTracksFromL2OIHit_*_*', 
        'keep *_hltL3TkTracksFromL2OIState_*_*', 
        'keep *_hltL3TkTracksFromL2_*_*', 
        'keep *_hltL3TrackCandidateFromL2IOHit_*_*', 
        'keep *_hltL3TrackCandidateFromL2OIHit_*_*', 
        'keep *_hltL3TrackCandidateFromL2OIState_*_*', 
        'keep *_hltL3TrajSeedIOHit_*_*', 
        'keep *_hltL3TrajSeedOIHit_*_*', 
        'keep *_hltL3TrajSeedOIState_*_*', 
        'keep *_hltL3TrajectorySeed_*_*', 
        'keep *_hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f18QL3crIsoRhoFiltered0p15_*_*', 
        'keep *_hltLeadingCentralJets30SingleTopEle27_*_*', 
        'keep *_hltLeadingCentralJets30SingleTopIsoMu24_*_*', 
        'keep *_hltMergedTracks_*_*', 
        'keep *_hltMet_*_*', 
        'keep *_hltMu8Ele17CaloIdTCaloIsoVLPixelMatchFilter_*_*', 
        'keep *_hltMuTrackJpsiCtfTrackCands_*_*', 
        'keep *_hltMuTrackJpsiPixelTrackCands_*_*', 
        'keep *_hltMuonCSCDigis_*_*', 
        'keep *_hltMuonCSCDigis_MuonCSCStripDigi_*', 
        'keep *_hltMuonCSCDigis_MuonCSCWireDigi_*', 
        'keep *_hltMuonDTDigis_*_*', 
        'keep *_hltMuonRPCDigis_*_*', 
        'keep *_hltOnlineBeamSpot_*_*', 
        'keep *_hltOnlinePrimaryVertices_*_*', 
        'keep *_hltOverlapFilterEle20LooseIsoPFTau20OldVersion_*_*', 
        'keep *_hltOverlapFilterIsoEle20MediumIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15IsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15MediumIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15TightIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu18LooseIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu18PFTau25TrackPt5Prong4_*_*', 
        'keep *_hltPFJetForBtag_*_*', 
        'keep *_hltPFMuonMerging_*_*', 
        'keep *_hltPFTau15TrackLooseIso_*_*', 
        'keep *_hltPFTau15Track_*_*', 
        'keep *_hltPFTau15_*_*', 
        'keep *_hltPFTau20IsoMuVertex_*_*', 
        'keep *_hltPFTau20TrackLooseIso_*_*', 
        'keep *_hltPFTau20Track_*_*', 
        'keep *_hltPFTau20_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolationProng4IsoMuVertex_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolationProng4_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolation_*_*', 
        'keep *_hltPFTau25TrackPt5_*_*', 
        'keep *_hltPFTau25_*_*', 
        'keep *_hltPFTau35TrackPt20LooseIsoProng2_*_*', 
        'keep *_hltPFTau35TrackPt20LooseIso_*_*', 
        'keep *_hltPFTau35TrackPt20_*_*', 
        'keep *_hltPFTau35Track_*_*', 
        'keep *_hltPFTau35_*_*', 
        'keep *_hltPFTauEleVertex20_*_*', 
        'keep *_hltPFTauJetTracksAssociator_*_*', 
        'keep *_hltPFTauMediumIso20TrackMediumIso_*_*', 
        'keep *_hltPFTauMediumIso20Track_*_*', 
        'keep *_hltPFTauMediumIso20_*_*', 
        'keep *_hltPFTauMediumIso35Track_*_*', 
        'keep *_hltPFTauMediumIso35_*_*', 
        'keep *_hltPFTauTagInfo_*_*', 
        'keep *_hltPFTauTightIso20TrackTightIso_*_*', 
        'keep *_hltPFTauTightIso20Track_*_*', 
        'keep *_hltPFTauTightIso20_*_*', 
        'keep *_hltPFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltParticleFlowForTaus_*_*', 
        'keep *_hltParticleFlow_*_*', 
        'keep *_hltPixelMatch3HitElectronsActivity_*_*', 
        'keep *_hltPixelMatch3HitElectronsL1Seeded_*_*', 
        'keep *_hltPixelMatchCleanElectronsL1Seeded_*_*', 
        'keep *_hltPixelMatchElectronsActivity_*_*', 
        'keep *_hltPixelMatchElectronsL1Iso_*_*', 
        'keep *_hltPixelMatchElectronsL1NonIso_*_*', 
        'keep *_hltPixelMatchElectronsL1Seeded_*_*', 
        'keep *_hltPixelTracksElectrons_*_*', 
        'keep *_hltPixelTracksMerged_*_*', 
        'keep *_hltPixelTracks_*_*', 
        'keep *_hltPixelVertices3DbbPhi_*_*', 
        'keep *_hltPixelVertices_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidateSC4_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidateSC5_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidate_*_*', 
        'keep *_hltRpcRecHits_*_*', 
        'keep *_hltSelector4CentralJetsL1FastJet_*_*', 
        'keep *_hltSelector8CentralJetsL1FastJet_*_*', 
        'keep *_hltSelectorJets20L1FastJet_*_*', 
        'keep *_hltSiPixelCluster_*_*', 
        'keep *_hltSiPixelClusters_*_*', 
        'keep *_hltSiStripClusters_*_*', 
        'keep *_hltSiStripRawToClustersFacility_*_*', 
        'keep *_hltSingleMu15L3Filtered15_*_*', 
        'keep *_hltSingleMuIsoL1s14L3IsoFiltered15eta2p1_*_*', 
        'keep *_hltSingleMuIsoL3IsoFiltered15_*_*', 
        'keep *_hltTowerMakerForAll_*_*', 
        'keep *_hltTowerMakerForMuons_*_*', 
        'keep *_hltTriggerSummaryAOD_*_*', 
        'keep *_hltTriggerSummaryRAW_*_*', 
        'keep *_hltTrimmedPixelVertices_*_*', 
        'keep *_hltVerticesL3_*_*', 
        'keep *_hltVerticesPFFilter_*_*', 
        'keep *_hltVerticesPFSelector_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_rawDataRepacker_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_virginRawDataRepacker_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep L2MuonTrajectorySeeds_hltL2MuonSeeds_*_*', 
        'keep L3MuonTrajectorySeeds_hltHIL3TrajSeedOIHit_*_*', 
        'keep L3MuonTrajectorySeeds_hltHIL3TrajectorySeed_*_*', 
        'keep L3MuonTrajectorySeeds_hltL3TrajSeedOIState_*_*', 
        'keep LumiScalerss_hltScalersRawToDigi_*_*', 
        'keep SiPixelClusteredmNewDetSetVector_hltSiPixelClusters_*_*', 
        'keep TrackCandidates_hltHIL3TrackCandidateFromL2OIHit_*_*', 
        'keep TrackCandidates_hltHIL3TrackCandidateFromL2OIState_*_*', 
        'keep TrackingRecHitsOwned_hltL3Muons_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep recoCaloJets_*_*_*', 
        'keep recoCaloMETs_*_*_*', 
        'keep recoCaloMETs_hltMet_*_*', 
        'keep recoCompositeCandidates_*_*_*', 
        'keep recoElectrons_*_*_*', 
        'keep recoIsolatedPixelTrackCandidates_*_*_*', 
        'keep recoMETs_*_*_*', 
        'keep recoPFJets_*_*_*', 
        'keep recoPFTaus_*_*_*', 
        'keep recoRecoChargedCandidates_*_*_*', 
        'keep recoRecoChargedCandidates_hltHIL3MuonCandidates_*_*', 
        'keep recoRecoChargedCandidates_hltL2MuonCandidates_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsoHLTClusterShape_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonEcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalForHE_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsoHLTClusterShape_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonEcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalForHE_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalIsol_*_*', 
        'keep recoRecoEcalCandidates_*_*_*', 
        'keep recoRecoEcalCandidates_hltL1IsoRecoEcalCandidate_*_*', 
        'keep recoRecoEcalCandidates_hltL1NonIsoRecoEcalCandidate_*_*', 
        'keep recoTrackExtras_hltHIL3MuonsOIHit_*_*', 
        'keep recoTrackExtras_hltHIL3MuonsOIState_*_*', 
        'keep recoTracks_hltHIL3MuonsOIHit_*_*', 
        'keep recoTracks_hltHIL3MuonsOIState_*_*', 
        'keep recoTracks_hltHIL3Muons_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2OIHit_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2OIState_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2_*_*', 
        'keep triggerTriggerEventWithRefs_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep triggerTriggerFilterObjectWithRefs_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'drop *_trackingtruthprod_*_*', 
        'drop *_electrontruth_*_*', 
        'keep *_mix_MergedTrackTruth_*', 
        'keep CrossingFramePlaybackInfoNew_*_*_*', 
        'keep PixelDigiSimLinkedmDetSetVector_simSiPixelDigis_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simSiStripDigis_*_*', 
        'drop *_mix_simSiPixelDigis*_*', 
        'drop *_mix_simSiStripDigis*_*', 
        'keep *_allTrackMCMatch_*_*' ) ),
    splitLevel = cms.untracked.int32(0)
)

process.RAWRECOEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring( ('drop *', 
        'drop *', 
        'keep DetIdedmEDCollection_siStripDigis_*_*', 
        'keep DetIdedmEDCollection_siPixelDigis_*_*', 
        'keep PixelFEDChanneledmNewDetSetVector_siPixelDigis_*_*', 
        'keep *_siPixelClusters_*_*', 
        'keep *_siStripClusters_*_*', 
        'keep ClusterSummary_clusterSummaryProducer_*_*', 
        'keep *_dt1DRecHits_*_*', 
        'keep *_dt4DSegments_*_*', 
        'keep *_dt1DCosmicRecHits_*_*', 
        'keep *_dt4DCosmicSegments_*_*', 
        'keep *_csc2DRecHits_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*', 
        'keep *_hbhereco_*_*', 
        'keep *_hbheprereco_*_*', 
        'keep *_hfprereco_*_*', 
        'keep *_hfreco_*_*', 
        'keep *_horeco_*_*', 
        'keep HBHERecHitsSorted_hbherecoMB_*_*', 
        'keep HORecHitsSorted_horecoMB_*_*', 
        'keep HFRecHitsSorted_hfrecoMB_*_*', 
        'keep ZDCDataFramesSorted_hcalDigis_*_*', 
        'keep ZDCDataFramesSorted_castorDigis_*_*', 
        'keep QIE10DataFrameHcalDataFrameContainer_hcalDigis_ZDC_*', 
        'keep ZDCRecHitsSorted_zdcreco_*_*', 
        'keep *_reducedHcalRecHits_*_*', 
        'keep *_castorreco_*_*', 
        'keep HcalUnpackerReport_castorDigis_*_*', 
        'keep HcalUnpackerReport_hcalDigiAlCaMB_*_*', 
        'keep HcalUnpackerReport_hcalDigis_*_*', 
        'keep *_ecalPreshowerRecHit_*_*', 
        'keep *_ecalRecHit_*_*', 
        'keep *_ecalCompactTrigPrim_*_*', 
        'keep *_ecalTPSkim_*_*', 
        'keep EBSrFlagsSorted_ecalDigis__*', 
        'keep EESrFlagsSorted_ecalDigis__*', 
        'keep *_selectDigi_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEE_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEB_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsES_*_*', 
        'keep *_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep *_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_*', 
        'keep *_particleFlowSuperClusterECAL_*_*', 
        'keep *_particleFlowSuperClusterOOTECAL_*_*', 
        'drop recoClusterShapes_*_*_*', 
        'drop recoBasicClustersToOnerecoClusterShapesAssociation_*_*_*', 
        'drop recoBasicClusters_multi5x5BasicClusters_multi5x5BarrelBasicClusters_*', 
        'drop recoSuperClusters_multi5x5SuperClusters_multi5x5BarrelSuperClusters_*', 
        'keep *_CkfElectronCandidates_*_*', 
        'keep *_GsfGlobalElectronTest_*_*', 
        'keep *_electronMergedSeeds_*_*', 
        'keep recoGsfTracks_electronGsfTracks_*_*', 
        'keep recoGsfTrackExtras_electronGsfTracks_*_*', 
        'keep recoTrackExtras_electronGsfTracks_*_*', 
        'keep TrackingRecHitsOwned_electronGsfTracks_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTrackExtras_generalTracks_*_*', 
        'keep TrackingRecHitsOwned_generalTracks_*_*', 
        'keep *_generalTracks_MVAValues_*', 
        'keep *_generalTracks_MVAVals_*', 
        'keep TrackingRecHitsOwned_extraFromSeeds_*_*', 
        'keep uints_extraFromSeeds_*_*', 
        'keep recoTracks_beamhaloTracks_*_*', 
        'keep recoTrackExtras_beamhaloTracks_*_*', 
        'keep TrackingRecHitsOwned_beamhaloTracks_*_*', 
        'keep recoTracks_conversionStepTracks_*_*', 
        'keep recoTrackExtras_conversionStepTracks_*_*', 
        'keep TrackingRecHitsOwned_conversionStepTracks_*_*', 
        'keep *_ctfPixelLess_*_*', 
        'keep *_dedxTruncated40_*_*', 
        'keep *_dedxHitInfo_*_*', 
        'keep *_dedxHarmonic2_*_*', 
        'keep *_dedxPixelHarmonic2_*_*', 
        'keep *_trackExtrapolator_*_*', 
        'keep recoTracks_cosmicDCTracks_*_*', 
        'keep recoTrackExtras_cosmicDCTracks_*_*', 
        'keep TrackingRecHitsOwned_cosmicDCTracks_*_*', 
        'keep *_ak4CaloJets_*_*', 
        'keep *_ak4PFJets_*_*', 
        'keep *_ak4PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHSSoftDrop_*_*', 
        'keep *_cmsTopTagPFJetsCHS_*_*', 
        'keep *_JetPlusTrackZSPCorJetAntiKt4_*_*', 
        'keep *_ak4TrackJets_*_*', 
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', 
        'keep *_caloTowers_*_*', 
        'keep *_towerMaker_*_*', 
        'keep *_CastorTowerReco_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertex_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertexPF_*_*', 
        'keep *_ak4JetTracksAssociatorAtCaloFace_*_*', 
        'keep *_ak4JetTracksAssociatorExplicit_*_*', 
        'keep *_ak4JetExtender_*_*', 
        'keep *_ak4JetID_*_*', 
        'keep *_ak5CastorJets_*_*', 
        'keep *_ak5CastorJetID_*_*', 
        'keep *_ak7CastorJets_*_*', 
        'keep *_ak7CastorJetID_*_*', 
        'keep *_fixedGridRhoAll_*_*', 
        'keep *_fixedGridRhoFastjetAll_*_*', 
        'keep *_fixedGridRhoFastjetAllTmp_*_*', 
        'keep *_fixedGridRhoFastjetAllCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentral_*_*', 
        'keep *_fixedGridRhoFastjetCentralCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralChargedPileUp_*_*', 
        'keep *_fixedGridRhoFastjetCentralNeutral_*_*', 
        'keep *_ak8PFJetsCHSSoftDropMass_*_*', 
        'keep recoCaloMETs_caloMet_*_*', 
        'keep recoCaloMETs_caloMetBE_*_*', 
        'keep recoCaloMETs_caloMetBEFO_*_*', 
        'keep recoCaloMETs_caloMetM_*_*', 
        'keep recoPFMETs_pfMet_*_*', 
        'keep recoPFMETs_pfChMet_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_*_*', 
        'keep recoHcalNoiseRBXs_hcalnoise_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoEcalHaloData_EcalHaloData_*_*', 
        'keep recoGlobalHaloData_GlobalHaloData_*_*', 
        'keep recoHcalHaloData_HcalHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep *_MuonSeed_*_*', 
        'keep *_ancientMuonSeed_*_*', 
        'keep *_displacedMuonSeeds_*_*', 
        'keep TrackingRecHitsOwned_globalMuons_*_*', 
        'keep TrackingRecHitsOwned_tevMuons_*_*', 
        'keep *_CosmicMuonSeed_*_*', 
        'keep recoTrackExtras_cosmicMuons_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons_*_*', 
        'keep recoTrackExtras_cosmicMuons1Leg_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons1Leg_*_*', 
        'keep recoTracks_cosmicsVetoTracks_*_*', 
        'keep *_muons_*_*', 
        'keep *_particleFlow_muons_*', 
        'drop *_muons_muons1stStep2muonsMap_*', 
        'drop recoIsoDepositedmValueMap_muons_*_*', 
        'drop doubleedmValueMap_muons_muPFIso*_*', 
        'keep recoTracks_standAloneMuons_*_*', 
        'keep recoTrackExtras_standAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_standAloneMuons_*_*', 
        'keep recoTracks_globalMuons_*_*', 
        'keep recoTrackExtras_globalMuons_*_*', 
        'keep recoTracks_tevMuons_*_*', 
        'keep recoTrackExtras_tevMuons_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_displacedTracks_*_*', 
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_*_*', 
        'keep recoTracks_displacedGlobalMuons_*_*', 
        'keep recoTrackExtras_displacedGlobalMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedGlobalMuons_*_*', 
        'keep recoTracks_cosmicMuons_*_*', 
        'keep recoMuons_muonsFromCosmics_*_*', 
        'keep recoTracks_cosmicMuons1Leg_*_*', 
        'keep recoMuons_muonsFromCosmics1Leg_*_*', 
        'keep recoTracks_refittedStandAloneMuons_*_*', 
        'keep recoTrackExtras_refittedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_refittedStandAloneMuons_*_*', 
        'keep recoTracks_displacedStandAloneMuons__*', 
        'keep recoTrackExtras_displacedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedStandAloneMuons_*_*', 
        'keep *_muIsoDepositTk_*_*', 
        'keep *_muIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muIsoDepositJets_*_*', 
        'keep *_muGlobalIsoDepositCtfTk_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muGlobalIsoDepositJets_*_*', 
        'keep *_softPFMuonsTagInfos_*_*', 
        'keep *_softPFElectronsTagInfos_*_*', 
        'keep *_softPFElectronBJetTags_*_*', 
        'keep *_softPFMuonBJetTags_*_*', 
        'keep *_pfImpactParameterTagInfos_*_*', 
        'keep *_pfTrackCountingHighEffBJetTags_*_*', 
        'keep *_pfJetProbabilityBJetTags_*_*', 
        'keep *_pfJetBProbabilityBJetTags_*_*', 
        'keep *_pfSecondaryVertexTagInfos_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderTagInfos_*_*', 
        'keep *_pfSimpleSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfCombinedSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfGhostTrackVertexTagInfos_*_*', 
        'keep *_pfGhostTrackBJetTags_*_*', 
        'keep *_pfCombinedMVAV2BJetTags_*_*', 
        'keep *_inclusiveCandidateSecondaryVertices_*_*', 
        'keep *_inclusiveCandidateSecondaryVerticesCvsL_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderCvsLTagInfos_*_*', 
        'keep *_pfCombinedCvsLJetTags_*_*', 
        'keep *_pfCombinedCvsBJetTags_*_*', 
        'keep *_pfChargeBJetTags_*_*', 
        'keep *_pfDeepCSVJetTags_*_*', 
        'keep *_pfDeepCMVAJetTags_*_*', 
        'keep recoRecoTauPiZeros_hpsPFTauProducer_pizeros_*', 
        'keep recoPFTaus_hpsPFTauProducer_*_*', 
        'keep *_hpsPFTauDiscriminationByDeadECALElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFinding_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingNewDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingOldDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseChargedIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByTightMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSum_*_*', 
        'keep *_hpsPFTauPUcorrPtSum_*_*', 
        'keep *_hpsPFTauChargedIsoPtSum_*_*', 
        'keep *_hpsPFTauTransverseImpactParameters_*_*', 
        'keep *_hpsPFTauFootprintCorrection_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeight_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6rawElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6LooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6MediumElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6TightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VTightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauChargedIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeightdR03_*_*', 
        'keep *_hpsPFTauFootprintCorrectiondR03_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalConedR03_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep  *_offlinePrimaryVertices__*', 
        'keep *_offlinePrimaryVerticesWithBS_*_*', 
        'keep *_offlinePrimaryVerticesFromCosmicTracks_*_*', 
        'keep *_nuclearInteractionMaker_*_*', 
        'keep *_generalV0Candidates_*_*', 
        'keep *_inclusiveSecondaryVertices_*_*', 
        'keep recoGsfElectronCores_gsfElectronCores_*_*', 
        'keep recoGsfElectronCores_gedGsfElectronCores_*_*', 
        'keep recoGsfElectrons_gsfElectrons_*_*', 
        'keep recoGsfElectrons_gedGsfElectrons_*_*', 
        'keep recoGsfElectronCores_uncleanedOnlyGsfElectronCores_*_*', 
        'keep recoGsfElectrons_uncleanedOnlyGsfElectrons_*_*', 
        'keep floatedmValueMap_eidRobustLoose_*_*', 
        'keep floatedmValueMap_eidRobustTight_*_*', 
        'keep floatedmValueMap_eidRobustHighEnergy_*_*', 
        'keep floatedmValueMap_eidLoose_*_*', 
        'keep floatedmValueMap_eidTight_*_*', 
        'keep *_egmGedGsfElectronPFIsolation_*_*', 
        'keep *_photonEcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonEcalPFClusterIsolationProducer_*_*', 
        'keep *_electronEcalPFClusterIsolationProducer_*_*', 
        'keep *_photonHcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonHcalPFClusterIsolationProducer_*_*', 
        'keep *_electronHcalPFClusterIsolationProducer_*_*', 
        'drop *_egmGsfElectronIDs_*_*', 
        'drop *_egmPhotonIDs_*_*', 
        'keep *_gedPhotonCore_*_*', 
        'keep *_gedPhotons_*_*', 
        'keep *_particleBasedIsolation_*_*', 
        'keep recoPhotons_mustachePhotons_*_*', 
        'keep recoPhotonCores_mustachePhotonCore_*_*', 
        'keep recoPhotons_ootPhotons_*_*', 
        'keep recoPhotonCores_ootPhotonCore_*_*', 
        'keep recoPhotons_photons_*_*', 
        'keep recoPhotonCores_photonCore_*_*', 
        'keep recoConversions_conversions_*_*', 
        'keep recoConversions_mustacheConversions_*_*', 
        'drop *_conversions_uncleanedConversions_*', 
        'drop *_gedPhotonsTmp_valMapPFEgammaCandToPhoton_*', 
        'keep recoConversions_allConversions_*_*', 
        'keep recoConversions_allConversionsOldEG_*_*', 
        'keep recoTracks_ckfOutInTracksFromConversions_*_*', 
        'keep recoTracks_ckfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfInOutTracksFromConversions_*_*', 
        'keep recoConversions_uncleanedOnlyAllConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep *_PhotonIDProdGED_*_*', 
        'keep recoRecoEcalCandidates_hfRecoEcalCandidate_*_*', 
        'keep *_hfEMClusters_*_*', 
        'keep *_gedGsfElectronCores_*_*', 
        'keep *_gedGsfElectrons_*_*', 
        'keep *_pixelTracks_*_*', 
        'keep *_pixelVertices_*_*', 
        'drop CaloTowersSorted_towerMakerPF_*_*', 
        'keep recoPFRecHits_particleFlowClusterECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHCAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterPS_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHBHE_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitPS_Cleaned_*', 
        'keep recoPFClusters_particleFlowClusterECAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHCAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHO_*_*', 
        'keep recoPFClusters_particleFlowClusterHF_*_*', 
        'keep recoPFClusters_particleFlowClusterPS_*_*', 
        'keep recoPFBlocks_particleFlowBlock_*_*', 
        'keep recoPFCandidates_particleFlowEGamma_*_*', 
        'keep recoCaloClusters_particleFlowEGamma_*_*', 
        'keep recoSuperClusters_particleFlowEGamma_*_*', 
        'keep recoConversions_particleFlowEGamma_*_*', 
        'keep recoPFCandidates_particleFlow_*_*', 
        'keep recoPFCandidates_particleFlowTmp_electrons_*', 
        'keep recoPFCandidates_particleFlowTmp_*_*', 
        'drop recoPFCandidates_particleFlowTmp__*', 
        'keep recoPFDisplacedVertexs_particleFlowDisplacedVertex_*_*', 
        'keep *_pfElectronTranslator_*_*', 
        'keep *_pfPhotonTranslator_*_*', 
        'keep *_particleFlow_electrons_*', 
        'keep *_particleFlow_photons_*', 
        'keep *_particleFlow_muons_*', 
        'keep *_trackerDrivenElectronSeeds_preid_*', 
        'keep *_particleFlowPtrs_*_*', 
        'keep *_particleFlowTmpPtrs_*_*', 
        'keep *_chargedHadronPFTrackIsolation_*_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep L1MuGMTReadoutCollection_gtDigis_*_*', 
        'keep L1GctEmCand*_gctDigis_*_*', 
        'keep L1GctJetCand*_gctDigis_*_*', 
        'keep L1GctEtHad*_gctDigis_*_*', 
        'keep L1GctEtMiss*_gctDigis_*_*', 
        'keep L1GctEtTotal*_gctDigis_*_*', 
        'keep L1GctHtMiss*_gctDigis_*_*', 
        'keep L1GctJetCounts*_gctDigis_*_*', 
        'keep L1GctHFRingEtSums*_gctDigis_*_*', 
        'keep L1GctHFBitCounts*_gctDigis_*_*', 
        'keep LumiDetails_lumiProducer_*_*', 
        'keep LumiSummary_lumiProducer_*_*', 
        'drop *_hlt*_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep L1AcceptBunchCrossings_scalersRawToDigi_*_*', 
        'keep L1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep Level1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep BeamSpotOnlines_scalersRawToDigi_*_*', 
        'keep DcsStatuss_scalersRawToDigi_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep *_tcdsDigis_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'keep *_pfIsolatedElectronsEI_*_*', 
        'keep *_pfIsolatedMuonsEI_*_*', 
        'keep recoPFJets_pfJetsEI_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTagsEI_*_*', 
        'keep recoPFTaus_pfTausEI_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByDecayModeFinding_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByIsolation_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*' ) ),
    splitLevel = cms.untracked.int32(0)
)

process.RAWRECOSIMHLTEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring( ('drop *', 
        'drop *', 
        'drop *', 
        'keep DetIdedmEDCollection_siStripDigis_*_*', 
        'keep DetIdedmEDCollection_siPixelDigis_*_*', 
        'keep PixelFEDChanneledmNewDetSetVector_siPixelDigis_*_*', 
        'keep *_siPixelClusters_*_*', 
        'keep *_siStripClusters_*_*', 
        'keep ClusterSummary_clusterSummaryProducer_*_*', 
        'keep *_dt1DRecHits_*_*', 
        'keep *_dt4DSegments_*_*', 
        'keep *_dt1DCosmicRecHits_*_*', 
        'keep *_dt4DCosmicSegments_*_*', 
        'keep *_csc2DRecHits_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*', 
        'keep *_hbhereco_*_*', 
        'keep *_hbheprereco_*_*', 
        'keep *_hfprereco_*_*', 
        'keep *_hfreco_*_*', 
        'keep *_horeco_*_*', 
        'keep HBHERecHitsSorted_hbherecoMB_*_*', 
        'keep HORecHitsSorted_horecoMB_*_*', 
        'keep HFRecHitsSorted_hfrecoMB_*_*', 
        'keep ZDCDataFramesSorted_hcalDigis_*_*', 
        'keep ZDCDataFramesSorted_castorDigis_*_*', 
        'keep QIE10DataFrameHcalDataFrameContainer_hcalDigis_ZDC_*', 
        'keep ZDCRecHitsSorted_zdcreco_*_*', 
        'keep *_reducedHcalRecHits_*_*', 
        'keep *_castorreco_*_*', 
        'keep HcalUnpackerReport_castorDigis_*_*', 
        'keep HcalUnpackerReport_hcalDigiAlCaMB_*_*', 
        'keep HcalUnpackerReport_hcalDigis_*_*', 
        'keep *_ecalPreshowerRecHit_*_*', 
        'keep *_ecalRecHit_*_*', 
        'keep *_ecalCompactTrigPrim_*_*', 
        'keep *_ecalTPSkim_*_*', 
        'keep EBSrFlagsSorted_ecalDigis__*', 
        'keep EESrFlagsSorted_ecalDigis__*', 
        'keep *_selectDigi_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEE_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEB_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsES_*_*', 
        'keep *_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep *_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_*', 
        'keep *_particleFlowSuperClusterECAL_*_*', 
        'keep *_particleFlowSuperClusterOOTECAL_*_*', 
        'drop recoClusterShapes_*_*_*', 
        'drop recoBasicClustersToOnerecoClusterShapesAssociation_*_*_*', 
        'drop recoBasicClusters_multi5x5BasicClusters_multi5x5BarrelBasicClusters_*', 
        'drop recoSuperClusters_multi5x5SuperClusters_multi5x5BarrelSuperClusters_*', 
        'keep *_CkfElectronCandidates_*_*', 
        'keep *_GsfGlobalElectronTest_*_*', 
        'keep *_electronMergedSeeds_*_*', 
        'keep recoGsfTracks_electronGsfTracks_*_*', 
        'keep recoGsfTrackExtras_electronGsfTracks_*_*', 
        'keep recoTrackExtras_electronGsfTracks_*_*', 
        'keep TrackingRecHitsOwned_electronGsfTracks_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTrackExtras_generalTracks_*_*', 
        'keep TrackingRecHitsOwned_generalTracks_*_*', 
        'keep *_generalTracks_MVAValues_*', 
        'keep *_generalTracks_MVAVals_*', 
        'keep TrackingRecHitsOwned_extraFromSeeds_*_*', 
        'keep uints_extraFromSeeds_*_*', 
        'keep recoTracks_beamhaloTracks_*_*', 
        'keep recoTrackExtras_beamhaloTracks_*_*', 
        'keep TrackingRecHitsOwned_beamhaloTracks_*_*', 
        'keep recoTracks_conversionStepTracks_*_*', 
        'keep recoTrackExtras_conversionStepTracks_*_*', 
        'keep TrackingRecHitsOwned_conversionStepTracks_*_*', 
        'keep *_ctfPixelLess_*_*', 
        'keep *_dedxTruncated40_*_*', 
        'keep *_dedxHitInfo_*_*', 
        'keep *_dedxHarmonic2_*_*', 
        'keep *_dedxPixelHarmonic2_*_*', 
        'keep *_trackExtrapolator_*_*', 
        'keep recoTracks_cosmicDCTracks_*_*', 
        'keep recoTrackExtras_cosmicDCTracks_*_*', 
        'keep TrackingRecHitsOwned_cosmicDCTracks_*_*', 
        'keep *_ak4CaloJets_*_*', 
        'keep *_ak4PFJets_*_*', 
        'keep *_ak4PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHSSoftDrop_*_*', 
        'keep *_cmsTopTagPFJetsCHS_*_*', 
        'keep *_JetPlusTrackZSPCorJetAntiKt4_*_*', 
        'keep *_ak4TrackJets_*_*', 
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', 
        'keep *_caloTowers_*_*', 
        'keep *_towerMaker_*_*', 
        'keep *_CastorTowerReco_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertex_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertexPF_*_*', 
        'keep *_ak4JetTracksAssociatorAtCaloFace_*_*', 
        'keep *_ak4JetTracksAssociatorExplicit_*_*', 
        'keep *_ak4JetExtender_*_*', 
        'keep *_ak4JetID_*_*', 
        'keep *_ak5CastorJets_*_*', 
        'keep *_ak5CastorJetID_*_*', 
        'keep *_ak7CastorJets_*_*', 
        'keep *_ak7CastorJetID_*_*', 
        'keep *_fixedGridRhoAll_*_*', 
        'keep *_fixedGridRhoFastjetAll_*_*', 
        'keep *_fixedGridRhoFastjetAllTmp_*_*', 
        'keep *_fixedGridRhoFastjetAllCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentral_*_*', 
        'keep *_fixedGridRhoFastjetCentralCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralChargedPileUp_*_*', 
        'keep *_fixedGridRhoFastjetCentralNeutral_*_*', 
        'keep *_ak8PFJetsCHSSoftDropMass_*_*', 
        'keep recoCaloMETs_caloMet_*_*', 
        'keep recoCaloMETs_caloMetBE_*_*', 
        'keep recoCaloMETs_caloMetBEFO_*_*', 
        'keep recoCaloMETs_caloMetM_*_*', 
        'keep recoPFMETs_pfMet_*_*', 
        'keep recoPFMETs_pfChMet_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_*_*', 
        'keep recoHcalNoiseRBXs_hcalnoise_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoEcalHaloData_EcalHaloData_*_*', 
        'keep recoGlobalHaloData_GlobalHaloData_*_*', 
        'keep recoHcalHaloData_HcalHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep *_MuonSeed_*_*', 
        'keep *_ancientMuonSeed_*_*', 
        'keep *_displacedMuonSeeds_*_*', 
        'keep TrackingRecHitsOwned_globalMuons_*_*', 
        'keep TrackingRecHitsOwned_tevMuons_*_*', 
        'keep *_CosmicMuonSeed_*_*', 
        'keep recoTrackExtras_cosmicMuons_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons_*_*', 
        'keep recoTrackExtras_cosmicMuons1Leg_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons1Leg_*_*', 
        'keep recoTracks_cosmicsVetoTracks_*_*', 
        'keep *_muons_*_*', 
        'keep *_particleFlow_muons_*', 
        'drop *_muons_muons1stStep2muonsMap_*', 
        'drop recoIsoDepositedmValueMap_muons_*_*', 
        'drop doubleedmValueMap_muons_muPFIso*_*', 
        'keep recoTracks_standAloneMuons_*_*', 
        'keep recoTrackExtras_standAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_standAloneMuons_*_*', 
        'keep recoTracks_globalMuons_*_*', 
        'keep recoTrackExtras_globalMuons_*_*', 
        'keep recoTracks_tevMuons_*_*', 
        'keep recoTrackExtras_tevMuons_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_displacedTracks_*_*', 
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_*_*', 
        'keep recoTracks_displacedGlobalMuons_*_*', 
        'keep recoTrackExtras_displacedGlobalMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedGlobalMuons_*_*', 
        'keep recoTracks_cosmicMuons_*_*', 
        'keep recoMuons_muonsFromCosmics_*_*', 
        'keep recoTracks_cosmicMuons1Leg_*_*', 
        'keep recoMuons_muonsFromCosmics1Leg_*_*', 
        'keep recoTracks_refittedStandAloneMuons_*_*', 
        'keep recoTrackExtras_refittedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_refittedStandAloneMuons_*_*', 
        'keep recoTracks_displacedStandAloneMuons__*', 
        'keep recoTrackExtras_displacedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedStandAloneMuons_*_*', 
        'keep *_muIsoDepositTk_*_*', 
        'keep *_muIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muIsoDepositJets_*_*', 
        'keep *_muGlobalIsoDepositCtfTk_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muGlobalIsoDepositJets_*_*', 
        'keep *_softPFMuonsTagInfos_*_*', 
        'keep *_softPFElectronsTagInfos_*_*', 
        'keep *_softPFElectronBJetTags_*_*', 
        'keep *_softPFMuonBJetTags_*_*', 
        'keep *_pfImpactParameterTagInfos_*_*', 
        'keep *_pfTrackCountingHighEffBJetTags_*_*', 
        'keep *_pfJetProbabilityBJetTags_*_*', 
        'keep *_pfJetBProbabilityBJetTags_*_*', 
        'keep *_pfSecondaryVertexTagInfos_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderTagInfos_*_*', 
        'keep *_pfSimpleSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfCombinedSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfGhostTrackVertexTagInfos_*_*', 
        'keep *_pfGhostTrackBJetTags_*_*', 
        'keep *_pfCombinedMVAV2BJetTags_*_*', 
        'keep *_inclusiveCandidateSecondaryVertices_*_*', 
        'keep *_inclusiveCandidateSecondaryVerticesCvsL_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderCvsLTagInfos_*_*', 
        'keep *_pfCombinedCvsLJetTags_*_*', 
        'keep *_pfCombinedCvsBJetTags_*_*', 
        'keep *_pfChargeBJetTags_*_*', 
        'keep *_pfDeepCSVJetTags_*_*', 
        'keep *_pfDeepCMVAJetTags_*_*', 
        'keep recoRecoTauPiZeros_hpsPFTauProducer_pizeros_*', 
        'keep recoPFTaus_hpsPFTauProducer_*_*', 
        'keep *_hpsPFTauDiscriminationByDeadECALElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFinding_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingNewDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingOldDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseChargedIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByTightMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSum_*_*', 
        'keep *_hpsPFTauPUcorrPtSum_*_*', 
        'keep *_hpsPFTauChargedIsoPtSum_*_*', 
        'keep *_hpsPFTauTransverseImpactParameters_*_*', 
        'keep *_hpsPFTauFootprintCorrection_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeight_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6rawElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6LooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6MediumElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6TightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VTightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauChargedIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeightdR03_*_*', 
        'keep *_hpsPFTauFootprintCorrectiondR03_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalConedR03_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep  *_offlinePrimaryVertices__*', 
        'keep *_offlinePrimaryVerticesWithBS_*_*', 
        'keep *_offlinePrimaryVerticesFromCosmicTracks_*_*', 
        'keep *_nuclearInteractionMaker_*_*', 
        'keep *_generalV0Candidates_*_*', 
        'keep *_inclusiveSecondaryVertices_*_*', 
        'keep recoGsfElectronCores_gsfElectronCores_*_*', 
        'keep recoGsfElectronCores_gedGsfElectronCores_*_*', 
        'keep recoGsfElectrons_gsfElectrons_*_*', 
        'keep recoGsfElectrons_gedGsfElectrons_*_*', 
        'keep recoGsfElectronCores_uncleanedOnlyGsfElectronCores_*_*', 
        'keep recoGsfElectrons_uncleanedOnlyGsfElectrons_*_*', 
        'keep floatedmValueMap_eidRobustLoose_*_*', 
        'keep floatedmValueMap_eidRobustTight_*_*', 
        'keep floatedmValueMap_eidRobustHighEnergy_*_*', 
        'keep floatedmValueMap_eidLoose_*_*', 
        'keep floatedmValueMap_eidTight_*_*', 
        'keep *_egmGedGsfElectronPFIsolation_*_*', 
        'keep *_photonEcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonEcalPFClusterIsolationProducer_*_*', 
        'keep *_electronEcalPFClusterIsolationProducer_*_*', 
        'keep *_photonHcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonHcalPFClusterIsolationProducer_*_*', 
        'keep *_electronHcalPFClusterIsolationProducer_*_*', 
        'drop *_egmGsfElectronIDs_*_*', 
        'drop *_egmPhotonIDs_*_*', 
        'keep *_gedPhotonCore_*_*', 
        'keep *_gedPhotons_*_*', 
        'keep *_particleBasedIsolation_*_*', 
        'keep recoPhotons_mustachePhotons_*_*', 
        'keep recoPhotonCores_mustachePhotonCore_*_*', 
        'keep recoPhotons_ootPhotons_*_*', 
        'keep recoPhotonCores_ootPhotonCore_*_*', 
        'keep recoPhotons_photons_*_*', 
        'keep recoPhotonCores_photonCore_*_*', 
        'keep recoConversions_conversions_*_*', 
        'keep recoConversions_mustacheConversions_*_*', 
        'drop *_conversions_uncleanedConversions_*', 
        'drop *_gedPhotonsTmp_valMapPFEgammaCandToPhoton_*', 
        'keep recoConversions_allConversions_*_*', 
        'keep recoConversions_allConversionsOldEG_*_*', 
        'keep recoTracks_ckfOutInTracksFromConversions_*_*', 
        'keep recoTracks_ckfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfInOutTracksFromConversions_*_*', 
        'keep recoConversions_uncleanedOnlyAllConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep *_PhotonIDProdGED_*_*', 
        'keep recoRecoEcalCandidates_hfRecoEcalCandidate_*_*', 
        'keep *_hfEMClusters_*_*', 
        'keep *_gedGsfElectronCores_*_*', 
        'keep *_gedGsfElectrons_*_*', 
        'keep *_pixelTracks_*_*', 
        'keep *_pixelVertices_*_*', 
        'drop CaloTowersSorted_towerMakerPF_*_*', 
        'keep recoPFRecHits_particleFlowClusterECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHCAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterPS_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHBHE_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitPS_Cleaned_*', 
        'keep recoPFClusters_particleFlowClusterECAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHCAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHO_*_*', 
        'keep recoPFClusters_particleFlowClusterHF_*_*', 
        'keep recoPFClusters_particleFlowClusterPS_*_*', 
        'keep recoPFBlocks_particleFlowBlock_*_*', 
        'keep recoPFCandidates_particleFlowEGamma_*_*', 
        'keep recoCaloClusters_particleFlowEGamma_*_*', 
        'keep recoSuperClusters_particleFlowEGamma_*_*', 
        'keep recoConversions_particleFlowEGamma_*_*', 
        'keep recoPFCandidates_particleFlow_*_*', 
        'keep recoPFCandidates_particleFlowTmp_electrons_*', 
        'keep recoPFCandidates_particleFlowTmp_*_*', 
        'drop recoPFCandidates_particleFlowTmp__*', 
        'keep recoPFDisplacedVertexs_particleFlowDisplacedVertex_*_*', 
        'keep *_pfElectronTranslator_*_*', 
        'keep *_pfPhotonTranslator_*_*', 
        'keep *_particleFlow_electrons_*', 
        'keep *_particleFlow_photons_*', 
        'keep *_particleFlow_muons_*', 
        'keep *_trackerDrivenElectronSeeds_preid_*', 
        'keep *_particleFlowPtrs_*_*', 
        'keep *_particleFlowTmpPtrs_*_*', 
        'keep *_chargedHadronPFTrackIsolation_*_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep L1MuGMTReadoutCollection_gtDigis_*_*', 
        'keep L1GctEmCand*_gctDigis_*_*', 
        'keep L1GctJetCand*_gctDigis_*_*', 
        'keep L1GctEtHad*_gctDigis_*_*', 
        'keep L1GctEtMiss*_gctDigis_*_*', 
        'keep L1GctEtTotal*_gctDigis_*_*', 
        'keep L1GctHtMiss*_gctDigis_*_*', 
        'keep L1GctJetCounts*_gctDigis_*_*', 
        'keep L1GctHFRingEtSums*_gctDigis_*_*', 
        'keep L1GctHFBitCounts*_gctDigis_*_*', 
        'keep LumiDetails_lumiProducer_*_*', 
        'keep LumiSummary_lumiProducer_*_*', 
        'drop *_hlt*_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep L1AcceptBunchCrossings_scalersRawToDigi_*_*', 
        'keep L1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep Level1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep BeamSpotOnlines_scalersRawToDigi_*_*', 
        'keep DcsStatuss_scalersRawToDigi_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep *_tcdsDigis_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'keep *_pfIsolatedElectronsEI_*_*', 
        'keep *_pfIsolatedMuonsEI_*_*', 
        'keep recoPFJets_pfJetsEI_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTagsEI_*_*', 
        'keep recoPFTaus_pfTausEI_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByDecayModeFinding_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByIsolation_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep *_ak4GenJets_*_*', 
        'keep *_ak8GenJets_*_*', 
        'keep *_ak4GenJetsNoNu_*_*', 
        'keep *_ak8GenJetsNoNu_*_*', 
        'keep *_genParticle_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep SimTracks_g4SimHits_*_*', 
        'keep SimVertexs_g4SimHits_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*', 
        'keep *_muonSimClassifier_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'drop *_hlt*_*_*', 
        'keep *_hltAK4CaloJetsCorrectedIDPassed_*_*', 
        'keep *_hltAK4CaloJetsIDPassed_*_*', 
        'keep *_hltAK4CaloJets_*_*', 
        'keep *_hltAK4PFJetsCorrected_*_*', 
        'keep *_hltAK4PFJetsForTaus_*_*', 
        'keep *_hltAK4PFJets_*_*', 
        'keep *_hltAlCaEtaEBRechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaEtaEBRechitsToDigis_*_*', 
        'keep *_hltAlCaEtaEERechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaEtaEERechitsToDigis_*_*', 
        'keep *_hltAlCaEtaRecHitsFilterEBonly_*_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegionalLowPU_etaEcalRecHitsES_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegional_etaEcalRecHitsES_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonly_*_*', 
        'keep *_hltAlCaEtaRecHitsFilter_*_*', 
        'keep *_hltAlCaPhiSymStream_*_*', 
        'keep *_hltAlCaPi0EBRechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaPi0EBRechitsToDigis_*_*', 
        'keep *_hltAlCaPi0EERechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaPi0EERechitsToDigis_*_*', 
        'keep *_hltAlCaPi0RecHitsFilterEBonly_*_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegionalLowPU_pi0EcalRecHitsES_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegional_pi0EcalRecHitsES_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonly_*_*', 
        'keep *_hltAlCaPi0RecHitsFilter_*_*', 
        'keep *_hltBLifetimeL25AssociatorbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL25BJetTagsbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL25TagInfosbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3AssociatorbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3BJetTagsbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3TagInfosbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeRegionalCtfWithMaterialTracksbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBSoftMuonDiJet110Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet110Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet20Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet20Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet40Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet40Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet70Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet70Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonMu5L3_*_*', 
        'keep *_hltCSVJetTagSingleTopEle27_*_*', 
        'keep *_hltCSVJetTagSingleTopIsoMu24_*_*', 
        'keep *_hltCaloJetCorrectedRegional_*_*', 
        'keep *_hltCaloJetCorrected_*_*', 
        'keep *_hltCaloJetL1FastJetCorrected_*_*', 
        'keep *_hltCaloStage2Digis_*_*', 
        'keep *_hltCleanedCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCleanedHiCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCombinedSecondaryVertexBJetTagsCalo_*_*', 
        'keep *_hltCombinedSecondaryVertexBJetTagsPF_*_*', 
        'keep *_hltConvPFTausTightIsoTrackFindingIsolation_*_*', 
        'keep *_hltConvPFTausTightIsoTrackFinding_*_*', 
        'keep *_hltConvPFTausTightIsoTrackPt5Isolation_*_*', 
        'keep *_hltConvPFTausTightIsoTrackPt5_*_*', 
        'keep *_hltConvPFTausTightIso_*_*', 
        'keep *_hltConvPFTausTrackFindingLooseIsolation_*_*', 
        'keep *_hltConvPFTausTrackFinding_*_*', 
        'keep *_hltConvPFTaus_*_*', 
        'keep *_hltCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCorrectedIslandEndcapSuperClustersHI_*_*', 
        'keep *_hltCsc2DRecHits_*_*', 
        'keep *_hltCscSegments_*_*', 
        'keep *_hltDisplacedhltIter4PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolationProng4L1HLTMatched_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolationProng4_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolation_*_*', 
        'keep *_hltDoublePFTau25TrackPt5_*_*', 
        'keep *_hltDoublePFTau25_*_*', 
        'keep *_hltDoublePFTauTightIso45Track5_*_*', 
        'keep *_hltDoublePFTauTightIso45Track_*_*', 
        'keep *_hltDt4DSegments_*_*', 
        'keep *_hltEcalPhiSymFilter_*_*', 
        'keep *_hltEcalRecHitAll_*_*', 
        'keep *_hltEcalRecHit_*_*', 
        'keep *_hltEgammaCandidates_*_*', 
        'keep *_hltEgammaGsfElectrons_*_*', 
        'keep *_hltEgammaGsfTracks_*_*', 
        'keep *_hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1SingleEG18orL1SingleEG20_*_*', 
        'keep *_hltEle20CaloIdVTTrkIdTDphiFilter_*_*', 
        'keep *_hltEle27WP85PixelMatchElectronsL1Seeded_*_*', 
        'keep *_hltElectronsVertex_*_*', 
        'keep *_hltFEDSelectorLumiPixels_*_*', 
        'keep *_hltFastPVPixelTracksMerger_*_*', 
        'keep *_hltFastPVPixelTracksRecover_*_*', 
        'keep *_hltFastPVPixelTracks_*_*', 
        'keep *_hltFastPVPixelVertices3D_*_*', 
        'keep *_hltFastPVPixelVertices_*_*', 
        'keep *_hltFastPixelBLifetimeL3TagInfos_*_*', 
        'keep *_hltFastPrimaryVertex_*_*', 
        'keep *_hltFilterDoubleIsoPFTau45Trk5LeadTrack5IsolationL1HLTMatched_*_*', 
        'keep *_hltFilterL2EtCutDoublePFIsoTau45Trk5_*_*', 
        'keep *_hltFilterL2EtCutSingleIsoPFTau35Trk20MET70_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20LeadTrackPt20_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20MET60LeadTrack20IsolationL1HLTMatched_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20MET70LeadTrack20IsolationL1HLTMatched_*_*', 
        'keep *_hltGmtStage2Digis_*_*', 
        'keep *_hltGoodOnlinePVs_*_*', 
        'keep *_hltGtStage2Digis_*_*', 
        'keep *_hltHICaloJetCorrected_*_*', 
        'keep *_hltHICaloJetIDPassed_*_*', 
        'keep *_hltHIGoodLooseTracks_*_*', 
        'keep *_hltHIPixel3PrimTracks_*_*', 
        'keep *_hltHISelectedVertex_*_*', 
        'keep *_hltHISiPixelClusters_*_*', 
        'keep *_hltHITIPTCorrectorHB_*_*', 
        'keep *_hltHITIPTCorrectorHE_*_*', 
        'keep *_hltHbhereco_*_*', 
        'keep *_hltHfreco_*_*', 
        'keep *_hltHiCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltHiCorrectedIslandEndcapSuperClustersHI_*_*', 
        'keep *_hltHiIslandSuperClustersHI_*_*', 
        'keep *_hltHoreco_*_*', 
        'keep *_hltIsolPixelTrackProdHB_*_*', 
        'keep *_hltIsolPixelTrackProdHE_*_*', 
        'keep *_hltIter0ElectronsTrackSelectionHighPurity_*_*', 
        'keep *_hltIter0HighPtTkMuPixelTracks_*_*', 
        'keep *_hltIter0HighPtTkMuTrackSelectionHighPurity_*_*', 
        'keep *_hltIter0PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter1Merged_*_*', 
        'keep *_hltIter1PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter2HighPtTkMuMerged_*_*', 
        'keep *_hltIter2HighPtTkMuTrackSelectionHighPurity_*_*', 
        'keep *_hltIter2MergedForBTag_*_*', 
        'keep *_hltIter2MergedForElectrons_*_*', 
        'keep *_hltIter2Merged_*_*', 
        'keep *_hltIter2PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter3Merged_*_*', 
        'keep *_hltIter4Merged_*_*', 
        'keep *_hltIterativeCone5PileupSubtractionCaloJets_*_*', 
        'keep *_hltL1HLTSingleIsoPFTau35Trk20Met60JetsMatch_*_*', 
        'keep *_hltL1IsoElectronTrackIsol_*_*', 
        'keep *_hltL1NonIsoElectronTrackIsol_*_*', 
        'keep *_hltL1SeededRecoEcalCandidate_*_*', 
        'keep *_hltL1sDoubleTauJet44Eta2p17orDoubleJet64Central_*_*', 
        'keep *_hltL1sDoubleTauJet44erorDoubleJetC64_*_*', 
        'keep *_hltL1sL1EG18er_*_*', 
        'keep *_hltL1sL1ETM36ORETM40_*_*', 
        'keep *_hltL1sL1Jet52ETM30_*_*', 
        'keep *_hltL1sL1SingleEG12_*_*', 
        'keep *_hltL1sL1SingleEG15_*_*', 
        'keep *_hltL1sL1SingleEG18orL1SingleEG20_*_*', 
        'keep *_hltL1sL1SingleMu10_*_*', 
        'keep *_hltL1sL1SingleMu14Eta2p1_*_*', 
        'keep *_hltL1sMu16Eta2p1_*_*', 
        'keep *_hltL2MuonCandidatesNoVtx_*_*', 
        'keep *_hltL2MuonCandidates_*_*', 
        'keep *_hltL2MuonSeeds_*_*', 
        'keep *_hltL2Muons_*_*', 
        'keep *_hltL2TauJets_*_*', 
        'keep *_hltL3MuonCandidates_*_*', 
        'keep *_hltL3MuonsIOHit_*_*', 
        'keep *_hltL3MuonsLinksCombination_*_*', 
        'keep *_hltL3MuonsOIHit_*_*', 
        'keep *_hltL3MuonsOIState_*_*', 
        'keep *_hltL3Muons_*_*', 
        'keep *_hltL3NoFiltersNoVtxMuonCandidates_*_*', 
        'keep *_hltL3NoFiltersNoVtxMuons_*_*', 
        'keep *_hltL3SecondaryVertexTagInfos_*_*', 
        'keep *_hltL3TkFromL2OICombination_*_*', 
        'keep *_hltL3TkTracksFromL2IOHit_*_*', 
        'keep *_hltL3TkTracksFromL2OIHit_*_*', 
        'keep *_hltL3TkTracksFromL2OIState_*_*', 
        'keep *_hltL3TkTracksFromL2_*_*', 
        'keep *_hltL3TrackCandidateFromL2IOHit_*_*', 
        'keep *_hltL3TrackCandidateFromL2OIHit_*_*', 
        'keep *_hltL3TrackCandidateFromL2OIState_*_*', 
        'keep *_hltL3TrajSeedIOHit_*_*', 
        'keep *_hltL3TrajSeedOIHit_*_*', 
        'keep *_hltL3TrajSeedOIState_*_*', 
        'keep *_hltL3TrajectorySeed_*_*', 
        'keep *_hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f18QL3crIsoRhoFiltered0p15_*_*', 
        'keep *_hltLeadingCentralJets30SingleTopEle27_*_*', 
        'keep *_hltLeadingCentralJets30SingleTopIsoMu24_*_*', 
        'keep *_hltMergedTracks_*_*', 
        'keep *_hltMet_*_*', 
        'keep *_hltMu8Ele17CaloIdTCaloIsoVLPixelMatchFilter_*_*', 
        'keep *_hltMuTrackJpsiCtfTrackCands_*_*', 
        'keep *_hltMuTrackJpsiPixelTrackCands_*_*', 
        'keep *_hltMuonCSCDigis_*_*', 
        'keep *_hltMuonCSCDigis_MuonCSCStripDigi_*', 
        'keep *_hltMuonCSCDigis_MuonCSCWireDigi_*', 
        'keep *_hltMuonDTDigis_*_*', 
        'keep *_hltMuonRPCDigis_*_*', 
        'keep *_hltOnlineBeamSpot_*_*', 
        'keep *_hltOnlinePrimaryVertices_*_*', 
        'keep *_hltOverlapFilterEle20LooseIsoPFTau20OldVersion_*_*', 
        'keep *_hltOverlapFilterIsoEle20MediumIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15IsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15MediumIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15TightIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu18LooseIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu18PFTau25TrackPt5Prong4_*_*', 
        'keep *_hltPFJetForBtag_*_*', 
        'keep *_hltPFMuonMerging_*_*', 
        'keep *_hltPFTau15TrackLooseIso_*_*', 
        'keep *_hltPFTau15Track_*_*', 
        'keep *_hltPFTau15_*_*', 
        'keep *_hltPFTau20IsoMuVertex_*_*', 
        'keep *_hltPFTau20TrackLooseIso_*_*', 
        'keep *_hltPFTau20Track_*_*', 
        'keep *_hltPFTau20_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolationProng4IsoMuVertex_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolationProng4_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolation_*_*', 
        'keep *_hltPFTau25TrackPt5_*_*', 
        'keep *_hltPFTau25_*_*', 
        'keep *_hltPFTau35TrackPt20LooseIsoProng2_*_*', 
        'keep *_hltPFTau35TrackPt20LooseIso_*_*', 
        'keep *_hltPFTau35TrackPt20_*_*', 
        'keep *_hltPFTau35Track_*_*', 
        'keep *_hltPFTau35_*_*', 
        'keep *_hltPFTauEleVertex20_*_*', 
        'keep *_hltPFTauJetTracksAssociator_*_*', 
        'keep *_hltPFTauMediumIso20TrackMediumIso_*_*', 
        'keep *_hltPFTauMediumIso20Track_*_*', 
        'keep *_hltPFTauMediumIso20_*_*', 
        'keep *_hltPFTauMediumIso35Track_*_*', 
        'keep *_hltPFTauMediumIso35_*_*', 
        'keep *_hltPFTauTagInfo_*_*', 
        'keep *_hltPFTauTightIso20TrackTightIso_*_*', 
        'keep *_hltPFTauTightIso20Track_*_*', 
        'keep *_hltPFTauTightIso20_*_*', 
        'keep *_hltPFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltParticleFlowForTaus_*_*', 
        'keep *_hltParticleFlow_*_*', 
        'keep *_hltPixelMatch3HitElectronsActivity_*_*', 
        'keep *_hltPixelMatch3HitElectronsL1Seeded_*_*', 
        'keep *_hltPixelMatchCleanElectronsL1Seeded_*_*', 
        'keep *_hltPixelMatchElectronsActivity_*_*', 
        'keep *_hltPixelMatchElectronsL1Iso_*_*', 
        'keep *_hltPixelMatchElectronsL1NonIso_*_*', 
        'keep *_hltPixelMatchElectronsL1Seeded_*_*', 
        'keep *_hltPixelTracksElectrons_*_*', 
        'keep *_hltPixelTracksMerged_*_*', 
        'keep *_hltPixelTracks_*_*', 
        'keep *_hltPixelVertices3DbbPhi_*_*', 
        'keep *_hltPixelVertices_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidateSC4_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidateSC5_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidate_*_*', 
        'keep *_hltRpcRecHits_*_*', 
        'keep *_hltSelector4CentralJetsL1FastJet_*_*', 
        'keep *_hltSelector8CentralJetsL1FastJet_*_*', 
        'keep *_hltSelectorJets20L1FastJet_*_*', 
        'keep *_hltSiPixelCluster_*_*', 
        'keep *_hltSiPixelClusters_*_*', 
        'keep *_hltSiStripClusters_*_*', 
        'keep *_hltSiStripRawToClustersFacility_*_*', 
        'keep *_hltSingleMu15L3Filtered15_*_*', 
        'keep *_hltSingleMuIsoL1s14L3IsoFiltered15eta2p1_*_*', 
        'keep *_hltSingleMuIsoL3IsoFiltered15_*_*', 
        'keep *_hltTowerMakerForAll_*_*', 
        'keep *_hltTowerMakerForMuons_*_*', 
        'keep *_hltTriggerSummaryAOD_*_*', 
        'keep *_hltTriggerSummaryRAW_*_*', 
        'keep *_hltTrimmedPixelVertices_*_*', 
        'keep *_hltVerticesL3_*_*', 
        'keep *_hltVerticesPFFilter_*_*', 
        'keep *_hltVerticesPFSelector_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_rawDataRepacker_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_virginRawDataRepacker_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep L2MuonTrajectorySeeds_hltL2MuonSeeds_*_*', 
        'keep L3MuonTrajectorySeeds_hltHIL3TrajSeedOIHit_*_*', 
        'keep L3MuonTrajectorySeeds_hltHIL3TrajectorySeed_*_*', 
        'keep L3MuonTrajectorySeeds_hltL3TrajSeedOIState_*_*', 
        'keep LumiScalerss_hltScalersRawToDigi_*_*', 
        'keep SiPixelClusteredmNewDetSetVector_hltSiPixelClusters_*_*', 
        'keep TrackCandidates_hltHIL3TrackCandidateFromL2OIHit_*_*', 
        'keep TrackCandidates_hltHIL3TrackCandidateFromL2OIState_*_*', 
        'keep TrackingRecHitsOwned_hltL3Muons_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep recoCaloJets_*_*_*', 
        'keep recoCaloMETs_*_*_*', 
        'keep recoCaloMETs_hltMet_*_*', 
        'keep recoCompositeCandidates_*_*_*', 
        'keep recoElectrons_*_*_*', 
        'keep recoIsolatedPixelTrackCandidates_*_*_*', 
        'keep recoMETs_*_*_*', 
        'keep recoPFJets_*_*_*', 
        'keep recoPFTaus_*_*_*', 
        'keep recoRecoChargedCandidates_*_*_*', 
        'keep recoRecoChargedCandidates_hltHIL3MuonCandidates_*_*', 
        'keep recoRecoChargedCandidates_hltL2MuonCandidates_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsoHLTClusterShape_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonEcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalForHE_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsoHLTClusterShape_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonEcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalForHE_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalIsol_*_*', 
        'keep recoRecoEcalCandidates_*_*_*', 
        'keep recoRecoEcalCandidates_hltL1IsoRecoEcalCandidate_*_*', 
        'keep recoRecoEcalCandidates_hltL1NonIsoRecoEcalCandidate_*_*', 
        'keep recoTrackExtras_hltHIL3MuonsOIHit_*_*', 
        'keep recoTrackExtras_hltHIL3MuonsOIState_*_*', 
        'keep recoTracks_hltHIL3MuonsOIHit_*_*', 
        'keep recoTracks_hltHIL3MuonsOIState_*_*', 
        'keep recoTracks_hltHIL3Muons_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2OIHit_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2OIState_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2_*_*', 
        'keep triggerTriggerEventWithRefs_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep triggerTriggerFilterObjectWithRefs_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*' ) ),
    splitLevel = cms.untracked.int32(0)
)

process.RAWSIMEventContent = cms.PSet(
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    eventAutoFlushCompressedSize = cms.untracked.int32(20971520),
    outputCommands = cms.untracked.vstring('drop *', 
        'drop *', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep *_g4SimHits_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep CSCDetIdCSCComparatorDigiMuonDigiCollection_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*', 
        'keep EBSrFlagsSorted_simEcalDigis_*_*', 
        'keep EESrFlagsSorted_simEcalDigis_*_*', 
        'keep CrossingFramePlaybackInfoNew_*_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int6stdbitsetstdpairs_*_AffectedAPVList_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'keep *_mix_MergedTrackTruth_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep recoGenJets_*_*_*', 
        'keep *_genParticle_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep *_MEtoEDMConverter_*_*', 
        'keep *_randomEngineStateProducer_*_*', 
        'keep *_logErrorHarvester_*_*'),
    splitLevel = cms.untracked.int32(0)
)

process.RAWSIMHLTEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring( ('drop *', 
        'drop *', 
        'drop *', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep *_g4SimHits_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep CSCDetIdCSCComparatorDigiMuonDigiCollection_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*', 
        'keep EBSrFlagsSorted_simEcalDigis_*_*', 
        'keep EESrFlagsSorted_simEcalDigis_*_*', 
        'keep CrossingFramePlaybackInfoNew_*_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int6stdbitsetstdpairs_*_AffectedAPVList_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'keep *_mix_MergedTrackTruth_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep recoGenJets_*_*_*', 
        'keep *_genParticle_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep *_MEtoEDMConverter_*_*', 
        'keep *_randomEngineStateProducer_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'drop *_hlt*_*_*', 
        'keep *_hltAK4CaloJetsCorrectedIDPassed_*_*', 
        'keep *_hltAK4CaloJetsIDPassed_*_*', 
        'keep *_hltAK4CaloJets_*_*', 
        'keep *_hltAK4PFJetsCorrected_*_*', 
        'keep *_hltAK4PFJetsForTaus_*_*', 
        'keep *_hltAK4PFJets_*_*', 
        'keep *_hltAlCaEtaEBRechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaEtaEBRechitsToDigis_*_*', 
        'keep *_hltAlCaEtaEERechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaEtaEERechitsToDigis_*_*', 
        'keep *_hltAlCaEtaRecHitsFilterEBonly_*_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegionalLowPU_etaEcalRecHitsES_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonlyRegional_etaEcalRecHitsES_*', 
        'keep *_hltAlCaEtaRecHitsFilterEEonly_*_*', 
        'keep *_hltAlCaEtaRecHitsFilter_*_*', 
        'keep *_hltAlCaPhiSymStream_*_*', 
        'keep *_hltAlCaPi0EBRechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaPi0EBRechitsToDigis_*_*', 
        'keep *_hltAlCaPi0EERechitsToDigisLowPU_*_*', 
        'keep *_hltAlCaPi0EERechitsToDigis_*_*', 
        'keep *_hltAlCaPi0RecHitsFilterEBonly_*_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegionalLowPU_pi0EcalRecHitsES_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonlyRegional_pi0EcalRecHitsES_*', 
        'keep *_hltAlCaPi0RecHitsFilterEEonly_*_*', 
        'keep *_hltAlCaPi0RecHitsFilter_*_*', 
        'keep *_hltBLifetimeL25AssociatorbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL25BJetTagsbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL25TagInfosbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3AssociatorbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3BJetTagsbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeL3TagInfosbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBLifetimeRegionalCtfWithMaterialTracksbbPhiL1FastJetFastPV_*_*', 
        'keep *_hltBSoftMuonDiJet110Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet110Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet20Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet20Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet40Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet40Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonDiJet70Mu5L3FilterByDR_*_*', 
        'keep *_hltBSoftMuonDiJet70Mu5SelL3BJetTagsByDR_*_*', 
        'keep *_hltBSoftMuonMu5L3_*_*', 
        'keep *_hltCSVJetTagSingleTopEle27_*_*', 
        'keep *_hltCSVJetTagSingleTopIsoMu24_*_*', 
        'keep *_hltCaloJetCorrectedRegional_*_*', 
        'keep *_hltCaloJetCorrected_*_*', 
        'keep *_hltCaloJetL1FastJetCorrected_*_*', 
        'keep *_hltCaloStage2Digis_*_*', 
        'keep *_hltCleanedCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCleanedHiCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCombinedSecondaryVertexBJetTagsCalo_*_*', 
        'keep *_hltCombinedSecondaryVertexBJetTagsPF_*_*', 
        'keep *_hltConvPFTausTightIsoTrackFindingIsolation_*_*', 
        'keep *_hltConvPFTausTightIsoTrackFinding_*_*', 
        'keep *_hltConvPFTausTightIsoTrackPt5Isolation_*_*', 
        'keep *_hltConvPFTausTightIsoTrackPt5_*_*', 
        'keep *_hltConvPFTausTightIso_*_*', 
        'keep *_hltConvPFTausTrackFindingLooseIsolation_*_*', 
        'keep *_hltConvPFTausTrackFinding_*_*', 
        'keep *_hltConvPFTaus_*_*', 
        'keep *_hltCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltCorrectedIslandEndcapSuperClustersHI_*_*', 
        'keep *_hltCsc2DRecHits_*_*', 
        'keep *_hltCscSegments_*_*', 
        'keep *_hltDisplacedhltIter4PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolationProng4L1HLTMatched_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolationProng4_*_*', 
        'keep *_hltDoublePFTau25TrackPt5MediumIsolation_*_*', 
        'keep *_hltDoublePFTau25TrackPt5_*_*', 
        'keep *_hltDoublePFTau25_*_*', 
        'keep *_hltDoublePFTauTightIso45Track5_*_*', 
        'keep *_hltDoublePFTauTightIso45Track_*_*', 
        'keep *_hltDt4DSegments_*_*', 
        'keep *_hltEcalPhiSymFilter_*_*', 
        'keep *_hltEcalRecHitAll_*_*', 
        'keep *_hltEcalRecHit_*_*', 
        'keep *_hltEgammaCandidates_*_*', 
        'keep *_hltEgammaGsfElectrons_*_*', 
        'keep *_hltEgammaGsfTracks_*_*', 
        'keep *_hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1SingleEG18orL1SingleEG20_*_*', 
        'keep *_hltEle20CaloIdVTTrkIdTDphiFilter_*_*', 
        'keep *_hltEle27WP85PixelMatchElectronsL1Seeded_*_*', 
        'keep *_hltElectronsVertex_*_*', 
        'keep *_hltFEDSelectorLumiPixels_*_*', 
        'keep *_hltFastPVPixelTracksMerger_*_*', 
        'keep *_hltFastPVPixelTracksRecover_*_*', 
        'keep *_hltFastPVPixelTracks_*_*', 
        'keep *_hltFastPVPixelVertices3D_*_*', 
        'keep *_hltFastPVPixelVertices_*_*', 
        'keep *_hltFastPixelBLifetimeL3TagInfos_*_*', 
        'keep *_hltFastPrimaryVertex_*_*', 
        'keep *_hltFilterDoubleIsoPFTau45Trk5LeadTrack5IsolationL1HLTMatched_*_*', 
        'keep *_hltFilterL2EtCutDoublePFIsoTau45Trk5_*_*', 
        'keep *_hltFilterL2EtCutSingleIsoPFTau35Trk20MET70_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20LeadTrackPt20_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20MET60LeadTrack20IsolationL1HLTMatched_*_*', 
        'keep *_hltFilterSingleIsoPFTau35Trk20MET70LeadTrack20IsolationL1HLTMatched_*_*', 
        'keep *_hltGmtStage2Digis_*_*', 
        'keep *_hltGoodOnlinePVs_*_*', 
        'keep *_hltGtStage2Digis_*_*', 
        'keep *_hltHICaloJetCorrected_*_*', 
        'keep *_hltHICaloJetIDPassed_*_*', 
        'keep *_hltHIGoodLooseTracks_*_*', 
        'keep *_hltHIPixel3PrimTracks_*_*', 
        'keep *_hltHISelectedVertex_*_*', 
        'keep *_hltHISiPixelClusters_*_*', 
        'keep *_hltHITIPTCorrectorHB_*_*', 
        'keep *_hltHITIPTCorrectorHE_*_*', 
        'keep *_hltHbhereco_*_*', 
        'keep *_hltHfreco_*_*', 
        'keep *_hltHiCorrectedIslandBarrelSuperClustersHI_*_*', 
        'keep *_hltHiCorrectedIslandEndcapSuperClustersHI_*_*', 
        'keep *_hltHiIslandSuperClustersHI_*_*', 
        'keep *_hltHoreco_*_*', 
        'keep *_hltIsolPixelTrackProdHB_*_*', 
        'keep *_hltIsolPixelTrackProdHE_*_*', 
        'keep *_hltIter0ElectronsTrackSelectionHighPurity_*_*', 
        'keep *_hltIter0HighPtTkMuPixelTracks_*_*', 
        'keep *_hltIter0HighPtTkMuTrackSelectionHighPurity_*_*', 
        'keep *_hltIter0PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter1Merged_*_*', 
        'keep *_hltIter1PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter2HighPtTkMuMerged_*_*', 
        'keep *_hltIter2HighPtTkMuTrackSelectionHighPurity_*_*', 
        'keep *_hltIter2MergedForBTag_*_*', 
        'keep *_hltIter2MergedForElectrons_*_*', 
        'keep *_hltIter2Merged_*_*', 
        'keep *_hltIter2PFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltIter3Merged_*_*', 
        'keep *_hltIter4Merged_*_*', 
        'keep *_hltIterativeCone5PileupSubtractionCaloJets_*_*', 
        'keep *_hltL1HLTSingleIsoPFTau35Trk20Met60JetsMatch_*_*', 
        'keep *_hltL1IsoElectronTrackIsol_*_*', 
        'keep *_hltL1NonIsoElectronTrackIsol_*_*', 
        'keep *_hltL1SeededRecoEcalCandidate_*_*', 
        'keep *_hltL1sDoubleTauJet44Eta2p17orDoubleJet64Central_*_*', 
        'keep *_hltL1sDoubleTauJet44erorDoubleJetC64_*_*', 
        'keep *_hltL1sL1EG18er_*_*', 
        'keep *_hltL1sL1ETM36ORETM40_*_*', 
        'keep *_hltL1sL1Jet52ETM30_*_*', 
        'keep *_hltL1sL1SingleEG12_*_*', 
        'keep *_hltL1sL1SingleEG15_*_*', 
        'keep *_hltL1sL1SingleEG18orL1SingleEG20_*_*', 
        'keep *_hltL1sL1SingleMu10_*_*', 
        'keep *_hltL1sL1SingleMu14Eta2p1_*_*', 
        'keep *_hltL1sMu16Eta2p1_*_*', 
        'keep *_hltL2MuonCandidatesNoVtx_*_*', 
        'keep *_hltL2MuonCandidates_*_*', 
        'keep *_hltL2MuonSeeds_*_*', 
        'keep *_hltL2Muons_*_*', 
        'keep *_hltL2TauJets_*_*', 
        'keep *_hltL3MuonCandidates_*_*', 
        'keep *_hltL3MuonsIOHit_*_*', 
        'keep *_hltL3MuonsLinksCombination_*_*', 
        'keep *_hltL3MuonsOIHit_*_*', 
        'keep *_hltL3MuonsOIState_*_*', 
        'keep *_hltL3Muons_*_*', 
        'keep *_hltL3NoFiltersNoVtxMuonCandidates_*_*', 
        'keep *_hltL3NoFiltersNoVtxMuons_*_*', 
        'keep *_hltL3SecondaryVertexTagInfos_*_*', 
        'keep *_hltL3TkFromL2OICombination_*_*', 
        'keep *_hltL3TkTracksFromL2IOHit_*_*', 
        'keep *_hltL3TkTracksFromL2OIHit_*_*', 
        'keep *_hltL3TkTracksFromL2OIState_*_*', 
        'keep *_hltL3TkTracksFromL2_*_*', 
        'keep *_hltL3TrackCandidateFromL2IOHit_*_*', 
        'keep *_hltL3TrackCandidateFromL2OIHit_*_*', 
        'keep *_hltL3TrackCandidateFromL2OIState_*_*', 
        'keep *_hltL3TrajSeedIOHit_*_*', 
        'keep *_hltL3TrajSeedOIHit_*_*', 
        'keep *_hltL3TrajSeedOIState_*_*', 
        'keep *_hltL3TrajectorySeed_*_*', 
        'keep *_hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f18QL3crIsoRhoFiltered0p15_*_*', 
        'keep *_hltLeadingCentralJets30SingleTopEle27_*_*', 
        'keep *_hltLeadingCentralJets30SingleTopIsoMu24_*_*', 
        'keep *_hltMergedTracks_*_*', 
        'keep *_hltMet_*_*', 
        'keep *_hltMu8Ele17CaloIdTCaloIsoVLPixelMatchFilter_*_*', 
        'keep *_hltMuTrackJpsiCtfTrackCands_*_*', 
        'keep *_hltMuTrackJpsiPixelTrackCands_*_*', 
        'keep *_hltMuonCSCDigis_*_*', 
        'keep *_hltMuonCSCDigis_MuonCSCStripDigi_*', 
        'keep *_hltMuonCSCDigis_MuonCSCWireDigi_*', 
        'keep *_hltMuonDTDigis_*_*', 
        'keep *_hltMuonRPCDigis_*_*', 
        'keep *_hltOnlineBeamSpot_*_*', 
        'keep *_hltOnlinePrimaryVertices_*_*', 
        'keep *_hltOverlapFilterEle20LooseIsoPFTau20OldVersion_*_*', 
        'keep *_hltOverlapFilterIsoEle20MediumIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15IsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15MediumIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu15TightIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu18LooseIsoPFTau20_*_*', 
        'keep *_hltOverlapFilterIsoMu18PFTau25TrackPt5Prong4_*_*', 
        'keep *_hltPFJetForBtag_*_*', 
        'keep *_hltPFMuonMerging_*_*', 
        'keep *_hltPFTau15TrackLooseIso_*_*', 
        'keep *_hltPFTau15Track_*_*', 
        'keep *_hltPFTau15_*_*', 
        'keep *_hltPFTau20IsoMuVertex_*_*', 
        'keep *_hltPFTau20TrackLooseIso_*_*', 
        'keep *_hltPFTau20Track_*_*', 
        'keep *_hltPFTau20_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolationProng4IsoMuVertex_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolationProng4_*_*', 
        'keep *_hltPFTau25TrackPt5MediumIsolation_*_*', 
        'keep *_hltPFTau25TrackPt5_*_*', 
        'keep *_hltPFTau25_*_*', 
        'keep *_hltPFTau35TrackPt20LooseIsoProng2_*_*', 
        'keep *_hltPFTau35TrackPt20LooseIso_*_*', 
        'keep *_hltPFTau35TrackPt20_*_*', 
        'keep *_hltPFTau35Track_*_*', 
        'keep *_hltPFTau35_*_*', 
        'keep *_hltPFTauEleVertex20_*_*', 
        'keep *_hltPFTauJetTracksAssociator_*_*', 
        'keep *_hltPFTauMediumIso20TrackMediumIso_*_*', 
        'keep *_hltPFTauMediumIso20Track_*_*', 
        'keep *_hltPFTauMediumIso20_*_*', 
        'keep *_hltPFTauMediumIso35Track_*_*', 
        'keep *_hltPFTauMediumIso35_*_*', 
        'keep *_hltPFTauTagInfo_*_*', 
        'keep *_hltPFTauTightIso20TrackTightIso_*_*', 
        'keep *_hltPFTauTightIso20Track_*_*', 
        'keep *_hltPFTauTightIso20_*_*', 
        'keep *_hltPFlowTrackSelectionHighPurity_*_*', 
        'keep *_hltParticleFlowForTaus_*_*', 
        'keep *_hltParticleFlow_*_*', 
        'keep *_hltPixelMatch3HitElectronsActivity_*_*', 
        'keep *_hltPixelMatch3HitElectronsL1Seeded_*_*', 
        'keep *_hltPixelMatchCleanElectronsL1Seeded_*_*', 
        'keep *_hltPixelMatchElectronsActivity_*_*', 
        'keep *_hltPixelMatchElectronsL1Iso_*_*', 
        'keep *_hltPixelMatchElectronsL1NonIso_*_*', 
        'keep *_hltPixelMatchElectronsL1Seeded_*_*', 
        'keep *_hltPixelTracksElectrons_*_*', 
        'keep *_hltPixelTracksMerged_*_*', 
        'keep *_hltPixelTracks_*_*', 
        'keep *_hltPixelVertices3DbbPhi_*_*', 
        'keep *_hltPixelVertices_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidateSC4_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidateSC5_*_*', 
        'keep *_hltRecoEcalSuperClusterActivityCandidate_*_*', 
        'keep *_hltRpcRecHits_*_*', 
        'keep *_hltSelector4CentralJetsL1FastJet_*_*', 
        'keep *_hltSelector8CentralJetsL1FastJet_*_*', 
        'keep *_hltSelectorJets20L1FastJet_*_*', 
        'keep *_hltSiPixelCluster_*_*', 
        'keep *_hltSiPixelClusters_*_*', 
        'keep *_hltSiStripClusters_*_*', 
        'keep *_hltSiStripRawToClustersFacility_*_*', 
        'keep *_hltSingleMu15L3Filtered15_*_*', 
        'keep *_hltSingleMuIsoL1s14L3IsoFiltered15eta2p1_*_*', 
        'keep *_hltSingleMuIsoL3IsoFiltered15_*_*', 
        'keep *_hltTowerMakerForAll_*_*', 
        'keep *_hltTowerMakerForMuons_*_*', 
        'keep *_hltTriggerSummaryAOD_*_*', 
        'keep *_hltTriggerSummaryRAW_*_*', 
        'keep *_hltTrimmedPixelVertices_*_*', 
        'keep *_hltVerticesL3_*_*', 
        'keep *_hltVerticesPFFilter_*_*', 
        'keep *_hltVerticesPFSelector_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_rawDataRepacker_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_virginRawDataRepacker_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep L2MuonTrajectorySeeds_hltL2MuonSeeds_*_*', 
        'keep L3MuonTrajectorySeeds_hltHIL3TrajSeedOIHit_*_*', 
        'keep L3MuonTrajectorySeeds_hltHIL3TrajectorySeed_*_*', 
        'keep L3MuonTrajectorySeeds_hltL3TrajSeedOIState_*_*', 
        'keep LumiScalerss_hltScalersRawToDigi_*_*', 
        'keep SiPixelClusteredmNewDetSetVector_hltSiPixelClusters_*_*', 
        'keep TrackCandidates_hltHIL3TrackCandidateFromL2OIHit_*_*', 
        'keep TrackCandidates_hltHIL3TrackCandidateFromL2OIState_*_*', 
        'keep TrackingRecHitsOwned_hltL3Muons_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep recoCaloJets_*_*_*', 
        'keep recoCaloMETs_*_*_*', 
        'keep recoCaloMETs_hltMet_*_*', 
        'keep recoCompositeCandidates_*_*_*', 
        'keep recoElectrons_*_*_*', 
        'keep recoIsolatedPixelTrackCandidates_*_*_*', 
        'keep recoMETs_*_*_*', 
        'keep recoPFJets_*_*_*', 
        'keep recoPFTaus_*_*_*', 
        'keep recoRecoChargedCandidates_*_*_*', 
        'keep recoRecoChargedCandidates_hltHIL3MuonCandidates_*_*', 
        'keep recoRecoChargedCandidates_hltL2MuonCandidates_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsoHLTClusterShape_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonEcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalForHE_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1IsolatedPhotonHcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsoHLTClusterShape_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonEcalIsol_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalForHE_*_*', 
        'keep recoRecoEcalCandidatesToValuefloatAssociation_hltL1NonIsolatedPhotonHcalIsol_*_*', 
        'keep recoRecoEcalCandidates_*_*_*', 
        'keep recoRecoEcalCandidates_hltL1IsoRecoEcalCandidate_*_*', 
        'keep recoRecoEcalCandidates_hltL1NonIsoRecoEcalCandidate_*_*', 
        'keep recoTrackExtras_hltHIL3MuonsOIHit_*_*', 
        'keep recoTrackExtras_hltHIL3MuonsOIState_*_*', 
        'keep recoTracks_hltHIL3MuonsOIHit_*_*', 
        'keep recoTracks_hltHIL3MuonsOIState_*_*', 
        'keep recoTracks_hltHIL3Muons_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2OIHit_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2OIState_*_*', 
        'keep recoTracks_hltHIL3TkTracksFromL2_*_*', 
        'keep triggerTriggerEventWithRefs_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep triggerTriggerFilterObjectWithRefs_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*' ) ),
    splitLevel = cms.untracked.int32(0)
)

process.RECODEBUGEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring( ('drop *', 
        'drop *', 
        'drop *', 
        'keep DetIdedmEDCollection_siStripDigis_*_*', 
        'keep DetIdedmEDCollection_siPixelDigis_*_*', 
        'keep PixelFEDChanneledmNewDetSetVector_siPixelDigis_*_*', 
        'keep *_siPixelClusters_*_*', 
        'keep *_siStripClusters_*_*', 
        'keep ClusterSummary_clusterSummaryProducer_*_*', 
        'keep *_dt1DRecHits_*_*', 
        'keep *_dt4DSegments_*_*', 
        'keep *_dt1DCosmicRecHits_*_*', 
        'keep *_dt4DCosmicSegments_*_*', 
        'keep *_csc2DRecHits_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*', 
        'keep *_hbhereco_*_*', 
        'keep *_hbheprereco_*_*', 
        'keep *_hfprereco_*_*', 
        'keep *_hfreco_*_*', 
        'keep *_horeco_*_*', 
        'keep HBHERecHitsSorted_hbherecoMB_*_*', 
        'keep HORecHitsSorted_horecoMB_*_*', 
        'keep HFRecHitsSorted_hfrecoMB_*_*', 
        'keep ZDCDataFramesSorted_hcalDigis_*_*', 
        'keep ZDCDataFramesSorted_castorDigis_*_*', 
        'keep QIE10DataFrameHcalDataFrameContainer_hcalDigis_ZDC_*', 
        'keep ZDCRecHitsSorted_zdcreco_*_*', 
        'keep *_reducedHcalRecHits_*_*', 
        'keep *_castorreco_*_*', 
        'keep HcalUnpackerReport_castorDigis_*_*', 
        'keep HcalUnpackerReport_hcalDigiAlCaMB_*_*', 
        'keep HcalUnpackerReport_hcalDigis_*_*', 
        'keep *_ecalPreshowerRecHit_*_*', 
        'keep *_ecalRecHit_*_*', 
        'keep *_ecalCompactTrigPrim_*_*', 
        'keep *_ecalTPSkim_*_*', 
        'keep EBSrFlagsSorted_ecalDigis__*', 
        'keep EESrFlagsSorted_ecalDigis__*', 
        'keep *_selectDigi_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEE_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEB_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsES_*_*', 
        'keep *_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep *_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_*', 
        'keep *_particleFlowSuperClusterECAL_*_*', 
        'keep *_particleFlowSuperClusterOOTECAL_*_*', 
        'drop recoClusterShapes_*_*_*', 
        'drop recoBasicClustersToOnerecoClusterShapesAssociation_*_*_*', 
        'drop recoBasicClusters_multi5x5BasicClusters_multi5x5BarrelBasicClusters_*', 
        'drop recoSuperClusters_multi5x5SuperClusters_multi5x5BarrelSuperClusters_*', 
        'keep *_CkfElectronCandidates_*_*', 
        'keep *_GsfGlobalElectronTest_*_*', 
        'keep *_electronMergedSeeds_*_*', 
        'keep recoGsfTracks_electronGsfTracks_*_*', 
        'keep recoGsfTrackExtras_electronGsfTracks_*_*', 
        'keep recoTrackExtras_electronGsfTracks_*_*', 
        'keep TrackingRecHitsOwned_electronGsfTracks_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTrackExtras_generalTracks_*_*', 
        'keep TrackingRecHitsOwned_generalTracks_*_*', 
        'keep *_generalTracks_MVAValues_*', 
        'keep *_generalTracks_MVAVals_*', 
        'keep TrackingRecHitsOwned_extraFromSeeds_*_*', 
        'keep uints_extraFromSeeds_*_*', 
        'keep recoTracks_beamhaloTracks_*_*', 
        'keep recoTrackExtras_beamhaloTracks_*_*', 
        'keep TrackingRecHitsOwned_beamhaloTracks_*_*', 
        'keep recoTracks_conversionStepTracks_*_*', 
        'keep recoTrackExtras_conversionStepTracks_*_*', 
        'keep TrackingRecHitsOwned_conversionStepTracks_*_*', 
        'keep *_ctfPixelLess_*_*', 
        'keep *_dedxTruncated40_*_*', 
        'keep *_dedxHitInfo_*_*', 
        'keep *_dedxHarmonic2_*_*', 
        'keep *_dedxPixelHarmonic2_*_*', 
        'keep *_trackExtrapolator_*_*', 
        'keep recoTracks_cosmicDCTracks_*_*', 
        'keep recoTrackExtras_cosmicDCTracks_*_*', 
        'keep TrackingRecHitsOwned_cosmicDCTracks_*_*', 
        'keep *_ak4CaloJets_*_*', 
        'keep *_ak4PFJets_*_*', 
        'keep *_ak4PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHSSoftDrop_*_*', 
        'keep *_cmsTopTagPFJetsCHS_*_*', 
        'keep *_JetPlusTrackZSPCorJetAntiKt4_*_*', 
        'keep *_ak4TrackJets_*_*', 
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', 
        'keep *_caloTowers_*_*', 
        'keep *_towerMaker_*_*', 
        'keep *_CastorTowerReco_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertex_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertexPF_*_*', 
        'keep *_ak4JetTracksAssociatorAtCaloFace_*_*', 
        'keep *_ak4JetTracksAssociatorExplicit_*_*', 
        'keep *_ak4JetExtender_*_*', 
        'keep *_ak4JetID_*_*', 
        'keep *_ak5CastorJets_*_*', 
        'keep *_ak5CastorJetID_*_*', 
        'keep *_ak7CastorJets_*_*', 
        'keep *_ak7CastorJetID_*_*', 
        'keep *_fixedGridRhoAll_*_*', 
        'keep *_fixedGridRhoFastjetAll_*_*', 
        'keep *_fixedGridRhoFastjetAllTmp_*_*', 
        'keep *_fixedGridRhoFastjetAllCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentral_*_*', 
        'keep *_fixedGridRhoFastjetCentralCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralChargedPileUp_*_*', 
        'keep *_fixedGridRhoFastjetCentralNeutral_*_*', 
        'keep *_ak8PFJetsCHSSoftDropMass_*_*', 
        'keep recoCaloMETs_caloMet_*_*', 
        'keep recoCaloMETs_caloMetBE_*_*', 
        'keep recoCaloMETs_caloMetBEFO_*_*', 
        'keep recoCaloMETs_caloMetM_*_*', 
        'keep recoPFMETs_pfMet_*_*', 
        'keep recoPFMETs_pfChMet_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_*_*', 
        'keep recoHcalNoiseRBXs_hcalnoise_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoEcalHaloData_EcalHaloData_*_*', 
        'keep recoGlobalHaloData_GlobalHaloData_*_*', 
        'keep recoHcalHaloData_HcalHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep *_MuonSeed_*_*', 
        'keep *_ancientMuonSeed_*_*', 
        'keep *_displacedMuonSeeds_*_*', 
        'keep TrackingRecHitsOwned_globalMuons_*_*', 
        'keep TrackingRecHitsOwned_tevMuons_*_*', 
        'keep *_CosmicMuonSeed_*_*', 
        'keep recoTrackExtras_cosmicMuons_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons_*_*', 
        'keep recoTrackExtras_cosmicMuons1Leg_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons1Leg_*_*', 
        'keep recoTracks_cosmicsVetoTracks_*_*', 
        'keep *_muons_*_*', 
        'keep *_particleFlow_muons_*', 
        'drop *_muons_muons1stStep2muonsMap_*', 
        'drop recoIsoDepositedmValueMap_muons_*_*', 
        'drop doubleedmValueMap_muons_muPFIso*_*', 
        'keep recoTracks_standAloneMuons_*_*', 
        'keep recoTrackExtras_standAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_standAloneMuons_*_*', 
        'keep recoTracks_globalMuons_*_*', 
        'keep recoTrackExtras_globalMuons_*_*', 
        'keep recoTracks_tevMuons_*_*', 
        'keep recoTrackExtras_tevMuons_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_displacedTracks_*_*', 
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_*_*', 
        'keep recoTracks_displacedGlobalMuons_*_*', 
        'keep recoTrackExtras_displacedGlobalMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedGlobalMuons_*_*', 
        'keep recoTracks_cosmicMuons_*_*', 
        'keep recoMuons_muonsFromCosmics_*_*', 
        'keep recoTracks_cosmicMuons1Leg_*_*', 
        'keep recoMuons_muonsFromCosmics1Leg_*_*', 
        'keep recoTracks_refittedStandAloneMuons_*_*', 
        'keep recoTrackExtras_refittedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_refittedStandAloneMuons_*_*', 
        'keep recoTracks_displacedStandAloneMuons__*', 
        'keep recoTrackExtras_displacedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedStandAloneMuons_*_*', 
        'keep *_muIsoDepositTk_*_*', 
        'keep *_muIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muIsoDepositJets_*_*', 
        'keep *_muGlobalIsoDepositCtfTk_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muGlobalIsoDepositJets_*_*', 
        'keep *_softPFMuonsTagInfos_*_*', 
        'keep *_softPFElectronsTagInfos_*_*', 
        'keep *_softPFElectronBJetTags_*_*', 
        'keep *_softPFMuonBJetTags_*_*', 
        'keep *_pfImpactParameterTagInfos_*_*', 
        'keep *_pfTrackCountingHighEffBJetTags_*_*', 
        'keep *_pfJetProbabilityBJetTags_*_*', 
        'keep *_pfJetBProbabilityBJetTags_*_*', 
        'keep *_pfSecondaryVertexTagInfos_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderTagInfos_*_*', 
        'keep *_pfSimpleSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfCombinedSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfGhostTrackVertexTagInfos_*_*', 
        'keep *_pfGhostTrackBJetTags_*_*', 
        'keep *_pfCombinedMVAV2BJetTags_*_*', 
        'keep *_inclusiveCandidateSecondaryVertices_*_*', 
        'keep *_inclusiveCandidateSecondaryVerticesCvsL_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderCvsLTagInfos_*_*', 
        'keep *_pfCombinedCvsLJetTags_*_*', 
        'keep *_pfCombinedCvsBJetTags_*_*', 
        'keep *_pfChargeBJetTags_*_*', 
        'keep *_pfDeepCSVJetTags_*_*', 
        'keep *_pfDeepCMVAJetTags_*_*', 
        'keep recoRecoTauPiZeros_hpsPFTauProducer_pizeros_*', 
        'keep recoPFTaus_hpsPFTauProducer_*_*', 
        'keep *_hpsPFTauDiscriminationByDeadECALElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFinding_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingNewDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingOldDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseChargedIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByTightMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSum_*_*', 
        'keep *_hpsPFTauPUcorrPtSum_*_*', 
        'keep *_hpsPFTauChargedIsoPtSum_*_*', 
        'keep *_hpsPFTauTransverseImpactParameters_*_*', 
        'keep *_hpsPFTauFootprintCorrection_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeight_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6rawElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6LooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6MediumElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6TightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VTightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauChargedIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeightdR03_*_*', 
        'keep *_hpsPFTauFootprintCorrectiondR03_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalConedR03_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep  *_offlinePrimaryVertices__*', 
        'keep *_offlinePrimaryVerticesWithBS_*_*', 
        'keep *_offlinePrimaryVerticesFromCosmicTracks_*_*', 
        'keep *_nuclearInteractionMaker_*_*', 
        'keep *_generalV0Candidates_*_*', 
        'keep *_inclusiveSecondaryVertices_*_*', 
        'keep recoGsfElectronCores_gsfElectronCores_*_*', 
        'keep recoGsfElectronCores_gedGsfElectronCores_*_*', 
        'keep recoGsfElectrons_gsfElectrons_*_*', 
        'keep recoGsfElectrons_gedGsfElectrons_*_*', 
        'keep recoGsfElectronCores_uncleanedOnlyGsfElectronCores_*_*', 
        'keep recoGsfElectrons_uncleanedOnlyGsfElectrons_*_*', 
        'keep floatedmValueMap_eidRobustLoose_*_*', 
        'keep floatedmValueMap_eidRobustTight_*_*', 
        'keep floatedmValueMap_eidRobustHighEnergy_*_*', 
        'keep floatedmValueMap_eidLoose_*_*', 
        'keep floatedmValueMap_eidTight_*_*', 
        'keep *_egmGedGsfElectronPFIsolation_*_*', 
        'keep *_photonEcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonEcalPFClusterIsolationProducer_*_*', 
        'keep *_electronEcalPFClusterIsolationProducer_*_*', 
        'keep *_photonHcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonHcalPFClusterIsolationProducer_*_*', 
        'keep *_electronHcalPFClusterIsolationProducer_*_*', 
        'drop *_egmGsfElectronIDs_*_*', 
        'drop *_egmPhotonIDs_*_*', 
        'keep *_gedPhotonCore_*_*', 
        'keep *_gedPhotons_*_*', 
        'keep *_particleBasedIsolation_*_*', 
        'keep recoPhotons_mustachePhotons_*_*', 
        'keep recoPhotonCores_mustachePhotonCore_*_*', 
        'keep recoPhotons_ootPhotons_*_*', 
        'keep recoPhotonCores_ootPhotonCore_*_*', 
        'keep recoPhotons_photons_*_*', 
        'keep recoPhotonCores_photonCore_*_*', 
        'keep recoConversions_conversions_*_*', 
        'keep recoConversions_mustacheConversions_*_*', 
        'drop *_conversions_uncleanedConversions_*', 
        'drop *_gedPhotonsTmp_valMapPFEgammaCandToPhoton_*', 
        'keep recoConversions_allConversions_*_*', 
        'keep recoConversions_allConversionsOldEG_*_*', 
        'keep recoTracks_ckfOutInTracksFromConversions_*_*', 
        'keep recoTracks_ckfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfInOutTracksFromConversions_*_*', 
        'keep recoConversions_uncleanedOnlyAllConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep *_PhotonIDProdGED_*_*', 
        'keep recoRecoEcalCandidates_hfRecoEcalCandidate_*_*', 
        'keep *_hfEMClusters_*_*', 
        'keep *_gedGsfElectronCores_*_*', 
        'keep *_gedGsfElectrons_*_*', 
        'keep *_pixelTracks_*_*', 
        'keep *_pixelVertices_*_*', 
        'drop CaloTowersSorted_towerMakerPF_*_*', 
        'keep recoPFRecHits_particleFlowClusterECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHCAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterPS_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHBHE_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitPS_Cleaned_*', 
        'keep recoPFClusters_particleFlowClusterECAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHCAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHO_*_*', 
        'keep recoPFClusters_particleFlowClusterHF_*_*', 
        'keep recoPFClusters_particleFlowClusterPS_*_*', 
        'keep recoPFBlocks_particleFlowBlock_*_*', 
        'keep recoPFCandidates_particleFlowEGamma_*_*', 
        'keep recoCaloClusters_particleFlowEGamma_*_*', 
        'keep recoSuperClusters_particleFlowEGamma_*_*', 
        'keep recoConversions_particleFlowEGamma_*_*', 
        'keep recoPFCandidates_particleFlow_*_*', 
        'keep recoPFCandidates_particleFlowTmp_electrons_*', 
        'keep recoPFCandidates_particleFlowTmp_*_*', 
        'drop recoPFCandidates_particleFlowTmp__*', 
        'keep recoPFDisplacedVertexs_particleFlowDisplacedVertex_*_*', 
        'keep *_pfElectronTranslator_*_*', 
        'keep *_pfPhotonTranslator_*_*', 
        'keep *_particleFlow_electrons_*', 
        'keep *_particleFlow_photons_*', 
        'keep *_particleFlow_muons_*', 
        'keep *_trackerDrivenElectronSeeds_preid_*', 
        'keep *_particleFlowPtrs_*_*', 
        'keep *_particleFlowTmpPtrs_*_*', 
        'keep *_chargedHadronPFTrackIsolation_*_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep L1MuGMTReadoutCollection_gtDigis_*_*', 
        'keep L1GctEmCand*_gctDigis_*_*', 
        'keep L1GctJetCand*_gctDigis_*_*', 
        'keep L1GctEtHad*_gctDigis_*_*', 
        'keep L1GctEtMiss*_gctDigis_*_*', 
        'keep L1GctEtTotal*_gctDigis_*_*', 
        'keep L1GctHtMiss*_gctDigis_*_*', 
        'keep L1GctJetCounts*_gctDigis_*_*', 
        'keep L1GctHFRingEtSums*_gctDigis_*_*', 
        'keep L1GctHFBitCounts*_gctDigis_*_*', 
        'keep LumiDetails_lumiProducer_*_*', 
        'keep LumiSummary_lumiProducer_*_*', 
        'drop *_hlt*_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep L1AcceptBunchCrossings_scalersRawToDigi_*_*', 
        'keep L1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep Level1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep BeamSpotOnlines_scalersRawToDigi_*_*', 
        'keep DcsStatuss_scalersRawToDigi_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep *_tcdsDigis_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'keep *_pfIsolatedElectronsEI_*_*', 
        'keep *_pfIsolatedMuonsEI_*_*', 
        'keep recoPFJets_pfJetsEI_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTagsEI_*_*', 
        'keep recoPFTaus_pfTausEI_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByDecayModeFinding_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByIsolation_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep *_ak4GenJets_*_*', 
        'keep *_ak8GenJets_*_*', 
        'keep *_ak4GenJetsNoNu_*_*', 
        'keep *_ak8GenJetsNoNu_*_*', 
        'keep *_genParticle_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep SimTracks_g4SimHits_*_*', 
        'keep SimVertexs_g4SimHits_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*', 
        'keep *_muonSimClassifier_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'drop *_trackingtruthprod_*_*', 
        'drop *_electrontruth_*_*', 
        'keep *_mix_MergedTrackTruth_*', 
        'keep CrossingFramePlaybackInfoNew_*_*_*', 
        'keep PixelDigiSimLinkedmDetSetVector_simSiPixelDigis_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simSiStripDigis_*_*', 
        'drop *_mix_simSiPixelDigis*_*', 
        'drop *_mix_simSiStripDigis*_*', 
        'keep *_allTrackMCMatch_*_*' ) ),
    splitLevel = cms.untracked.int32(0)
)

process.RECOEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring( ('drop *', 
        'keep DetIdedmEDCollection_siStripDigis_*_*', 
        'keep DetIdedmEDCollection_siPixelDigis_*_*', 
        'keep PixelFEDChanneledmNewDetSetVector_siPixelDigis_*_*', 
        'keep *_siPixelClusters_*_*', 
        'keep *_siStripClusters_*_*', 
        'keep ClusterSummary_clusterSummaryProducer_*_*', 
        'keep *_dt1DRecHits_*_*', 
        'keep *_dt4DSegments_*_*', 
        'keep *_dt1DCosmicRecHits_*_*', 
        'keep *_dt4DCosmicSegments_*_*', 
        'keep *_csc2DRecHits_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*', 
        'keep *_hbhereco_*_*', 
        'keep *_hbheprereco_*_*', 
        'keep *_hfprereco_*_*', 
        'keep *_hfreco_*_*', 
        'keep *_horeco_*_*', 
        'keep HBHERecHitsSorted_hbherecoMB_*_*', 
        'keep HORecHitsSorted_horecoMB_*_*', 
        'keep HFRecHitsSorted_hfrecoMB_*_*', 
        'keep ZDCDataFramesSorted_hcalDigis_*_*', 
        'keep ZDCDataFramesSorted_castorDigis_*_*', 
        'keep QIE10DataFrameHcalDataFrameContainer_hcalDigis_ZDC_*', 
        'keep ZDCRecHitsSorted_zdcreco_*_*', 
        'keep *_reducedHcalRecHits_*_*', 
        'keep *_castorreco_*_*', 
        'keep HcalUnpackerReport_castorDigis_*_*', 
        'keep HcalUnpackerReport_hcalDigiAlCaMB_*_*', 
        'keep HcalUnpackerReport_hcalDigis_*_*', 
        'keep *_ecalPreshowerRecHit_*_*', 
        'keep *_ecalRecHit_*_*', 
        'keep *_ecalCompactTrigPrim_*_*', 
        'keep *_ecalTPSkim_*_*', 
        'keep EBSrFlagsSorted_ecalDigis__*', 
        'keep EESrFlagsSorted_ecalDigis__*', 
        'keep *_selectDigi_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEE_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEB_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsES_*_*', 
        'keep *_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep *_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_*', 
        'keep *_particleFlowSuperClusterECAL_*_*', 
        'keep *_particleFlowSuperClusterOOTECAL_*_*', 
        'drop recoClusterShapes_*_*_*', 
        'drop recoBasicClustersToOnerecoClusterShapesAssociation_*_*_*', 
        'drop recoBasicClusters_multi5x5BasicClusters_multi5x5BarrelBasicClusters_*', 
        'drop recoSuperClusters_multi5x5SuperClusters_multi5x5BarrelSuperClusters_*', 
        'keep *_CkfElectronCandidates_*_*', 
        'keep *_GsfGlobalElectronTest_*_*', 
        'keep *_electronMergedSeeds_*_*', 
        'keep recoGsfTracks_electronGsfTracks_*_*', 
        'keep recoGsfTrackExtras_electronGsfTracks_*_*', 
        'keep recoTrackExtras_electronGsfTracks_*_*', 
        'keep TrackingRecHitsOwned_electronGsfTracks_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTrackExtras_generalTracks_*_*', 
        'keep TrackingRecHitsOwned_generalTracks_*_*', 
        'keep *_generalTracks_MVAValues_*', 
        'keep *_generalTracks_MVAVals_*', 
        'keep TrackingRecHitsOwned_extraFromSeeds_*_*', 
        'keep uints_extraFromSeeds_*_*', 
        'keep recoTracks_beamhaloTracks_*_*', 
        'keep recoTrackExtras_beamhaloTracks_*_*', 
        'keep TrackingRecHitsOwned_beamhaloTracks_*_*', 
        'keep recoTracks_conversionStepTracks_*_*', 
        'keep recoTrackExtras_conversionStepTracks_*_*', 
        'keep TrackingRecHitsOwned_conversionStepTracks_*_*', 
        'keep *_ctfPixelLess_*_*', 
        'keep *_dedxTruncated40_*_*', 
        'keep *_dedxHitInfo_*_*', 
        'keep *_dedxHarmonic2_*_*', 
        'keep *_dedxPixelHarmonic2_*_*', 
        'keep *_trackExtrapolator_*_*', 
        'keep recoTracks_cosmicDCTracks_*_*', 
        'keep recoTrackExtras_cosmicDCTracks_*_*', 
        'keep TrackingRecHitsOwned_cosmicDCTracks_*_*', 
        'keep *_ak4CaloJets_*_*', 
        'keep *_ak4PFJets_*_*', 
        'keep *_ak4PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHSSoftDrop_*_*', 
        'keep *_cmsTopTagPFJetsCHS_*_*', 
        'keep *_JetPlusTrackZSPCorJetAntiKt4_*_*', 
        'keep *_ak4TrackJets_*_*', 
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', 
        'keep *_caloTowers_*_*', 
        'keep *_towerMaker_*_*', 
        'keep *_CastorTowerReco_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertex_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertexPF_*_*', 
        'keep *_ak4JetTracksAssociatorAtCaloFace_*_*', 
        'keep *_ak4JetTracksAssociatorExplicit_*_*', 
        'keep *_ak4JetExtender_*_*', 
        'keep *_ak4JetID_*_*', 
        'keep *_ak5CastorJets_*_*', 
        'keep *_ak5CastorJetID_*_*', 
        'keep *_ak7CastorJets_*_*', 
        'keep *_ak7CastorJetID_*_*', 
        'keep *_fixedGridRhoAll_*_*', 
        'keep *_fixedGridRhoFastjetAll_*_*', 
        'keep *_fixedGridRhoFastjetAllTmp_*_*', 
        'keep *_fixedGridRhoFastjetAllCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentral_*_*', 
        'keep *_fixedGridRhoFastjetCentralCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralChargedPileUp_*_*', 
        'keep *_fixedGridRhoFastjetCentralNeutral_*_*', 
        'keep *_ak8PFJetsCHSSoftDropMass_*_*', 
        'keep recoCaloMETs_caloMet_*_*', 
        'keep recoCaloMETs_caloMetBE_*_*', 
        'keep recoCaloMETs_caloMetBEFO_*_*', 
        'keep recoCaloMETs_caloMetM_*_*', 
        'keep recoPFMETs_pfMet_*_*', 
        'keep recoPFMETs_pfChMet_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_*_*', 
        'keep recoHcalNoiseRBXs_hcalnoise_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoEcalHaloData_EcalHaloData_*_*', 
        'keep recoGlobalHaloData_GlobalHaloData_*_*', 
        'keep recoHcalHaloData_HcalHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep *_MuonSeed_*_*', 
        'keep *_ancientMuonSeed_*_*', 
        'keep *_displacedMuonSeeds_*_*', 
        'keep TrackingRecHitsOwned_globalMuons_*_*', 
        'keep TrackingRecHitsOwned_tevMuons_*_*', 
        'keep *_CosmicMuonSeed_*_*', 
        'keep recoTrackExtras_cosmicMuons_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons_*_*', 
        'keep recoTrackExtras_cosmicMuons1Leg_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons1Leg_*_*', 
        'keep recoTracks_cosmicsVetoTracks_*_*', 
        'keep *_muons_*_*', 
        'keep *_particleFlow_muons_*', 
        'drop *_muons_muons1stStep2muonsMap_*', 
        'drop recoIsoDepositedmValueMap_muons_*_*', 
        'drop doubleedmValueMap_muons_muPFIso*_*', 
        'keep recoTracks_standAloneMuons_*_*', 
        'keep recoTrackExtras_standAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_standAloneMuons_*_*', 
        'keep recoTracks_globalMuons_*_*', 
        'keep recoTrackExtras_globalMuons_*_*', 
        'keep recoTracks_tevMuons_*_*', 
        'keep recoTrackExtras_tevMuons_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_displacedTracks_*_*', 
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_*_*', 
        'keep recoTracks_displacedGlobalMuons_*_*', 
        'keep recoTrackExtras_displacedGlobalMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedGlobalMuons_*_*', 
        'keep recoTracks_cosmicMuons_*_*', 
        'keep recoMuons_muonsFromCosmics_*_*', 
        'keep recoTracks_cosmicMuons1Leg_*_*', 
        'keep recoMuons_muonsFromCosmics1Leg_*_*', 
        'keep recoTracks_refittedStandAloneMuons_*_*', 
        'keep recoTrackExtras_refittedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_refittedStandAloneMuons_*_*', 
        'keep recoTracks_displacedStandAloneMuons__*', 
        'keep recoTrackExtras_displacedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedStandAloneMuons_*_*', 
        'keep *_muIsoDepositTk_*_*', 
        'keep *_muIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muIsoDepositJets_*_*', 
        'keep *_muGlobalIsoDepositCtfTk_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muGlobalIsoDepositJets_*_*', 
        'keep *_softPFMuonsTagInfos_*_*', 
        'keep *_softPFElectronsTagInfos_*_*', 
        'keep *_softPFElectronBJetTags_*_*', 
        'keep *_softPFMuonBJetTags_*_*', 
        'keep *_pfImpactParameterTagInfos_*_*', 
        'keep *_pfTrackCountingHighEffBJetTags_*_*', 
        'keep *_pfJetProbabilityBJetTags_*_*', 
        'keep *_pfJetBProbabilityBJetTags_*_*', 
        'keep *_pfSecondaryVertexTagInfos_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderTagInfos_*_*', 
        'keep *_pfSimpleSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfCombinedSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfGhostTrackVertexTagInfos_*_*', 
        'keep *_pfGhostTrackBJetTags_*_*', 
        'keep *_pfCombinedMVAV2BJetTags_*_*', 
        'keep *_inclusiveCandidateSecondaryVertices_*_*', 
        'keep *_inclusiveCandidateSecondaryVerticesCvsL_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderCvsLTagInfos_*_*', 
        'keep *_pfCombinedCvsLJetTags_*_*', 
        'keep *_pfCombinedCvsBJetTags_*_*', 
        'keep *_pfChargeBJetTags_*_*', 
        'keep *_pfDeepCSVJetTags_*_*', 
        'keep *_pfDeepCMVAJetTags_*_*', 
        'keep recoRecoTauPiZeros_hpsPFTauProducer_pizeros_*', 
        'keep recoPFTaus_hpsPFTauProducer_*_*', 
        'keep *_hpsPFTauDiscriminationByDeadECALElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFinding_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingNewDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingOldDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseChargedIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByTightMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSum_*_*', 
        'keep *_hpsPFTauPUcorrPtSum_*_*', 
        'keep *_hpsPFTauChargedIsoPtSum_*_*', 
        'keep *_hpsPFTauTransverseImpactParameters_*_*', 
        'keep *_hpsPFTauFootprintCorrection_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeight_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6rawElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6LooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6MediumElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6TightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VTightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauChargedIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeightdR03_*_*', 
        'keep *_hpsPFTauFootprintCorrectiondR03_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalConedR03_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep  *_offlinePrimaryVertices__*', 
        'keep *_offlinePrimaryVerticesWithBS_*_*', 
        'keep *_offlinePrimaryVerticesFromCosmicTracks_*_*', 
        'keep *_nuclearInteractionMaker_*_*', 
        'keep *_generalV0Candidates_*_*', 
        'keep *_inclusiveSecondaryVertices_*_*', 
        'keep recoGsfElectronCores_gsfElectronCores_*_*', 
        'keep recoGsfElectronCores_gedGsfElectronCores_*_*', 
        'keep recoGsfElectrons_gsfElectrons_*_*', 
        'keep recoGsfElectrons_gedGsfElectrons_*_*', 
        'keep recoGsfElectronCores_uncleanedOnlyGsfElectronCores_*_*', 
        'keep recoGsfElectrons_uncleanedOnlyGsfElectrons_*_*', 
        'keep floatedmValueMap_eidRobustLoose_*_*', 
        'keep floatedmValueMap_eidRobustTight_*_*', 
        'keep floatedmValueMap_eidRobustHighEnergy_*_*', 
        'keep floatedmValueMap_eidLoose_*_*', 
        'keep floatedmValueMap_eidTight_*_*', 
        'keep *_egmGedGsfElectronPFIsolation_*_*', 
        'keep *_photonEcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonEcalPFClusterIsolationProducer_*_*', 
        'keep *_electronEcalPFClusterIsolationProducer_*_*', 
        'keep *_photonHcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonHcalPFClusterIsolationProducer_*_*', 
        'keep *_electronHcalPFClusterIsolationProducer_*_*', 
        'drop *_egmGsfElectronIDs_*_*', 
        'drop *_egmPhotonIDs_*_*', 
        'keep *_gedPhotonCore_*_*', 
        'keep *_gedPhotons_*_*', 
        'keep *_particleBasedIsolation_*_*', 
        'keep recoPhotons_mustachePhotons_*_*', 
        'keep recoPhotonCores_mustachePhotonCore_*_*', 
        'keep recoPhotons_ootPhotons_*_*', 
        'keep recoPhotonCores_ootPhotonCore_*_*', 
        'keep recoPhotons_photons_*_*', 
        'keep recoPhotonCores_photonCore_*_*', 
        'keep recoConversions_conversions_*_*', 
        'keep recoConversions_mustacheConversions_*_*', 
        'drop *_conversions_uncleanedConversions_*', 
        'drop *_gedPhotonsTmp_valMapPFEgammaCandToPhoton_*', 
        'keep recoConversions_allConversions_*_*', 
        'keep recoConversions_allConversionsOldEG_*_*', 
        'keep recoTracks_ckfOutInTracksFromConversions_*_*', 
        'keep recoTracks_ckfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfInOutTracksFromConversions_*_*', 
        'keep recoConversions_uncleanedOnlyAllConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep *_PhotonIDProdGED_*_*', 
        'keep recoRecoEcalCandidates_hfRecoEcalCandidate_*_*', 
        'keep *_hfEMClusters_*_*', 
        'keep *_gedGsfElectronCores_*_*', 
        'keep *_gedGsfElectrons_*_*', 
        'keep *_pixelTracks_*_*', 
        'keep *_pixelVertices_*_*', 
        'drop CaloTowersSorted_towerMakerPF_*_*', 
        'keep recoPFRecHits_particleFlowClusterECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHCAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterPS_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHBHE_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitPS_Cleaned_*', 
        'keep recoPFClusters_particleFlowClusterECAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHCAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHO_*_*', 
        'keep recoPFClusters_particleFlowClusterHF_*_*', 
        'keep recoPFClusters_particleFlowClusterPS_*_*', 
        'keep recoPFBlocks_particleFlowBlock_*_*', 
        'keep recoPFCandidates_particleFlowEGamma_*_*', 
        'keep recoCaloClusters_particleFlowEGamma_*_*', 
        'keep recoSuperClusters_particleFlowEGamma_*_*', 
        'keep recoConversions_particleFlowEGamma_*_*', 
        'keep recoPFCandidates_particleFlow_*_*', 
        'keep recoPFCandidates_particleFlowTmp_electrons_*', 
        'keep recoPFCandidates_particleFlowTmp_*_*', 
        'drop recoPFCandidates_particleFlowTmp__*', 
        'keep recoPFDisplacedVertexs_particleFlowDisplacedVertex_*_*', 
        'keep *_pfElectronTranslator_*_*', 
        'keep *_pfPhotonTranslator_*_*', 
        'keep *_particleFlow_electrons_*', 
        'keep *_particleFlow_photons_*', 
        'keep *_particleFlow_muons_*', 
        'keep *_trackerDrivenElectronSeeds_preid_*', 
        'keep *_particleFlowPtrs_*_*', 
        'keep *_particleFlowTmpPtrs_*_*', 
        'keep *_chargedHadronPFTrackIsolation_*_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep L1MuGMTReadoutCollection_gtDigis_*_*', 
        'keep L1GctEmCand*_gctDigis_*_*', 
        'keep L1GctJetCand*_gctDigis_*_*', 
        'keep L1GctEtHad*_gctDigis_*_*', 
        'keep L1GctEtMiss*_gctDigis_*_*', 
        'keep L1GctEtTotal*_gctDigis_*_*', 
        'keep L1GctHtMiss*_gctDigis_*_*', 
        'keep L1GctJetCounts*_gctDigis_*_*', 
        'keep L1GctHFRingEtSums*_gctDigis_*_*', 
        'keep L1GctHFBitCounts*_gctDigis_*_*', 
        'keep LumiDetails_lumiProducer_*_*', 
        'keep LumiSummary_lumiProducer_*_*', 
        'drop *_hlt*_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep L1AcceptBunchCrossings_scalersRawToDigi_*_*', 
        'keep L1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep Level1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep BeamSpotOnlines_scalersRawToDigi_*_*', 
        'keep DcsStatuss_scalersRawToDigi_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep *_tcdsDigis_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'keep *_pfIsolatedElectronsEI_*_*', 
        'keep *_pfIsolatedMuonsEI_*_*', 
        'keep recoPFJets_pfJetsEI_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTagsEI_*_*', 
        'keep recoPFTaus_pfTausEI_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByDecayModeFinding_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByIsolation_*_*', 
        'keep recoPFMETs_pfMetEI_*_*' ) ),
    splitLevel = cms.untracked.int32(0)
)

process.RECOSIMEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring( ('drop *', 
        'drop *', 
        'keep DetIdedmEDCollection_siStripDigis_*_*', 
        'keep DetIdedmEDCollection_siPixelDigis_*_*', 
        'keep PixelFEDChanneledmNewDetSetVector_siPixelDigis_*_*', 
        'keep *_siPixelClusters_*_*', 
        'keep *_siStripClusters_*_*', 
        'keep ClusterSummary_clusterSummaryProducer_*_*', 
        'keep *_dt1DRecHits_*_*', 
        'keep *_dt4DSegments_*_*', 
        'keep *_dt1DCosmicRecHits_*_*', 
        'keep *_dt4DCosmicSegments_*_*', 
        'keep *_csc2DRecHits_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*', 
        'keep *_hbhereco_*_*', 
        'keep *_hbheprereco_*_*', 
        'keep *_hfprereco_*_*', 
        'keep *_hfreco_*_*', 
        'keep *_horeco_*_*', 
        'keep HBHERecHitsSorted_hbherecoMB_*_*', 
        'keep HORecHitsSorted_horecoMB_*_*', 
        'keep HFRecHitsSorted_hfrecoMB_*_*', 
        'keep ZDCDataFramesSorted_hcalDigis_*_*', 
        'keep ZDCDataFramesSorted_castorDigis_*_*', 
        'keep QIE10DataFrameHcalDataFrameContainer_hcalDigis_ZDC_*', 
        'keep ZDCRecHitsSorted_zdcreco_*_*', 
        'keep *_reducedHcalRecHits_*_*', 
        'keep *_castorreco_*_*', 
        'keep HcalUnpackerReport_castorDigis_*_*', 
        'keep HcalUnpackerReport_hcalDigiAlCaMB_*_*', 
        'keep HcalUnpackerReport_hcalDigis_*_*', 
        'keep *_ecalPreshowerRecHit_*_*', 
        'keep *_ecalRecHit_*_*', 
        'keep *_ecalCompactTrigPrim_*_*', 
        'keep *_ecalTPSkim_*_*', 
        'keep EBSrFlagsSorted_ecalDigis__*', 
        'keep EESrFlagsSorted_ecalDigis__*', 
        'keep *_selectDigi_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEE_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEB_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsES_*_*', 
        'keep *_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep *_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_*', 
        'keep *_particleFlowSuperClusterECAL_*_*', 
        'keep *_particleFlowSuperClusterOOTECAL_*_*', 
        'drop recoClusterShapes_*_*_*', 
        'drop recoBasicClustersToOnerecoClusterShapesAssociation_*_*_*', 
        'drop recoBasicClusters_multi5x5BasicClusters_multi5x5BarrelBasicClusters_*', 
        'drop recoSuperClusters_multi5x5SuperClusters_multi5x5BarrelSuperClusters_*', 
        'keep *_CkfElectronCandidates_*_*', 
        'keep *_GsfGlobalElectronTest_*_*', 
        'keep *_electronMergedSeeds_*_*', 
        'keep recoGsfTracks_electronGsfTracks_*_*', 
        'keep recoGsfTrackExtras_electronGsfTracks_*_*', 
        'keep recoTrackExtras_electronGsfTracks_*_*', 
        'keep TrackingRecHitsOwned_electronGsfTracks_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTrackExtras_generalTracks_*_*', 
        'keep TrackingRecHitsOwned_generalTracks_*_*', 
        'keep *_generalTracks_MVAValues_*', 
        'keep *_generalTracks_MVAVals_*', 
        'keep TrackingRecHitsOwned_extraFromSeeds_*_*', 
        'keep uints_extraFromSeeds_*_*', 
        'keep recoTracks_beamhaloTracks_*_*', 
        'keep recoTrackExtras_beamhaloTracks_*_*', 
        'keep TrackingRecHitsOwned_beamhaloTracks_*_*', 
        'keep recoTracks_conversionStepTracks_*_*', 
        'keep recoTrackExtras_conversionStepTracks_*_*', 
        'keep TrackingRecHitsOwned_conversionStepTracks_*_*', 
        'keep *_ctfPixelLess_*_*', 
        'keep *_dedxTruncated40_*_*', 
        'keep *_dedxHitInfo_*_*', 
        'keep *_dedxHarmonic2_*_*', 
        'keep *_dedxPixelHarmonic2_*_*', 
        'keep *_trackExtrapolator_*_*', 
        'keep recoTracks_cosmicDCTracks_*_*', 
        'keep recoTrackExtras_cosmicDCTracks_*_*', 
        'keep TrackingRecHitsOwned_cosmicDCTracks_*_*', 
        'keep *_ak4CaloJets_*_*', 
        'keep *_ak4PFJets_*_*', 
        'keep *_ak4PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHSSoftDrop_*_*', 
        'keep *_cmsTopTagPFJetsCHS_*_*', 
        'keep *_JetPlusTrackZSPCorJetAntiKt4_*_*', 
        'keep *_ak4TrackJets_*_*', 
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', 
        'keep *_caloTowers_*_*', 
        'keep *_towerMaker_*_*', 
        'keep *_CastorTowerReco_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertex_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertexPF_*_*', 
        'keep *_ak4JetTracksAssociatorAtCaloFace_*_*', 
        'keep *_ak4JetTracksAssociatorExplicit_*_*', 
        'keep *_ak4JetExtender_*_*', 
        'keep *_ak4JetID_*_*', 
        'keep *_ak5CastorJets_*_*', 
        'keep *_ak5CastorJetID_*_*', 
        'keep *_ak7CastorJets_*_*', 
        'keep *_ak7CastorJetID_*_*', 
        'keep *_fixedGridRhoAll_*_*', 
        'keep *_fixedGridRhoFastjetAll_*_*', 
        'keep *_fixedGridRhoFastjetAllTmp_*_*', 
        'keep *_fixedGridRhoFastjetAllCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentral_*_*', 
        'keep *_fixedGridRhoFastjetCentralCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralChargedPileUp_*_*', 
        'keep *_fixedGridRhoFastjetCentralNeutral_*_*', 
        'keep *_ak8PFJetsCHSSoftDropMass_*_*', 
        'keep recoCaloMETs_caloMet_*_*', 
        'keep recoCaloMETs_caloMetBE_*_*', 
        'keep recoCaloMETs_caloMetBEFO_*_*', 
        'keep recoCaloMETs_caloMetM_*_*', 
        'keep recoPFMETs_pfMet_*_*', 
        'keep recoPFMETs_pfChMet_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_*_*', 
        'keep recoHcalNoiseRBXs_hcalnoise_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoEcalHaloData_EcalHaloData_*_*', 
        'keep recoGlobalHaloData_GlobalHaloData_*_*', 
        'keep recoHcalHaloData_HcalHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*', 
        'keep *_MuonSeed_*_*', 
        'keep *_ancientMuonSeed_*_*', 
        'keep *_displacedMuonSeeds_*_*', 
        'keep TrackingRecHitsOwned_globalMuons_*_*', 
        'keep TrackingRecHitsOwned_tevMuons_*_*', 
        'keep *_CosmicMuonSeed_*_*', 
        'keep recoTrackExtras_cosmicMuons_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons_*_*', 
        'keep recoTrackExtras_cosmicMuons1Leg_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons1Leg_*_*', 
        'keep recoTracks_cosmicsVetoTracks_*_*', 
        'keep *_muons_*_*', 
        'keep *_particleFlow_muons_*', 
        'drop *_muons_muons1stStep2muonsMap_*', 
        'drop recoIsoDepositedmValueMap_muons_*_*', 
        'drop doubleedmValueMap_muons_muPFIso*_*', 
        'keep recoTracks_standAloneMuons_*_*', 
        'keep recoTrackExtras_standAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_standAloneMuons_*_*', 
        'keep recoTracks_globalMuons_*_*', 
        'keep recoTrackExtras_globalMuons_*_*', 
        'keep recoTracks_tevMuons_*_*', 
        'keep recoTrackExtras_tevMuons_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_displacedTracks_*_*', 
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_*_*', 
        'keep recoTracks_displacedGlobalMuons_*_*', 
        'keep recoTrackExtras_displacedGlobalMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedGlobalMuons_*_*', 
        'keep recoTracks_cosmicMuons_*_*', 
        'keep recoMuons_muonsFromCosmics_*_*', 
        'keep recoTracks_cosmicMuons1Leg_*_*', 
        'keep recoMuons_muonsFromCosmics1Leg_*_*', 
        'keep recoTracks_refittedStandAloneMuons_*_*', 
        'keep recoTrackExtras_refittedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_refittedStandAloneMuons_*_*', 
        'keep recoTracks_displacedStandAloneMuons__*', 
        'keep recoTrackExtras_displacedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedStandAloneMuons_*_*', 
        'keep *_muIsoDepositTk_*_*', 
        'keep *_muIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muIsoDepositJets_*_*', 
        'keep *_muGlobalIsoDepositCtfTk_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muGlobalIsoDepositJets_*_*', 
        'keep *_softPFMuonsTagInfos_*_*', 
        'keep *_softPFElectronsTagInfos_*_*', 
        'keep *_softPFElectronBJetTags_*_*', 
        'keep *_softPFMuonBJetTags_*_*', 
        'keep *_pfImpactParameterTagInfos_*_*', 
        'keep *_pfTrackCountingHighEffBJetTags_*_*', 
        'keep *_pfJetProbabilityBJetTags_*_*', 
        'keep *_pfJetBProbabilityBJetTags_*_*', 
        'keep *_pfSecondaryVertexTagInfos_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderTagInfos_*_*', 
        'keep *_pfSimpleSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfCombinedSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfGhostTrackVertexTagInfos_*_*', 
        'keep *_pfGhostTrackBJetTags_*_*', 
        'keep *_pfCombinedMVAV2BJetTags_*_*', 
        'keep *_inclusiveCandidateSecondaryVertices_*_*', 
        'keep *_inclusiveCandidateSecondaryVerticesCvsL_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderCvsLTagInfos_*_*', 
        'keep *_pfCombinedCvsLJetTags_*_*', 
        'keep *_pfCombinedCvsBJetTags_*_*', 
        'keep *_pfChargeBJetTags_*_*', 
        'keep *_pfDeepCSVJetTags_*_*', 
        'keep *_pfDeepCMVAJetTags_*_*', 
        'keep recoRecoTauPiZeros_hpsPFTauProducer_pizeros_*', 
        'keep recoPFTaus_hpsPFTauProducer_*_*', 
        'keep *_hpsPFTauDiscriminationByDeadECALElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFinding_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingNewDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingOldDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseChargedIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByTightMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSum_*_*', 
        'keep *_hpsPFTauPUcorrPtSum_*_*', 
        'keep *_hpsPFTauChargedIsoPtSum_*_*', 
        'keep *_hpsPFTauTransverseImpactParameters_*_*', 
        'keep *_hpsPFTauFootprintCorrection_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeight_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6rawElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6LooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6MediumElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6TightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VTightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauChargedIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeightdR03_*_*', 
        'keep *_hpsPFTauFootprintCorrectiondR03_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalConedR03_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep  *_offlinePrimaryVertices__*', 
        'keep *_offlinePrimaryVerticesWithBS_*_*', 
        'keep *_offlinePrimaryVerticesFromCosmicTracks_*_*', 
        'keep *_nuclearInteractionMaker_*_*', 
        'keep *_generalV0Candidates_*_*', 
        'keep *_inclusiveSecondaryVertices_*_*', 
        'keep recoGsfElectronCores_gsfElectronCores_*_*', 
        'keep recoGsfElectronCores_gedGsfElectronCores_*_*', 
        'keep recoGsfElectrons_gsfElectrons_*_*', 
        'keep recoGsfElectrons_gedGsfElectrons_*_*', 
        'keep recoGsfElectronCores_uncleanedOnlyGsfElectronCores_*_*', 
        'keep recoGsfElectrons_uncleanedOnlyGsfElectrons_*_*', 
        'keep floatedmValueMap_eidRobustLoose_*_*', 
        'keep floatedmValueMap_eidRobustTight_*_*', 
        'keep floatedmValueMap_eidRobustHighEnergy_*_*', 
        'keep floatedmValueMap_eidLoose_*_*', 
        'keep floatedmValueMap_eidTight_*_*', 
        'keep *_egmGedGsfElectronPFIsolation_*_*', 
        'keep *_photonEcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonEcalPFClusterIsolationProducer_*_*', 
        'keep *_electronEcalPFClusterIsolationProducer_*_*', 
        'keep *_photonHcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonHcalPFClusterIsolationProducer_*_*', 
        'keep *_electronHcalPFClusterIsolationProducer_*_*', 
        'drop *_egmGsfElectronIDs_*_*', 
        'drop *_egmPhotonIDs_*_*', 
        'keep *_gedPhotonCore_*_*', 
        'keep *_gedPhotons_*_*', 
        'keep *_particleBasedIsolation_*_*', 
        'keep recoPhotons_mustachePhotons_*_*', 
        'keep recoPhotonCores_mustachePhotonCore_*_*', 
        'keep recoPhotons_ootPhotons_*_*', 
        'keep recoPhotonCores_ootPhotonCore_*_*', 
        'keep recoPhotons_photons_*_*', 
        'keep recoPhotonCores_photonCore_*_*', 
        'keep recoConversions_conversions_*_*', 
        'keep recoConversions_mustacheConversions_*_*', 
        'drop *_conversions_uncleanedConversions_*', 
        'drop *_gedPhotonsTmp_valMapPFEgammaCandToPhoton_*', 
        'keep recoConversions_allConversions_*_*', 
        'keep recoConversions_allConversionsOldEG_*_*', 
        'keep recoTracks_ckfOutInTracksFromConversions_*_*', 
        'keep recoTracks_ckfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfInOutTracksFromConversions_*_*', 
        'keep recoConversions_uncleanedOnlyAllConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep *_PhotonIDProdGED_*_*', 
        'keep recoRecoEcalCandidates_hfRecoEcalCandidate_*_*', 
        'keep *_hfEMClusters_*_*', 
        'keep *_gedGsfElectronCores_*_*', 
        'keep *_gedGsfElectrons_*_*', 
        'keep *_pixelTracks_*_*', 
        'keep *_pixelVertices_*_*', 
        'drop CaloTowersSorted_towerMakerPF_*_*', 
        'keep recoPFRecHits_particleFlowClusterECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHCAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterPS_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHBHE_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitPS_Cleaned_*', 
        'keep recoPFClusters_particleFlowClusterECAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHCAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHO_*_*', 
        'keep recoPFClusters_particleFlowClusterHF_*_*', 
        'keep recoPFClusters_particleFlowClusterPS_*_*', 
        'keep recoPFBlocks_particleFlowBlock_*_*', 
        'keep recoPFCandidates_particleFlowEGamma_*_*', 
        'keep recoCaloClusters_particleFlowEGamma_*_*', 
        'keep recoSuperClusters_particleFlowEGamma_*_*', 
        'keep recoConversions_particleFlowEGamma_*_*', 
        'keep recoPFCandidates_particleFlow_*_*', 
        'keep recoPFCandidates_particleFlowTmp_electrons_*', 
        'keep recoPFCandidates_particleFlowTmp_*_*', 
        'drop recoPFCandidates_particleFlowTmp__*', 
        'keep recoPFDisplacedVertexs_particleFlowDisplacedVertex_*_*', 
        'keep *_pfElectronTranslator_*_*', 
        'keep *_pfPhotonTranslator_*_*', 
        'keep *_particleFlow_electrons_*', 
        'keep *_particleFlow_photons_*', 
        'keep *_particleFlow_muons_*', 
        'keep *_trackerDrivenElectronSeeds_preid_*', 
        'keep *_particleFlowPtrs_*_*', 
        'keep *_particleFlowTmpPtrs_*_*', 
        'keep *_chargedHadronPFTrackIsolation_*_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1GtRecord_*_*', 
        'keep *_l1GtTriggerMenuLite_*_*', 
        'keep *_conditionsInEdm_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep L1MuGMTReadoutCollection_gtDigis_*_*', 
        'keep L1GctEmCand*_gctDigis_*_*', 
        'keep L1GctJetCand*_gctDigis_*_*', 
        'keep L1GctEtHad*_gctDigis_*_*', 
        'keep L1GctEtMiss*_gctDigis_*_*', 
        'keep L1GctEtTotal*_gctDigis_*_*', 
        'keep L1GctHtMiss*_gctDigis_*_*', 
        'keep L1GctJetCounts*_gctDigis_*_*', 
        'keep L1GctHFRingEtSums*_gctDigis_*_*', 
        'keep L1GctHFBitCounts*_gctDigis_*_*', 
        'keep LumiDetails_lumiProducer_*_*', 
        'keep LumiSummary_lumiProducer_*_*', 
        'drop *_hlt*_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep L1AcceptBunchCrossings_scalersRawToDigi_*_*', 
        'keep L1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep Level1TriggerScalerss_scalersRawToDigi_*_*', 
        'keep LumiScalerss_scalersRawToDigi_*_*', 
        'keep BeamSpotOnlines_scalersRawToDigi_*_*', 
        'keep DcsStatuss_scalersRawToDigi_*_*', 
        'keep DcsStatuss_hltScalersRawToDigi_*_*', 
        'keep *_tcdsDigis_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'keep *_pfIsolatedElectronsEI_*_*', 
        'keep *_pfIsolatedMuonsEI_*_*', 
        'keep recoPFJets_pfJetsEI_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTagsEI_*_*', 
        'keep recoPFTaus_pfTausEI_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByDecayModeFinding_*_*', 
        'keep recoPFTauDiscriminator_pfTausDiscriminationByIsolation_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep *_ak4GenJets_*_*', 
        'keep *_ak8GenJets_*_*', 
        'keep *_ak4GenJetsNoNu_*_*', 
        'keep *_ak8GenJetsNoNu_*_*', 
        'keep *_genParticle_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep SimTracks_g4SimHits_*_*', 
        'keep SimVertexs_g4SimHits_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*', 
        'keep *_muonSimClassifier_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*' ) ),
    splitLevel = cms.untracked.int32(0)
)

process.REDIGIEventContent = cms.PSet(
    inputCommands = cms.untracked.vstring('drop *', 
        'keep *_g4SimHits_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep *_randomEngineStateProducer_*_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'drop *_randomEngineStateProducer_*_*')
)

process.REGENEventContent = cms.PSet(
    inputCommands = cms.untracked.vstring('keep *', 
        'drop *_genParticles_*_*', 
        'drop *_genParticlesForJets_*_*', 
        'drop *_kt4GenJets_*_*', 
        'drop *_kt6GenJets_*_*', 
        'drop *_iterativeCone5GenJets_*_*', 
        'drop *_ak4GenJets_*_*', 
        'drop *_ak7GenJets_*_*', 
        'drop *_ak8GenJets_*_*', 
        'drop *_ak4GenJetsNoNu_*_*', 
        'drop *_ak8GenJetsNoNu_*_*', 
        'drop *_genCandidatesForMET_*_*', 
        'drop *_genParticlesForMETAllVisible_*_*', 
        'drop *_genMetCalo_*_*', 
        'drop *_genMetCaloAndNonPrompt_*_*', 
        'drop *_genMetTrue_*_*', 
        'drop *_genMetIC5GenJs_*_*')
)

process.REPACKRAWEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *', 
        'drop FEDRawDataCollection_*_*_*', 
        'keep FEDRawDataCollection_rawDataRepacker_*_*', 
        'keep FEDRawDataCollection_virginRawDataRepacker_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'drop FEDRawDataCollection_source_*_*', 
        'drop FEDRawDataCollection_rawDataCollector_*_*'),
    splitLevel = cms.untracked.int32(0)
)

process.REPACKRAWSIMEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *', 
        'drop FEDRawDataCollection_*_*_*', 
        'keep FEDRawDataCollection_rawDataRepacker_*_*', 
        'keep FEDRawDataCollection_virginRawDataRepacker_*_*', 
        'keep  FEDRawDataCollection_rawDataCollector_*_*', 
        'keep  FEDRawDataCollection_source_*_*', 
        'drop *_hlt*_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep GlobalObjectMapRecord_hltGtStage2ObjectMap_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep triggerTriggerEvent_*_*_*', 
        'keep *_hltFEDSelectorL1_*_*', 
        'keep *_hltScoutingCaloPacker_*_*', 
        'keep *_hltScoutingEgammaPacker_*_*', 
        'keep *_hltScoutingMuonPackerCalo_*_*', 
        'keep *_hltScoutingMuonPacker_*_*', 
        'keep *_hltScoutingPFPacker_*_*', 
        'keep *_hltScoutingPrimaryVertexPackerCaloMuon_*_*', 
        'keep *_hltScoutingPrimaryVertexPacker_*_*', 
        'keep edmTriggerResults_*_*_*', 
        'keep *_g4SimHits_*_*', 
        'keep edmHepMCProduct_source_*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep CSCDetIdCSCComparatorDigiMuonDigiCollection_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*', 
        'keep EBSrFlagsSorted_simEcalDigis_*_*', 
        'keep EESrFlagsSorted_simEcalDigis_*_*', 
        'keep CrossingFramePlaybackInfoNew_*_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int6stdbitsetstdpairs_*_AffectedAPVList_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'keep *_mix_MergedTrackTruth_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*', 
        'keep recoGenJets_*_*_*', 
        'keep *_genParticle_*_*', 
        'keep recoGenMETs_*_*_*', 
        'keep FEDRawDataCollection_source_*_*', 
        'keep FEDRawDataCollection_rawDataCollector_*_*', 
        'keep *_MEtoEDMConverter_*_*', 
        'keep *_randomEngineStateProducer_*_*', 
        'keep *_logErrorHarvester_*_*', 
        'drop FEDRawDataCollection_source_*_*', 
        'drop FEDRawDataCollection_rawDataCollector_*_*'),
    splitLevel = cms.untracked.int32(0)
)

process.RESIMEventContent = cms.PSet(
    inputCommands = cms.untracked.vstring('drop *', 
        'keep *_randomEngineStateProducer_*_*', 
        'keep LHERunInfoProduct_*_*_*', 
        'keep LHEEventProduct_*_*_*', 
        'keep GenRunInfoProduct_generator_*_*', 
        'keep GenLumiInfoHeader_generator_*_*', 
        'keep GenLumiInfoProduct_generator_*_*', 
        'keep GenEventInfoProduct_generator_*_*', 
        'keep edmHepMCProduct_generatorSmeared_*_*', 
        'keep GenFilterInfo_*_*_*', 
        'keep *_genParticles_*_*')
)

process.RecoBTagAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_softPFElectronBJetTags_*_*', 
        'keep *_softPFMuonBJetTags_*_*', 
        'keep *_pfTrackCountingHighEffBJetTags_*_*', 
        'keep *_pfJetProbabilityBJetTags_*_*', 
        'keep *_pfJetBProbabilityBJetTags_*_*', 
        'keep *_pfSimpleSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfCombinedSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfGhostTrackBJetTags_*_*', 
        'keep *_pfCombinedMVAV2BJetTags_*_*', 
        'keep *_inclusiveCandidateSecondaryVertices_*_*', 
        'keep *_inclusiveCandidateSecondaryVerticesCvsL_*_*', 
        'keep *_pfCombinedCvsLJetTags_*_*', 
        'keep *_pfCombinedCvsBJetTags_*_*', 
        'keep *_pfChargeBJetTags_*_*', 
        'keep *_pfDeepCSVJetTags_*_*', 
        'keep *_pfDeepCMVAJetTags_*_*')
)

process.RecoBTagFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_softPFMuonsTagInfos_*_*', 
        'keep *_softPFElectronsTagInfos_*_*', 
        'keep *_softPFElectronBJetTags_*_*', 
        'keep *_softPFMuonBJetTags_*_*', 
        'keep *_pfImpactParameterTagInfos_*_*', 
        'keep *_pfTrackCountingHighEffBJetTags_*_*', 
        'keep *_pfJetProbabilityBJetTags_*_*', 
        'keep *_pfJetBProbabilityBJetTags_*_*', 
        'keep *_pfSecondaryVertexTagInfos_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderTagInfos_*_*', 
        'keep *_pfSimpleSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfCombinedSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfGhostTrackVertexTagInfos_*_*', 
        'keep *_pfGhostTrackBJetTags_*_*', 
        'keep *_pfCombinedMVAV2BJetTags_*_*', 
        'keep *_inclusiveCandidateSecondaryVertices_*_*', 
        'keep *_inclusiveCandidateSecondaryVerticesCvsL_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderCvsLTagInfos_*_*', 
        'keep *_pfCombinedCvsLJetTags_*_*', 
        'keep *_pfCombinedCvsBJetTags_*_*', 
        'keep *_pfChargeBJetTags_*_*', 
        'keep *_pfDeepCSVJetTags_*_*', 
        'keep *_pfDeepCMVAJetTags_*_*')
)

process.RecoBTagRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_softPFMuonsTagInfos_*_*', 
        'keep *_softPFElectronsTagInfos_*_*', 
        'keep *_softPFElectronBJetTags_*_*', 
        'keep *_softPFMuonBJetTags_*_*', 
        'keep *_pfImpactParameterTagInfos_*_*', 
        'keep *_pfTrackCountingHighEffBJetTags_*_*', 
        'keep *_pfJetProbabilityBJetTags_*_*', 
        'keep *_pfJetBProbabilityBJetTags_*_*', 
        'keep *_pfSecondaryVertexTagInfos_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderTagInfos_*_*', 
        'keep *_pfSimpleSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfSimpleInclusiveSecondaryVertexHighEffBJetTags_*_*', 
        'keep *_pfCombinedSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfCombinedInclusiveSecondaryVertexV2BJetTags_*_*', 
        'keep *_pfGhostTrackVertexTagInfos_*_*', 
        'keep *_pfGhostTrackBJetTags_*_*', 
        'keep *_pfCombinedMVAV2BJetTags_*_*', 
        'keep *_inclusiveCandidateSecondaryVertices_*_*', 
        'keep *_inclusiveCandidateSecondaryVerticesCvsL_*_*', 
        'keep *_pfInclusiveSecondaryVertexFinderCvsLTagInfos_*_*', 
        'keep *_pfCombinedCvsLJetTags_*_*', 
        'keep *_pfCombinedCvsBJetTags_*_*', 
        'keep *_pfChargeBJetTags_*_*', 
        'keep *_pfDeepCSVJetTags_*_*', 
        'keep *_pfDeepCMVAJetTags_*_*')
)

process.RecoBTauAOD = cms.PSet(
    outputCommands = cms.untracked.vstring()
)

process.RecoBTauFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring()
)

process.RecoBTauRECO = cms.PSet(
    outputCommands = cms.untracked.vstring()
)

process.RecoCTPPSAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep TotemTriggerCounters_totemTriggerRawToDigi_*_*', 
        'keep TotemFEDInfos_totemRPRawToDigi_*_*', 
        'keep TotemRPDigiedmDetSetVector_totemRPRawToDigi_*_*', 
        'keep TotemVFATStatusedmDetSetVector_totemRPRawToDigi_*_*', 
        'keep TotemRPClusteredmDetSetVector_totemRPClusterProducer_*_*', 
        'keep TotemRPRecHitedmDetSetVector_totemRPRecHitProducer_*_*', 
        'keep TotemRPUVPatternedmDetSetVector_totemRPUVPatternFinder_*_*', 
        'keep TotemRPLocalTrackedmDetSetVector_totemRPLocalTrackFitter_*_*', 
        'keep TotemFEDInfos_ctppsDiamondRawToDigi_*_*', 
        'keep CTPPSDiamondDigiedmDetSetVector_ctppsDiamondRawToDigi_*_*', 
        'keep TotemVFATStatusedmDetSetVector_ctppsDiamondRawToDigi_*_*', 
        'keep CTPPSDiamondRecHitedmDetSetVector_ctppsDiamondRecHits_*_*', 
        'keep CTPPSDiamondLocalTrackedmDetSetVector_ctppsDiamondLocalTracks_*_*', 
        'keep CTPPSPixelDigiedmDetSetVector_ctppsPixelDigis_*_*', 
        'keep CTPPSPixelClusteredmDetSetVector_ctppsPixelClusters_*_*', 
        'keep CTPPSPixelRecHitedmDetSetVector_ctppsPixelRecHits_*_*', 
        'keep CTPPSPixelLocalTrackedmDetSetVector_ctppsPixelLocalTracks_*_*', 
        'keep CTPPSLocalTrackLites_ctppsLocalTrackLiteProducer_*_*')
)

process.RecoCTPPSFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep TotemTriggerCounters_totemTriggerRawToDigi_*_*', 
        'keep TotemFEDInfos_totemRPRawToDigi_*_*', 
        'keep TotemRPDigiedmDetSetVector_totemRPRawToDigi_*_*', 
        'keep TotemVFATStatusedmDetSetVector_totemRPRawToDigi_*_*', 
        'keep TotemRPClusteredmDetSetVector_totemRPClusterProducer_*_*', 
        'keep TotemRPRecHitedmDetSetVector_totemRPRecHitProducer_*_*', 
        'keep TotemRPUVPatternedmDetSetVector_totemRPUVPatternFinder_*_*', 
        'keep TotemRPLocalTrackedmDetSetVector_totemRPLocalTrackFitter_*_*', 
        'keep TotemFEDInfos_ctppsDiamondRawToDigi_*_*', 
        'keep CTPPSDiamondDigiedmDetSetVector_ctppsDiamondRawToDigi_*_*', 
        'keep TotemVFATStatusedmDetSetVector_ctppsDiamondRawToDigi_*_*', 
        'keep CTPPSDiamondRecHitedmDetSetVector_ctppsDiamondRecHits_*_*', 
        'keep CTPPSDiamondLocalTrackedmDetSetVector_ctppsDiamondLocalTracks_*_*', 
        'keep CTPPSPixelDigiedmDetSetVector_ctppsPixelDigis_*_*', 
        'keep CTPPSPixelClusteredmDetSetVector_ctppsPixelClusters_*_*', 
        'keep CTPPSPixelRecHitedmDetSetVector_ctppsPixelRecHits_*_*', 
        'keep CTPPSPixelLocalTrackedmDetSetVector_ctppsPixelLocalTracks_*_*', 
        'keep CTPPSLocalTrackLites_ctppsLocalTrackLiteProducer_*_*')
)

process.RecoCTPPSRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep TotemTriggerCounters_totemTriggerRawToDigi_*_*', 
        'keep TotemFEDInfos_totemRPRawToDigi_*_*', 
        'keep TotemRPDigiedmDetSetVector_totemRPRawToDigi_*_*', 
        'keep TotemVFATStatusedmDetSetVector_totemRPRawToDigi_*_*', 
        'keep TotemRPClusteredmDetSetVector_totemRPClusterProducer_*_*', 
        'keep TotemRPRecHitedmDetSetVector_totemRPRecHitProducer_*_*', 
        'keep TotemRPUVPatternedmDetSetVector_totemRPUVPatternFinder_*_*', 
        'keep TotemRPLocalTrackedmDetSetVector_totemRPLocalTrackFitter_*_*', 
        'keep TotemFEDInfos_ctppsDiamondRawToDigi_*_*', 
        'keep CTPPSDiamondDigiedmDetSetVector_ctppsDiamondRawToDigi_*_*', 
        'keep TotemVFATStatusedmDetSetVector_ctppsDiamondRawToDigi_*_*', 
        'keep CTPPSDiamondRecHitedmDetSetVector_ctppsDiamondRecHits_*_*', 
        'keep CTPPSDiamondLocalTrackedmDetSetVector_ctppsDiamondLocalTracks_*_*', 
        'keep CTPPSPixelDigiedmDetSetVector_ctppsPixelDigis_*_*', 
        'keep CTPPSPixelClusteredmDetSetVector_ctppsPixelClusters_*_*', 
        'keep CTPPSPixelRecHitedmDetSetVector_ctppsPixelRecHits_*_*', 
        'keep CTPPSPixelLocalTrackedmDetSetVector_ctppsPixelLocalTracks_*_*', 
        'keep CTPPSLocalTrackLites_ctppsLocalTrackLiteProducer_*_*')
)

process.RecoEcalAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_selectDigi_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEB_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEE_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsES_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep recoCaloClusters_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_hybridSuperClusters_uncleanOnlyHybridSuperClusters_*', 
        'keep recoCaloClusters_multi5x5SuperClusters_multi5x5EndcapBasicClusters_*', 
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_*', 
        'keep recoSuperClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoCaloClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoSuperClusters_particleFlowSuperClusterOOTECAL_*_*', 
        'keep recoCaloClusters_particleFlowSuperClusterOOTECAL_*_*')
)

process.RecoEcalFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_selectDigi_*_*', 
        'keep *_reducedEcalRecHitsEB_*_*', 
        'keep *_reducedEcalRecHitsEE_*_*', 
        'keep *_reducedEcalRecHitsES_*_*', 
        'keep *_interestingEcalDetId*_*_*', 
        'keep *_ecalWeightUncalibRecHit_*_*', 
        'keep *_ecalPreshowerRecHit_*_*', 
        'keep *_hybridSuperClusters_*_*', 
        'keep *_correctedHybridSuperClusters_*_*', 
        'keep *_multi5x5*_*_*', 
        'keep *_correctedMulti5x5*_*_*', 
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_*', 
        'keep *_particleFlowSuperClusterECAL_*_*', 
        'keep *_particleFlowSuperClusterOOTECAL_*_*', 
        'drop recoBasicClusters_multi5x5BasicClusters_multi5x5BarrelBasicClusters_*', 
        'drop recoSuperClusters_multi5x5SuperClusters_multi5x5BarrelSuperClusters_*')
)

process.RecoEcalRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_selectDigi_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEE_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsEB_*_*', 
        'keep EcalRecHitsSorted_reducedEcalRecHitsES_*_*', 
        'keep *_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep *_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerClusterShapes_multi5x5PreshowerClusterShape_*_*', 
        'keep *_particleFlowSuperClusterECAL_*_*', 
        'keep *_particleFlowSuperClusterOOTECAL_*_*', 
        'drop recoClusterShapes_*_*_*', 
        'drop recoBasicClustersToOnerecoClusterShapesAssociation_*_*_*', 
        'drop recoBasicClusters_multi5x5BasicClusters_multi5x5BarrelBasicClusters_*', 
        'drop recoSuperClusters_multi5x5SuperClusters_multi5x5BarrelSuperClusters_*')
)

process.RecoEgammaAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoGsfElectronCores_gsfElectronCores_*_*', 
        'keep recoGsfElectronCores_gedGsfElectronCores_*_*', 
        'keep recoGsfElectrons_gsfElectrons_*_*', 
        'keep recoGsfElectrons_gedGsfElectrons_*_*', 
        'keep recoGsfElectronCores_uncleanedOnlyGsfElectronCores_*_*', 
        'keep recoGsfElectrons_uncleanedOnlyGsfElectrons_*_*', 
        'keep floatedmValueMap_eidRobustLoose_*_*', 
        'keep floatedmValueMap_eidRobustTight_*_*', 
        'keep floatedmValueMap_eidRobustHighEnergy_*_*', 
        'keep floatedmValueMap_eidLoose_*_*', 
        'keep floatedmValueMap_eidTight_*_*', 
        'keep *_egmGedGsfElectronPFIsolation_*_*', 
        'keep *_photonEcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonEcalPFClusterIsolationProducer_*_*', 
        'keep *_electronEcalPFClusterIsolationProducer_*_*', 
        'keep *_photonHcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonHcalPFClusterIsolationProducer_*_*', 
        'keep *_electronHcalPFClusterIsolationProducer_*_*', 
        'drop *_egmGsfElectronIDs_*_*', 
        'drop *_egmPhotonIDs_*_*', 
        'keep recoPhotonCores_gedPhotonCore_*_*', 
        'keep recoPhotons_gedPhotons_*_*', 
        'keep *_particleBasedIsolation_*_*', 
        'drop *_gedPhotons_valMapPFEgammaCandToPhoton_*', 
        'keep recoPhotonCores_photonCore_*_*', 
        'keep recoPhotons_photons_*_*', 
        'keep recoPhotonCores_ootPhotonCore_*_*', 
        'keep recoPhotons_ootPhotons_*_*', 
        'keep recoConversions_conversions_*_*', 
        'keep recoConversions_mustacheConversions_*_*', 
        'drop *_conversions_uncleanedConversions_*', 
        'keep recoConversions_allConversions_*_*', 
        'keep recoConversions_allConversionsOldEG_*_*', 
        'keep recoTracks_ckfOutInTracksFromConversions_*_*', 
        'keep recoTracks_ckfInOutTracksFromConversions_*_*', 
        'keep recoConversions_uncleanedOnlyAllConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep *_PhotonIDProdGED_*_*', 
        'keep *_hfRecoEcalCandidate_*_*', 
        'keep recoRecoEcalCandidates_hfRecoEcalCandidate_*_*', 
        'keep *_hfEMClusters_*_*', 
        'keep *_gedGsfElectronCores_*_*', 
        'keep *_gedGsfElectrons_*_*')
)

process.RecoEgammaFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_gsfElectronCores_*_*', 
        'keep *_gsfElectrons_*_*', 
        'keep *_uncleanedOnlyGsfElectronCores_*_*', 
        'keep *_uncleanedOnlyGsfElectrons_*_*', 
        'keep *_eidRobustLoose_*_*', 
        'keep *_eidRobustTight_*_*', 
        'keep *_eidRobustHighEnergy_*_*', 
        'keep *_eidLoose_*_*', 
        'keep *_eidTight_*_*', 
        'keep *_egmGedGsfElectronPF*Isolation_*_*', 
        'keep *_egmGsfElectronIDs_*_*', 
        'keep *_egmPhotonIDs_*_*', 
        'keep *_photonEcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonEcalPFClusterIsolationProducer_*_*', 
        'keep *_electronEcalPFClusterIsolationProducer_*_*', 
        'keep *_photonHcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonHcalPFClusterIsolationProducer_*_*', 
        'keep *_electronHcalPFClusterIsolationProducer_*_*', 
        'keep *_conversions_*_*', 
        'keep *_mustacheConversions_*_*', 
        'drop *_conversions_uncleanedConversions_*', 
        'keep *_gedPhotonCore_*_*', 
        'keep *_gedPhotonsTmp_*_*', 
        'keep *_gedPhotons_*_*', 
        'keep *_particleBasedIsolation_*_*', 
        'keep *_photonCore_*_*', 
        'keep *_photons_*_*', 
        'keep *_mustachePhotonCore_*_*', 
        'keep *_mustachePhotons_*_*', 
        'keep *_ootPhotonCore_*_*', 
        'keep *_ootPhotons_*_*', 
        'keep *_allConversions_*_*', 
        'keep *_allConversionsOldEG_*_*', 
        'keep *_ckfOutInTracksFrom*Conversions_*_*', 
        'keep *_ckfInOutTracksFrom*Conversions_*_*', 
        'keep *_uncleanedOnlyAllConversions_*_*', 
        'keep *_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep *_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep *_PhotonIDProdGED_*_*', 
        'keep *_hfRecoEcalCandidate_*_*', 
        'keep *_hfEMClusters_*_*', 
        'keep *_gedGsfElectronCores_*_*', 
        'keep *_gedGsfElectrons_*_*')
)

process.RecoEgammaRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoGsfElectronCores_gsfElectronCores_*_*', 
        'keep recoGsfElectronCores_gedGsfElectronCores_*_*', 
        'keep recoGsfElectrons_gsfElectrons_*_*', 
        'keep recoGsfElectrons_gedGsfElectrons_*_*', 
        'keep recoGsfElectronCores_uncleanedOnlyGsfElectronCores_*_*', 
        'keep recoGsfElectrons_uncleanedOnlyGsfElectrons_*_*', 
        'keep floatedmValueMap_eidRobustLoose_*_*', 
        'keep floatedmValueMap_eidRobustTight_*_*', 
        'keep floatedmValueMap_eidRobustHighEnergy_*_*', 
        'keep floatedmValueMap_eidLoose_*_*', 
        'keep floatedmValueMap_eidTight_*_*', 
        'keep *_egmGedGsfElectronPFIsolation_*_*', 
        'keep *_photonEcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonEcalPFClusterIsolationProducer_*_*', 
        'keep *_electronEcalPFClusterIsolationProducer_*_*', 
        'keep *_photonHcalPFClusterIsolationProducer_*_*', 
        'keep *_ootPhotonHcalPFClusterIsolationProducer_*_*', 
        'keep *_electronHcalPFClusterIsolationProducer_*_*', 
        'drop *_egmGsfElectronIDs_*_*', 
        'drop *_egmPhotonIDs_*_*', 
        'keep *_gedPhotonCore_*_*', 
        'keep *_gedPhotons_*_*', 
        'keep *_particleBasedIsolation_*_*', 
        'keep recoPhotons_mustachePhotons_*_*', 
        'keep recoPhotonCores_mustachePhotonCore_*_*', 
        'keep recoPhotons_ootPhotons_*_*', 
        'keep recoPhotonCores_ootPhotonCore_*_*', 
        'keep recoPhotons_photons_*_*', 
        'keep recoPhotonCores_photonCore_*_*', 
        'keep recoConversions_conversions_*_*', 
        'keep recoConversions_mustacheConversions_*_*', 
        'drop *_conversions_uncleanedConversions_*', 
        'drop *_gedPhotonsTmp_valMapPFEgammaCandToPhoton_*', 
        'keep recoConversions_allConversions_*_*', 
        'keep recoConversions_allConversionsOldEG_*_*', 
        'keep recoTracks_ckfOutInTracksFromConversions_*_*', 
        'keep recoTracks_ckfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_ckfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_ckfInOutTracksFromConversions_*_*', 
        'keep recoConversions_uncleanedOnlyAllConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTracks_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep recoTrackExtras_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfOutInTracksFromConversions_*_*', 
        'keep TrackingRecHitsOwned_uncleanedOnlyCkfInOutTracksFromConversions_*_*', 
        'keep *_PhotonIDProd_*_*', 
        'keep *_PhotonIDProdGED_*_*', 
        'keep recoRecoEcalCandidates_hfRecoEcalCandidate_*_*', 
        'keep *_hfEMClusters_*_*', 
        'keep *_gedGsfElectronCores_*_*', 
        'keep *_gedGsfElectrons_*_*')
)

process.RecoGenJetsAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_ak4GenJets_*_*', 
        'keep *_ak8GenJets_*_*', 
        'keep *_ak4GenJetsNoNu_*_*', 
        'keep *_ak8GenJetsNoNu_*_*', 
        'keep *_genParticle_*_*')
)

process.RecoGenJetsFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoGenJets_*_*_*', 
        'keep *_genParticle_*_*')
)

process.RecoGenJetsRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_ak4GenJets_*_*', 
        'keep *_ak8GenJets_*_*', 
        'keep *_ak4GenJetsNoNu_*_*', 
        'keep *_ak8GenJetsNoNu_*_*', 
        'keep *_genParticle_*_*')
)

process.RecoGenMETAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoGenMETs_*_*_*')
)

process.RecoGenMETFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoGenMETs_*_*_*')
)

process.RecoGenMETRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoGenMETs_*_*_*')
)

process.RecoHcalNoiseAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('drop recoHcalNoiseRBXs_hcalnoise_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*')
)

process.RecoHcalNoiseFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoHcalNoiseRBXs_hcalnoise_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*')
)

process.RecoHcalNoiseRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoHcalNoiseRBXs_hcalnoise_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*')
)

process.RecoJetsAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_ak4CaloJets_*_*', 
        'keep *_ak4PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHSSoftDrop_*_*', 
        'keep *_cmsTopTagPFJetsCHS_*_*', 
        'keep *_ak4PFJets_*_*', 
        'keep *_JetPlusTrackZSPCorJetAntiKt4_*_*', 
        'keep *_ak4TrackJets_*_*', 
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', 
        'keep *_caloTowers_*_*', 
        'keep *_CastorTowerReco_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertex_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertexPF_*_*', 
        'keep *_ak4JetTracksAssociatorExplicit_*_*', 
        'keep *_ak4JetExtender_*_*', 
        'keep *_ak4JetID_*_*', 
        'keep *_ak5CastorJets_*_*', 
        'keep *_ak5CastorJetID_*_*', 
        'keep *_ak7CastorJets_*_*', 
        'keep *_ak7CastorJetID_*_*', 
        'keep *_fixedGridRhoAll_*_*', 
        'keep *_fixedGridRhoFastjetAll_*_*', 
        'keep *_fixedGridRhoFastjetAllTmp_*_*', 
        'keep *_fixedGridRhoFastjetCentral_*_*', 
        'keep *_fixedGridRhoFastjetAllCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralChargedPileUp_*_*', 
        'keep *_fixedGridRhoFastjetCentralNeutral_*_*', 
        'drop doubles_*Jets_rhos_*', 
        'drop doubles_*Jets_sigmas_*', 
        'keep *_ak8PFJetsCHSSoftDropMass_*_*')
)

process.RecoJetsFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoCaloJets_*_*_*', 
        'keep recoPFJets_*_*_*', 
        'keep recoTrackJets_*_*_*', 
        'keep recoJPTJets_*_*_*', 
        'keep recoBasicJets_*_*_*', 
        'keep *_caloTowers_*_*', 
        'keep *_towerMaker_*_*', 
        'keep *_CastorTowerReco_*_*', 
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', 
        'keep *_kt4JetTracksAssociatorAtVertex_*_*', 
        'keep *_kt4JetTracksAssociatorAtCaloFace_*_*', 
        'keep *_kt4JetExtender_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertex*_*_*', 
        'keep *_ak4JetTracksAssociatorAtCaloFace*_*_*', 
        'keep *_ak4JetExtender_*_*', 
        'keep *_ak4JetTracksAssociatorExplicit_*_*', 
        'keep *_ak7JetTracksAssociatorAtVertex*_*_*', 
        'keep *_ak7JetTracksAssociatorAtCaloFace*_*_*', 
        'keep *_ak7JetExtender_*_*', 
        'keep *_*JetID_*_*', 
        'keep *_kt4CaloJets_*_*', 
        'keep *_kt6CaloJets_*_*', 
        'keep *_ak4CaloJets_*_*', 
        'keep *_ak5CaloJets_*_*', 
        'keep *_ak7CaloJets_*_*', 
        'keep *_kt4PFJets_*_*', 
        'keep *_kt6PFJets_*_*', 
        'keep *_ak4PFJets_*_*', 
        'keep *_ak5PFJets_*_*', 
        'keep *_ak7PFJets_*_*', 
        'keep *_JetPlusTrackZSPCorJetAntiKt4_*_*', 
        'keep *_ak4TrackJets_*_*', 
        'keep *_kt4TrackJets_*_*', 
        'keep *_ak5CastorJets_*_*', 
        'keep *_ak5CastorJetID_*_*', 
        'keep *_ak7CastorJets_*_*', 
        'keep *_ak7CastorJetID_*_*', 
        'keep *_fixedGridRho*_*_*', 
        'keep *_ca*Mass_*_*', 
        'keep *_ak*Mass_*_*')
)

process.RecoJetsRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_ak4CaloJets_*_*', 
        'keep *_ak4PFJets_*_*', 
        'keep *_ak4PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHS_*_*', 
        'keep *_ak8PFJetsCHSSoftDrop_*_*', 
        'keep *_cmsTopTagPFJetsCHS_*_*', 
        'keep *_JetPlusTrackZSPCorJetAntiKt4_*_*', 
        'keep *_ak4TrackJets_*_*', 
        'keep recoRecoChargedRefCandidates_trackRefsForJets_*_*', 
        'keep *_caloTowers_*_*', 
        'keep *_towerMaker_*_*', 
        'keep *_CastorTowerReco_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertex_*_*', 
        'keep *_ak4JetTracksAssociatorAtVertexPF_*_*', 
        'keep *_ak4JetTracksAssociatorAtCaloFace_*_*', 
        'keep *_ak4JetTracksAssociatorExplicit_*_*', 
        'keep *_ak4JetExtender_*_*', 
        'keep *_ak4JetID_*_*', 
        'keep *_ak5CastorJets_*_*', 
        'keep *_ak5CastorJetID_*_*', 
        'keep *_ak7CastorJets_*_*', 
        'keep *_ak7CastorJetID_*_*', 
        'keep *_fixedGridRhoAll_*_*', 
        'keep *_fixedGridRhoFastjetAll_*_*', 
        'keep *_fixedGridRhoFastjetAllTmp_*_*', 
        'keep *_fixedGridRhoFastjetAllCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentral_*_*', 
        'keep *_fixedGridRhoFastjetCentralCalo_*_*', 
        'keep *_fixedGridRhoFastjetCentralChargedPileUp_*_*', 
        'keep *_fixedGridRhoFastjetCentralNeutral_*_*', 
        'keep *_ak8PFJetsCHSSoftDropMass_*_*')
)

process.RecoLocalCaloAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_castorreco_*_*', 
        'keep *_reducedHcalRecHits_*_*', 
        'keep HcalUnpackerReport_castorDigis_*_*', 
        'keep HcalUnpackerReport_hcalDigiAlCaMB_*_*', 
        'keep HcalUnpackerReport_hcalDigis_*_*')
)

process.RecoLocalCaloFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_hbhereco_*_*', 
        'keep *_hbheprereco_*_*', 
        'keep *_hfprereco_*_*', 
        'keep *_hfreco_*_*', 
        'keep *_horeco_*_*', 
        'keep HBHERecHitsSorted_hbherecoMB_*_*', 
        'keep HBHERecHitsSorted_hbheprerecoMB_*_*', 
        'keep HORecHitsSorted_horecoMB_*_*', 
        'keep HFRecHitsSorted_hfrecoMB_*_*', 
        'keep ZDCDataFramesSorted_*Digis_*_*', 
        'keep ZDCRecHitsSorted_*_*_*', 
        'keep QIE10DataFrameHcalDataFrameContainer_hcalDigis_ZDC_*', 
        'keep *_reducedHcalRecHits_*_*', 
        'keep *_castorreco_*_*', 
        'keep HcalUnpackerReport_*_*_*', 
        'keep *_ecalMultiFitUncalibRecHit_*_*', 
        'keep *_ecalPreshowerRecHit_*_*', 
        'keep *_ecalRecHit_*_*', 
        'keep EBSrFlagsSorted_ecalDigis__*', 
        'keep EESrFlagsSorted_ecalDigis__*')
)

process.RecoLocalCaloRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_hbhereco_*_*', 
        'keep *_hbheprereco_*_*', 
        'keep *_hfprereco_*_*', 
        'keep *_hfreco_*_*', 
        'keep *_horeco_*_*', 
        'keep HBHERecHitsSorted_hbherecoMB_*_*', 
        'keep HORecHitsSorted_horecoMB_*_*', 
        'keep HFRecHitsSorted_hfrecoMB_*_*', 
        'keep ZDCDataFramesSorted_hcalDigis_*_*', 
        'keep ZDCDataFramesSorted_castorDigis_*_*', 
        'keep QIE10DataFrameHcalDataFrameContainer_hcalDigis_ZDC_*', 
        'keep ZDCRecHitsSorted_zdcreco_*_*', 
        'keep *_reducedHcalRecHits_*_*', 
        'keep *_castorreco_*_*', 
        'keep HcalUnpackerReport_castorDigis_*_*', 
        'keep HcalUnpackerReport_hcalDigiAlCaMB_*_*', 
        'keep HcalUnpackerReport_hcalDigis_*_*', 
        'keep *_ecalPreshowerRecHit_*_*', 
        'keep *_ecalRecHit_*_*', 
        'keep *_ecalCompactTrigPrim_*_*', 
        'keep *_ecalTPSkim_*_*', 
        'keep EBSrFlagsSorted_ecalDigis__*', 
        'keep EESrFlagsSorted_ecalDigis__*')
)

process.RecoLocalFastTimeAOD = cms.PSet(
    outputCommands = cms.untracked.vstring()
)

process.RecoLocalFastTimeFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_ftlUncalibratedRecHits_*_*', 
        'keep *_ftlRecHits_*_*')
)

process.RecoLocalFastTimeRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_ftlRecHits_*_*')
)

process.RecoLocalMuonAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_dt4DSegments_*_*', 
        'keep *_dt4DCosmicSegments_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*')
)

process.RecoLocalMuonFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_dt1DRecHits_*_*', 
        'keep *_dt4DSegments_*_*', 
        'keep *_csc2DRecHits_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*')
)

process.RecoLocalMuonRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_dt1DRecHits_*_*', 
        'keep *_dt4DSegments_*_*', 
        'keep *_dt1DCosmicRecHits_*_*', 
        'keep *_dt4DCosmicSegments_*_*', 
        'keep *_csc2DRecHits_*_*', 
        'keep *_cscSegments_*_*', 
        'keep *_rpcRecHits_*_*')
)

process.RecoLocalTrackerAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep ClusterSummary_clusterSummaryProducer_*_*')
)

process.RecoLocalTrackerFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep DetIdedmEDCollection_siStripDigis_*_*', 
        'keep DetIdedmEDCollection_siPixelDigis_*_*', 
        'keep PixelFEDChanneledmNewDetSetVector_siPixelDigis_*_*', 
        'keep *_siPixelClusters_*_*', 
        'keep *_siStripClusters_*_*', 
        'keep *_clusterSummaryProducer_*_*')
)

process.RecoLocalTrackerRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep DetIdedmEDCollection_siStripDigis_*_*', 
        'keep DetIdedmEDCollection_siPixelDigis_*_*', 
        'keep PixelFEDChanneledmNewDetSetVector_siPixelDigis_*_*', 
        'keep *_siPixelClusters_*_*', 
        'keep *_siStripClusters_*_*', 
        'keep ClusterSummary_clusterSummaryProducer_*_*')
)

process.RecoMETAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoCaloMETs_caloMet_*_*', 
        'keep recoCaloMETs_caloMetBE_*_*', 
        'keep recoCaloMETs_caloMetBEFO_*_*', 
        'keep recoCaloMETs_caloMetM_*_*', 
        'keep recoPFMETs_pfMet_*_*', 
        'keep recoPFMETs_pfChMet_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_*_*', 
        'drop recoHcalNoiseRBXs_*_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*', 
        'keep recoGlobalHaloData_GlobalHaloData_*_*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*')
)

process.RecoMETFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoCaloMETs_caloMet_*_*', 
        'keep recoCaloMETs_caloMetBE_*_*', 
        'keep recoCaloMETs_caloMetBEFO_*_*', 
        'keep recoCaloMETs_caloMetM_*_*', 
        'keep recoPFMETs_pfMet_*_*', 
        'keep recoPFMETs_pfChMet_*_*', 
        'keep recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_*_*', 
        'keep recoHcalNoiseRBXs_hcalnoise_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*', 
        'keep *HaloData_*_*_*', 
        'keep *BeamHaloSummary_BeamHaloSummary_*_*')
)

process.RecoMETRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoCaloMETs_caloMet_*_*', 
        'keep recoCaloMETs_caloMetBE_*_*', 
        'keep recoCaloMETs_caloMetBEFO_*_*', 
        'keep recoCaloMETs_caloMetM_*_*', 
        'keep recoPFMETs_pfMet_*_*', 
        'keep recoPFMETs_pfChMet_*_*', 
        'keep recoPFMETs_pfMetEI_*_*', 
        'keep recoMuonMETCorrectionDataedmValueMap_muonMETValueMapProducer_*_*', 
        'keep recoHcalNoiseRBXs_hcalnoise_*_*', 
        'keep HcalNoiseSummary_hcalnoise_*_*', 
        'keep recoCSCHaloData_CSCHaloData_*_*', 
        'keep recoEcalHaloData_EcalHaloData_*_*', 
        'keep recoGlobalHaloData_GlobalHaloData_*_*', 
        'keep recoHcalHaloData_HcalHaloData_*_*', 
        'keep recoBeamHaloSummary_BeamHaloSummary_*_*')
)

process.RecoMuonAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_muons_*_*', 
        'keep *_particleFlow_muons_*', 
        'drop *_muons_muons1stStep2muonsMap_*', 
        'drop recoIsoDepositedmValueMap_muons_*_*', 
        'drop doubleedmValueMap_muons_muPFIso*_*', 
        'keep recoTracks_standAloneMuons_*_*', 
        'keep recoTrackExtras_standAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_standAloneMuons_*_*', 
        'keep recoTracks_globalMuons_*_*', 
        'keep recoTrackExtras_globalMuons_*_*', 
        'keep recoTracks_tevMuons_*_*', 
        'keep recoTrackExtras_tevMuons_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_displacedTracks_*_*', 
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_*_*', 
        'keep recoTracks_displacedGlobalMuons_*_*', 
        'keep recoTrackExtras_displacedGlobalMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedGlobalMuons_*_*', 
        'keep recoTracks_cosmicMuons_*_*', 
        'keep recoMuons_muonsFromCosmics_*_*', 
        'keep recoTracks_cosmicMuons1Leg_*_*', 
        'keep recoMuons_muonsFromCosmics1Leg_*_*', 
        'keep recoTracks_refittedStandAloneMuons_*_*', 
        'keep recoTrackExtras_refittedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_refittedStandAloneMuons_*_*', 
        'keep recoTracks_displacedStandAloneMuons__*', 
        'keep recoTrackExtras_displacedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedStandAloneMuons_*_*')
)

process.RecoMuonFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_MuonSeed_*_*', 
        'keep *_ancientMuonSeed_*_*', 
        'keep *_displacedMuonSeeds_*_*', 
        'keep TrackingRecHitsOwned_globalMuons_*_*', 
        'keep TrackingRecHitsOwned_tevMuons_*_*', 
        'keep *_CosmicMuonSeed_*_*', 
        'keep recoTrackExtras_cosmicMuons_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons_*_*', 
        'keep recoTrackExtras_cosmicMuons1Leg_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons1Leg_*_*', 
        'keep recoTracks_cosmicsVetoTracks_*_*', 
        'keep *_muons_*_*', 
        'keep *_particleFlow_muons_*', 
        'drop *_muons_muons1stStep2muonsMap_*', 
        'drop recoIsoDepositedmValueMap_muons_*_*', 
        'drop doubleedmValueMap_muons_muPFIso*_*', 
        'keep recoTracks_standAloneMuons_*_*', 
        'keep recoTrackExtras_standAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_standAloneMuons_*_*', 
        'keep recoTracks_globalMuons_*_*', 
        'keep recoTrackExtras_globalMuons_*_*', 
        'keep recoTracks_tevMuons_*_*', 
        'keep recoTrackExtras_tevMuons_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_displacedTracks_*_*', 
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_*_*', 
        'keep recoTracks_displacedGlobalMuons_*_*', 
        'keep recoTrackExtras_displacedGlobalMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedGlobalMuons_*_*', 
        'keep recoTracks_cosmicMuons_*_*', 
        'keep recoMuons_muonsFromCosmics_*_*', 
        'keep recoTracks_cosmicMuons1Leg_*_*', 
        'keep recoMuons_muonsFromCosmics1Leg_*_*', 
        'keep recoTracks_refittedStandAloneMuons_*_*', 
        'keep recoTrackExtras_refittedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_refittedStandAloneMuons_*_*', 
        'keep recoTracks_displacedStandAloneMuons__*', 
        'keep recoTrackExtras_displacedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedStandAloneMuons_*_*', 
        'keep *_muIsoDepositTk_*_*', 
        'keep *_muIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muIsoDepositJets_*_*', 
        'keep *_muGlobalIsoDepositCtfTk_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muGlobalIsoDepositJets_*_*')
)

process.RecoMuonIsolationAOD = cms.PSet(
    outputCommands = cms.untracked.vstring()
)

process.RecoMuonIsolationFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_muIsoDepositTk_*_*', 
        'keep *_muIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muIsoDepositJets_*_*', 
        'keep *_muGlobalIsoDepositCtfTk_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muGlobalIsoDepositJets_*_*')
)

process.RecoMuonIsolationParamGlobal = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_muParamGlobalIsoDepositGsTk_*_*', 
        'keep *_muParamGlobalIsoDepositCalEcal_*_*', 
        'keep *_muParamGlobalIsoDepositCalHcal_*_*', 
        'keep *_muParamGlobalIsoDepositCtfTk_*_*', 
        'keep *_muParamGlobalIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muParamGlobalIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muParamGlobalIsoDepositJets_*_*')
)

process.RecoMuonIsolationRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_muIsoDepositTk_*_*', 
        'keep *_muIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muIsoDepositJets_*_*', 
        'keep *_muGlobalIsoDepositCtfTk_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muGlobalIsoDepositJets_*_*')
)

process.RecoMuonRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_MuonSeed_*_*', 
        'keep *_ancientMuonSeed_*_*', 
        'keep *_displacedMuonSeeds_*_*', 
        'keep TrackingRecHitsOwned_globalMuons_*_*', 
        'keep TrackingRecHitsOwned_tevMuons_*_*', 
        'keep *_CosmicMuonSeed_*_*', 
        'keep recoTrackExtras_cosmicMuons_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons_*_*', 
        'keep recoTrackExtras_cosmicMuons1Leg_*_*', 
        'keep TrackingRecHitsOwned_cosmicMuons1Leg_*_*', 
        'keep recoTracks_cosmicsVetoTracks_*_*', 
        'keep *_muons_*_*', 
        'keep *_particleFlow_muons_*', 
        'drop *_muons_muons1stStep2muonsMap_*', 
        'drop recoIsoDepositedmValueMap_muons_*_*', 
        'drop doubleedmValueMap_muons_muPFIso*_*', 
        'keep recoTracks_standAloneMuons_*_*', 
        'keep recoTrackExtras_standAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_standAloneMuons_*_*', 
        'keep recoTracks_globalMuons_*_*', 
        'keep recoTrackExtras_globalMuons_*_*', 
        'keep recoTracks_tevMuons_*_*', 
        'keep recoTrackExtras_tevMuons_*_*', 
        'keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_displacedTracks_*_*', 
        'keep recoTracksToOnerecoTracksAssociation_tevMuons_*_*', 
        'keep recoTracks_displacedGlobalMuons_*_*', 
        'keep recoTrackExtras_displacedGlobalMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedGlobalMuons_*_*', 
        'keep recoTracks_cosmicMuons_*_*', 
        'keep recoMuons_muonsFromCosmics_*_*', 
        'keep recoTracks_cosmicMuons1Leg_*_*', 
        'keep recoMuons_muonsFromCosmics1Leg_*_*', 
        'keep recoTracks_refittedStandAloneMuons_*_*', 
        'keep recoTrackExtras_refittedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_refittedStandAloneMuons_*_*', 
        'keep recoTracks_displacedStandAloneMuons__*', 
        'keep recoTrackExtras_displacedStandAloneMuons_*_*', 
        'keep TrackingRecHitsOwned_displacedStandAloneMuons_*_*', 
        'keep *_muIsoDepositTk_*_*', 
        'keep *_muIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muIsoDepositJets_*_*', 
        'keep *_muGlobalIsoDepositCtfTk_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorTowers_*_*', 
        'keep *_muGlobalIsoDepositCalByAssociatorHits_*_*', 
        'keep *_muGlobalIsoDepositJets_*_*')
)

process.RecoParticleFlowAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('drop CaloTowersSorted_towerMakerPF_*_*', 
        'drop *_pfElectronTranslator_*_*', 
        'keep recoPFRecHits_particleFlowClusterECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHCAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterPS_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHBHE_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitPS_Cleaned_*', 
        'keep recoCaloClusters_particleFlowEGamma_*_*', 
        'keep recoSuperClusters_particleFlowEGamma_*_*', 
        'keep recoCaloClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoSuperClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoConversions_particleFlowEGamma_*_*', 
        'keep recoPFCandidates_particleFlow_*_*', 
        'keep recoPFCandidates_particleFlowTmp_*_*', 
        'drop recoPFCandidates_particleFlowTmp__*', 
        'keep *_particleFlow_electrons_*', 
        'keep *_particleFlow_photons_*', 
        'keep *_particleFlow_muons_*', 
        'keep recoCaloClusters_pfElectronTranslator_*_*', 
        'keep recoPreshowerClusters_pfElectronTranslator_*_*', 
        'keep recoSuperClusters_pfElectronTranslator_*_*', 
        'keep recoCaloClusters_pfPhotonTranslator_*_*', 
        'keep recoPreshowerClusters_pfPhotonTranslator_*_*', 
        'keep recoSuperClusters_pfPhotonTranslator_*_*', 
        'keep recoPhotons_pfPhotonTranslator_*_*', 
        'keep recoPhotonCores_pfPhotonTranslator_*_*', 
        'keep recoConversions_pfPhotonTranslator_*_*', 
        'keep *_particleFlowPtrs_*_*', 
        'keep *_particleFlowTmpPtrs_*_*', 
        'keep *_chargedHadronPFTrackIsolation_*_*')
)

process.RecoParticleFlowFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('drop CaloTowersSorted_towerMakerPF_*_*', 
        'keep recoPFRecHits_particleFlowClusterECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHCAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterPS_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHBHE_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitPS_Cleaned_*', 
        'keep recoPFClusters_particleFlowClusterECAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHCAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHO_*_*', 
        'keep recoPFClusters_particleFlowClusterHF_*_*', 
        'keep recoPFClusters_particleFlowClusterPS_*_*', 
        'keep recoPFBlocks_particleFlowBlock_*_*', 
        'keep recoPFCandidates_particleFlowEGamma_*_*', 
        'keep recoCaloClusters_particleFlowEGamma_*_*', 
        'keep recoSuperClusters_particleFlowEGamma_*_*', 
        'keep recoConversions_particleFlowEGamma_*_*', 
        'keep recoPFCandidates_particleFlow_*_*', 
        'keep recoPFCandidates_particleFlowTmp_*_*', 
        'drop recoPFCandidates_particleFlowTmp__*', 
        'keep recoPFDisplacedVertexs_particleFlowDisplacedVertex_*_*', 
        'keep *_pfElectronTranslator_*_*', 
        'keep *_pfPhotonTranslator_*_*', 
        'keep *_particleFlow_electrons_*', 
        'keep *_particleFlow_photons_*', 
        'keep *_trackerDrivenElectronSeeds_preid_*', 
        'keep *_particleFlowPtrs_*_*', 
        'keep *_particleFlowTmpPtrs_*_*', 
        'keep *_chargedHadronPFTrackIsolation_*_*')
)

process.RecoParticleFlowRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('drop CaloTowersSorted_towerMakerPF_*_*', 
        'keep recoPFRecHits_particleFlowClusterECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHCAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowClusterPS_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitECAL_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHBHE_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHF_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitHO_Cleaned_*', 
        'keep recoPFRecHits_particleFlowRecHitPS_Cleaned_*', 
        'keep recoPFClusters_particleFlowClusterECAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHCAL_*_*', 
        'keep recoPFClusters_particleFlowClusterHO_*_*', 
        'keep recoPFClusters_particleFlowClusterHF_*_*', 
        'keep recoPFClusters_particleFlowClusterPS_*_*', 
        'keep recoPFBlocks_particleFlowBlock_*_*', 
        'keep recoPFCandidates_particleFlowEGamma_*_*', 
        'keep recoCaloClusters_particleFlowEGamma_*_*', 
        'keep recoSuperClusters_particleFlowEGamma_*_*', 
        'keep recoConversions_particleFlowEGamma_*_*', 
        'keep recoPFCandidates_particleFlow_*_*', 
        'keep recoPFCandidates_particleFlowTmp_electrons_*', 
        'keep recoPFCandidates_particleFlowTmp_*_*', 
        'drop recoPFCandidates_particleFlowTmp__*', 
        'keep recoPFDisplacedVertexs_particleFlowDisplacedVertex_*_*', 
        'keep *_pfElectronTranslator_*_*', 
        'keep *_pfPhotonTranslator_*_*', 
        'keep *_particleFlow_electrons_*', 
        'keep *_particleFlow_photons_*', 
        'keep *_particleFlow_muons_*', 
        'keep *_trackerDrivenElectronSeeds_preid_*', 
        'keep *_particleFlowPtrs_*_*', 
        'keep *_particleFlowTmpPtrs_*_*', 
        'keep *_chargedHadronPFTrackIsolation_*_*')
)

process.RecoPixelVertexingFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_pixelTracks_*_*', 
        'keep *_pixelVertices_*_*')
)

process.RecoPixelVertexingRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_pixelTracks_*_*', 
        'keep *_pixelVertices_*_*')
)

process.RecoTauTagAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoRecoTauPiZeros_hpsPFTauProducer_pizeros_*', 
        'keep recoPFTaus_hpsPFTauProducer_*_*', 
        'keep *_hpsPFTauDiscriminationByDeadECALElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFinding_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingNewDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingOldDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseChargedIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByTightMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSum_*_*', 
        'keep *_hpsPFTauPUcorrPtSum_*_*', 
        'keep *_hpsPFTauChargedIsoPtSum_*_*', 
        'keep *_hpsPFTauTransverseImpactParameters_*_*', 
        'keep *_hpsPFTauFootprintCorrection_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeight_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6rawElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6LooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6MediumElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6TightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VTightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauChargedIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeightdR03_*_*', 
        'keep *_hpsPFTauFootprintCorrectiondR03_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalConedR03_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*')
)

process.RecoTauTagFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_ak4PFJetsRecoTauPiZeros_*_*', 
        'keep *_hpsPFTauProducer_*_*', 
        'keep *_hpsPFTauDiscrimination*_*_*', 
        'keep *_hpsPFTau*PtSum_*_*', 
        'keep *_hpsPFTauTransverseImpactParameters_*_*')
)

process.RecoTauTagRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoRecoTauPiZeros_hpsPFTauProducer_pizeros_*', 
        'keep recoPFTaus_hpsPFTauProducer_*_*', 
        'keep *_hpsPFTauDiscriminationByDeadECALElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFinding_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingNewDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByDecayModeFindingOldDMs_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseChargedIsolation_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByRawCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3Hits_*_*', 
        'keep *_hpsPFTauDiscriminationByTightCombinedIsolationDBSumPtCorr3HitsdR03_*_*', 
        'keep *_hpsPFTauDiscriminationByTightMuonRejection3_*_*', 
        'keep *_hpsPFTauDiscriminationByPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSum_*_*', 
        'keep *_hpsPFTauPUcorrPtSum_*_*', 
        'keep *_hpsPFTauChargedIsoPtSum_*_*', 
        'keep *_hpsPFTauTransverseImpactParameters_*_*', 
        'keep *_hpsPFTauFootprintCorrection_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeight_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalCone_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6rawElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VLooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6LooseElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6MediumElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6TightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByMVA6VTightElectronRejection_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWoldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWoldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWnewDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWnewDMwLT_*_*', 
        'keep *_hpsPFTauChargedIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumdR03_*_*', 
        'keep *_hpsPFTauNeutralIsoPtSumWeightdR03_*_*', 
        'keep *_hpsPFTauFootprintCorrectiondR03_*_*', 
        'keep *_hpsPFTauPhotonPtSumOutsideSignalConedR03_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1DBdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByIsolationMVArun2v1PWdR03oldDMwLTraw_*_*', 
        'keep *_hpsPFTauDiscriminationByVLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByLooseIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByMediumIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*', 
        'keep *_hpsPFTauDiscriminationByVVTightIsolationMVArun2v1PWdR03oldDMwLT_*_*')
)

process.RecoTrackerAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoTracks_generalTracks_*_*', 
        'keep recoTracks_conversionStepTracks_*_*', 
        'keep recoTracks_beamhaloTracks_*_*', 
        'keep recoTracks_ctfPixelLess_*_*', 
        'keep *_dedxHarmonic2_*_*', 
        'keep *_dedxPixelHarmonic2_*_*', 
        'keep *_dedxHitInfo_*_*', 
        'keep *_trackExtrapolator_*_*', 
        'keep *_generalTracks_MVAValues_*', 
        'keep *_generalTracks_MVAVals_*')
)

process.RecoTrackerFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoTracks_generalTracks_*_*', 
        'keep recoTrackExtras_generalTracks_*_*', 
        'keep TrackingRecHitsOwned_generalTracks_*_*', 
        'keep *_generalTracks_MVAValues_*', 
        'keep TrackingRecHitsOwned_extraFromSeeds_*_*', 
        'keep uints_extraFromSeeds_*_*', 
        'keep recoTracks_beamhaloTracks_*_*', 
        'keep recoTrackExtras_beamhaloTracks_*_*', 
        'keep TrackingRecHitsOwned_beamhaloTracks_*_*', 
        'keep recoTracks_conversionStepTracks_*_*', 
        'keep recoTrackExtras_conversionStepTracks_*_*', 
        'keep TrackingRecHitsOwned_conversionStepTracks_*_*', 
        'keep *_ctfPixelLess_*_*', 
        'keep *_dedxTruncated40_*_*', 
        'keep *_dedxHitInfo_*_*', 
        'keep *_dedxHarmonic2_*_*', 
        'keep *_dedxPixelHarmonic2_*_*', 
        'keep *_trackExtrapolator_*_*', 
        'keep recoTracks_cosmicDCTracks_*_*', 
        'keep recoTrackExtras_cosmicDCTracks_*_*', 
        'keep TrackingRecHitsOwned_cosmicDCTracks_*_*')
)

process.RecoTrackerRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoTracks_generalTracks_*_*', 
        'keep recoTrackExtras_generalTracks_*_*', 
        'keep TrackingRecHitsOwned_generalTracks_*_*', 
        'keep *_generalTracks_MVAValues_*', 
        'keep *_generalTracks_MVAVals_*', 
        'keep TrackingRecHitsOwned_extraFromSeeds_*_*', 
        'keep uints_extraFromSeeds_*_*', 
        'keep recoTracks_beamhaloTracks_*_*', 
        'keep recoTrackExtras_beamhaloTracks_*_*', 
        'keep TrackingRecHitsOwned_beamhaloTracks_*_*', 
        'keep recoTracks_conversionStepTracks_*_*', 
        'keep recoTrackExtras_conversionStepTracks_*_*', 
        'keep TrackingRecHitsOwned_conversionStepTracks_*_*', 
        'keep *_ctfPixelLess_*_*', 
        'keep *_dedxTruncated40_*_*', 
        'keep *_dedxHitInfo_*_*', 
        'keep *_dedxHarmonic2_*_*', 
        'keep *_dedxPixelHarmonic2_*_*', 
        'keep *_trackExtrapolator_*_*', 
        'keep recoTracks_cosmicDCTracks_*_*', 
        'keep recoTrackExtras_cosmicDCTracks_*_*', 
        'keep TrackingRecHitsOwned_cosmicDCTracks_*_*')
)

process.RecoVertexAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep  *_offlinePrimaryVertices__*', 
        'keep *_offlinePrimaryVerticesWithBS_*_*', 
        'keep *_offlinePrimaryVerticesFromCosmicTracks_*_*', 
        'keep *_nuclearInteractionMaker_*_*', 
        'keep *_generalV0Candidates_*_*', 
        'keep *_inclusiveSecondaryVertices_*_*')
)

process.RecoVertexFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep  *_offlinePrimaryVertices__*', 
        'keep *_offlinePrimaryVerticesWithBS_*_*', 
        'keep *_offlinePrimaryVerticesFromCosmicTracks_*_*', 
        'keep *_nuclearInteractionMaker_*_*', 
        'keep *_generalV0Candidates_*_*', 
        'keep *_inclusiveSecondaryVertices_*_*')
)

process.RecoVertexRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep  *_offlinePrimaryVertices__*', 
        'keep *_offlinePrimaryVerticesWithBS_*_*', 
        'keep *_offlinePrimaryVerticesFromCosmicTracks_*_*', 
        'keep *_nuclearInteractionMaker_*_*', 
        'keep *_generalV0Candidates_*_*', 
        'keep *_inclusiveSecondaryVertices_*_*')
)

process.SimCalorimetryAOD = cms.PSet(
    outputCommands = cms.untracked.vstring()
)

process.SimCalorimetryFEVTDEBUG = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_simEcalDigis_*_*', 
        'keep *_simEcalPreshowerDigis_*_*', 
        'keep *_simEcalTriggerPrimitiveDigis_*_*', 
        'keep *_simEcalEBTriggerPrimitiveDigis_*_*', 
        'keep *_simHcalDigis_*_*', 
        'keep ZDCDataFramesSorted_simHcalUnsuppressedDigis_*_*', 
        'drop ZDCDataFramesSorted_mix_simHcalUnsuppressedDigis*_*', 
        'keep *_simHcalTriggerPrimitiveDigis_*_*', 
        'keep *_mix_HcalSamples_*', 
        'keep *_mixData_HcalSamples_*', 
        'keep *_mix_HcalHits_*', 
        'keep *_mixData_HcalHits_*', 
        'keep *_DMHcalTriggerPrimitiveDigis_*_*')
)

process.SimCalorimetryRAW = cms.PSet(
    outputCommands = cms.untracked.vstring('keep EBSrFlagsSorted_simEcalDigis_*_*', 
        'keep EESrFlagsSorted_simEcalDigis_*_*')
)

process.SimCalorimetryRECO = cms.PSet(
    outputCommands = cms.untracked.vstring()
)

process.SimG4CoreAOD = cms.PSet(
    outputCommands = cms.untracked.vstring()
)

process.SimG4CoreHLTAODSIM = cms.PSet(
    outputCommands = cms.untracked.vstring('keep SimVertexs_g4SimHits_*_*')
)

process.SimG4CoreRAW = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_g4SimHits_*_*', 
        'keep edmHepMCProduct_source_*_*')
)

process.SimG4CoreRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep edmHepMCProduct_source_*_*', 
        'keep SimTracks_g4SimHits_*_*', 
        'keep SimVertexs_g4SimHits_*_*')
)

process.SimGeneralAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep PileupSummaryInfos_*_*_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*')
)

process.SimGeneralFEVTDEBUG = cms.PSet(
    outputCommands = cms.untracked.vstring('drop *_trackingtruthprod_*_*', 
        'drop *_electrontruth_*_*', 
        'keep *_mix_MergedTrackTruth_*', 
        'keep CrossingFramePlaybackInfoNew_*_*_*')
)

process.SimGeneralRAW = cms.PSet(
    outputCommands = cms.untracked.vstring('keep CrossingFramePlaybackInfoNew_*_*_*', 
        'keep PileupSummaryInfos_*_*_*', 
        'keep int6stdbitsetstdpairs_*_AffectedAPVList_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*', 
        'keep *_mix_MergedTrackTruth_*')
)

process.SimGeneralRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep PileupSummaryInfos_*_*_*', 
        'keep int_*_bunchSpacing_*', 
        'keep *_genPUProtons_*_*')
)

process.SimMuonAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_muonSimClassifier_*_*')
)

process.SimMuonFEVTDEBUG = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_simMuonCSCDigis_*_*', 
        'keep *_simMuonDTDigis_*_*', 
        'keep *_simMuonRPCDigis_*_*')
)

process.SimMuonRAW = cms.PSet(
    outputCommands = cms.untracked.vstring('keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep CSCDetIdCSCComparatorDigiMuonDigiCollection_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*')
)

process.SimMuonRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep StripDigiSimLinkedmDetSetVector_simMuonCSCDigis_*_*', 
        'keep DTLayerIdDTDigiSimLinkMuonDigiCollection_simMuonDTDigis_*_*', 
        'keep RPCDigiSimLinkedmDetSetVector_simMuonRPCDigis_*_*', 
        'keep *_muonSimClassifier_*_*')
)

process.SimTrackerAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_allTrackMCMatch_*_*')
)

process.SimTrackerDEBUG = cms.PSet(
    outputCommands = cms.untracked.vstring('keep PixelDigiSimLinkedmDetSetVector_simSiPixelDigis_*_*', 
        'keep StripDigiSimLinkedmDetSetVector_simSiStripDigis_*_*', 
        'drop *_mix_simSiPixelDigis*_*', 
        'drop *_mix_simSiStripDigis*_*', 
        'keep *_allTrackMCMatch_*_*')
)

process.SimTrackerFEVTDEBUG = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_simSiPixelDigis_*_*', 
        'keep *_simSiStripDigis_*_*', 
        'drop *_mix_simSiPixelDigis*_*', 
        'drop *_mix_simSiStripDigis*_*', 
        'keep *_allTrackMCMatch_*_*', 
        'keep *_trackingParticleRecoTrackAsssociation_*_*', 
        'keep *_assoc2secStepTk_*_*', 
        'keep *_assoc2thStepTk_*_*', 
        'keep *_assoc2GsfTracks_*_*', 
        'keep *_assocOutInConversionTracks_*_*', 
        'keep *_assocInOutConversionTracks_*_*')
)

process.SimTrackerRAW = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_allTrackMCMatch_*_*')
)

process.SimTrackerRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_allTrackMCMatch_*_*')
)

process.TcdsEventContent = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_tcdsDigis_*_*')
)

process.TrackingToolsAOD = cms.PSet(
    outputCommands = cms.untracked.vstring('keep recoTracks_GsfGlobalElectronTest_*_*', 
        'keep recoGsfTracks_electronGsfTracks_*_*')
)

process.TrackingToolsFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_CkfElectronCandidates_*_*', 
        'keep *_GsfGlobalElectronTest_*_*', 
        'keep *_electronMergedSeeds_*_*', 
        'keep *_electronGsfTracks_*_*')
)

process.TrackingToolsRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_CkfElectronCandidates_*_*', 
        'keep *_GsfGlobalElectronTest_*_*', 
        'keep *_electronMergedSeeds_*_*', 
        'keep recoGsfTracks_electronGsfTracks_*_*', 
        'keep recoGsfTrackExtras_electronGsfTracks_*_*', 
        'keep recoTrackExtras_electronGsfTracks_*_*', 
        'keep TrackingRecHitsOwned_electronGsfTracks_*_*')
)

process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('CMS3'),
    name = cms.untracked.string('CMS3 test configuration'),
    version = cms.untracked.string('$Revision: 1.11 $')
)

process.ecalLocalRecoAOD = cms.PSet(
    outputCommands = cms.untracked.vstring()
)

process.ecalLocalRecoFEVT = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_ecalMultiFitUncalibRecHit_*_*', 
        'keep *_ecalPreshowerRecHit_*_*', 
        'keep *_ecalRecHit_*_*', 
        'keep EBSrFlagsSorted_ecalDigis__*', 
        'keep EESrFlagsSorted_ecalDigis__*')
)

process.ecalLocalRecoRECO = cms.PSet(
    outputCommands = cms.untracked.vstring('keep *_ecalPreshowerRecHit_*_*', 
        'keep *_ecalRecHit_*_*', 
        'keep *_ecalCompactTrigPrim_*_*', 
        'keep *_ecalTPSkim_*_*', 
        'keep EBSrFlagsSorted_ecalDigis__*', 
        'keep EESrFlagsSorted_ecalDigis__*')
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.mvaEleID_Fall17_iso_V1_producer_config = cms.PSet(
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clipLower = cms.vdouble(-inf, -inf, -1.0, -inf, -inf, 
        -inf, -inf, -inf, -inf, -inf, 
        -1.0, -inf, -inf, -inf, -inf, 
        -inf, -inf, -0.06, -0.6, -0.2, 
        -inf, -inf, -inf, -inf, -inf),
    clipUpper = cms.vdouble(inf, inf, 2.0, 5.0, inf, 
        inf, inf, inf, 10.0, 200.0, 
        inf, inf, inf, inf, 20.0, 
        20.0, inf, 0.06, 0.6, 0.2, 
        inf, inf, inf, inf, inf),
    conversionsAOD = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    ebSplit = cms.double(0.8),
    ebeeSplit = cms.double(1.479),
    mvaName = cms.string('ElectronMVAEstimatorRun2Fall17Iso'),
    mvaTag = cms.string('V1'),
    ptSplit = cms.double(10.0),
    varNames = cms.vstring('ele_oldsigmaietaieta', 
        'ele_oldsigmaiphiiphi', 
        'ele_oldcircularity', 
        'ele_oldr9', 
        'ele_scletawidth', 
        'ele_sclphiwidth', 
        'ele_oldhe', 
        'ele_kfhits', 
        'ele_kfchi2', 
        'ele_gsfchi2', 
        'ele_fbrem', 
        'ele_gsfhits', 
        'ele_expected_inner_hits', 
        'ele_conversionVertexFitProbability', 
        'ele_ep', 
        'ele_eelepout', 
        'ele_IoEmIop', 
        'ele_deltaetain', 
        'ele_deltaphiin', 
        'ele_deltaetaseed', 
        'ele_pfPhotonIso', 
        'ele_pfChargedHadIso', 
        'ele_pfNeutralHadIso', 
        'rho', 
        'ele_psEoverEraw'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_5_2017_puinfo_iso_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_5_2017_puinfo_iso_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_5_2017_puinfo_iso_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_10_2017_puinfo_iso_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_10_2017_puinfo_iso_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_10_2017_puinfo_iso_BDT.weights.xml')
)

process.mvaEleID_Fall17_noIso_V1_producer_config = cms.PSet(
    beamSpot = cms.InputTag("offlineBeamSpot"),
    clipLower = cms.vdouble(-inf, -inf, -1.0, -inf, -inf, 
        -inf, -inf, -inf, -inf, -inf, 
        -1.0, -inf, -inf, -inf, -inf, 
        -inf, -inf, -0.06, -0.6, -0.2, 
        -inf, -inf),
    clipUpper = cms.vdouble(inf, inf, 2.0, 5.0, inf, 
        inf, inf, inf, 10.0, 200.0, 
        inf, inf, inf, inf, 20.0, 
        20.0, inf, 0.06, 0.6, 0.2, 
        inf, inf),
    conversionsAOD = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    ebSplit = cms.double(0.8),
    ebeeSplit = cms.double(1.479),
    mvaName = cms.string('ElectronMVAEstimatorRun2Fall17NoIso'),
    mvaTag = cms.string('V1'),
    ptSplit = cms.double(10.0),
    varNames = cms.vstring('ele_oldsigmaietaieta', 
        'ele_oldsigmaiphiiphi', 
        'ele_oldcircularity', 
        'ele_oldr9', 
        'ele_scletawidth', 
        'ele_sclphiwidth', 
        'ele_oldhe', 
        'ele_kfhits', 
        'ele_kfchi2', 
        'ele_gsfchi2', 
        'ele_fbrem', 
        'ele_gsfhits', 
        'ele_expected_inner_hits', 
        'ele_conversionVertexFitProbability', 
        'ele_ep', 
        'ele_eelepout', 
        'ele_IoEmIop', 
        'ele_deltaetain', 
        'ele_deltaphiin', 
        'ele_deltaetaseed', 
        'rho', 
        'ele_psEoverEraw'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_5_2017_puinfo_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_5_2017_puinfo_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_5_2017_puinfo_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_10_2017_puinfo_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_10_2017_puinfo_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_10_2017_puinfo_BDT.weights.xml')
)

process.mvaEleID_PHYS14_PU20bx25_nonTrig_V1_producer_config = cms.PSet(
    mvaName = cms.string('ElectronMVAEstimatorRun2Phys14NonTrig'),
    mvaTag = cms.string('25nsV1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB1_5_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB2_5_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EE_5_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB1_10_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB2_10_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EE_10_oldscenario2phys14_BDT.weights.xml')
)

process.mvaEleID_PHYS14_PU20bx25_nonTrig_V1_wp80 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Phys14NonTrig25nsV1Categories"),
        mvaCuts = cms.vdouble(-0.253, 0.081, -0.081, 0.965, 0.917, 
            0.683),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Phys14NonTrig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp80'),
    isPOGApproved = cms.untracked.bool(False)
)

process.mvaEleID_PHYS14_PU20bx25_nonTrig_V1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Phys14NonTrig25nsV1Categories"),
        mvaCuts = cms.vdouble(-0.483, -0.267, -0.323, 0.933, 0.825, 
            0.337),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Phys14NonTrig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-PHYS14-PU20bx25-nonTrig-V1-wp90'),
    isPOGApproved = cms.untracked.bool(False)
)

process.mvaEleID_Spring15_25ns_Trig_V1_producer_config = cms.PSet(
    beamSpot = cms.InputTag("offlineBeamSpot"),
    conversionsAOD = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    mvaName = cms.string('ElectronMVAEstimatorRun2Spring15Trig'),
    mvaTag = cms.string('25nsV1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml')
)

process.mvaEleID_Spring15_25ns_Trig_V1_wp80 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
        mvaCuts = cms.vdouble(0.988153, 0.96791, 0.841729),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-25ns-Trig-V1-wp80'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring15_25ns_Trig_V1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
        mvaCuts = cms.vdouble(0.972153, 0.922126, 0.610764),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-25ns-Trig-V1-wp90'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring15_25ns_nonTrig_V1_producer_config = cms.PSet(
    beamSpot = cms.InputTag("offlineBeamSpot"),
    conversionsAOD = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    mvaName = cms.string('ElectronMVAEstimatorRun2Spring15NonTrig'),
    mvaTag = cms.string('25nsV1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml')
)

process.mvaEleID_Spring15_25ns_nonTrig_V1_wp80 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
        mvaCuts = cms.vdouble(0.287435, 0.221846, -0.303263, 0.967083, 0.929117, 
            0.726311),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-25ns-nonTrig-V1-wp80'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring15_25ns_nonTrig_V1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
        mvaCuts = cms.vdouble(-0.083313, -0.235222, -0.67099, 0.913286, 0.805013, 
            0.358969),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-25ns-nonTrig-V1-wp90'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring15_25ns_nonTrig_V1_wpLoose = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
        mvaCuts = cms.vdouble(-0.265, -0.556, -0.551, -0.072, -0.286, 
            -0.267),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-25ns-nonTrig-V1-wpLoose'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring15_50ns_Trig_V1_producer_config = cms.PSet(
    beamSpot = cms.InputTag("offlineBeamSpot"),
    conversionsAOD = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    mvaName = cms.string('ElectronMVAEstimatorRun2Spring15Trig'),
    mvaTag = cms.string('50nsV1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml')
)

process.mvaEleID_Spring15_50ns_Trig_V1_wp80 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig50nsV1Categories"),
        mvaCuts = cms.vdouble(0.981841, 0.946762, 0.79704),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig50nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-50ns-Trig-V1-wp80'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring15_50ns_Trig_V1_wp90 = cms.PSet(
    cutFlow = cms.VPSet(cms.PSet(
        cutName = cms.string('GsfEleMVACut'),
        isIgnored = cms.bool(False),
        mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig50nsV1Categories"),
        mvaCuts = cms.vdouble(0.953843, 0.849994, 0.514118),
        mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig50nsV1Values"),
        needsAdditionalProducts = cms.bool(True)
    )),
    idName = cms.string('mvaEleID-Spring15-50ns-Trig-V1-wp90'),
    isPOGApproved = cms.untracked.bool(True)
)

process.mvaEleID_Spring16_GeneralPurpose_V1_producer_config = cms.PSet(
    beamSpot = cms.InputTag("offlineBeamSpot"),
    conversionsAOD = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    mvaName = cms.string('ElectronMVAEstimatorRun2Spring16GeneralPurpose'),
    mvaTag = cms.string('V1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB1_10.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB2_10.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EE_10.weights.xml')
)

process.mvaEleID_Spring16_HZZ_V1_producer_config = cms.PSet(
    beamSpot = cms.InputTag("offlineBeamSpot"),
    conversionsAOD = cms.InputTag("allConversions"),
    conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
    mvaName = cms.string('ElectronMVAEstimatorRun2Spring16HZZ'),
    mvaTag = cms.string('V1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_5.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_5.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_5.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_10.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_10.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_10.weights.xml')
)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(False)
)

process.pset = cms.PSet(
    etaMax = cms.double(-2.901376),
    etaMin = cms.double(-5.2),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('egammaHFMinus'),
    px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
    py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
    type = cms.int32(7),
    varType = cms.int32(0)
)

process.mvaConfigsForEleProducer = cms.VPSet(cms.PSet(
    mvaName = cms.string('ElectronMVAEstimatorRun2Phys14NonTrig'),
    mvaTag = cms.string('25nsV1'),
    weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB1_5_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB2_5_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EE_5_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB1_10_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB2_10_oldscenario2phys14_BDT.weights.xml', 
        'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EE_10_oldscenario2phys14_BDT.weights.xml')
), 
    cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        conversionsAOD = cms.InputTag("allConversions"),
        conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
        mvaName = cms.string('ElectronMVAEstimatorRun2Spring15NonTrig'),
        mvaTag = cms.string('25nsV1'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml')
    ), 
    cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        conversionsAOD = cms.InputTag("allConversions"),
        conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
        mvaName = cms.string('ElectronMVAEstimatorRun2Spring15Trig'),
        mvaTag = cms.string('50nsV1'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml')
    ), 
    cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        conversionsAOD = cms.InputTag("allConversions"),
        conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
        mvaName = cms.string('ElectronMVAEstimatorRun2Spring15Trig'),
        mvaTag = cms.string('25nsV1'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml')
    ), 
    cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        conversionsAOD = cms.InputTag("allConversions"),
        conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
        mvaName = cms.string('ElectronMVAEstimatorRun2Spring16HZZ'),
        mvaTag = cms.string('V1'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_5.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_5.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_5.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_10.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_10.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_10.weights.xml')
    ), 
    cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        conversionsAOD = cms.InputTag("allConversions"),
        conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
        mvaName = cms.string('ElectronMVAEstimatorRun2Spring16GeneralPurpose'),
        mvaTag = cms.string('V1'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB1_10.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB2_10.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EE_10.weights.xml')
    ), 
    cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        clipLower = cms.vdouble(-inf, -inf, -1.0, -inf, -inf, 
            -inf, -inf, -inf, -inf, -inf, 
            -1.0, -inf, -inf, -inf, -inf, 
            -inf, -inf, -0.06, -0.6, -0.2, 
            -inf, -inf),
        clipUpper = cms.vdouble(inf, inf, 2.0, 5.0, inf, 
            inf, inf, inf, 10.0, 200.0, 
            inf, inf, inf, inf, 20.0, 
            20.0, inf, 0.06, 0.6, 0.2, 
            inf, inf),
        conversionsAOD = cms.InputTag("allConversions"),
        conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
        ebSplit = cms.double(0.8),
        ebeeSplit = cms.double(1.479),
        mvaName = cms.string('ElectronMVAEstimatorRun2Fall17NoIso'),
        mvaTag = cms.string('V1'),
        ptSplit = cms.double(10.0),
        varNames = cms.vstring('ele_oldsigmaietaieta', 
            'ele_oldsigmaiphiiphi', 
            'ele_oldcircularity', 
            'ele_oldr9', 
            'ele_scletawidth', 
            'ele_sclphiwidth', 
            'ele_oldhe', 
            'ele_kfhits', 
            'ele_kfchi2', 
            'ele_gsfchi2', 
            'ele_fbrem', 
            'ele_gsfhits', 
            'ele_expected_inner_hits', 
            'ele_conversionVertexFitProbability', 
            'ele_ep', 
            'ele_eelepout', 
            'ele_IoEmIop', 
            'ele_deltaetain', 
            'ele_deltaphiin', 
            'ele_deltaetaseed', 
            'rho', 
            'ele_psEoverEraw'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_5_2017_puinfo_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_5_2017_puinfo_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_5_2017_puinfo_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_10_2017_puinfo_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_10_2017_puinfo_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_10_2017_puinfo_BDT.weights.xml')
    ), 
    cms.PSet(
        beamSpot = cms.InputTag("offlineBeamSpot"),
        clipLower = cms.vdouble(-inf, -inf, -1.0, -inf, -inf, 
            -inf, -inf, -inf, -inf, -inf, 
            -1.0, -inf, -inf, -inf, -inf, 
            -inf, -inf, -0.06, -0.6, -0.2, 
            -inf, -inf, -inf, -inf, -inf),
        clipUpper = cms.vdouble(inf, inf, 2.0, 5.0, inf, 
            inf, inf, inf, 10.0, 200.0, 
            inf, inf, inf, inf, 20.0, 
            20.0, inf, 0.06, 0.6, 0.2, 
            inf, inf, inf, inf, inf),
        conversionsAOD = cms.InputTag("allConversions"),
        conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
        ebSplit = cms.double(0.8),
        ebeeSplit = cms.double(1.479),
        mvaName = cms.string('ElectronMVAEstimatorRun2Fall17Iso'),
        mvaTag = cms.string('V1'),
        ptSplit = cms.double(10.0),
        varNames = cms.vstring('ele_oldsigmaietaieta', 
            'ele_oldsigmaiphiiphi', 
            'ele_oldcircularity', 
            'ele_oldr9', 
            'ele_scletawidth', 
            'ele_sclphiwidth', 
            'ele_oldhe', 
            'ele_kfhits', 
            'ele_kfchi2', 
            'ele_gsfchi2', 
            'ele_fbrem', 
            'ele_gsfhits', 
            'ele_expected_inner_hits', 
            'ele_conversionVertexFitProbability', 
            'ele_ep', 
            'ele_eelepout', 
            'ele_IoEmIop', 
            'ele_deltaetain', 
            'ele_deltaphiin', 
            'ele_deltaetaseed', 
            'ele_pfPhotonIso', 
            'ele_pfChargedHadIso', 
            'ele_pfNeutralHadIso', 
            'rho', 
            'ele_psEoverEraw'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_5_2017_puinfo_iso_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_5_2017_puinfo_iso_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_5_2017_puinfo_iso_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_10_2017_puinfo_iso_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_10_2017_puinfo_iso_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_10_2017_puinfo_iso_BDT.weights.xml')
    ))

process.patMultPhiCorrParams_T0pcT1SmearTxy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T0pcT1SmearTxy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T0pcT1T2SmearTxy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T0pcT1T2SmearTxy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T0pcT1T2Txy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T0pcT1T2Txy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T0pcT1Txy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T0pcT1Txy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T0pcTxy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T0pcTxy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T1SmearTxy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T1SmearTxy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T1T2SmearTxy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T1T2SmearTxy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T1T2Txy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T1T2Txy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T1Txy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_T1Txy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_Txy_25ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
    py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
        py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
        py = cms.vdouble(0.00798098092474, -0.000103998219585),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00305719113962, -0.00032676418359),
        py = cms.vdouble(-0.00345131507897, 0.000164816815994),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.000159031461755, 0.00012231873804),
        py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
        py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
        py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
        py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
        py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
        py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
        py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
        py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.patMultPhiCorrParams_Txy_50ns = cms.VPSet(cms.PSet(
    etaMax = cms.double(2.7),
    etaMin = cms.double(0),
    fx = cms.string('(x*[0])+(sq(x)*[1])'),
    fy = cms.string('(x*[0])+(sq(x)*[1])'),
    name = cms.string('hEtaPlus'),
    px = cms.vdouble(-0.00220049396857, 4.86017686051e-07),
    py = cms.vdouble(0.000301784350668, -2.55951949068e-07),
    type = cms.int32(1),
    varType = cms.int32(0)
), 
    cms.PSet(
        etaMax = cms.double(0),
        etaMin = cms.double(-2.7),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaMinus'),
        px = cms.vdouble(-0.000217969078412, 3.0200051094e-07),
        py = cms.vdouble(-0.0014606200538, -2.29508676725e-06),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.392),
        etaMin = cms.double(-1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0Barrel'),
        px = cms.vdouble(-0.0135587323577, 5.55593286464e-05),
        py = cms.vdouble(0.00848783774079, -0.00022596699093),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3),
        etaMin = cms.double(1.392),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapPlus'),
        px = cms.vdouble(-0.00285895832031, -6.08161900014e-05),
        py = cms.vdouble(-0.00934018266651, 0.000259105827172),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.392),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('h0EndcapMinus'),
        px = cms.vdouble(-0.00537876208774, 0.000209817129512),
        py = cms.vdouble(0.011148063877, -4.44149746767e-06),
        type = cms.int32(5),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(1.479),
        etaMin = cms.double(-1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaBarrel'),
        px = cms.vdouble(-0.00192842680623, 2.61152485314e-06),
        py = cms.vdouble(-0.000507607323037, 4.48832037695e-06),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(3.0),
        etaMin = cms.double(1.479),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapPlus'),
        px = cms.vdouble(-0.000519297328533, -2.0682880001e-05),
        py = cms.vdouble(0.00282867507264, 6.66930895313e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-1.479),
        etaMin = cms.double(-3.0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('gammaEndcapMinus'),
        px = cms.vdouble(-0.00103112559755, 1.33699817646e-05),
        py = cms.vdouble(-0.00209888421545, -3.30667819828e-05),
        type = cms.int32(4),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFPlus'),
        px = cms.vdouble(-0.000392672935556, -9.65693700264e-07),
        py = cms.vdouble(0.000114612488388, -3.44552389568e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hHFMinus'),
        px = cms.vdouble(-0.00093227965176, 7.74599924874e-07),
        py = cms.vdouble(-2.95036363418e-05, -7.98830257983e-07),
        type = cms.int32(6),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(5.2),
        etaMin = cms.double(2.901376),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFPlus'),
        px = cms.vdouble(0.00275218993341, -1.69184089548e-05),
        py = cms.vdouble(-0.00113061539637, 6.05994897808e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ), 
    cms.PSet(
        etaMax = cms.double(-2.901376),
        etaMin = cms.double(-5.2),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('egammaHFMinus'),
        px = cms.vdouble(0.00136623849956, -5.55451851761e-06),
        py = cms.vdouble(0.00117549065237, -6.54719096891e-06),
        type = cms.int32(7),
        varType = cms.int32(0)
    ))

process.ak10PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak10PFL1FastL2L3'),
    src = cms.InputTag("ak10PFJets")
)


process.ak10PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak10PFL1FastL2L3Residual'),
    src = cms.InputTag("ak10PFJets")
)


process.ak10PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak10PFL1L2L3'),
    src = cms.InputTag("ak10PFJets")
)


process.ak10PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak10PFL1L2L3Residual'),
    src = cms.InputTag("ak10PFJets")
)


process.ak10PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak10PFL2L3'),
    src = cms.InputTag("ak10PFJets")
)


process.ak10PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak10PFL2L3Residual'),
    src = cms.InputTag("ak10PFJets")
)


process.ak1PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak1PFL1FastL2L3'),
    src = cms.InputTag("ak1PFJets")
)


process.ak1PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak1PFL1FastL2L3Residual'),
    src = cms.InputTag("ak1PFJets")
)


process.ak1PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak1PFL1L2L3'),
    src = cms.InputTag("ak1PFJets")
)


process.ak1PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak1PFL1L2L3Residual'),
    src = cms.InputTag("ak1PFJets")
)


process.ak1PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak1PFL2L3'),
    src = cms.InputTag("ak1PFJets")
)


process.ak1PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak1PFL2L3Residual'),
    src = cms.InputTag("ak1PFJets")
)


process.ak2PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak2PFL1FastL2L3'),
    src = cms.InputTag("ak2PFJets")
)


process.ak2PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak2PFL1FastL2L3Residual'),
    src = cms.InputTag("ak2PFJets")
)


process.ak2PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak2PFL1L2L3'),
    src = cms.InputTag("ak2PFJets")
)


process.ak2PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak2PFL1L2L3Residual'),
    src = cms.InputTag("ak2PFJets")
)


process.ak2PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak2PFL2L3'),
    src = cms.InputTag("ak2PFJets")
)


process.ak2PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak2PFL2L3Residual'),
    src = cms.InputTag("ak2PFJets")
)


process.ak3PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak3PFL1FastL2L3'),
    src = cms.InputTag("ak3PFJets")
)


process.ak3PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak3PFL1FastL2L3Residual'),
    src = cms.InputTag("ak3PFJets")
)


process.ak3PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak3PFL1L2L3'),
    src = cms.InputTag("ak3PFJets")
)


process.ak3PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak3PFL1L2L3Residual'),
    src = cms.InputTag("ak3PFJets")
)


process.ak3PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak3PFL2L3'),
    src = cms.InputTag("ak3PFJets")
)


process.ak3PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak3PFL2L3Residual'),
    src = cms.InputTag("ak3PFJets")
)


process.ak4CaloJetsL1FastL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak4CaloL1FastL2L3'),
    src = cms.InputTag("ak4CaloJets")
)


process.ak4CaloJetsL1FastL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak4CaloL1FastL2L3Residual'),
    src = cms.InputTag("ak4CaloJets")
)


process.ak4CaloJetsL1L2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak4CaloL1L2L3'),
    src = cms.InputTag("ak4CaloJets")
)


process.ak4CaloJetsL1L2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak4CaloL1L2L3Residual'),
    src = cms.InputTag("ak4CaloJets")
)


process.ak4CaloJetsL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak4CaloL2L3'),
    src = cms.InputTag("ak4CaloJets")
)


process.ak4CaloJetsL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak4CaloL2L3Residual'),
    src = cms.InputTag("ak4CaloJets")
)


process.ak4CaloL1FastL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1FastjetCorrector", "ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector")
)


process.ak4CaloL1FastL2L3L6Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1FastjetCorrector", "ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector", "ak4CaloL6SLBCorrector")
)


process.ak4CaloL1FastL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1FastjetCorrector", "ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector", "ak4CaloResidualCorrector")
)


process.ak4CaloL1FastjetCorrector = cms.EDProducer("L1FastjetCorrectorProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.ak4CaloL1L2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1OffsetCorrector", "ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector")
)


process.ak4CaloL1L2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1OffsetCorrector", "ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector", "ak4CaloResidualCorrector")
)


process.ak4CaloL1OffsetCorrector = cms.EDProducer("L1OffsetCorrectorProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.InputTag("offlinePrimaryVertices")
)


process.ak4CaloL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector")
)


process.ak4CaloL2L3L6Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector", "ak4CaloL6SLBCorrector")
)


process.ak4CaloL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL2RelativeCorrector", "ak4CaloL3AbsoluteCorrector", "ak4CaloResidualCorrector")
)


process.ak4CaloL2RelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L2Relative')
)


process.ak4CaloL3AbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L3Absolute')
)


process.ak4CaloL6SLBCorrector = cms.EDProducer("L6SLBCorrectorProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak4CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak4CaloJetsSoftMuonTagInfos")
)


process.ak4CaloResidualCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L2L3Residual')
)


process.ak4JPTJetsL1FastL2L3 = cms.EDProducer("JPTJetCorrectionProducer",
    correctors = cms.vstring('ak4JPTL1FastL2L3'),
    src = cms.InputTag("JetPlusTrackZSPCorJetAntiKt4")
)


process.ak4JPTJetsL1FastL2L3Residual = cms.EDProducer("JPTJetCorrectionProducer",
    correctors = cms.vstring('ak4JPTL1FastL2L3Residual'),
    src = cms.InputTag("JetPlusTrackZSPCorJetAntiKt4")
)


process.ak4JPTJetsL1L2L3 = cms.EDProducer("JPTJetCorrectionProducer",
    correctors = cms.vstring('ak4JPTL1L2L3'),
    src = cms.InputTag("JetPlusTrackZSPCorJetAntiKt4")
)


process.ak4JPTJetsL1L2L3Residual = cms.EDProducer("JPTJetCorrectionProducer",
    correctors = cms.vstring('ak4JPTL1L2L3Residual'),
    src = cms.InputTag("JetPlusTrackZSPCorJetAntiKt4")
)


process.ak4JPTJetsL2L3 = cms.EDProducer("JPTJetCorrectionProducer",
    correctors = cms.vstring('ak4JPTL2L3'),
    src = cms.InputTag("JetPlusTrackZSPCorJetAntiKt4")
)


process.ak4JPTJetsL2L3Residual = cms.EDProducer("JPTJetCorrectionProducer",
    correctors = cms.vstring('ak4JPTL2L3Residual'),
    src = cms.InputTag("JetPlusTrackZSPCorJetAntiKt4")
)


process.ak4JPTL1FastL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1FastjetCorrector", "ak4L1JPTFastjetCorrector", "ak4JPTL2RelativeCorrector", "ak4JPTL3AbsoluteCorrector")
)


process.ak4JPTL1FastL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1FastjetCorrector", "ak4L1JPTFastjetCorrector", "ak4JPTL2RelativeCorrector", "ak4JPTL3AbsoluteCorrector", "ak4JPTResidualCorrector")
)


process.ak4JPTL1L2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1OffsetCorrector", "ak4L1JPTOffsetCorrector", "ak4JPTL2RelativeCorrector", "ak4JPTL3AbsoluteCorrector")
)


process.ak4JPTL1L2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1OffsetCorrector", "ak4L1JPTOffsetCorrector", "ak4JPTL2RelativeCorrector", "ak4JPTL3AbsoluteCorrector", "ak4JPTResidualCorrector")
)


process.ak4JPTL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1OffsetCorrector", "ak4L1JPTOffsetCorrector", "ak4JPTL2RelativeCorrector", "ak4JPTL3AbsoluteCorrector")
)


process.ak4JPTL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4CaloL1OffsetCorrector", "ak4L1JPTOffsetCorrector", "ak4JPTL2RelativeCorrector", "ak4JPTL3AbsoluteCorrector", "ak4JPTResidualCorrector")
)


process.ak4JPTL2RelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L2Relative')
)


process.ak4JPTL3AbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L3Absolute')
)


process.ak4JPTResidualCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L2L3Residual')
)


process.ak4L1JPTFastjetCorrector = cms.EDProducer("L1JPTOffsetCorrectorProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L1JPTOffset'),
    offsetService = cms.InputTag("ak4CaloL1FastjetCorrector")
)


process.ak4L1JPTOffsetCorrector = cms.EDProducer("L1JPTOffsetCorrectorProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L1JPTOffset'),
    offsetService = cms.InputTag("ak4CaloL1OffsetCorrector")
)


process.ak4PFCHSL1FastL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFCHSL1FastjetCorrector", "ak4PFCHSL2RelativeCorrector", "ak4PFCHSL3AbsoluteCorrector")
)


process.ak4PFCHSL1FastL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFCHSL1FastjetCorrector", "ak4PFCHSL2RelativeCorrector", "ak4PFCHSL3AbsoluteCorrector", "ak4PFCHSResidualCorrector")
)


process.ak4PFCHSL1FastjetCorrector = cms.EDProducer("L1FastjetCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak4PFCHSL1L2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFCHSL1OffsetCorrector", "ak4PFCHSL2RelativeCorrector", "ak4PFCHSL3AbsoluteCorrector")
)


process.ak4PFCHSL1L2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFCHSL1OffsetCorrector", "ak4PFCHSL2RelativeCorrector", "ak4PFCHSL3AbsoluteCorrector", "ak4PFCHSResidualCorrector")
)


process.ak4PFCHSL1OffsetCorrector = cms.EDProducer("L1OffsetCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.InputTag("offlinePrimaryVertices")
)


process.ak4PFCHSL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFCHSL2RelativeCorrector", "ak4PFCHSL3AbsoluteCorrector")
)


process.ak4PFCHSL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFCHSL2RelativeCorrector", "ak4PFCHSL3AbsoluteCorrector", "ak4PFCHSResidualCorrector")
)


process.ak4PFCHSL2RelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2Relative')
)


process.ak4PFCHSL3AbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L3Absolute')
)


process.ak4PFCHSResidualCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak4PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak4PFL1FastL2L3'),
    src = cms.InputTag("ak4PFJets")
)


process.ak4PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak4PFL1FastL2L3Residual'),
    src = cms.InputTag("ak4PFJets")
)


process.ak4PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak4PFL1L2L3'),
    src = cms.InputTag("ak4PFJets")
)


process.ak4PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak4PFL1L2L3Residual'),
    src = cms.InputTag("ak4PFJets")
)


process.ak4PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak4PFL2L3'),
    src = cms.InputTag("ak4PFJets")
)


process.ak4PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak4PFL2L3Residual'),
    src = cms.InputTag("ak4PFJets")
)


process.ak4PFL1FastL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL1FastjetCorrector", "ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector")
)


process.ak4PFL1FastL2L3L6Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL1FastjetCorrector", "ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector", "ak4PFL6SLBCorrector")
)


process.ak4PFL1FastL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL1FastjetCorrector", "ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector", "ak4PFResidualCorrector")
)


process.ak4PFL1FastjetCorrector = cms.EDProducer("L1FastjetCorrectorProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak4PFL1L2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL1OffsetCorrector", "ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector")
)


process.ak4PFL1L2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL1OffsetCorrector", "ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector", "ak4PFResidualCorrector")
)


process.ak4PFL1OffsetCorrector = cms.EDProducer("L1OffsetCorrectorProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.InputTag("offlinePrimaryVertices")
)


process.ak4PFL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector")
)


process.ak4PFL2L3L6Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector", "ak4PFL6SLBCorrector")
)


process.ak4PFL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFL2RelativeCorrector", "ak4PFL3AbsoluteCorrector", "ak4PFResidualCorrector")
)


process.ak4PFL2RelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L2Relative')
)


process.ak4PFL3AbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L3Absolute')
)


process.ak4PFL6SLBCorrector = cms.EDProducer("L6SLBCorrectorProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak4PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak4PFJetsSoftMuonTagInfos")
)


process.ak4PFPuppiL1FastL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFPuppiL1FastjetCorrector", "ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector")
)


process.ak4PFPuppiL1FastL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFPuppiL1FastjetCorrector", "ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector", "ak4PFPuppiResidualCorrector")
)


process.ak4PFPuppiL1FastjetCorrector = cms.EDProducer("L1FastjetCorrectorProducer",
    algorithm = cms.string('AK4PFPuppi'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak4PFPuppiL1L2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFPuppiL1OffsetCorrector", "ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector")
)


process.ak4PFPuppiL1L2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFPuppiL1OffsetCorrector", "ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector", "ak4PFPuppiResidualCorrector")
)


process.ak4PFPuppiL1OffsetCorrector = cms.EDProducer("L1OffsetCorrectorProducer",
    algorithm = cms.string('AK4PFPuppi'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.InputTag("offlinePrimaryVertices")
)


process.ak4PFPuppiL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector")
)


process.ak4PFPuppiL2L3ResidualCorrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4PFPuppiL2RelativeCorrector", "ak4PFPuppiL3AbsoluteCorrector", "ak4PFPuppiResidualCorrector")
)


process.ak4PFPuppiL2RelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFPuppi'),
    level = cms.string('L2Relative')
)


process.ak4PFPuppiL3AbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFPuppi'),
    level = cms.string('L3Absolute')
)


process.ak4PFPuppiResidualCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PFPuppi'),
    level = cms.string('L2L3Residual')
)


process.ak4PFResidualCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L2L3Residual')
)


process.ak4TrackJetsL2L3 = cms.EDProducer("TrackJetCorrectionProducer",
    correctors = cms.vstring('ak4TrackL2L3'),
    src = cms.InputTag("ak4TrackJets")
)


process.ak4TrackL2L3Corrector = cms.EDProducer("ChainedJetCorrectorProducer",
    correctors = cms.VInputTag("ak4TrackL2RelativeCorrector", "ak4TrackL3AbsoluteCorrector")
)


process.ak4TrackL2RelativeCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4TRK'),
    level = cms.string('L2Relative')
)


process.ak4TrackL3AbsoluteCorrector = cms.EDProducer("LXXXCorrectorProducer",
    algorithm = cms.string('AK4TRK'),
    level = cms.string('L3Absolute')
)


process.ak5PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak5PFL1FastL2L3'),
    src = cms.InputTag("ak5PFJets")
)


process.ak5PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak5PFL1FastL2L3Residual'),
    src = cms.InputTag("ak5PFJets")
)


process.ak5PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak5PFL1L2L3'),
    src = cms.InputTag("ak5PFJets")
)


process.ak5PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak5PFL1L2L3Residual'),
    src = cms.InputTag("ak5PFJets")
)


process.ak5PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak5PFL2L3'),
    src = cms.InputTag("ak5PFJets")
)


process.ak5PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak5PFL2L3Residual'),
    src = cms.InputTag("ak5PFJets")
)


process.ak6PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak6PFL1FastL2L3'),
    src = cms.InputTag("ak6PFJets")
)


process.ak6PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak6PFL1FastL2L3Residual'),
    src = cms.InputTag("ak6PFJets")
)


process.ak6PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak6PFL1L2L3'),
    src = cms.InputTag("ak6PFJets")
)


process.ak6PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak6PFL1L2L3Residual'),
    src = cms.InputTag("ak6PFJets")
)


process.ak6PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak6PFL2L3'),
    src = cms.InputTag("ak6PFJets")
)


process.ak6PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak6PFL2L3Residual'),
    src = cms.InputTag("ak6PFJets")
)


process.ak7CaloJetsL1FastL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak7CaloL1FastL2L3'),
    src = cms.InputTag("ak7CaloJets")
)


process.ak7CaloJetsL1FastL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak7CaloL1FastL2L3Residual'),
    src = cms.InputTag("ak7CaloJets")
)


process.ak7CaloJetsL1L2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak7CaloL1L2L3'),
    src = cms.InputTag("ak7CaloJets")
)


process.ak7CaloJetsL1L2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak7CaloL1L2L3Residual'),
    src = cms.InputTag("ak7CaloJets")
)


process.ak7CaloJetsL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak7CaloL2L3'),
    src = cms.InputTag("ak7CaloJets")
)


process.ak7CaloJetsL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ak7CaloL2L3Residual'),
    src = cms.InputTag("ak7CaloJets")
)


process.ak7PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak7PFL1FastL2L3'),
    src = cms.InputTag("ak7PFJets")
)


process.ak7PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak7PFL1FastL2L3Residual'),
    src = cms.InputTag("ak7PFJets")
)


process.ak7PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak7PFL1L2L3'),
    src = cms.InputTag("ak7PFJets")
)


process.ak7PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak7PFL1L2L3Residual'),
    src = cms.InputTag("ak7PFJets")
)


process.ak7PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak7PFL2L3'),
    src = cms.InputTag("ak7PFJets")
)


process.ak7PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak7PFL2L3Residual'),
    src = cms.InputTag("ak7PFJets")
)


process.ak8PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak8PFL1FastL2L3'),
    src = cms.InputTag("ak8PFJets")
)


process.ak8PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak8PFL1FastL2L3Residual'),
    src = cms.InputTag("ak8PFJets")
)


process.ak8PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak8PFL1L2L3'),
    src = cms.InputTag("ak8PFJets")
)


process.ak8PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak8PFL1L2L3Residual'),
    src = cms.InputTag("ak8PFJets")
)


process.ak8PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak8PFL2L3'),
    src = cms.InputTag("ak8PFJets")
)


process.ak8PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak8PFL2L3Residual'),
    src = cms.InputTag("ak8PFJets")
)


process.ak9PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak9PFL1FastL2L3'),
    src = cms.InputTag("ak9PFJets")
)


process.ak9PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak9PFL1FastL2L3Residual'),
    src = cms.InputTag("ak9PFJets")
)


process.ak9PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak9PFL1L2L3'),
    src = cms.InputTag("ak9PFJets")
)


process.ak9PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak9PFL1L2L3Residual'),
    src = cms.InputTag("ak9PFJets")
)


process.ak9PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak9PFL2L3'),
    src = cms.InputTag("ak9PFJets")
)


process.ak9PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ak9PFL2L3Residual'),
    src = cms.InputTag("ak9PFJets")
)


process.basicJetsForMet = cms.EDProducer("PATJetCleanerForType1MET",
    jetCorrEtaMax = cms.double(9.9),
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L2L3Residual"),
    offsetCorrLabel = cms.InputTag("L1FastJet"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("patJetsReapplyJEC"),
    type1JetPtThreshold = cms.double(15.0)
)


process.candToGenAssMaker = cms.EDProducer("CandToGenAssMaker",
    ak8JetsInputTag = cms.InputTag("subJetMaker","ak8jetsp4"),
    electronsInputTag = cms.InputTag("electronMaker","elsp4"),
    genJetsInputTag = cms.InputTag("slimmedGenJets"),
    genParticlesInputTagPacked = cms.InputTag("packedGenParticles"),
    genParticlesInputTagPruned = cms.InputTag("prunedGenParticles"),
    jetsInputTag = cms.InputTag("jetMaker","jetsp4"),
    muonsInputTag = cms.InputTag("muonMaker","musp4"),
    ntuplePackedGenParticles = cms.bool(False),
    pfJetsInputTag = cms.InputTag("pfJetMaker","pfjetsp4"),
    photonsInputTag = cms.InputTag("photonMaker","photonsp4"),
    tracksInputTag = cms.InputTag("trackMaker","trkstrkp4"),
    vPIDsToExclude = cms.untracked.vint32(12, 14, 16, 18, 1000022)
)


process.cleanedPatJets = cms.EDProducer("PATJetCleaner",
    checkOverlaps = cms.PSet(
        electrons = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False),
            src = cms.InputTag("slimmedElectrons")
        ),
        muons = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.5),
            pairCut = cms.string(''),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(False),
            src = cms.InputTag("slimmedMuons")
        )
    ),
    finalCut = cms.string(''),
    preselection = cms.string(''),
    src = cms.InputTag("jetSelectorForMet")
)


process.corrPfMetType1 = cms.EDProducer("PFJetMETcorrInputProducer",
    jetCorrEtaMax = cms.double(9.9),
    jetCorrLabel = cms.InputTag("ak4PFCHSL1FastL2L3Corrector"),
    jetCorrLabelRes = cms.InputTag("ak4PFCHSL1FastL2L3ResidualCorrector"),
    offsetCorrLabel = cms.InputTag("ak4PFCHSL1FastjetCorrector"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("ak4PFJetsCHS"),
    type1JetPtThreshold = cms.double(15.0)
)


process.corrPfMetType2 = cms.EDProducer("Type2CorrectionProducer",
    srcUnclEnergySums = cms.VInputTag(cms.InputTag("corrPfMetType1","type2"), cms.InputTag("corrPfMetType1","offset"), cms.InputTag("pfCandMETcorr")),
    type2CorrFormula = cms.string('A'),
    type2CorrParameter = cms.PSet(
        A = cms.double(1.4)
    )
)


process.egmGsfElectronIDs = cms.EDProducer("VersionedGsfElectronIdProducer",
    physicsObjectIDs = cms.VPSet(cms.PSet(
        idDefinition = cms.PSet(
            cutFlow = cms.VPSet(cms.PSet(
                cutName = cms.string('GsfEleMVACut'),
                isIgnored = cms.bool(False),
                mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
                mvaCuts = cms.vdouble(0.287435, 0.221846, -0.303263, 0.967083, 0.929117, 
                    0.726311),
                mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                needsAdditionalProducts = cms.bool(True)
            )),
            idName = cms.string('mvaEleID-Spring15-25ns-nonTrig-V1-wp80'),
            isPOGApproved = cms.untracked.bool(True)
        ),
        idMD5 = cms.string('113c47ceaea0fa687b8bd6d880eb4957'),
        isPOGApproved = cms.untracked.bool(True)
    ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVACut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
                    mvaCuts = cms.vdouble(-0.083313, -0.235222, -0.67099, 0.913286, 0.805013, 
                        0.358969),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Spring15-25ns-nonTrig-V1-wp90'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string('ac4fdc160eefe9eae7338601c02ed4bb'),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVACut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
                    mvaCuts = cms.vdouble(-0.265, -0.556, -0.551, -0.072, -0.286, 
                        -0.267),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Spring15-25ns-nonTrig-V1-wpLoose'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string('99ff36834c4342110d84ea2350a1229c'),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVACut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
                    mvaCuts = cms.vdouble(0.988153, 0.96791, 0.841729),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Spring15-25ns-Trig-V1-wp80'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string('81046ab478185af337be1be9b30948ae'),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVACut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
                    mvaCuts = cms.vdouble(0.972153, 0.922126, 0.610764),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Spring15-25ns-Trig-V1-wp90'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string('bb430b638bf3a4d970627021b0da63ae'),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVACut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),
                    mvaCuts = cms.vdouble(0.940962684155, 0.899208843708, 0.758484721184),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Spring16-GeneralPurpose-V1-wp80'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string('b490bc0b0af2d5f3e9efea562370af2a'),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVACut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),
                    mvaCuts = cms.vdouble(0.836695742607, 0.715337944031, 0.356799721718),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Spring16-GeneralPurpose-V1-wp90'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string('14c153aaf3c207deb3ad4932586647a7'),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVACut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring16HZZV1Categories"),
                    mvaCuts = cms.vdouble(-0.211, -0.396, -0.215, -0.87, -0.838, 
                        -0.763),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring16HZZV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Spring16-HZZ-V1-wpLoose'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string('1797cc03eb62387e10266fca72ea10cd'),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVAExpoScalingCut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17IsoV1Categories"),
                    mvaCuts = cms.vdouble(0.972550955975, 2.97659326151, 0.26538587364, 0.95080381416, 2.66335005587, 
                        0.235582049926, 0.93650371676, 1.57654423239, 3.06701528922, 0.989656208772, 
                        10.342490512, 0.402041564174, 0.981923265653, 9.05548836482, 0.772674931169, 
                        0.962509820174, 8.42589315557, 2.29161526151),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17IsoV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Fall17-iso-V1-wp80'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string(''),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVAExpoScalingCut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17IsoV1Categories"),
                    mvaCuts = cms.vdouble(0.93870703961, 2.65255852282, 0.822264716415, 0.894880292568, 2.76456703588, 
                        0.41233812187, -1830.85836611, -36578.1105538, -1831.20835781, 0.971767483761, 
                        8.9128509851, 1.97124149404, 0.945874502327, 8.83104420393, 2.40849932041, 
                        0.897911201209, 9.81408214417, 4.17158169489),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17IsoV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Fall17-iso-V1-wp90'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string(''),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVACut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17IsoV1Categories"),
                    mvaCuts = cms.vdouble(-0.0956408614642, -0.282299169819, -0.0546668229696, -0.833466688584, -0.767700024757, 
                        -0.691730599565),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17IsoV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Fall17-iso-V1-wpLoose'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string(''),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVAExpoScalingCut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17NoIsoV1Categories"),
                    mvaCuts = cms.vdouble(0.953024095656, 2.7591425841, 0.466964471855, 0.933656476396, 2.70927628427, 
                        0.335122865992, 0.931313368837, 1.58219348007, 3.88894626197, 0.982526856494, 
                        8.70260145586, 1.19748615966, 0.972750945793, 8.17952563102, 1.71117550947, 
                        0.956261953954, 8.10984536628, 3.01392769913),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Fall17-noIso-V1-wp80'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string(''),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVAExpoScalingCut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17NoIsoV1Categories"),
                    mvaCuts = cms.vdouble(0.916511282697, 2.73817035551, 1.03549199648, 0.865573832222, 2.40279446526, 
                        0.797561561328, -3016.03505523, -52140.6185633, -3016.30293872, 0.961654281613, 
                        8.75794383789, 3.13902003216, 0.931925801143, 8.84605743257, 3.59850637933, 
                        0.8899260781, 10.1242341159, 4.35279125072),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Fall17-noIso-V1-wp90'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string(''),
            isPOGApproved = cms.untracked.bool(True)
        ), 
        cms.PSet(
            idDefinition = cms.PSet(
                cutFlow = cms.VPSet(cms.PSet(
                    cutName = cms.string('GsfEleMVACut'),
                    isIgnored = cms.bool(False),
                    mvaCategoriesMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17NoIsoV1Categories"),
                    mvaCuts = cms.vdouble(-0.132858672938, -0.317653009588, -0.0799205914719, -0.856871961305, -0.810764214158, 
                        -0.717926593302),
                    mvaValueMapName = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
                    needsAdditionalProducts = cms.bool(True)
                )),
                idName = cms.string('mvaEleID-Fall17-noIso-V1-wpLoose'),
                isPOGApproved = cms.untracked.bool(True)
            ),
            idMD5 = cms.string(''),
            isPOGApproved = cms.untracked.bool(True)
        )),
    physicsObjectSrc = cms.InputTag("slimmedElectrons","","PAT")
)


process.electronMVAValueMapProducer = cms.EDProducer("ElectronMVAValueMapProducer",
    mvaConfigurations = cms.VPSet(cms.PSet(
        mvaName = cms.string('ElectronMVAEstimatorRun2Phys14NonTrig'),
        mvaTag = cms.string('25nsV1'),
        weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB1_5_oldscenario2phys14_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB2_5_oldscenario2phys14_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EE_5_oldscenario2phys14_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB1_10_oldscenario2phys14_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EB2_10_oldscenario2phys14_BDT.weights.xml', 
            'RecoEgamma/ElectronIdentification/data/PHYS14/EIDmva_EE_10_oldscenario2phys14_BDT.weights.xml')
    ), 
        cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            conversionsAOD = cms.InputTag("allConversions"),
            conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
            mvaName = cms.string('ElectronMVAEstimatorRun2Spring15NonTrig'),
            mvaTag = cms.string('25nsV1'),
            weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_5_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldNonTrigSpring15_ConvVarCwoBoolean_TMVA412_FullStatLowPt_PairNegWeightsGlobal_BDT.weights.xml')
        ), 
        cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            conversionsAOD = cms.InputTag("allConversions"),
            conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
            mvaName = cms.string('ElectronMVAEstimatorRun2Spring15Trig'),
            mvaTag = cms.string('50nsV1'),
            weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldTrigSpring15_50ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml')
        ), 
        cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            conversionsAOD = cms.InputTag("allConversions"),
            conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
            mvaName = cms.string('ElectronMVAEstimatorRun2Spring15Trig'),
            mvaTag = cms.string('25nsV1'),
            weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB1_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EB2_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring15/EIDmva_EE_10_oldTrigSpring15_25ns_data_1_VarD_TMVA412_Sig6BkgAll_MG_noSpec_BDT.weights.xml')
        ), 
        cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            conversionsAOD = cms.InputTag("allConversions"),
            conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
            mvaName = cms.string('ElectronMVAEstimatorRun2Spring16HZZ'),
            mvaTag = cms.string('V1'),
            weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_5.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_5.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_5.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB1_10.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EB2_10.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_HZZ_V1/electronID_mva_Spring16_HZZ_V1_EE_10.weights.xml')
        ), 
        cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            conversionsAOD = cms.InputTag("allConversions"),
            conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
            mvaName = cms.string('ElectronMVAEstimatorRun2Spring16GeneralPurpose'),
            mvaTag = cms.string('V1'),
            weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB1_10.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EB2_10.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Spring16_GeneralPurpose_V1/electronID_mva_Spring16_GeneralPurpose_V1_EE_10.weights.xml')
        ), 
        cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            clipLower = cms.vdouble(-inf, -inf, -1.0, -inf, -inf, 
                -inf, -inf, -inf, -inf, -inf, 
                -1.0, -inf, -inf, -inf, -inf, 
                -inf, -inf, -0.06, -0.6, -0.2, 
                -inf, -inf),
            clipUpper = cms.vdouble(inf, inf, 2.0, 5.0, inf, 
                inf, inf, inf, 10.0, 200.0, 
                inf, inf, inf, inf, 20.0, 
                20.0, inf, 0.06, 0.6, 0.2, 
                inf, inf),
            conversionsAOD = cms.InputTag("allConversions"),
            conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
            ebSplit = cms.double(0.8),
            ebeeSplit = cms.double(1.479),
            mvaName = cms.string('ElectronMVAEstimatorRun2Fall17NoIso'),
            mvaTag = cms.string('V1'),
            ptSplit = cms.double(10.0),
            varNames = cms.vstring('ele_oldsigmaietaieta', 
                'ele_oldsigmaiphiiphi', 
                'ele_oldcircularity', 
                'ele_oldr9', 
                'ele_scletawidth', 
                'ele_sclphiwidth', 
                'ele_oldhe', 
                'ele_kfhits', 
                'ele_kfchi2', 
                'ele_gsfchi2', 
                'ele_fbrem', 
                'ele_gsfhits', 
                'ele_expected_inner_hits', 
                'ele_conversionVertexFitProbability', 
                'ele_ep', 
                'ele_eelepout', 
                'ele_IoEmIop', 
                'ele_deltaetain', 
                'ele_deltaphiin', 
                'ele_deltaetaseed', 
                'rho', 
                'ele_psEoverEraw'),
            weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_5_2017_puinfo_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_5_2017_puinfo_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_5_2017_puinfo_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_10_2017_puinfo_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_10_2017_puinfo_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_10_2017_puinfo_BDT.weights.xml')
        ), 
        cms.PSet(
            beamSpot = cms.InputTag("offlineBeamSpot"),
            clipLower = cms.vdouble(-inf, -inf, -1.0, -inf, -inf, 
                -inf, -inf, -inf, -inf, -inf, 
                -1.0, -inf, -inf, -inf, -inf, 
                -inf, -inf, -0.06, -0.6, -0.2, 
                -inf, -inf, -inf, -inf, -inf),
            clipUpper = cms.vdouble(inf, inf, 2.0, 5.0, inf, 
                inf, inf, inf, 10.0, 200.0, 
                inf, inf, inf, inf, 20.0, 
                20.0, inf, 0.06, 0.6, 0.2, 
                inf, inf, inf, inf, inf),
            conversionsAOD = cms.InputTag("allConversions"),
            conversionsMiniAOD = cms.InputTag("reducedEgamma","reducedConversions"),
            ebSplit = cms.double(0.8),
            ebeeSplit = cms.double(1.479),
            mvaName = cms.string('ElectronMVAEstimatorRun2Fall17Iso'),
            mvaTag = cms.string('V1'),
            ptSplit = cms.double(10.0),
            varNames = cms.vstring('ele_oldsigmaietaieta', 
                'ele_oldsigmaiphiiphi', 
                'ele_oldcircularity', 
                'ele_oldr9', 
                'ele_scletawidth', 
                'ele_sclphiwidth', 
                'ele_oldhe', 
                'ele_kfhits', 
                'ele_kfchi2', 
                'ele_gsfchi2', 
                'ele_fbrem', 
                'ele_gsfhits', 
                'ele_expected_inner_hits', 
                'ele_conversionVertexFitProbability', 
                'ele_ep', 
                'ele_eelepout', 
                'ele_IoEmIop', 
                'ele_deltaetain', 
                'ele_deltaphiin', 
                'ele_deltaetaseed', 
                'ele_pfPhotonIso', 
                'ele_pfChargedHadIso', 
                'ele_pfNeutralHadIso', 
                'rho', 
                'ele_psEoverEraw'),
            weightFileNames = cms.vstring('RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_5_2017_puinfo_iso_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_5_2017_puinfo_iso_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_5_2017_puinfo_iso_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB1_10_2017_puinfo_iso_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EB2_10_2017_puinfo_iso_BDT.weights.xml', 
                'RecoEgamma/ElectronIdentification/data/Fall17/EIDmva_EE_10_2017_puinfo_iso_BDT.weights.xml')
        )),
    src = cms.InputTag("gedGsfElectrons"),
    srcMiniAOD = cms.InputTag("slimmedElectrons","","PAT")
)


process.electronMaker = cms.EDProducer("ElectronMaker",
    aliasPrefix = cms.untracked.string('els'),
    bFieldInputTag = cms.InputTag("eventMaker","evtbField"),
    beamSpotInputTag = cms.InputTag("beamSpotMaker","evtbsp4"),
    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    cms2scsseeddetidInputTag = cms.InputTag("scMaker"),
    ebReducedRecHitCollectionTag = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    eeReducedRecHitCollectionTag = cms.InputTag("reducedEgamma","reducedEERecHits"),
    eidLHTag = cms.InputTag("egammaIDLikelihood"),
    electronHEEPIdMap = cms.InputTag("egmGsfElectronIDs","heepElectronID-HEEPV60"),
    electronLooseIdMap = cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
    electronMediumIdMap = cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    electronTightIdMap = cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
    electronVIDFall17IsoMvaCatMap = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17IsoV1Categories"),
    electronVIDFall17IsoMvaValueMap = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17IsoV1Values"),
    electronVIDFall17NoIsoMvaCatMap = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17NoIsoV1Categories"),
    electronVIDFall17NoIsoMvaValueMap = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Fall17NoIsoV1Values"),
    electronVIDNonTrigMvaCatMap = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Categories"),
    electronVIDNonTrigMvaValueMap = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values"),
    electronVIDNonTrigMvaWP80IdMap = cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring15-25ns-nonTrig-V1-wp80"),
    electronVIDNonTrigMvaWP90IdMap = cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring15-25ns-nonTrig-V1-wp90"),
    electronVIDSpring16GPMvaCatMap = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Categories"),
    electronVIDSpring16GPMvaValueMap = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring16GeneralPurposeV1Values"),
    electronVIDSpring16HZZMvaCatMap = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring16HZZV1Categories"),
    electronVIDSpring16HZZMvaValueMap = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring16HZZV1Values"),
    electronVIDTrigMvaCatMap = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories"),
    electronVIDTrigMvaValueMap = cms.InputTag("electronMVAValueMapProducer","ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
    electronVIDTrigMvaWP80IdMap = cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring15-25ns-Trig-V1-wp80"),
    electronVIDTrigMvaWP90IdMap = cms.InputTag("egmGsfElectronIDs","mvaEleID-Spring15-25ns-Trig-V1-wp90"),
    electronVetoIdMap = cms.InputTag("egmGsfElectronIDs","cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
    electronsInputTag = cms.InputTag("slimmedElectrons"),
    esReducedRecHitCollectionTag = cms.InputTag("reducedEgamma","reducedESRecHits"),
    gsftracksInputTag = cms.InputTag("electronGsfTracks"),
    minAbsDcot = cms.double(0.02),
    minAbsDist = cms.double(0.02),
    minSharedFractionOfHits = cms.double(0.45),
    miniIsoAllValueMap = cms.InputTag("isoForEle","miniIsoAll"),
    miniIsoChgValueMap = cms.InputTag("isoForEle","miniIsoChg"),
    pfCandsInputTag = cms.InputTag("packedPFCandidates"),
    pfJetsInputTag = cms.InputTag("slimmedJets","","PAT"),
    recoConversionInputTag = cms.InputTag("reducedEgamma","reducedConversions"),
    rhoInputTag = cms.InputTag("fastJetMaker","evtrho"),
    trksInputTag = cms.InputTag("generalTracks"),
    vtxInputTag = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.eventMaker = cms.EDProducer("EventMaker",
    CMS3tag = cms.string('V04'),
    aliasPrefix = cms.untracked.string('evt'),
    datasetName = cms.string('/DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'),
    dcsTag = cms.InputTag("scalersRawToDigi"),
    isData = cms.bool(False),
    scalersTag = cms.InputTag("scalersRawToDigi")
)


process.fixedGridRhoAllMaker = cms.EDProducer("EnergyDensityMaker",
    alias = cms.untracked.string('evt_fixgrid_all_rho'),
    input = cms.InputTag("fixedGridRhoAll","","RECO")
)


process.fixedGridRhoFastJetAllCaloMaker = cms.EDProducer("EnergyDensityMaker",
    alias = cms.untracked.string('evt_fixgridfastjet_allcalo_rho'),
    input = cms.InputTag("fixedGridRhoFastjetAllCalo","","RECO")
)


process.fixedGridRhoFastJetAllCentralMaker = cms.EDProducer("EnergyDensityMaker",
    alias = cms.untracked.string('evt_fixgridfastjet_central_rho'),
    input = cms.InputTag("fixedGridRhoFastjetCentral","","RECO")
)


process.fixedGridRhoFastJetAllMaker = cms.EDProducer("EnergyDensityMaker",
    alias = cms.untracked.string('evt_fixgridfastjet_all_rho'),
    input = cms.InputTag("fixedGridRhoFastjetAll","","RECO")
)


process.fixedGridRhoFastJetCentralCaloMaker = cms.EDProducer("EnergyDensityMaker",
    alias = cms.untracked.string('evt_fixgridfastjet_centralcalo_rho'),
    input = cms.InputTag("fixedGridRhoFastjetCentralCalo","","RECO")
)


process.fixedGridRhoFastJetCentralChargedPileUpMaker = cms.EDProducer("EnergyDensityMaker",
    alias = cms.untracked.string('evt_fixgridfastjet_centralchargedpileup_rho'),
    input = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp","","RECO")
)


process.fixedGridRhoFastJetCentralNeutralMaker = cms.EDProducer("EnergyDensityMaker",
    alias = cms.untracked.string('evt_fixgridfastjet_centralneutral_rho'),
    input = cms.InputTag("fixedGridRhoFastjetCentralNeutral","","RECO")
)


process.genJetMaker = cms.EDProducer("GenJetMaker",
    aliasPostfix = cms.untracked.string('NoMuNoNu'),
    genJetMinPtCut = cms.double(10.0),
    genJetsInputTag = cms.InputTag("slimmedGenJets")
)


process.genMaker = cms.EDProducer("GenMaker",
    LHEInputTag = cms.InputTag("externalLHEProducer"),
    aliasPrefix = cms.untracked.string('genps'),
    exclusiveCrossSection = cms.untracked.double(0.0),
    genEvtInfoInputTag = cms.InputTag("generator"),
    genParticlesInputTag = cms.InputTag("prunedGenParticles"),
    inclusiveCrossSection = cms.untracked.double(0.0),
    kfactor = cms.untracked.double(1.0),
    ntupleDaughters = cms.bool(True),
    ntupleOnlyStatus3 = cms.bool(False),
    ntuplePackedGenParticles = cms.bool(False),
    packedGenParticlesInputTag = cms.InputTag("packedGenParticles"),
    vmetPIDs = cms.untracked.vint32(1000022)
)


process.genMetExtractor = cms.EDProducer("GenMETExtractor",
    metSource = cms.InputTag("slimmedMETs","","@skipCurrentProcess")
)


process.hltMaker = cms.EDProducer("HLTMaker",
    aliasPrefix = cms.untracked.string('hlt'),
    fillTriggerObjects = cms.untracked.bool(False),
    processName = cms.untracked.string('HLT'),
    prunedTriggerNames = cms.untracked.vstring('HLT*Mu*', 
        'HLT*Ele*', 
        'HLT*EG*', 
        'HLT*Photon*', 
        'HLT*Jet*', 
        'HLT*MET*'),
    triggerObjectsName = cms.untracked.string('selectedPatTrigger'),
    triggerPrescaleInputTag = cms.untracked.string('patTrigger')
)


process.hypDilepMaker = cms.EDProducer("HypDilepMaker",
    LooseLepton_PtCut = cms.double(10.0),
    TightLepton_PtCut = cms.double(10.0),
    aliasPrefix = cms.untracked.string('hyp'),
    electronsInputTag = cms.InputTag("electronMaker"),
    elsChargeInputTag = cms.InputTag("electronMaker","elscharge"),
    elsTypeInputTag = cms.InputTag("electronMaker","elstype"),
    elsp4InputTag = cms.InputTag("electronMaker","elsp4"),
    hypJetMaxEtaCut = cms.double(5.0),
    hypJetMinPtCut = cms.double(30.0),
    jetsInputTag = cms.InputTag("pfJetMaker","pfjetsp4"),
    metInputTag = cms.InputTag("pfmetMaker"),
    muToGenInputTag = cms.InputTag("muToGenAssMaker"),
    musChargeInputTag = cms.InputTag("muonMaker","muscharge"),
    musTypeInputTag = cms.InputTag("muonMaker","mustype"),
    musp4InputTag = cms.InputTag("muonMaker","musp4"),
    useSTAMuons = cms.bool(False)
)


process.ic5CaloJetsL1FastL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ic5CaloL1FastL2L3'),
    src = cms.InputTag("iterativeCone5CaloJets")
)


process.ic5CaloJetsL1FastL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ic5CaloL1FastL2L3Residual'),
    src = cms.InputTag("iterativeCone5CaloJets")
)


process.ic5CaloJetsL1L2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ic5CaloL1L2L3'),
    src = cms.InputTag("iterativeCone5CaloJets")
)


process.ic5CaloJetsL1L2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ic5CaloL1L2L3Residual'),
    src = cms.InputTag("iterativeCone5CaloJets")
)


process.ic5CaloJetsL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ic5CaloL2L3'),
    src = cms.InputTag("iterativeCone5CaloJets")
)


process.ic5CaloJetsL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('ic5CaloL2L3Residual'),
    src = cms.InputTag("iterativeCone5CaloJets")
)


process.ic5PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ic5PFL1FastL2L3'),
    src = cms.InputTag("iterativeCone5PFJets")
)


process.ic5PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ic5PFL1FastL2L3Residual'),
    src = cms.InputTag("iterativeCone5PFJets")
)


process.ic5PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ic5PFL1L2L3'),
    src = cms.InputTag("iterativeCone5PFJets")
)


process.ic5PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ic5PFL1L2L3Residual'),
    src = cms.InputTag("iterativeCone5PFJets")
)


process.ic5PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ic5PFL2L3'),
    src = cms.InputTag("iterativeCone5PFJets")
)


process.ic5PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('ic5PFL2L3Residual'),
    src = cms.InputTag("iterativeCone5PFJets")
)


process.isoForEle = cms.EDProducer("EleIsoValueMapProducer",
    EAFile_MiniIso = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Spring15/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_25ns.txt'),
    EAFile_PFIso = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Summer16/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_80X.txt'),
    relative = cms.bool(False),
    rho_MiniIso = cms.InputTag("fixedGridRhoFastjetAll"),
    rho_PFIso = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("slimmedElectrons")
)


process.isoForMu = cms.EDProducer("MuonIsoValueMapProducer",
    EAFile_MiniIso = cms.FileInPath('PhysicsTools/NanoAOD/data/effAreaMuons_cone03_pfNeuHadronsAndPhotons_80X.txt'),
    relative = cms.bool(False),
    rho_MiniIso = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("slimmedMuons")
)


process.isoTrackMaker = cms.EDProducer("IsoTrackMaker",
    isoTracksTag = cms.InputTag("isolatedTracks","","PAT"),
    lostTracksTag = cms.InputTag("lostTracks","","PAT"),
    pT_cut = cms.double(5.0),
    pT_cut_noIso = cms.double(20.0),
    pfCandidatesTag = cms.InputTag("packedPFCandidates","","PAT")
)


process.kt4CaloJetsL1FastL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt4CaloL1FastL2L3'),
    src = cms.InputTag("kt4CaloJets")
)


process.kt4CaloJetsL1FastL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt4CaloL1FastL2L3Residual'),
    src = cms.InputTag("kt4CaloJets")
)


process.kt4CaloJetsL1L2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt4CaloL1L2L3'),
    src = cms.InputTag("kt4CaloJets")
)


process.kt4CaloJetsL1L2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt4CaloL1L2L3Residual'),
    src = cms.InputTag("kt4CaloJets")
)


process.kt4CaloJetsL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt4CaloL2L3'),
    src = cms.InputTag("kt4CaloJets")
)


process.kt4CaloJetsL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt4CaloL2L3Residual'),
    src = cms.InputTag("kt4CaloJets")
)


process.kt4PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt4PFL1FastL2L3'),
    src = cms.InputTag("kt4PFJets")
)


process.kt4PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt4PFL1FastL2L3Residual'),
    src = cms.InputTag("kt4PFJets")
)


process.kt4PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt4PFL1L2L3'),
    src = cms.InputTag("kt4PFJets")
)


process.kt4PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt4PFL1L2L3Residual'),
    src = cms.InputTag("kt4PFJets")
)


process.kt4PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt4PFL2L3'),
    src = cms.InputTag("kt4PFJets")
)


process.kt4PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt4PFL2L3Residual'),
    src = cms.InputTag("kt4PFJets")
)


process.kt6CaloJetsL1FastL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt6CaloL1FastL2L3'),
    src = cms.InputTag("kt6CaloJets")
)


process.kt6CaloJetsL1FastL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt6CaloL1FastL2L3Residual'),
    src = cms.InputTag("kt6CaloJets")
)


process.kt6CaloJetsL1L2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt6CaloL1L2L3'),
    src = cms.InputTag("kt6CaloJets")
)


process.kt6CaloJetsL1L2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt6CaloL1L2L3Residual'),
    src = cms.InputTag("kt6CaloJets")
)


process.kt6CaloJetsL2L3 = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt6CaloL2L3'),
    src = cms.InputTag("kt6CaloJets")
)


process.kt6CaloJetsL2L3Residual = cms.EDProducer("CaloJetCorrectionProducer",
    correctors = cms.vstring('kt6CaloL2L3Residual'),
    src = cms.InputTag("kt6CaloJets")
)


process.kt6PFJetsL1FastL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt6PFL1FastL2L3'),
    src = cms.InputTag("kt6PFJets")
)


process.kt6PFJetsL1FastL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt6PFL1FastL2L3Residual'),
    src = cms.InputTag("kt6PFJets")
)


process.kt6PFJetsL1L2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt6PFL1L2L3'),
    src = cms.InputTag("kt6PFJets")
)


process.kt6PFJetsL1L2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt6PFL1L2L3Residual'),
    src = cms.InputTag("kt6PFJets")
)


process.kt6PFJetsL2L3 = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt6PFL2L3'),
    src = cms.InputTag("kt6PFJets")
)


process.kt6PFJetsL2L3Residual = cms.EDProducer("PFJetCorrectionProducer",
    correctors = cms.vstring('kt6PFL2L3Residual'),
    src = cms.InputTag("kt6PFJets")
)


process.metFilterMaker = cms.EDProducer("MetFilterMaker",
    aliasPrefix = cms.untracked.string('filt'),
    filtersInputTag = cms.InputTag("TriggerResults"),
    processName = cms.untracked.string('PAT')
)


process.metrawCalo = cms.EDProducer("RecoMETExtractor",
    correctionLevel = cms.string('rawCalo'),
    metSource = cms.InputTag("slimmedMETs","","@skipCurrentProcess")
)


process.muonMaker = cms.EDProducer("MuonMaker",
    aliasPrefix = cms.untracked.string('mus'),
    cosmicCompat = cms.InputTag("muons","cosmicsVeto"),
    miniIsoAllValueMap = cms.InputTag("isoForMu","miniIsoAll"),
    miniIsoChgValueMap = cms.InputTag("isoForMu","miniIsoChg"),
    muonsInputTag = cms.InputTag("slimmedMuons"),
    pfCandsInputTag = cms.InputTag("packedPFCandidates"),
    pfJetsInputTag = cms.InputTag("slimmedJets","","PAT"),
    pfNoPileUpInputTag_ = cms.InputTag("pfNoPileUp"),
    tevMuonsName = cms.string('tevMuons'),
    vtxInputTag = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.particleFlowDisplacedVertex = cms.EDProducer("PFDisplacedVertexProducer",
    avfParameters = cms.PSet(
        Tini = cms.double(256.0),
        ratio = cms.double(0.25),
        sigmacut = cms.double(6.0)
    ),
    debug = cms.untracked.bool(False),
    longSize = cms.double(5),
    mainVertexLabel = cms.InputTag("offlinePrimaryVertices"),
    minAdaptWeight = cms.double(0.5),
    offlineBeamSpotLabel = cms.InputTag("offlineBeamSpot"),
    primaryVertexCut = cms.double(1.8),
    switchOff2TrackVertex = cms.untracked.bool(True),
    tecCut = cms.double(220),
    tobCut = cms.double(100),
    tracksSelectorParameters = cms.PSet(
        bSelectTracks = cms.bool(True),
        dxy_min = cms.double(0.2),
        nChi2_max = cms.double(5.0),
        nChi2_min = cms.double(0.5),
        nHits_min = cms.int32(6),
        nOuterHits_max = cms.int32(9),
        pt_min = cms.double(0.2),
        quality = cms.string('HighPurity')
    ),
    transvSize = cms.double(1.0),
    verbose = cms.untracked.bool(False),
    vertexCandidatesLabel = cms.InputTag("particleFlowDisplacedVertexCandidate"),
    vertexIdentifierParameters = cms.PSet(
        angles = cms.vdouble(15, 15),
        bIdentifyVertices = cms.bool(True),
        logPrimSec_min = cms.double(0.0),
        looper_eta_max = cms.double(0.1),
        masses = cms.vdouble(0.05, 0.485, 0.515, 0.48, 0.52, 
            1.107, 1.125, 0.2),
        pt_kink_min = cms.double(3.0),
        pt_min = cms.double(0.5)
    )
)


process.particleFlowPtrs = cms.EDProducer("PFCandidateFwdPtrProducer",
    src = cms.InputTag("particleFlow")
)


process.patCHSMet = cms.EDProducer("PATMETProducer",
    addEfficiencies = cms.bool(False),
    addGenMET = cms.bool(False),
    addMuonCorrections = cms.bool(False),
    addResolutions = cms.bool(False),
    computeMETSignificance = cms.bool(False),
    computeMETSignificant = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    genMETSource = cms.InputTag("genMetTrue"),
    metSource = cms.InputTag("pfMetCHS"),
    muonSource = cms.InputTag("muons"),
    parameters = cms.PSet(
        dRMatch = cms.double(0.4),
        jetThreshold = cms.double(15),
        jeta = cms.vdouble(0.8, 1.3, 1.9, 2.5),
        jpar = cms.vdouble(1.29, 1.19, 1.07, 1.13, 1.12),
        pjpar = cms.vdouble(-0.04, 0.6504)
    ),
    resolutions = cms.PSet(

    ),
    srcJetResPhi = cms.string('AK4PFchs_phi'),
    srcJetResPt = cms.string('AK4PFchs_pt'),
    srcJetSF = cms.string('AK4PFchs'),
    srcJets = cms.InputTag("cleanedPatJets"),
    srcLeptons = cms.VInputTag("selectedPatElectrons", "selectedPatMuons", "selectedPatPhotons"),
    srcPFCands = cms.InputTag("particleFlow"),
    srcRho = cms.InputTag("fixedGridRhoAll"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.patCaloMet = cms.EDProducer("PATMETProducer",
    addEfficiencies = cms.bool(False),
    addGenMET = cms.bool(False),
    addMuonCorrections = cms.bool(False),
    addResolutions = cms.bool(False),
    computeMETSignificance = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    genMETSource = cms.InputTag("genMetTrue"),
    metSource = cms.InputTag("metrawCalo"),
    muonSource = cms.InputTag("muons"),
    parameters = cms.PSet(
        dRMatch = cms.double(0.4),
        jetThreshold = cms.double(15),
        jeta = cms.vdouble(0.8, 1.3, 1.9, 2.5),
        jpar = cms.vdouble(1.29, 1.19, 1.07, 1.13, 1.12),
        pjpar = cms.vdouble(-0.04, 0.6504)
    ),
    resolutions = cms.PSet(

    ),
    srcJetResPhi = cms.string('AK4PFchs_phi'),
    srcJetResPt = cms.string('AK4PFchs_pt'),
    srcJetSF = cms.string('AK4PFchs'),
    srcJets = cms.InputTag("cleanedPatJets"),
    srcLeptons = cms.VInputTag("selectedPatElectrons", "selectedPatMuons", "selectedPatPhotons"),
    srcPFCands = cms.InputTag("particleFlow"),
    srcRho = cms.InputTag("fixedGridRhoAll"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.patJetCorrFactorsReapplyJEC = cms.EDProducer("JetCorrFactorsProducer",
    emf = cms.bool(False),
    extraJPTOffset = cms.string('L1FastJet'),
    flavorType = cms.string('J'),
    levels = cms.vstring('L1FastJet', 
        'L2Relative', 
        'L3Absolute'),
    payload = cms.string('AK4PFchs'),
    primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("slimmedJets"),
    useNPV = cms.bool(True),
    useRho = cms.bool(True)
)


process.patJetsReapplyJEC = cms.EDProducer("PATJetUpdater",
    addBTagInfo = cms.bool(True),
    addDiscriminators = cms.bool(True),
    addJetCorrFactors = cms.bool(True),
    addTagInfos = cms.bool(False),
    discriminatorSources = cms.VInputTag(),
    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC")),
    jetSource = cms.InputTag("slimmedJets"),
    printWarning = cms.bool(True),
    tagInfoSources = cms.VInputTag(),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.patMETs = cms.EDProducer("PATMETProducer",
    addEfficiencies = cms.bool(False),
    addGenMET = cms.bool(True),
    addMuonCorrections = cms.bool(False),
    addResolutions = cms.bool(False),
    computeMETSignificance = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    genMETSource = cms.InputTag("genMetTrue"),
    metSource = cms.InputTag("pfMetT1"),
    muonSource = cms.InputTag("muons"),
    parameters = cms.PSet(
        dRMatch = cms.double(0.4),
        jetThreshold = cms.double(15),
        jeta = cms.vdouble(0.8, 1.3, 1.9, 2.5),
        jpar = cms.vdouble(1.29, 1.19, 1.07, 1.13, 1.12),
        pjpar = cms.vdouble(-0.04, 0.6504)
    ),
    resolutions = cms.PSet(

    ),
    srcJetResPhi = cms.string('AK4PFchs_phi'),
    srcJetResPt = cms.string('AK4PFchs_pt'),
    srcJetSF = cms.string('AK4PFchs'),
    srcJets = cms.InputTag("cleanedPatJets"),
    srcLeptons = cms.VInputTag("selectedPatElectrons", "selectedPatMuons", "selectedPatPhotons"),
    srcPFCands = cms.InputTag("particleFlow"),
    srcRho = cms.InputTag("fixedGridRhoAll"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.patPFMet = cms.EDProducer("PATMETProducer",
    addEfficiencies = cms.bool(False),
    addGenMET = cms.bool(True),
    addMuonCorrections = cms.bool(False),
    addResolutions = cms.bool(False),
    computeMETSignificance = cms.bool(True),
    efficiencies = cms.PSet(

    ),
    genMETSource = cms.InputTag("genMetExtractor"),
    metSource = cms.InputTag("pfMet"),
    muonSource = cms.InputTag("muons"),
    parameters = cms.PSet(
        dRMatch = cms.double(0.4),
        jetThreshold = cms.double(15),
        jeta = cms.vdouble(0.8, 1.3, 1.9, 2.5),
        jpar = cms.vdouble(1.29, 1.19, 1.07, 1.13, 1.12),
        pjpar = cms.vdouble(-0.04, 0.6504)
    ),
    resolutions = cms.PSet(

    ),
    srcJetResPhi = cms.string('AK4PFchs_phi'),
    srcJetResPt = cms.string('AK4PFchs_pt'),
    srcJetSF = cms.string('AK4PFchs'),
    srcJets = cms.InputTag("cleanedPatJets"),
    srcLeptons = cms.VInputTag(cms.InputTag("slimmedElectrons"), cms.InputTag("slimmedMuons"), cms.InputTag("slimmedPhotons")),
    srcPFCands = cms.InputTag("packedPFCandidates"),
    srcRho = cms.InputTag("fixedGridRhoAll"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.patPFMetT0Corr = cms.EDProducer("Type0PFMETcorrInputProducer",
    correction = cms.PSet(
        formula = cms.string('(x<35)?(-( [0]+x*[1]+pow(x, 2)*[2]+pow(x, 3)*[3] )):(-( [0]+35*[1]+pow(35, 2)*[2]+pow(35, 3)*[3] ))'),
        par0 = cms.double(-0.181414),
        par1 = cms.double(-0.476934),
        par2 = cms.double(0.00863564),
        par3 = cms.double(-4.94181e-05)
    ),
    minDz = cms.double(0.2),
    srcHardScatterVertex = cms.InputTag("selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0"),
    srcPFCandidateToVertexAssociations = cms.InputTag("pfCandidateToVertexAssociation")
)


process.patPFMetT0pcT1 = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetT0Corr"))
)


process.patPFMetT0pcT1Smear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetT0Corr"))
)


process.patPFMetT0pcT1T2 = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetT2Corr","type2"), cms.InputTag("patPFMetT0Corr"))
)


process.patPFMetT0pcT1T2Smear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetT2SmearCorr","type2"), cms.InputTag("patPFMetT0Corr"))
)


process.patPFMetT0pcT1T2Txy = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetT2Corr","type2"), cms.InputTag("patPFMetT0Corr"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT0pcT1T2TxySmear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetT2SmearCorr","type2"), cms.InputTag("patPFMetT0Corr"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT0pcT1Txy = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetT0Corr"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT0pcT1TxySmear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetT0Corr"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT1 = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"))
)


process.patPFMetT1ElectronEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrElectronEnDown"))
)


process.patPFMetT1ElectronEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrElectronEnUp"))
)


process.patPFMetT1JetEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetEnDown"))
)


process.patPFMetT1JetEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetEnUp"))
)


process.patPFMetT1JetResDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetResDown"))
)


process.patPFMetT1JetResUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetResUp"))
)


process.patPFMetT1MuonEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrMuonEnDown"))
)


process.patPFMetT1MuonEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrMuonEnUp"))
)


process.patPFMetT1PhotonEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrPhotonEnDown"))
)


process.patPFMetT1PhotonEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrPhotonEnUp"))
)


process.patPFMetT1Smear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"))
)


process.patPFMetT1SmearElectronEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrElectronEnDown"))
)


process.patPFMetT1SmearElectronEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrElectronEnUp"))
)


process.patPFMetT1SmearJetEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetEnDown"))
)


process.patPFMetT1SmearJetEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrJetEnUp"))
)


process.patPFMetT1SmearJetResDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrSmearedJetResDown"))
)


process.patPFMetT1SmearJetResUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrSmearedJetResUp"))
)


process.patPFMetT1SmearMuonEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrMuonEnDown"))
)


process.patPFMetT1SmearMuonEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrMuonEnUp"))
)


process.patPFMetT1SmearPhotonEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrPhotonEnDown"))
)


process.patPFMetT1SmearPhotonEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrPhotonEnUp"))
)


process.patPFMetT1SmearTauEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrTauEnDown"))
)


process.patPFMetT1SmearTauEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrTauEnUp"))
)


process.patPFMetT1SmearUnclusteredEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrUnclusteredEnDown"))
)


process.patPFMetT1SmearUnclusteredEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1Smear"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrUnclusteredEnUp"))
)


process.patPFMetT1T2 = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetT2Corr","type2"))
)


process.patPFMetT1T2Corr = cms.EDProducer("PATPFJetMETcorrInputProducer",
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L3Absolute"),
    offsetCorrLabel = cms.InputTag("L1FastJet"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("cleanedPatJets"),
    type1JetPtThreshold = cms.double(15.0)
)


process.patPFMetT1T2Smear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetT2SmearCorr","type2"))
)


process.patPFMetT1T2SmearCorr = cms.EDProducer("PATPFJetMETcorrInputProducer",
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L3Absolute"),
    offsetCorrLabel = cms.InputTag("L1FastJet"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("selectedPatJetsForMetT1T2SmearCorr"),
    type1JetPtThreshold = cms.double(15.0)
)


process.patPFMetT1T2Txy = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetT2Corr","type2"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT1T2TxySmear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetT2SmearCorr","type2"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT1TauEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrTauEnDown"))
)


process.patPFMetT1TauEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrTauEnUp"))
)


process.patPFMetT1Txy = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2Corr","type1"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT1TxySmear = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetT1T2SmearCorr","type1"), cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetT1UnclusteredEnDown = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrUnclusteredEnDown"))
)


process.patPFMetT1UnclusteredEnUp = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMetT1"),
    srcCorrections = cms.VInputTag(cms.InputTag("shiftedPatMETCorrUnclusteredEnUp"))
)


process.patPFMetT2Corr = cms.EDProducer("PATPFJetMETcorrInputProducer",
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L3Absolute"),
    offsetCorrLabel = cms.InputTag("L1FastJet"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("cleanedPatJets"),
    type1JetPtThreshold = cms.double(15.0)
)


process.patPFMetT2SmearCorr = cms.EDProducer("PATPFJetMETcorrInputProducer",
    jetCorrLabel = cms.InputTag("L3Absolute"),
    jetCorrLabelRes = cms.InputTag("L3Absolute"),
    offsetCorrLabel = cms.InputTag("L1FastJet"),
    skipEM = cms.bool(True),
    skipEMfractionThreshold = cms.double(0.9),
    skipMuonSelection = cms.string('isGlobalMuon | isStandAloneMuon'),
    skipMuons = cms.bool(True),
    src = cms.InputTag("selectedPatJetsForMetT2SmearCorr"),
    type1JetPtThreshold = cms.double(15.0)
)


process.patPFMetTxy = cms.EDProducer("CorrectedPATMETProducer",
    src = cms.InputTag("patPFMet"),
    srcCorrections = cms.VInputTag(cms.InputTag("patPFMetTxyCorr"))
)


process.patPFMetTxyCorr = cms.EDProducer("MultShiftMETcorrInputProducer",
    parameters = cms.VPSet(cms.PSet(
        etaMax = cms.double(2.7),
        etaMin = cms.double(0),
        fx = cms.string('(x*[0])+(sq(x)*[1])'),
        fy = cms.string('(x*[0])+(sq(x)*[1])'),
        name = cms.string('hEtaPlus'),
        px = cms.vdouble(-0.00229295500096, 3.15487850373e-07),
        py = cms.vdouble(0.000114282381437, -1.58467325852e-08),
        type = cms.int32(1),
        varType = cms.int32(0)
    ), 
        cms.PSet(
            etaMax = cms.double(0),
            etaMin = cms.double(-2.7),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('hEtaMinus'),
            px = cms.vdouble(-0.000198571488347, -1.94054852726e-07),
            py = cms.vdouble(-0.00137832489313, -2.02238617742e-06),
            type = cms.int32(1),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(1.392),
            etaMin = cms.double(-1.392),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('h0Barrel'),
            px = cms.vdouble(-0.0153652906396, -3.80210366974e-05),
            py = cms.vdouble(0.00798098092474, -0.000103998219585),
            type = cms.int32(5),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(3),
            etaMin = cms.double(1.392),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('h0EndcapPlus'),
            px = cms.vdouble(-0.00305719113962, -0.00032676418359),
            py = cms.vdouble(-0.00345131507897, 0.000164816815994),
            type = cms.int32(5),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-1.392),
            etaMin = cms.double(-3.0),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('h0EndcapMinus'),
            px = cms.vdouble(-0.000159031461755, 0.00012231873804),
            py = cms.vdouble(0.0260436390996, -8.17994745657e-05),
            type = cms.int32(5),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(1.479),
            etaMin = cms.double(-1.479),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('gammaBarrel'),
            px = cms.vdouble(-0.00163144589987, 3.17557692226e-06),
            py = cms.vdouble(-0.000710945802217, 6.45810884842e-06),
            type = cms.int32(4),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(3.0),
            etaMin = cms.double(1.479),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('gammaEndcapPlus'),
            px = cms.vdouble(-0.00108893779312, -2.53584544941e-05),
            py = cms.vdouble(0.00188026342884, 8.15028097381e-05),
            type = cms.int32(4),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-1.479),
            etaMin = cms.double(-3.0),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('gammaEndcapMinus'),
            px = cms.vdouble(-0.00130486432072, 1.72313009972e-05),
            py = cms.vdouble(-0.00367119684052, -1.63143116342e-05),
            type = cms.int32(4),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(5.2),
            etaMin = cms.double(2.901376),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('hHFPlus'),
            px = cms.vdouble(-0.000218928792083, -1.0492437382e-06),
            py = cms.vdouble(2.7982430778e-05, -6.87804028426e-08),
            type = cms.int32(6),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-2.901376),
            etaMin = cms.double(-5.2),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('hHFMinus'),
            px = cms.vdouble(-0.000851170798547, 3.18768998961e-07),
            py = cms.vdouble(6.10447368609e-05, -5.92655106387e-07),
            type = cms.int32(6),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(5.2),
            etaMin = cms.double(2.901376),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('egammaHFPlus'),
            px = cms.vdouble(0.00138084425101, -6.39459000901e-06),
            py = cms.vdouble(-0.000532336534523, 2.21305870813e-06),
            type = cms.int32(7),
            varType = cms.int32(0)
        ), 
        cms.PSet(
            etaMax = cms.double(-2.901376),
            etaMin = cms.double(-5.2),
            fx = cms.string('(x*[0])+(sq(x)*[1])'),
            fy = cms.string('(x*[0])+(sq(x)*[1])'),
            name = cms.string('egammaHFMinus'),
            px = cms.vdouble(0.00102598393499, -3.37284909389e-06),
            py = cms.vdouble(0.000439449053802, -2.3750891943e-06),
            type = cms.int32(7),
            varType = cms.int32(0)
        )),
    srcPFlow = cms.InputTag("packedPFCandidates"),
    vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.patSmearedJets = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFchs'),
    algopt = cms.string('AK4PFchs_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJets"),
    useDeterministicSeed = cms.bool(True),
    variation = cms.int32(0)
)


process.patTrkMet = cms.EDProducer("PATMETProducer",
    addEfficiencies = cms.bool(False),
    addGenMET = cms.bool(False),
    addMuonCorrections = cms.bool(False),
    addResolutions = cms.bool(False),
    computeMETSignificance = cms.bool(False),
    computeMETSignificant = cms.bool(False),
    efficiencies = cms.PSet(

    ),
    genMETSource = cms.InputTag("genMetTrue"),
    metSource = cms.InputTag("pfMetTrk"),
    muonSource = cms.InputTag("muons"),
    parameters = cms.PSet(
        dRMatch = cms.double(0.4),
        jetThreshold = cms.double(15),
        jeta = cms.vdouble(0.8, 1.3, 1.9, 2.5),
        jpar = cms.vdouble(1.29, 1.19, 1.07, 1.13, 1.12),
        pjpar = cms.vdouble(-0.04, 0.6504)
    ),
    resolutions = cms.PSet(

    ),
    srcJetResPhi = cms.string('AK4PFchs_phi'),
    srcJetResPt = cms.string('AK4PFchs_pt'),
    srcJetSF = cms.string('AK4PFchs'),
    srcJets = cms.InputTag("cleanedPatJets"),
    srcLeptons = cms.VInputTag("selectedPatElectrons", "selectedPatMuons", "selectedPatPhotons"),
    srcPFCands = cms.InputTag("particleFlow"),
    srcRho = cms.InputTag("fixedGridRhoAll"),
    userData = cms.PSet(
        userCands = cms.PSet(
            src = cms.VInputTag("")
        ),
        userClasses = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFloats = cms.PSet(
            src = cms.VInputTag("")
        ),
        userFunctionLabels = cms.vstring(),
        userFunctions = cms.vstring(),
        userInts = cms.PSet(
            src = cms.VInputTag("")
        )
    )
)


process.pdfinfoMaker = cms.EDProducer("PDFInfoMaker",
    aliasPrefix = cms.untracked.string('pdfinfo'),
    genEventInfoInputTag = cms.string('generator'),
    hepmcHandle = cms.string('generator')
)


process.pfCandMETcorr = cms.EDProducer("PFCandMETcorrInputProducer",
    src = cms.InputTag("pfCandsNotInJetsForMetCorr")
)


process.pfCandidateMaker = cms.EDProducer("PFCandidateMaker",
    aliasPrefix = cms.untracked.string('pfcands'),
    minPt = cms.double(5.0),
    pfCandidatesTag = cms.InputTag("packedPFCandidates","","PAT")
)


process.pfCandidateToVertexAssociation = cms.EDProducer("PFCand_AssoMap",
    AssociationType = cms.InputTag("Both"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    ConversionsCollection = cms.InputTag("allConversions"),
    FinalAssociation = cms.untracked.int32(1),
    GetCleanedCollections = cms.bool(False),
    MaxNumberOfAssociations = cms.int32(1),
    NIVertexCollection = cms.InputTag("particleFlowDisplacedVertex"),
    PFCandidateCollection = cms.InputTag("particleFlow"),
    UseBeamSpotCompatibility = cms.untracked.bool(True),
    V0KshortCollection = cms.InputTag("generalV0Candidates","Kshort"),
    V0LambdaCollection = cms.InputTag("generalV0Candidates","Lambda"),
    VertexCollection = cms.InputTag("offlinePrimaryVertices"),
    doReassociation = cms.bool(True),
    ignoreMissingCollection = cms.bool(True),
    nTrackWeight = cms.double(0.001)
)


process.pfCandsForUnclusteredUnc = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("pfCandsNoJetsNoEleNoMuNoTau"),
    veto = cms.InputTag("slimmedPhotons")
)


process.pfCandsNoJets = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("packedPFCandidates"),
    veto = cms.InputTag("cleanedPatJets")
)


process.pfCandsNoJetsNoEle = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("pfCandsNoJets"),
    veto = cms.InputTag("slimmedElectrons")
)


process.pfCandsNoJetsNoEleNoMu = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("pfCandsNoJetsNoEle"),
    veto = cms.InputTag("slimmedMuons")
)


process.pfCandsNoJetsNoEleNoMuNoTau = cms.EDProducer("CandPtrProjector",
    src = cms.InputTag("pfCandsNoJetsNoEleNoMu"),
    veto = cms.InputTag("slimmedTaus")
)


process.pfCandsNotInJetsForMetCorr = cms.EDProducer("PFCandidateFromFwdPtrProducer",
    src = cms.InputTag("pfCandsNotInJetsPtrForMetCorr")
)


process.pfCandsNotInJetsPtrForMetCorr = cms.EDProducer("TPPFJetsOnPFCandidates",
    bottomCollection = cms.InputTag("particleFlowPtrs"),
    enable = cms.bool(True),
    name = cms.untracked.string('noJet'),
    topCollection = cms.InputTag("pfJetsPtrForMetCorr"),
    verbose = cms.untracked.bool(False)
)


process.pfJetMaker = cms.EDProducer("PFJetMaker",
    aliasPrefix = cms.untracked.string('pfjets'),
    pfCandidatesTag = cms.InputTag("packedPFCandidates","","PAT"),
    pfJetPtCut = cms.double(5.0),
    pfJetsInputTag = cms.InputTag("slimmedJets","","PAT")
)


process.pfJetPUPPIMaker = cms.EDProducer("PFJetMaker",
    aliasPrefix = cms.untracked.string('pfjets_puppi'),
    pfCandidatesTag = cms.InputTag("packedPFCandidates","","PAT"),
    pfJetPtCut = cms.double(20.0),
    pfJetsInputTag = cms.InputTag("slimmedJetsPuppi","","PAT")
)


process.pfJetsPtrForMetCorr = cms.EDProducer("PFJetFwdPtrProducer",
    src = cms.InputTag("ak4PFJets")
)


process.pfMETcorrType0 = cms.EDProducer("Type0PFMETcorrInputProducer",
    correction = cms.PSet(
        formula = cms.string('(x<35)?(-( [0]+x*[1]+pow(x, 2)*[2]+pow(x, 3)*[3] )):(-( [0]+35*[1]+pow(35, 2)*[2]+pow(35, 3)*[3] ))'),
        par0 = cms.double(-0.181414),
        par1 = cms.double(-0.476934),
        par2 = cms.double(0.00863564),
        par3 = cms.double(-4.94181e-05)
    ),
    minDz = cms.double(0.2),
    srcHardScatterVertex = cms.InputTag("selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0"),
    srcPFCandidateToVertexAssociations = cms.InputTag("pfCandidateToVertexAssociation")
)


process.pfMet = cms.EDProducer("RecoMETExtractor",
    correctionLevel = cms.string('raw'),
    metSource = cms.InputTag("slimmedMETs","","@skipCurrentProcess")
)


process.pfMetCHS = cms.EDProducer("PFMETProducer",
    alias = cms.string('pfMet'),
    calculateSignificance = cms.bool(False),
    globalThreshold = cms.double(0.0),
    src = cms.InputTag("pfCHS")
)


process.pfMetTrk = cms.EDProducer("PFMETProducer",
    alias = cms.string('pfMet'),
    calculateSignificance = cms.bool(False),
    globalThreshold = cms.double(0.0),
    src = cms.InputTag("pfTrk")
)


process.pfmetMaker = cms.EDProducer("PFMETMaker",
    aliasPrefix = cms.untracked.string('evt'),
    doUncertainties = cms.bool(True),
    onlySaveTwoVector = cms.bool(False),
    pfMetInputTag_ = cms.InputTag("slimmedMETs","","PAT")
)


process.pfmetpuppiMaker = cms.EDProducer("PFMETMaker",
    aliasPrefix = cms.untracked.string('evt_puppi'),
    doUncertainties = cms.bool(True),
    onlySaveTwoVector = cms.bool(False),
    pfMetInputTag_ = cms.InputTag("slimmedMETsPuppi","","PAT")
)


process.pftauMaker = cms.EDProducer("PFTauMaker",
    aliasPrefix = cms.untracked.string('taus_pf'),
    pftausInputTag = cms.InputTag("slimmedTaus")
)


process.photonMaker = cms.EDProducer("PhotonMaker",
    aliasPrefix = cms.untracked.string('photons'),
    minEt = cms.double(10.0),
    photonsInputTag = cms.InputTag("slimmedPhotons")
)


process.puSummaryInfoMaker = cms.EDProducer("PUSummaryInfoMaker",
    PUInfoInputTag = cms.InputTag("slimmedAddPileupInfo"),
    aliasPrefix = cms.untracked.string('puInfo')
)


process.randomEngineStateProducer = cms.EDProducer("RandomEngineStateProducer")


process.secondaryVertexMaker = cms.EDProducer("SecondaryVertexMaker",
    aliasPrefix = cms.untracked.string('svs'),
    inclusiveVertexInputTag = cms.InputTag("slimmedSecondaryVertices","","PAT"),
    primaryVertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.shiftedPatElectronEnDown = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("slimmedElectrons"),
    uncertainty = cms.string('((abs(y)<1.479)?(0.006+0*x):(0.015+0*x))')
)


process.shiftedPatElectronEnUp = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(1.0),
    src = cms.InputTag("slimmedElectrons"),
    uncertainty = cms.string('((abs(y)<1.479)?(0.006+0*x):(0.015+0*x))')
)


process.shiftedPatJetEnDown = cms.EDProducer("ShiftedPATJetProducer",
    addResidualJES = cms.bool(True),
    jetCorrLabelUpToL3 = cms.InputTag("ak4PFCHSL1FastL2L3Corrector"),
    jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector"),
    jetCorrPayloadName = cms.string('AK4PFchs'),
    jetCorrUncertaintyTag = cms.string('Uncertainty'),
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("cleanedPatJets")
)


process.shiftedPatJetEnUp = cms.EDProducer("ShiftedPATJetProducer",
    addResidualJES = cms.bool(True),
    jetCorrLabelUpToL3 = cms.InputTag("ak4PFCHSL1FastL2L3Corrector"),
    jetCorrLabelUpToL3Res = cms.InputTag("ak4PFCHSL1FastL2L3Corrector"),
    jetCorrPayloadName = cms.string('AK4PFchs'),
    jetCorrUncertaintyTag = cms.string('Uncertainty'),
    shiftBy = cms.double(1.0),
    src = cms.InputTag("cleanedPatJets")
)


process.shiftedPatJetResDown = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFchs'),
    algopt = cms.string('AK4PFchs_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJets"),
    useDeterministicSeed = cms.bool(True),
    variation = cms.int32(-101)
)


process.shiftedPatJetResUp = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFchs'),
    algopt = cms.string('AK4PFchs_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJets"),
    useDeterministicSeed = cms.bool(True),
    variation = cms.int32(101)
)


process.shiftedPatMETCorrElectronEnDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("slimmedElectrons"),
    srcShifted = cms.InputTag("shiftedPatElectronEnDown")
)


process.shiftedPatMETCorrElectronEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("slimmedElectrons"),
    srcShifted = cms.InputTag("shiftedPatElectronEnUp")
)


process.shiftedPatMETCorrJetEnDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("cleanedPatJets"),
    srcShifted = cms.InputTag("shiftedPatJetEnDown")
)


process.shiftedPatMETCorrJetEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("cleanedPatJets"),
    srcShifted = cms.InputTag("shiftedPatJetEnUp")
)


process.shiftedPatMETCorrJetResDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("cleanedPatJets"),
    srcShifted = cms.InputTag("shiftedPatJetResDown")
)


process.shiftedPatMETCorrJetResUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("cleanedPatJets"),
    srcShifted = cms.InputTag("shiftedPatJetResUp")
)


process.shiftedPatMETCorrMuonEnDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("slimmedMuons"),
    srcShifted = cms.InputTag("shiftedPatMuonEnDown")
)


process.shiftedPatMETCorrMuonEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("slimmedMuons"),
    srcShifted = cms.InputTag("shiftedPatMuonEnUp")
)


process.shiftedPatMETCorrPhotonEnDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("slimmedPhotons"),
    srcShifted = cms.InputTag("shiftedPatPhotonEnDown")
)


process.shiftedPatMETCorrPhotonEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("slimmedPhotons"),
    srcShifted = cms.InputTag("shiftedPatPhotonEnUp")
)


process.shiftedPatMETCorrSmearedJetResDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("selectedPatJetsForMetT1T2SmearCorr"),
    srcShifted = cms.InputTag("shiftedPatSmearedJetResDown")
)


process.shiftedPatMETCorrSmearedJetResUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("selectedPatJetsForMetT1T2SmearCorr"),
    srcShifted = cms.InputTag("shiftedPatSmearedJetResUp")
)


process.shiftedPatMETCorrTauEnDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("slimmedTaus"),
    srcShifted = cms.InputTag("shiftedPatTauEnDown")
)


process.shiftedPatMETCorrTauEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("slimmedTaus"),
    srcShifted = cms.InputTag("shiftedPatTauEnUp")
)


process.shiftedPatMETCorrUnclusteredEnDown = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfCandsForUnclusteredUnc"),
    srcShifted = cms.InputTag("shiftedPatUnclusteredEnDown")
)


process.shiftedPatMETCorrUnclusteredEnUp = cms.EDProducer("ShiftedParticleMETcorrInputProducer",
    srcOriginal = cms.InputTag("pfCandsForUnclusteredUnc"),
    srcShifted = cms.InputTag("shiftedPatUnclusteredEnUp")
)


process.shiftedPatMuonEnDown = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("slimmedMuons"),
    uncertainty = cms.string('((x<100)?(0.002+0*y):(0.05+0*y))')
)


process.shiftedPatMuonEnUp = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(1.0),
    src = cms.InputTag("slimmedMuons"),
    uncertainty = cms.string('((x<100)?(0.002+0*y):(0.05+0*y))')
)


process.shiftedPatPhotonEnDown = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("slimmedPhotons"),
    uncertainty = cms.string('((abs(y)<1.479)?(0.01+0*x):(0.025+0*x))')
)


process.shiftedPatPhotonEnUp = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(1.0),
    src = cms.InputTag("slimmedPhotons"),
    uncertainty = cms.string('((abs(y)<1.479)?(0.01+0*x):(0.025+0*x))')
)


process.shiftedPatSmearedJetResDown = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFchs'),
    algopt = cms.string('AK4PFchs_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJets"),
    useDeterministicSeed = cms.bool(True),
    variation = cms.int32(-1)
)


process.shiftedPatSmearedJetResUp = cms.EDProducer("SmearedPATJetProducer",
    algo = cms.string('AK4PFchs'),
    algopt = cms.string('AK4PFchs_pt'),
    dPtMaxFactor = cms.double(3),
    dRMax = cms.double(0.2),
    debug = cms.untracked.bool(False),
    enabled = cms.bool(True),
    genJets = cms.InputTag("slimmedGenJets"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    seed = cms.uint32(37428479),
    skipGenMatching = cms.bool(False),
    src = cms.InputTag("cleanedPatJets"),
    useDeterministicSeed = cms.bool(True),
    variation = cms.int32(1)
)


process.shiftedPatTauEnDown = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("slimmedTaus"),
    uncertainty = cms.string('0.03+0*x*y')
)


process.shiftedPatTauEnUp = cms.EDProducer("ShiftedParticleProducer",
    shiftBy = cms.double(1.0),
    src = cms.InputTag("slimmedTaus"),
    uncertainty = cms.string('0.03+0*x*y')
)


process.shiftedPatUnclusteredEnDown = cms.EDProducer("ShiftedParticleProducer",
    binning = cms.VPSet(cms.PSet(
        binSelection = cms.string('charge!=0'),
        binUncertainty = cms.string('sqrt(pow(0.00009*x,2)+pow(0.0085/sqrt(sin(2*atan(exp(-y)))),2))')
    ), 
        cms.PSet(
            binSelection = cms.string('pdgId==130'),
            binUncertainty = cms.string('((abs(y)<1.3)?(min(0.25,sqrt(0.64/x+0.0025))):(min(0.30,sqrt(1.0/x+0.0016))))'),
            energyDependency = cms.bool(True)
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==22'),
            binUncertainty = cms.string('sqrt(0.0009/x+0.000001)+0*y'),
            energyDependency = cms.bool(True)
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==1 || pdgId==2'),
            binUncertainty = cms.string('sqrt(1./x+0.0025)+0*y'),
            energyDependency = cms.bool(True)
        )),
    shiftBy = cms.double(-1.0),
    src = cms.InputTag("pfCandsForUnclusteredUnc")
)


process.shiftedPatUnclusteredEnUp = cms.EDProducer("ShiftedParticleProducer",
    binning = cms.VPSet(cms.PSet(
        binSelection = cms.string('charge!=0'),
        binUncertainty = cms.string('sqrt(pow(0.00009*x,2)+pow(0.0085/sqrt(sin(2*atan(exp(-y)))),2))')
    ), 
        cms.PSet(
            binSelection = cms.string('pdgId==130'),
            binUncertainty = cms.string('((abs(y)<1.3)?(min(0.25,sqrt(0.64/x+0.0025))):(min(0.30,sqrt(1.0/x+0.0016))))'),
            energyDependency = cms.bool(True)
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==22'),
            binUncertainty = cms.string('sqrt(0.0009/x+0.000001)+0*y'),
            energyDependency = cms.bool(True)
        ), 
        cms.PSet(
            binSelection = cms.string('pdgId==1 || pdgId==2'),
            binUncertainty = cms.string('sqrt(1./x+0.0025)+0*y'),
            energyDependency = cms.bool(True)
        )),
    shiftBy = cms.double(1.0),
    src = cms.InputTag("pfCandsForUnclusteredUnc")
)


process.slimmedMETs = cms.EDProducer("PATMETSlimmer",
    caloMET = cms.InputTag("patCaloMet"),
    chsMET = cms.InputTag("patCHSMet"),
    rawVariation = cms.InputTag("patPFMet"),
    runningOnMiniAOD = cms.bool(True),
    src = cms.InputTag("patPFMetT1"),
    t01Variation = cms.InputTag("slimmedMETs","","@skipCurrentProcess"),
    t1SmearedVarsAndUncs = cms.InputTag("patPFMetT1Smear%s"),
    t1Uncertainties = cms.InputTag("patPFMetT1%s"),
    tXYUncForT1 = cms.InputTag("patPFMetT1Txy"),
    trkMET = cms.InputTag("patTrkMet")
)


process.subJetMaker = cms.EDProducer("SubJetMaker",
    lessBranches = cms.bool(False),
    pfCandidatesTag = cms.InputTag("packedPFCandidates"),
    pfJetPtCut = cms.double(200),
    pfJetsInputTag = cms.InputTag("slimmedJetsAK8")
)


process.unpackedTracksAndVertices = cms.EDProducer("PATTrackAndVertexUnpacker",
    additionalTracks = cms.InputTag("lostTracks"),
    packedCandidates = cms.InputTag("packedPFCandidates"),
    slimmedSecondaryVertices = cms.InputTag("slimmedSecondaryVertices"),
    slimmedVertices = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.vertexMaker = cms.EDProducer("VertexMaker",
    aliasPrefix = cms.untracked.string('vtxs'),
    primaryVertexInputTag = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.jetSelectorForMet = cms.EDFilter("PATJetSelector",
    cut = cms.string('pt>15 && abs(eta)<9.9'),
    cutLoose = cms.string(''),
    nLoose = cms.uint32(0),
    src = cms.InputTag("basicJetsForMet")
)


process.pfCHS = cms.EDFilter("CandPtrSelector",
    cut = cms.string('fromPV(0)>0'),
    src = cms.InputTag("packedPFCandidates")
)


process.pfTrk = cms.EDFilter("CandPtrSelector",
    cut = cms.string('fromPV(0) > 0 && charge()!=0'),
    src = cms.InputTag("packedPFCandidates")
)


process.selectedPatJetsForMetT1T2Corr = cms.EDFilter("PATJetSelector",
    cut = cms.string('abs(eta) < 9.9'),
    filter = cms.bool(False),
    src = cms.InputTag("patJets")
)


process.selectedPatJetsForMetT1T2SmearCorr = cms.EDFilter("PATJetSelector",
    cut = cms.string('abs(eta) < 9.9'),
    filter = cms.bool(False),
    src = cms.InputTag("patSmearedJets")
)


process.selectedPatJetsForMetT2Corr = cms.EDFilter("PATJetSelector",
    cut = cms.string('abs(eta) > 9.9'),
    filter = cms.bool(False),
    src = cms.InputTag("patJets")
)


process.selectedPatJetsForMetT2SmearCorr = cms.EDFilter("PATJetSelector",
    cut = cms.string('abs(eta) > 9.9'),
    filter = cms.bool(False),
    src = cms.InputTag("patSmearedJets")
)


process.selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0 = cms.EDFilter("PATSingleVertexSelector",
    filter = cms.bool(False),
    mode = cms.string('firstVertex'),
    vertices = cms.InputTag("selectedVerticesForPFMEtCorrType0")
)


process.selectedVerticesForPFMEtCorrType0 = cms.EDFilter("VertexSelector",
    cut = cms.string('isValid & ndof >= 4 & chi2 > 0 & tracksSize > 0 & abs(z) < 24 & abs(position.Rho) < 2.'),
    filter = cms.bool(False),
    src = cms.InputTag("offlinePrimaryVertices")
)


process.out = cms.OutputModule("PoolOutputModule",
    basketSize = cms.untracked.int32(376832),
    dropMetaData = cms.untracked.string('ALL'),
    fileName = cms.untracked.string('ntuple.root'),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_*Maker*_*_CMS3*')
)


process.DQMStore = cms.Service("DQMStore",
    LSbasedMode = cms.untracked.bool(False),
    collateHistograms = cms.untracked.bool(False),
    enableMultiThread = cms.untracked.bool(False),
    forceResetOnBeginLumi = cms.untracked.bool(False),
    referenceFileName = cms.untracked.string(''),
    verbose = cms.untracked.int32(0),
    verboseQT = cms.untracked.int32(0)
)


process.MessageLogger = cms.Service("MessageLogger",
    FrameworkJobReport = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        optionalPSet = cms.untracked.bool(True)
    ),
    categories = cms.untracked.vstring('FwkJob', 
        'FwkReport', 
        'FwkSummary', 
        'Root_NoDictionary'),
    cerr = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        FwkReport = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(100)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1)
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        noTimeStamps = cms.untracked.bool(False),
        optionalPSet = cms.untracked.bool(True),
        threshold = cms.untracked.string('')
    ),
    cerr_stats = cms.untracked.PSet(
        optionalPSet = cms.untracked.bool(True),
        output = cms.untracked.string('cerr'),
        threshold = cms.untracked.string('WARNING')
    ),
    cout = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    debugModules = cms.untracked.vstring(),
    debugs = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    default = cms.untracked.PSet(

    ),
    destinations = cms.untracked.vstring('warnings', 
        'errors', 
        'infos', 
        'debugs', 
        'cout', 
        'cerr'),
    errors = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    fwkJobReports = cms.untracked.vstring('FrameworkJobReport'),
    infos = cms.untracked.PSet(
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        optionalPSet = cms.untracked.bool(True),
        placeholder = cms.untracked.bool(True)
    ),
    statistics = cms.untracked.vstring('cerr_stats'),
    suppressDebug = cms.untracked.vstring(),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring('ecalLaserCorrFilter', 
        'manystripclus53X', 
        'toomanystripclus53X'),
    warnings = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    )
)


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
    CTPPSFastRecHits = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(1357987)
    ),
    LHCTransport = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(87654321)
    ),
    MuonSimHits = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(987346)
    ),
    VtxSmeared = cms.PSet(
        engineName = cms.untracked.string('HepJamesRandom'),
        initialSeed = cms.untracked.uint32(98765432)
    ),
    ecalPreshowerRecHit = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(6541321)
    ),
    ecalRecHit = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(654321)
    ),
    externalLHEProducer = cms.PSet(
        engineName = cms.untracked.string('HepJamesRandom'),
        initialSeed = cms.untracked.uint32(234567)
    ),
    famosPileUp = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(918273)
    ),
    famosSimHits = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(13579)
    ),
    fastTrackerRecHits = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(24680)
    ),
    g4SimHits = cms.PSet(
        engineName = cms.untracked.string('HepJamesRandom'),
        initialSeed = cms.untracked.uint32(11)
    ),
    generator = cms.PSet(
        engineName = cms.untracked.string('HepJamesRandom'),
        initialSeed = cms.untracked.uint32(123456789)
    ),
    hbhereco = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(541321)
    ),
    hfreco = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(541321)
    ),
    hiSignal = cms.PSet(
        engineName = cms.untracked.string('HepJamesRandom'),
        initialSeed = cms.untracked.uint32(123456789)
    ),
    hiSignalG4SimHits = cms.PSet(
        engineName = cms.untracked.string('HepJamesRandom'),
        initialSeed = cms.untracked.uint32(11)
    ),
    hiSignalLHCTransport = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(88776655)
    ),
    horeco = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(541321)
    ),
    l1ParamMuons = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(6453209)
    ),
    mix = cms.PSet(
        engineName = cms.untracked.string('HepJamesRandom'),
        initialSeed = cms.untracked.uint32(12345)
    ),
    mixData = cms.PSet(
        engineName = cms.untracked.string('HepJamesRandom'),
        initialSeed = cms.untracked.uint32(12345)
    ),
    mixGenPU = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(918273)
    ),
    mixRecoTracks = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(918273)
    ),
    mixSimCaloHits = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(918273)
    ),
    paramMuons = cms.PSet(
        engineName = cms.untracked.string('TRandom3'),
        initialSeed = cms.untracked.uint32(54525)
    ),
    saveFileName = cms.untracked.string(''),
    simBeamSpotFilter = cms.PSet(
        engineName = cms.untracked.string('HepJamesRandom'),
        initialSeed = cms.untracked.uint32(87654321)
    ),
    simMuonCSCDigis = cms.PSet(
        engineName = cms.untracked.string('HepJamesRandom'),
        initialSeed = cms.untracked.uint32(11223344)
    ),
    simMuonDTDigis = cms.PSet(
        engineName = cms.untracked.string('HepJamesRandom'),
        initialSeed = cms.untracked.uint32(1234567)
    ),
    simMuonRPCDigis = cms.PSet(
        engineName = cms.untracked.string('HepJamesRandom'),
        initialSeed = cms.untracked.uint32(1234567)
    ),
    simSiStripDigiSimLink = cms.PSet(
        engineName = cms.untracked.string('HepJamesRandom'),
        initialSeed = cms.untracked.uint32(1234567)
    )
)


process.Timing = cms.Service("Timing",
    summaryOnly = cms.untracked.bool(True)
)


process.CSCGeometryESModule = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    debugV = cms.untracked.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useDDD = cms.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.CaloGeometryBuilder = cms.ESProducer("CaloGeometryBuilder",
    SelectedCalos = cms.vstring('HCAL', 
        'ZDC', 
        'CASTOR', 
        'EcalBarrel', 
        'EcalEndcap', 
        'EcalPreshower', 
        'TOWER')
)


process.CaloTopologyBuilder = cms.ESProducer("CaloTopologyBuilder")


process.CaloTowerGeometryFromDBEP = cms.ESProducer("CaloTowerGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.CaloTowerTopologyEP = cms.ESProducer("CaloTowerTopologyEP")


process.CastorDbProducer = cms.ESProducer("CastorDbProducer",
    appendToDataLabel = cms.string('')
)


process.CastorGeometryFromDBEP = cms.ESProducer("CastorGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.DTGeometryESModule = cms.ESProducer("DTGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.EcalBarrelGeometryFromDBEP = cms.ESProducer("EcalBarrelGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalElectronicsMappingBuilder = cms.ESProducer("EcalElectronicsMappingBuilder")


process.EcalEndcapGeometryFromDBEP = cms.ESProducer("EcalEndcapGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalLaserCorrectionService = cms.ESProducer("EcalLaserCorrectionService")


process.EcalPreshowerGeometryFromDBEP = cms.ESProducer("EcalPreshowerGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalTrigTowerConstituentsMapBuilder = cms.ESProducer("EcalTrigTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/EcalMapping/data/EndCap_TTMap.txt')
)


process.GlobalTrackingGeometryESProducer = cms.ESProducer("GlobalTrackingGeometryESProducer")


process.HcalAlignmentEP = cms.ESProducer("HcalAlignmentEP")


process.HcalGeometryFromDBEP = cms.ESProducer("HcalGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.MuonDetLayerGeometryESProducer = cms.ESProducer("MuonDetLayerGeometryESProducer")


process.ParabolicParametrizedMagneticFieldProducer = cms.ESProducer("AutoParametrizedMagneticFieldProducer",
    label = cms.untracked.string('ParabolicMf'),
    valueOverride = cms.int32(-1),
    version = cms.string('Parabolic')
)


process.RPCGeometryESModule = cms.ESProducer("RPCGeometryESModule",
    compatibiltyWith11 = cms.untracked.bool(True),
    useDDD = cms.untracked.bool(False)
)


process.SiStripRecHitMatcherESProducer = cms.ESProducer("SiStripRecHitMatcherESProducer",
    ComponentName = cms.string('StandardMatcher'),
    NSigmaInside = cms.double(3.0),
    PreFilter = cms.bool(False)
)


process.StripCPEfromTrackAngleESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('StripCPEfromTrackAngle'),
    ComponentType = cms.string('StripCPEfromTrackAngle'),
    parameters = cms.PSet(
        mLC_P0 = cms.double(-0.326),
        mLC_P1 = cms.double(0.618),
        mLC_P2 = cms.double(0.3),
        mTEC_P0 = cms.double(-1.885),
        mTEC_P1 = cms.double(0.471),
        mTIB_P0 = cms.double(-0.742),
        mTIB_P1 = cms.double(0.202),
        mTID_P0 = cms.double(-1.427),
        mTID_P1 = cms.double(0.433),
        mTOB_P0 = cms.double(-1.026),
        mTOB_P1 = cms.double(0.253),
        maxChgOneMIP = cms.double(6000.0),
        useLegacyError = cms.bool(False)
    )
)


process.TrackerRecoGeometryESProducer = cms.ESProducer("TrackerRecoGeometryESProducer")


process.TransientTrackBuilderESProducer = cms.ESProducer("TransientTrackBuilderESProducer",
    ComponentName = cms.string('TransientTrackBuilder')
)


process.VolumeBasedMagneticFieldESProducer = cms.ESProducer("VolumeBasedMagneticFieldESProducerFromDB",
    debugBuilder = cms.untracked.bool(False),
    label = cms.untracked.string(''),
    valueOverride = cms.int32(-1)
)


process.ZdcGeometryFromDBEP = cms.ESProducer("ZdcGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.ak10PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak10PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak10PFCHSL1Fastjet', 
        'ak10PFCHSL2Relative', 
        'ak10PFCHSL3Absolute', 
        'ak10PFCHSResidual')
)


process.ak10PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak10PFCHSL1Offset', 
        'ak10PFCHSL2Relative', 
        'ak10PFCHSL3Absolute', 
        'ak10PFCHSResidual')
)


process.ak10PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak10PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak10PFCHSL2Relative', 
        'ak10PFCHSL3Absolute')
)


process.ak10PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak10PFCHSL2Relative', 
        'ak10PFCHSL3Absolute', 
        'ak10PFCHSResidual')
)


process.ak10PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L2Relative')
)


process.ak10PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L3Absolute')
)


process.ak10PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak10PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak10PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak10PFL1Fastjet', 
        'ak10PFL2Relative', 
        'ak10PFL3Absolute', 
        'ak10PFResidual')
)


process.ak10PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak10PFL1Offset', 
        'ak10PFL2Relative', 
        'ak10PFL3Absolute', 
        'ak10PFResidual')
)


process.ak10PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak10PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak10PFL2Relative', 
        'ak10PFL3Absolute')
)


process.ak10PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak10PFL2Relative', 
        'ak10PFL3Absolute', 
        'ak10PFResidual')
)


process.ak10PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L2Relative')
)


process.ak10PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L3Absolute')
)


process.ak10PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK10PF'),
    level = cms.string('L2L3Residual')
)


process.ak1PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak1PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak1PFCHSL1Fastjet', 
        'ak1PFCHSL2Relative', 
        'ak1PFCHSL3Absolute', 
        'ak1PFCHSResidual')
)


process.ak1PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak1PFCHSL1Offset', 
        'ak1PFCHSL2Relative', 
        'ak1PFCHSL3Absolute', 
        'ak1PFCHSResidual')
)


process.ak1PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak1PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak1PFCHSL2Relative', 
        'ak1PFCHSL3Absolute')
)


process.ak1PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak1PFCHSL2Relative', 
        'ak1PFCHSL3Absolute', 
        'ak1PFCHSResidual')
)


process.ak1PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L2Relative')
)


process.ak1PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L3Absolute')
)


process.ak1PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak1PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak1PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak1PFL1Fastjet', 
        'ak1PFL2Relative', 
        'ak1PFL3Absolute', 
        'ak1PFResidual')
)


process.ak1PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak1PFL1Offset', 
        'ak1PFL2Relative', 
        'ak1PFL3Absolute', 
        'ak1PFResidual')
)


process.ak1PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak1PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak1PFL2Relative', 
        'ak1PFL3Absolute')
)


process.ak1PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak1PFL2Relative', 
        'ak1PFL3Absolute', 
        'ak1PFResidual')
)


process.ak1PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L2Relative')
)


process.ak1PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L3Absolute')
)


process.ak1PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK1PF'),
    level = cms.string('L2L3Residual')
)


process.ak2PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak2PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak2PFCHSL1Fastjet', 
        'ak2PFCHSL2Relative', 
        'ak2PFCHSL3Absolute', 
        'ak2PFCHSResidual')
)


process.ak2PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak2PFCHSL1Offset', 
        'ak2PFCHSL2Relative', 
        'ak2PFCHSL3Absolute', 
        'ak2PFCHSResidual')
)


process.ak2PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak2PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak2PFCHSL2Relative', 
        'ak2PFCHSL3Absolute')
)


process.ak2PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak2PFCHSL2Relative', 
        'ak2PFCHSL3Absolute', 
        'ak2PFCHSResidual')
)


process.ak2PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L2Relative')
)


process.ak2PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L3Absolute')
)


process.ak2PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak2PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak2PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak2PFL1Fastjet', 
        'ak2PFL2Relative', 
        'ak2PFL3Absolute', 
        'ak2PFResidual')
)


process.ak2PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak2PFL1Offset', 
        'ak2PFL2Relative', 
        'ak2PFL3Absolute', 
        'ak2PFResidual')
)


process.ak2PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak2PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak2PFL2Relative', 
        'ak2PFL3Absolute')
)


process.ak2PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak2PFL2Relative', 
        'ak2PFL3Absolute', 
        'ak2PFResidual')
)


process.ak2PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L2Relative')
)


process.ak2PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L3Absolute')
)


process.ak2PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK2PF'),
    level = cms.string('L2L3Residual')
)


process.ak3PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak3PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak3PFCHSL1Fastjet', 
        'ak3PFCHSL2Relative', 
        'ak3PFCHSL3Absolute', 
        'ak3PFCHSResidual')
)


process.ak3PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak3PFCHSL1Offset', 
        'ak3PFCHSL2Relative', 
        'ak3PFCHSL3Absolute', 
        'ak3PFCHSResidual')
)


process.ak3PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak3PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak3PFCHSL2Relative', 
        'ak3PFCHSL3Absolute')
)


process.ak3PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak3PFCHSL2Relative', 
        'ak3PFCHSL3Absolute', 
        'ak3PFCHSResidual')
)


process.ak3PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L2Relative')
)


process.ak3PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L3Absolute')
)


process.ak3PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak3PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak3PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak3PFL1Fastjet', 
        'ak3PFL2Relative', 
        'ak3PFL3Absolute', 
        'ak3PFResidual')
)


process.ak3PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak3PFL1Offset', 
        'ak3PFL2Relative', 
        'ak3PFL3Absolute', 
        'ak3PFResidual')
)


process.ak3PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak3PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak3PFL2Relative', 
        'ak3PFL3Absolute')
)


process.ak3PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak3PFL2Relative', 
        'ak3PFL3Absolute', 
        'ak3PFResidual')
)


process.ak3PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L2Relative')
)


process.ak3PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L3Absolute')
)


process.ak3PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK3PF'),
    level = cms.string('L2L3Residual')
)


process.ak4CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Fastjet', 
        'ak4CaloL2Relative', 
        'ak4CaloL3Absolute')
)


process.ak4CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Fastjet', 
        'ak4CaloL2Relative', 
        'ak4CaloL3Absolute', 
        'ak4CaloL6SLB')
)


process.ak4CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Fastjet', 
        'ak4CaloL2Relative', 
        'ak4CaloL3Absolute', 
        'ak4CaloResidual')
)


process.ak4CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.ak4CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Offset', 
        'ak4CaloL2Relative', 
        'ak4CaloL3Absolute')
)


process.ak4CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Offset', 
        'ak4CaloL2Relative', 
        'ak4CaloL3Absolute', 
        'ak4CaloResidual')
)


process.ak4CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak4CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL2Relative', 
        'ak4CaloL3Absolute')
)


process.ak4CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL2Relative', 
        'ak4CaloL3Absolute', 
        'ak4CaloL6SLB')
)


process.ak4CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL2Relative', 
        'ak4CaloL3Absolute', 
        'ak4CaloResidual')
)


process.ak4CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L2Relative')
)


process.ak4CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L3Absolute')
)


process.ak4CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak4CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak4CaloJetsSoftMuonTagInfos")
)


process.ak4CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5Calo'),
    level = cms.string('L2L3Residual')
)


process.ak4JPTL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4L1JPTFastjet', 
        'ak4JPTL2Relative', 
        'ak4JPTL3Absolute')
)


process.ak4JPTL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4L1JPTFastjet', 
        'ak4JPTL2Relative', 
        'ak4JPTL3Absolute', 
        'ak4JPTResidual')
)


process.ak4JPTL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4L1JPTOffset', 
        'ak4JPTL2Relative', 
        'ak4JPTL3Absolute')
)


process.ak4JPTL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4L1JPTOffset', 
        'ak4JPTL2Relative', 
        'ak4JPTL3Absolute', 
        'ak4JPTResidual')
)


process.ak4JPTL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak4JPTL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4L1JPTOffset', 
        'ak4JPTL2Relative', 
        'ak4JPTL3Absolute')
)


process.ak4JPTL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4L1JPTOffset', 
        'ak4JPTL2Relative', 
        'ak4JPTL3Absolute', 
        'ak4JPTResidual')
)


process.ak4JPTL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L2Relative')
)


process.ak4JPTL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L3Absolute')
)


process.ak4JPTResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L2L3Residual')
)


process.ak4L1JPTFastjet = cms.ESProducer("L1JPTOffsetCorrectionESProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L1JPTOffset'),
    offsetService = cms.string('ak4CaloL1Fastjet')
)


process.ak4L1JPTOffset = cms.ESProducer("L1JPTOffsetCorrectionESProducer",
    algorithm = cms.string('AK4JPT'),
    level = cms.string('L1JPTOffset'),
    offsetService = cms.string('ak4CaloL1Offset')
)


process.ak4PFCHSL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFCHSL1Fastjet', 
        'ak4PFCHSL2Relative', 
        'ak4PFCHSL3Absolute')
)


process.ak4PFCHSL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFCHSL1Fastjet', 
        'ak4PFCHSL2Relative', 
        'ak4PFCHSL3Absolute', 
        'ak4PFCHSResidual')
)


process.ak4PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak4PFCHSL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFCHSL1Offset', 
        'ak4PFCHSL2Relative', 
        'ak4PFCHSL3Absolute')
)


process.ak4PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFCHSL1Offset', 
        'ak4PFCHSL2Relative', 
        'ak4PFCHSL3Absolute', 
        'ak4PFCHSResidual')
)


process.ak4PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak4PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFCHSL2Relative', 
        'ak4PFCHSL3Absolute')
)


process.ak4PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFCHSL2Relative', 
        'ak4PFCHSL3Absolute', 
        'ak4PFCHSResidual')
)


process.ak4PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2Relative')
)


process.ak4PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L3Absolute')
)


process.ak4PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak4PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'ak4PFL2Relative', 
        'ak4PFL3Absolute')
)


process.ak4PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'ak4PFL2Relative', 
        'ak4PFL3Absolute', 
        'ak4PFL6SLB')
)


process.ak4PFL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'ak4PFL2Relative', 
        'ak4PFL3Absolute', 
        'ak4PFResidual')
)


process.ak4PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak4PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Offset', 
        'ak4PFL2Relative', 
        'ak4PFL3Absolute')
)


process.ak4PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Offset', 
        'ak4PFL2Relative', 
        'ak4PFL3Absolute', 
        'ak4PFResidual')
)


process.ak4PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak4PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL2Relative', 
        'ak4PFL3Absolute')
)


process.ak4PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL2Relative', 
        'ak4PFL3Absolute', 
        'ak4PFL6SLB')
)


process.ak4PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL2Relative', 
        'ak4PFL3Absolute', 
        'ak4PFResidual')
)


process.ak4PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L2Relative')
)


process.ak4PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L3Absolute')
)


process.ak4PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak4PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak4PFJetsSoftMuonTagInfos")
)


process.ak4PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK4PF'),
    level = cms.string('L2L3Residual')
)


process.ak4TrackL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Fastjet', 
        'ak4TrackL2Relative', 
        'ak4TrackL3Absolute')
)


process.ak4TrackL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4TrackL2Relative', 
        'ak4TrackL3Absolute')
)


process.ak4TrackL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5TRK'),
    level = cms.string('L2Relative')
)


process.ak4TrackL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5TRK'),
    level = cms.string('L3Absolute')
)


process.ak5PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak5PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak5PFCHSL1Fastjet', 
        'ak5PFCHSL2Relative', 
        'ak5PFCHSL3Absolute', 
        'ak5PFCHSResidual')
)


process.ak5PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak5PFCHSL1Offset', 
        'ak5PFCHSL2Relative', 
        'ak5PFCHSL3Absolute', 
        'ak5PFCHSResidual')
)


process.ak5PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak5PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak5PFCHSL2Relative', 
        'ak5PFCHSL3Absolute')
)


process.ak5PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak5PFCHSL2Relative', 
        'ak5PFCHSL3Absolute', 
        'ak5PFCHSResidual')
)


process.ak5PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L2Relative')
)


process.ak5PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L3Absolute')
)


process.ak5PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak5PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak5PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak5PFL1Fastjet', 
        'ak5PFL2Relative', 
        'ak5PFL3Absolute', 
        'ak5PFResidual')
)


process.ak5PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak5PFL1Offset', 
        'ak5PFL2Relative', 
        'ak5PFL3Absolute', 
        'ak5PFResidual')
)


process.ak5PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak5PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak5PFL2Relative', 
        'ak5PFL3Absolute')
)


process.ak5PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak5PFL2Relative', 
        'ak5PFL3Absolute', 
        'ak5PFResidual')
)


process.ak5PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L2Relative')
)


process.ak5PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L3Absolute')
)


process.ak5PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK5PF'),
    level = cms.string('L2L3Residual')
)


process.ak6PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak6PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak6PFCHSL1Fastjet', 
        'ak6PFCHSL2Relative', 
        'ak6PFCHSL3Absolute', 
        'ak6PFCHSResidual')
)


process.ak6PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak6PFCHSL1Offset', 
        'ak6PFCHSL2Relative', 
        'ak6PFCHSL3Absolute', 
        'ak6PFCHSResidual')
)


process.ak6PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak6PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak6PFCHSL2Relative', 
        'ak6PFCHSL3Absolute')
)


process.ak6PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak6PFCHSL2Relative', 
        'ak6PFCHSL3Absolute', 
        'ak6PFCHSResidual')
)


process.ak6PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L2Relative')
)


process.ak6PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L3Absolute')
)


process.ak6PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak6PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak6PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak6PFL1Fastjet', 
        'ak6PFL2Relative', 
        'ak6PFL3Absolute', 
        'ak6PFResidual')
)


process.ak6PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak6PFL1Offset', 
        'ak6PFL2Relative', 
        'ak6PFL3Absolute', 
        'ak6PFResidual')
)


process.ak6PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak6PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak6PFL2Relative', 
        'ak6PFL3Absolute')
)


process.ak6PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak6PFL2Relative', 
        'ak6PFL3Absolute', 
        'ak6PFResidual')
)


process.ak6PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L2Relative')
)


process.ak6PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L3Absolute')
)


process.ak6PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK6PF'),
    level = cms.string('L2L3Residual')
)


process.ak7CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Fastjet', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute')
)


process.ak7CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7CaloL1Offset', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloL6SLB')
)


process.ak7CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7CaloL1Fastjet', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloResidual')
)


process.ak7CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.ak7CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7CaloL1Offset', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute')
)


process.ak7CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7CaloL1Offset', 
        'ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloResidual')
)


process.ak7CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak7CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7CaloL2Relative', 
        'ak7CaloL3Absolute')
)


process.ak7CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloL6SLB')
)


process.ak7CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7CaloL2Relative', 
        'ak7CaloL3Absolute', 
        'ak7CaloResidual')
)


process.ak7CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L2Relative')
)


process.ak7CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L3Absolute')
)


process.ak7CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak7CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak7CaloJetsSoftMuonTagInfos")
)


process.ak7CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7Calo'),
    level = cms.string('L2L3Residual')
)


process.ak7JPTL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7L1JPTFastjet', 
        'ak7JPTL2Relative', 
        'ak7JPTL3Absolute')
)


process.ak7JPTL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7L1JPTFastjet', 
        'ak7JPTL2Relative', 
        'ak7JPTL3Absolute', 
        'ak7JPTResidual')
)


process.ak7JPTL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7L1JPTOffset', 
        'ak7JPTL2Relative', 
        'ak7JPTL3Absolute')
)


process.ak7JPTL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7L1JPTOffset', 
        'ak7JPTL2Relative', 
        'ak7JPTL3Absolute', 
        'ak7JPTResidual')
)


process.ak7JPTL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7L1JPTOffset', 
        'ak7JPTL2Relative', 
        'ak7JPTL3Absolute')
)


process.ak7L1JPTFastjet = cms.ESProducer("L1JPTOffsetCorrectionESProducer",
    algorithm = cms.string('AK7JPT'),
    level = cms.string('L1JPTOffset'),
    offsetService = cms.string('ak7CaloL1Fastjet')
)


process.ak7L1JPTOffset = cms.ESProducer("L1JPTOffsetCorrectionESProducer",
    algorithm = cms.string('AK7JPT'),
    level = cms.string('L1JPTOffset'),
    offsetService = cms.string('ak7CaloL1Offset')
)


process.ak7PFCHSL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFCHSL1Fastjet', 
        'ak7PFCHSL2Relative', 
        'ak7PFCHSL3Absolute')
)


process.ak7PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak7PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFCHSL1Fastjet', 
        'ak7PFCHSL2Relative', 
        'ak7PFCHSL3Absolute', 
        'ak7PFCHSResidual')
)


process.ak7PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFCHSL1Offset', 
        'ak7PFCHSL2Relative', 
        'ak7PFCHSL3Absolute', 
        'ak7PFCHSResidual')
)


process.ak7PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak7PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFCHSL2Relative', 
        'ak7PFCHSL3Absolute')
)


process.ak7PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFCHSL2Relative', 
        'ak7PFCHSL3Absolute', 
        'ak7PFCHSResidual')
)


process.ak7PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L2Relative')
)


process.ak7PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L3Absolute')
)


process.ak7PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak7PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute')
)


process.ak7PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFL6SLB')
)


process.ak7PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak7PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFL1Fastjet', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFResidual')
)


process.ak7PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFL1Offset', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute')
)


process.ak7PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFL1Offset', 
        'ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFResidual')
)


process.ak7PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak7PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFL2Relative', 
        'ak7PFL3Absolute')
)


process.ak7PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFL6SLB')
)


process.ak7PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak7PFL2Relative', 
        'ak7PFL3Absolute', 
        'ak7PFResidual')
)


process.ak7PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L2Relative')
)


process.ak7PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L3Absolute')
)


process.ak7PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ak7PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ak7PFJetsSoftMuonTagInfos")
)


process.ak7PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK7PF'),
    level = cms.string('L2L3Residual')
)


process.ak8PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak8PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak8PFCHSL1Fastjet', 
        'ak8PFCHSL2Relative', 
        'ak8PFCHSL3Absolute', 
        'ak8PFCHSResidual')
)


process.ak8PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak8PFCHSL1Offset', 
        'ak8PFCHSL2Relative', 
        'ak8PFCHSL3Absolute', 
        'ak8PFCHSResidual')
)


process.ak8PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak8PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak8PFCHSL2Relative', 
        'ak8PFCHSL3Absolute')
)


process.ak8PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak8PFCHSL2Relative', 
        'ak8PFCHSL3Absolute', 
        'ak8PFCHSResidual')
)


process.ak8PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L2Relative')
)


process.ak8PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L3Absolute')
)


process.ak8PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak8PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak8PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak8PFL1Fastjet', 
        'ak8PFL2Relative', 
        'ak8PFL3Absolute', 
        'ak8PFResidual')
)


process.ak8PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak8PFL1Offset', 
        'ak8PFL2Relative', 
        'ak8PFL3Absolute', 
        'ak8PFResidual')
)


process.ak8PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak8PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak8PFL2Relative', 
        'ak8PFL3Absolute')
)


process.ak8PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak8PFL2Relative', 
        'ak8PFL3Absolute', 
        'ak8PFResidual')
)


process.ak8PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L2Relative')
)


process.ak8PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L3Absolute')
)


process.ak8PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK8PF'),
    level = cms.string('L2L3Residual')
)


process.ak9PFCHSL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak9PFCHSL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak9PFCHSL1Fastjet', 
        'ak9PFCHSL2Relative', 
        'ak9PFCHSL3Absolute', 
        'ak9PFCHSResidual')
)


process.ak9PFCHSL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak9PFCHSL1Offset', 
        'ak9PFCHSL2Relative', 
        'ak9PFCHSL3Absolute', 
        'ak9PFCHSResidual')
)


process.ak9PFCHSL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak9PFCHSL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak9PFCHSL2Relative', 
        'ak9PFCHSL3Absolute')
)


process.ak9PFCHSL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak9PFCHSL2Relative', 
        'ak9PFCHSL3Absolute', 
        'ak9PFCHSResidual')
)


process.ak9PFCHSL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L2Relative')
)


process.ak9PFCHSL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L3Absolute')
)


process.ak9PFCHSResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PFchs'),
    level = cms.string('L2L3Residual')
)


process.ak9PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ak9PFL1FastjetL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak9PFL1Fastjet', 
        'ak9PFL2Relative', 
        'ak9PFL3Absolute', 
        'ak9PFResidual')
)


process.ak9PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak9PFL1Offset', 
        'ak9PFL2Relative', 
        'ak9PFL3Absolute', 
        'ak9PFResidual')
)


process.ak9PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ak9PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak9PFL2Relative', 
        'ak9PFL3Absolute')
)


process.ak9PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak9PFL2Relative', 
        'ak9PFL3Absolute', 
        'ak9PFResidual')
)


process.ak9PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L2Relative')
)


process.ak9PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L3Absolute')
)


process.ak9PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('AK9PF'),
    level = cms.string('L2L3Residual')
)


process.fakeForIdealAlignment = cms.ESProducer("FakeAlignmentProducer",
    appendToDataLabel = cms.string('fakeForIdeal')
)


process.hcalDDDRecConstants = cms.ESProducer("HcalDDDRecConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalDDDSimConstants = cms.ESProducer("HcalDDDSimConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalTopologyIdeal = cms.ESProducer("HcalTopologyIdealEP",
    Exclude = cms.untracked.string(''),
    MergePosition = cms.untracked.bool(True),
    appendToDataLabel = cms.string('')
)


process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
    dump = cms.untracked.vstring(''),
    file = cms.untracked.string('')
)


process.ic5CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Fastjet', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute')
)


process.ic5CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5CaloL1Offset', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloL6SLB')
)


process.ic5CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5CaloL1Fastjet', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloResidual')
)


process.ic5CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.ic5CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5CaloL1Offset', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute')
)


process.ic5CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5CaloL1Offset', 
        'ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloResidual')
)


process.ic5CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ic5CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5CaloL2Relative', 
        'ic5CaloL3Absolute')
)


process.ic5CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloL6SLB')
)


process.ic5CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5CaloL2Relative', 
        'ic5CaloL3Absolute', 
        'ic5CaloResidual')
)


process.ic5CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L2Relative')
)


process.ic5CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L3Absolute')
)


process.ic5CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ic5CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ic5CaloJetsSoftMuonTagInfos")
)


process.ic5CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5Calo'),
    level = cms.string('L2L3Residual')
)


process.ic5PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute')
)


process.ic5PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFL6SLB')
)


process.ic5PFL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5PFL1Fastjet', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFResidual')
)


process.ic5PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.ic5PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5PFL1Offset', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute')
)


process.ic5PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5PFL1Offset', 
        'ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFResidual')
)


process.ic5PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.ic5PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5PFL2Relative', 
        'ic5PFL3Absolute')
)


process.ic5PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFL6SLB')
)


process.ic5PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ic5PFL2Relative', 
        'ic5PFL3Absolute', 
        'ic5PFResidual')
)


process.ic5PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L2Relative')
)


process.ic5PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L3Absolute')
)


process.ic5PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("ic5PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("ic5PFJetsSoftMuonTagInfos")
)


process.ic5PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('IC5PF'),
    level = cms.string('L2L3Residual')
)


process.idealForDigiCSCGeometry = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    debugV = cms.untracked.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useDDD = cms.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.idealForDigiDTGeometry = cms.ESProducer("DTGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.idealForDigiTrackerGeometry = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.kt4CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Fastjet', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute')
)


process.kt4CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4CaloL1Offset', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloL6SLB')
)


process.kt4CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4CaloL1Fastjet', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloResidual')
)


process.kt4CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.kt4CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4CaloL1Offset', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute')
)


process.kt4CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4CaloL1Offset', 
        'kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloResidual')
)


process.kt4CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.kt4CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4CaloL2Relative', 
        'kt4CaloL3Absolute')
)


process.kt4CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloL6SLB')
)


process.kt4CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4CaloL2Relative', 
        'kt4CaloL3Absolute', 
        'kt4CaloResidual')
)


process.kt4CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L2Relative')
)


process.kt4CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L3Absolute')
)


process.kt4CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("kt4CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("kt4CaloJetsSoftMuonTagInfos")
)


process.kt4CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4Calo'),
    level = cms.string('L2L3Residual')
)


process.kt4PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute')
)


process.kt4PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFL6SLB')
)


process.kt4PFL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4PFL1Fastjet', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFResidual')
)


process.kt4PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.kt4PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4PFL1Offset', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute')
)


process.kt4PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4PFL1Offset', 
        'kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFResidual')
)


process.kt4PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.kt4PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4PFL2Relative', 
        'kt4PFL3Absolute')
)


process.kt4PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFL6SLB')
)


process.kt4PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt4PFL2Relative', 
        'kt4PFL3Absolute', 
        'kt4PFResidual')
)


process.kt4PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L2Relative')
)


process.kt4PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L3Absolute')
)


process.kt4PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("kt4PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("kt4PFJetsSoftMuonTagInfos")
)


process.kt4PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT4PF'),
    level = cms.string('L2L3Residual')
)


process.kt6CaloL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4CaloL1Fastjet', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute')
)


process.kt6CaloL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6CaloL1Offset', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloL6SLB')
)


process.kt6CaloL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6CaloL1Fastjet', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloResidual')
)


process.kt6CaloL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAllCalo")
)


process.kt6CaloL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6CaloL1Offset', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute')
)


process.kt6CaloL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6CaloL1Offset', 
        'kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloResidual')
)


process.kt6CaloL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.kt6CaloL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6CaloL2Relative', 
        'kt6CaloL3Absolute')
)


process.kt6CaloL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloL6SLB')
)


process.kt6CaloL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6CaloL2Relative', 
        'kt6CaloL3Absolute', 
        'kt6CaloResidual')
)


process.kt6CaloL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L2Relative')
)


process.kt6CaloL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L3Absolute')
)


process.kt6CaloL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(True),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("kt6CaloJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("kt6CaloJetsSoftMuonTagInfos")
)


process.kt6CaloResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6Calo'),
    level = cms.string('L2L3Residual')
)


process.kt6PFL1FastL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute')
)


process.kt6PFL1FastL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('ak4PFL1Fastjet', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFL6SLB')
)


process.kt6PFL1FastL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6PFL1Fastjet', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFResidual')
)


process.kt6PFL1Fastjet = cms.ESProducer("L1FastjetCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L1FastJet'),
    srcRho = cms.InputTag("fixedGridRhoFastjetAll")
)


process.kt6PFL1L2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6PFL1Offset', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute')
)


process.kt6PFL1L2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6PFL1Offset', 
        'kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFResidual')
)


process.kt6PFL1Offset = cms.ESProducer("L1OffsetCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L1Offset'),
    minVtxNdof = cms.int32(4),
    vertexCollection = cms.string('offlinePrimaryVertices')
)


process.kt6PFL2L3 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6PFL2Relative', 
        'kt6PFL3Absolute')
)


process.kt6PFL2L3L6 = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFL6SLB')
)


process.kt6PFL2L3Residual = cms.ESProducer("JetCorrectionESChain",
    correctors = cms.vstring('kt6PFL2Relative', 
        'kt6PFL3Absolute', 
        'kt6PFResidual')
)


process.kt6PFL2Relative = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L2Relative')
)


process.kt6PFL3Absolute = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L3Absolute')
)


process.kt6PFL6SLB = cms.ESProducer("L6SLBCorrectionESProducer",
    addMuonToJet = cms.bool(False),
    algorithm = cms.string(''),
    level = cms.string('L6SLB'),
    srcBTagInfoElectron = cms.InputTag("kt6PFJetsSoftElectronTagInfos"),
    srcBTagInfoMuon = cms.InputTag("kt6PFJetsSoftMuonTagInfos")
)


process.kt6PFResidual = cms.ESProducer("LXXXCorrectionESProducer",
    algorithm = cms.string('KT6PF'),
    level = cms.string('L2L3Residual')
)


process.siPixelQualityESProducer = cms.ESProducer("SiPixelQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(cms.PSet(
        record = cms.string('SiPixelQualityFromDbRcd'),
        tag = cms.string('')
    ), 
        cms.PSet(
            record = cms.string('SiPixelDetVOffRcd'),
            tag = cms.string('')
        ))
)


process.siStripBackPlaneCorrectionDepESProducer = cms.ESProducer("SiStripBackPlaneCorrectionDepESProducer",
    BackPlaneCorrectionDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    BackPlaneCorrectionPeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    )
)


process.siStripGainESProducer = cms.ESProducer("SiStripGainESProducer",
    APVGain = cms.VPSet(cms.PSet(
        Label = cms.untracked.string(''),
        NormalizationFactor = cms.untracked.double(1.0),
        Record = cms.string('SiStripApvGainRcd')
    ), 
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGain2Rcd')
        )),
    AutomaticNormalization = cms.bool(False),
    appendToDataLabel = cms.string(''),
    printDebug = cms.untracked.bool(False)
)


process.siStripLorentzAngleDepESProducer = cms.ESProducer("SiStripLorentzAngleDepESProducer",
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    ),
    LorentzAngleDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripLorentzAngleRcd')
    ),
    LorentzAnglePeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripLorentzAngleRcd')
    )
)


process.siStripQualityESProducer = cms.ESProducer("SiStripQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(cms.PSet(
        record = cms.string('SiStripDetVOffRcd'),
        tag = cms.string('')
    ), 
        cms.PSet(
            record = cms.string('SiStripDetCablingRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('RunInfoRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadChannelRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadFiberRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadModuleRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadStripRcd'),
            tag = cms.string('')
        )),
    PrintDebugOutput = cms.bool(False),
    ReduceGranularity = cms.bool(False),
    ThresholdForReducedGranularity = cms.double(0.3),
    UseEmptyRunInfo = cms.bool(False),
    appendToDataLabel = cms.string('')
)


process.sistripconn = cms.ESProducer("SiStripConnectivity")


process.stripCPEESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('stripCPE'),
    ComponentType = cms.string('SimpleStripCPE'),
    parameters = cms.PSet(

    )
)


process.trackerGeometryDB = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.trackerNumberingGeometryDB = cms.ESProducer("TrackerGeometricDetESModule",
    appendToDataLabel = cms.string(''),
    fromDDD = cms.bool(False)
)


process.trackerTopology = cms.ESProducer("TrackerTopologyEP",
    appendToDataLabel = cms.string('')
)


process.GlobalTag = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    DumpStat = cms.untracked.bool(False),
    ReconnectEachRun = cms.untracked.bool(False),
    RefreshAlways = cms.untracked.bool(False),
    RefreshEachRun = cms.untracked.bool(False),
    RefreshOpenIOVs = cms.untracked.bool(False),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    globaltag = cms.string('94X_mc2017_realistic_v14'),
    pfnPostfix = cms.untracked.string(''),
    pfnPrefix = cms.untracked.string(''),
    snapshotTime = cms.string(''),
    toGet = cms.VPSet()
)


process.HepPDTESSource = cms.ESSource("HepPDTESSource",
    pdtFileName = cms.FileInPath('SimGeneral/HepPDTESSource/data/pythiaparticle.tbl')
)


process.eegeom = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('EcalMappingRcd')
)


process.es_hardcode = cms.ESSource("HcalHardcodeCalibrations",
    GainWidthsForTrigPrims = cms.bool(False),
    HBRecalibration = cms.bool(False),
    HBmeanenergies = cms.FileInPath('CalibCalorimetry/HcalPlugins/data/meanenergiesHB.txt'),
    HBreCalibCutoff = cms.double(20.0),
    HERecalibration = cms.bool(False),
    HEmeanenergies = cms.FileInPath('CalibCalorimetry/HcalPlugins/data/meanenergiesHE.txt'),
    HEreCalibCutoff = cms.double(20.0),
    HFRecalParameterBlock = cms.PSet(
        HFdepthOneParameterA = cms.vdouble(0.004123, 0.00602, 0.008201, 0.010489, 0.013379, 
            0.016997, 0.021464, 0.027371, 0.034195, 0.044807, 
            0.058939, 0.125497),
        HFdepthOneParameterB = cms.vdouble(-4e-06, -2e-06, 0.0, 4e-06, 1.5e-05, 
            2.6e-05, 6.3e-05, 8.4e-05, 0.00016, 0.000107, 
            0.000425, 0.000209),
        HFdepthTwoParameterA = cms.vdouble(0.002861, 0.004168, 0.0064, 0.008388, 0.011601, 
            0.014425, 0.018633, 0.023232, 0.028274, 0.035447, 
            0.051579, 0.086593),
        HFdepthTwoParameterB = cms.vdouble(-2e-06, -0.0, -7e-06, -6e-06, -2e-06, 
            1e-06, 1.9e-05, 3.1e-05, 6.7e-05, 1.2e-05, 
            0.000157, -3e-06)
    ),
    HFRecalibration = cms.bool(False),
    SiPMCharacteristics = cms.VPSet(cms.PSet(
        crosstalk = cms.double(0.0),
        nonlin1 = cms.double(1.0),
        nonlin2 = cms.double(0.0),
        nonlin3 = cms.double(0.0),
        pixels = cms.int32(36000)
    ), 
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(2500)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.17),
            nonlin1 = cms.double(1.00985),
            nonlin2 = cms.double(7.84089e-06),
            nonlin3 = cms.double(2.86282e-10),
            pixels = cms.int32(27370)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.196),
            nonlin1 = cms.double(1.00546),
            nonlin2 = cms.double(6.40239e-06),
            nonlin3 = cms.double(1.27011e-10),
            pixels = cms.int32(38018)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.17),
            nonlin1 = cms.double(1.00985),
            nonlin2 = cms.double(7.84089e-06),
            nonlin3 = cms.double(2.86282e-10),
            pixels = cms.int32(27370)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.196),
            nonlin1 = cms.double(1.00546),
            nonlin2 = cms.double(6.40239e-06),
            nonlin3 = cms.double(1.27011e-10),
            pixels = cms.int32(38018)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(0)
        )),
    hb = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.19),
        gainWidth = cms.vdouble(0.0),
        mcShape = cms.int32(125),
        pedestal = cms.double(3.285),
        pedestalWidth = cms.double(0.809),
        photoelectronsToAnalog = cms.double(0.3305),
        qieOffset = cms.vdouble(-0.49, 1.8, 7.2, 37.9),
        qieSlope = cms.vdouble(0.912, 0.917, 0.922, 0.923),
        qieType = cms.int32(0),
        recoShape = cms.int32(105)
    ),
    hbUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.01, 0.015),
        doRadiationDamage = cms.bool(True),
        gain = cms.vdouble(0.000439367311072),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(203),
        pedestal = cms.double(17.3),
        pedestalWidth = cms.double(1.5),
        photoelectronsToAnalog = cms.double(44.0),
        qieOffset = cms.vdouble(0.0, 0.0, 0.0, 0.0),
        qieSlope = cms.vdouble(0.05376, 0.05376, 0.05376, 0.05376),
        qieType = cms.int32(2),
        radiationDamage = cms.PSet(
            depVsNeutrons = cms.vdouble(5.69e-11, 7.9e-11),
            depVsTemp = cms.double(0.0631),
            intlumiOffset = cms.double(150),
            intlumiToNeutrons = cms.double(367000000.0),
            temperatureBase = cms.double(20),
            temperatureNew = cms.double(-5)
        ),
        recoShape = cms.int32(203)
    ),
    he = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.23),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(125),
        pedestal = cms.double(3.163),
        pedestalWidth = cms.double(0.9698),
        photoelectronsToAnalog = cms.double(0.3305),
        qieOffset = cms.vdouble(-0.38, 2.0, 7.6, 39.6),
        qieSlope = cms.vdouble(0.912, 0.916, 0.92, 0.922),
        qieType = cms.int32(0),
        recoShape = cms.int32(105)
    ),
    heUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.01, 0.015),
        doRadiationDamage = cms.bool(True),
        gain = cms.vdouble(0.000439367311072),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(203),
        pedestal = cms.double(17.3),
        pedestalWidth = cms.double(1.5),
        photoelectronsToAnalog = cms.double(44.0),
        qieOffset = cms.vdouble(0.0, 0.0, 0.0, 0.0),
        qieSlope = cms.vdouble(0.05376, 0.05376, 0.05376, 0.05376),
        qieType = cms.int32(2),
        radiationDamage = cms.PSet(
            depVsNeutrons = cms.vdouble(5.69e-11, 7.9e-11),
            depVsTemp = cms.double(0.0631),
            intlumiOffset = cms.double(75),
            intlumiToNeutrons = cms.double(29200000.0),
            temperatureBase = cms.double(20),
            temperatureNew = cms.double(5)
        ),
        recoShape = cms.int32(203)
    ),
    hf = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.14, 0.135),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(301),
        pedestal = cms.double(9.354),
        pedestalWidth = cms.double(2.516),
        photoelectronsToAnalog = cms.double(0.0),
        qieOffset = cms.vdouble(-0.87, 1.4, 7.8, -29.6),
        qieSlope = cms.vdouble(0.359, 0.358, 0.36, 0.367),
        qieType = cms.int32(0),
        recoShape = cms.int32(301)
    ),
    hfUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.14, 0.135),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(301),
        pedestal = cms.double(13.33),
        pedestalWidth = cms.double(3.33),
        photoelectronsToAnalog = cms.double(0.0),
        qieOffset = cms.vdouble(0.0697, -0.7405, 12.38, -671.9),
        qieSlope = cms.vdouble(0.297, 0.298, 0.298, 0.313),
        qieType = cms.int32(1),
        recoShape = cms.int32(301)
    ),
    ho = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.006, 0.0087),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(201),
        pedestal = cms.double(12.06),
        pedestalWidth = cms.double(0.6285),
        photoelectronsToAnalog = cms.double(4.0),
        qieOffset = cms.vdouble(-0.44, 1.4, 7.1, 38.5),
        qieSlope = cms.vdouble(0.907, 0.915, 0.92, 0.921),
        qieType = cms.int32(0),
        recoShape = cms.int32(201)
    ),
    iLumi = cms.double(-1.0),
    killHE = cms.bool(False),
    testHEPlan1 = cms.bool(False),
    testHFQIE10 = cms.bool(False),
    toGet = cms.untracked.vstring('GainWidths'),
    useHBUpgrade = cms.bool(False),
    useHEUpgrade = cms.bool(False),
    useHFUpgrade = cms.bool(False),
    useHOUpgrade = cms.bool(True),
    useLayer0Weight = cms.bool(False)
)


process.l1tPS = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    toGet = cms.VPSet(cms.PSet(
        record = cms.string('L1TGlobalPrescalesVetosRcd'),
        tag = cms.string('L1TGlobalPrescalesVetos_passThrough_mc')
    ))
)


process.prefer("es_hardcode")

process.prefer("l1tPS")

process.type0PFMEtCorrectionPFCandToVertexAssociationTask = cms.Task(process.particleFlowDisplacedVertex, process.pfCandidateToVertexAssociation, process.selectedPrimaryVertexHighestPtTrackSumForPFMEtCorrType0, process.selectedVerticesForPFMEtCorrType0)


process.ak4CaloL2L3ResidualCorrectorTask = cms.Task(process.ak4CaloL2L3ResidualCorrector, process.ak4CaloL2RelativeCorrector, process.ak4CaloL3AbsoluteCorrector, process.ak4CaloResidualCorrector)


process.ak4PFCHSL2L3ResidualCorrectorTask = cms.Task(process.ak4PFCHSL2L3ResidualCorrector, process.ak4PFCHSL2RelativeCorrector, process.ak4PFCHSL3AbsoluteCorrector, process.ak4PFCHSResidualCorrector)


process.ak4CaloL1L2L3ResidualCorrectorTask = cms.Task(process.ak4CaloL1L2L3ResidualCorrector, process.ak4CaloL1OffsetCorrector, process.ak4CaloL2RelativeCorrector, process.ak4CaloL3AbsoluteCorrector, process.ak4CaloResidualCorrector)


process.ak4PFPuppiL1L2L3ResidualCorrectorTask = cms.Task(process.ak4PFPuppiL1L2L3ResidualCorrector, process.ak4PFPuppiL1OffsetCorrector, process.ak4PFPuppiL2RelativeCorrector, process.ak4PFPuppiL3AbsoluteCorrector, process.ak4PFPuppiResidualCorrector)


process.ak4L1JPTOffsetCorrectorTask = cms.Task(process.ak4CaloL1OffsetCorrector, process.ak4L1JPTOffsetCorrector)


process.ak4JPTL1L2L3ResidualCorrectorTask = cms.Task(process.ak4JPTL1L2L3ResidualCorrector, process.ak4JPTL2RelativeCorrector, process.ak4JPTL3AbsoluteCorrector, process.ak4JPTResidualCorrector, process.ak4L1JPTOffsetCorrectorTask)


process.patPFMetTxyCorrTask = cms.Task(process.patPFMetTxyCorr)


process.producePatPFMETCorrectionsTask = cms.Task(process.patPFMet, process.patPFMetT0Corr, process.patPFMetT0pcT1, process.patPFMetT0pcT1T2, process.patPFMetT1, process.patPFMetT1T2, process.patPFMetT1T2Corr, process.patPFMetT2Corr, process.pfCandMETcorr, process.pfCandsNotInJetsForMetCorr, process.selectedPatJetsForMetT1T2Corr, process.selectedPatJetsForMetT2Corr, process.type0PFMEtCorrectionPFCandToVertexAssociationTask)


process.ak4PFL1L2L3CorrectorTask = cms.Task(process.ak4PFL1L2L3Corrector, process.ak4PFL1OffsetCorrector, process.ak4PFL2RelativeCorrector, process.ak4PFL3AbsoluteCorrector)


process.ak4CaloL1L2L3CorrectorTask = cms.Task(process.ak4CaloL1L2L3Corrector, process.ak4CaloL1OffsetCorrector, process.ak4CaloL2RelativeCorrector, process.ak4CaloL3AbsoluteCorrector)


process.ak4CaloL1FastL2L3CorrectorTask = cms.Task(process.ak4CaloL1FastL2L3Corrector, process.ak4CaloL1FastjetCorrector, process.ak4CaloL2RelativeCorrector, process.ak4CaloL3AbsoluteCorrector)


process.ak4PFPuppiL1FastL2L3ResidualCorrectorTask = cms.Task(process.ak4PFPuppiL1FastL2L3ResidualCorrector, process.ak4PFPuppiL1FastjetCorrector, process.ak4PFPuppiL2RelativeCorrector, process.ak4PFPuppiL3AbsoluteCorrector, process.ak4PFPuppiResidualCorrector)


process.ak4CaloL2L3L6CorrectorTask = cms.Task(process.ak4CaloL2L3L6Corrector, process.ak4CaloL2RelativeCorrector, process.ak4CaloL3AbsoluteCorrector, process.ak4CaloL6SLBCorrector)


process.ak4L1JPTFastjetCorrectorTask = cms.Task(process.ak4CaloL1FastjetCorrector, process.ak4L1JPTFastjetCorrector)


process.patPFMetT2SmearCorrTask = cms.Task(process.patPFMetT1T2SmearCorr, process.patPFMetT2SmearCorr, process.patSmearedJets, process.selectedPatJetsForMetT1T2SmearCorr, process.selectedPatJetsForMetT2SmearCorr)


process.ak4PFL1FastL2L3CorrectorTask = cms.Task(process.ak4PFL1FastL2L3Corrector, process.ak4PFL1FastjetCorrector, process.ak4PFL2RelativeCorrector, process.ak4PFL3AbsoluteCorrector)


process.ak4PFL2L3L6CorrectorTask = cms.Task(process.ak4PFL2L3L6Corrector, process.ak4PFL2RelativeCorrector, process.ak4PFL3AbsoluteCorrector, process.ak4PFL6SLBCorrector)


process.patPFMetT2CorrTask = cms.Task(process.patPFMetT2Corr)


process.producePatPFMETCorrectionsUncTask = cms.Task(process.patPFMet, process.patPFMetT0Corr, process.patPFMetT1T2Corr, process.patPFMetT2Corr, process.pfCandMETcorr, process.pfCandsNotInJetsForMetCorr, process.selectedPatJetsForMetT1T2Corr, process.selectedPatJetsForMetT2Corr, process.type0PFMEtCorrectionPFCandToVertexAssociationTask)


process.ak4JPTL2L3ResidualCorrectorTask = cms.Task(process.ak4JPTL2L3ResidualCorrector, process.ak4JPTL2RelativeCorrector, process.ak4JPTL3AbsoluteCorrector, process.ak4JPTResidualCorrector, process.ak4L1JPTOffsetCorrectorTask)


process.ak4JPTL2L3CorrectorTask = cms.Task(process.ak4JPTL2L3Corrector, process.ak4JPTL2RelativeCorrector, process.ak4JPTL3AbsoluteCorrector, process.ak4L1JPTOffsetCorrectorTask)


process.ak4PFPuppiL1FastL2L3CorrectorTask = cms.Task(process.ak4PFPuppiL1FastL2L3Corrector, process.ak4PFPuppiL1FastjetCorrector, process.ak4PFPuppiL2RelativeCorrector, process.ak4PFPuppiL3AbsoluteCorrector)


process.ak4PFPuppiL1L2L3CorrectorTask = cms.Task(process.ak4PFPuppiL1L2L3Corrector, process.ak4PFPuppiL1OffsetCorrector, process.ak4PFPuppiL2RelativeCorrector, process.ak4PFPuppiL3AbsoluteCorrector)


process.ak4TrackL2L3CorrectorTask = cms.Task(process.ak4TrackL2L3Corrector, process.ak4TrackL2RelativeCorrector, process.ak4TrackL3AbsoluteCorrector)


process.patPFMetT0CorrTask = cms.Task(process.patPFMetT0Corr, process.type0PFMEtCorrectionPFCandToVertexAssociationTask)


process.type0PFMEtCorrectionTask = cms.Task(process.pfMETcorrType0, process.type0PFMEtCorrectionPFCandToVertexAssociationTask)


process.ak4PFL1FastL2L3ResidualCorrectorTask = cms.Task(process.ak4PFL1FastL2L3ResidualCorrector, process.ak4PFL1FastjetCorrector, process.ak4PFL2RelativeCorrector, process.ak4PFL3AbsoluteCorrector, process.ak4PFResidualCorrector)


process.patPFMetSmearCorrTask = cms.Task(process.patPFMetT1T2SmearCorr, process.patSmearedJets, process.selectedPatJetsForMetT1T2SmearCorr)


process.ak4CaloL1FastL2L3ResidualCorrectorTask = cms.Task(process.ak4CaloL1FastL2L3ResidualCorrector, process.ak4CaloL1FastjetCorrector, process.ak4CaloL2RelativeCorrector, process.ak4CaloL3AbsoluteCorrector, process.ak4CaloResidualCorrector)


process.ak4PFCHSL1FastL2L3CorrectorTask = cms.Task(process.ak4PFCHSL1FastL2L3Corrector, process.ak4PFCHSL1FastjetCorrector, process.ak4PFCHSL2RelativeCorrector, process.ak4PFCHSL3AbsoluteCorrector)


process.ak4PFL1L2L3ResidualCorrectorTask = cms.Task(process.ak4PFL1L2L3ResidualCorrector, process.ak4PFL1OffsetCorrector, process.ak4PFL2RelativeCorrector, process.ak4PFL3AbsoluteCorrector, process.ak4PFResidualCorrector)


process.ak4PFL2L3ResidualCorrectorTask = cms.Task(process.ak4PFL2L3ResidualCorrector, process.ak4PFL2RelativeCorrector, process.ak4PFL3AbsoluteCorrector, process.ak4PFResidualCorrector)


process.ak4JPTL1L2L3CorrectorTask = cms.Task(process.ak4JPTL1L2L3Corrector, process.ak4JPTL2RelativeCorrector, process.ak4JPTL3AbsoluteCorrector, process.ak4L1JPTOffsetCorrectorTask)


process.ak4PFCHSL2L3CorrectorTask = cms.Task(process.ak4PFCHSL2L3Corrector, process.ak4PFCHSL2RelativeCorrector, process.ak4PFCHSL3AbsoluteCorrector)


process.ak4CaloL2L3CorrectorTask = cms.Task(process.ak4CaloL2L3Corrector, process.ak4CaloL2RelativeCorrector, process.ak4CaloL3AbsoluteCorrector)


process.patPFMetT1T2CorrTask = cms.Task(process.patPFMetT1T2Corr, process.selectedPatJetsForMetT1T2Corr)


process.ak4PFL1FastL2L3L6CorrectorTask = cms.Task(process.ak4PFL1FastL2L3L6Corrector, process.ak4PFL1FastjetCorrector, process.ak4PFL2RelativeCorrector, process.ak4PFL3AbsoluteCorrector, process.ak4PFL6SLBCorrector)


process.ak4PFCHSL1L2L3CorrectorTask = cms.Task(process.ak4PFCHSL1L2L3Corrector, process.ak4PFCHSL1OffsetCorrector, process.ak4PFCHSL2RelativeCorrector, process.ak4PFCHSL3AbsoluteCorrector)


process.ak4PFCHSL1L2L3ResidualCorrectorTask = cms.Task(process.ak4PFCHSL1L2L3ResidualCorrector, process.ak4PFCHSL1OffsetCorrector, process.ak4PFCHSL2RelativeCorrector, process.ak4PFCHSL3AbsoluteCorrector, process.ak4PFCHSResidualCorrector)


process.ak4PFPuppiL2L3ResidualCorrectorTask = cms.Task(process.ak4PFPuppiL2L3ResidualCorrector, process.ak4PFPuppiL2RelativeCorrector, process.ak4PFPuppiL3AbsoluteCorrector, process.ak4PFPuppiResidualCorrector)


process.ak4PFPuppiL2L3CorrectorTask = cms.Task(process.ak4PFPuppiL2L3Corrector, process.ak4PFPuppiL2RelativeCorrector, process.ak4PFPuppiL3AbsoluteCorrector)


process.ak4JPTL1FastL2L3CorrectorTask = cms.Task(process.ak4JPTL1FastL2L3Corrector, process.ak4JPTL2RelativeCorrector, process.ak4JPTL3AbsoluteCorrector, process.ak4L1JPTFastjetCorrectorTask)


process.ak4PFCHSL1FastL2L3ResidualCorrectorTask = cms.Task(process.ak4PFCHSL1FastL2L3ResidualCorrector, process.ak4PFCHSL1FastjetCorrector, process.ak4PFCHSL2RelativeCorrector, process.ak4PFCHSL3AbsoluteCorrector, process.ak4PFCHSResidualCorrector)


process.correctionTermsPfMetType1Type2Task = cms.Task(process.ak4PFCHSL1FastL2L3CorrectorTask, process.ak4PFCHSL1FastL2L3ResidualCorrectorTask, process.corrPfMetType1, process.corrPfMetType2, process.particleFlowPtrs, process.pfCandMETcorr, process.pfCandsNotInJetsForMetCorr, process.pfCandsNotInJetsPtrForMetCorr, process.pfJetsPtrForMetCorr)


process.ak4PFL2L3CorrectorTask = cms.Task(process.ak4PFL2L3Corrector, process.ak4PFL2RelativeCorrector, process.ak4PFL3AbsoluteCorrector)


process.ak4CaloL1FastL2L3L6CorrectorTask = cms.Task(process.ak4CaloL1FastL2L3L6Corrector, process.ak4CaloL1FastjetCorrector, process.ak4CaloL2RelativeCorrector, process.ak4CaloL3AbsoluteCorrector, process.ak4CaloL6SLBCorrector)


process.ak4JPTL1FastL2L3ResidualCorrectorTask = cms.Task(process.ak4JPTL1FastL2L3ResidualCorrector, process.ak4JPTL2RelativeCorrector, process.ak4JPTL3AbsoluteCorrector, process.ak4JPTResidualCorrector, process.ak4L1JPTFastjetCorrectorTask)


process.jetCorrectorsTask = cms.Task(process.ak4CaloL1FastL2L3CorrectorTask, process.ak4CaloL1FastL2L3L6CorrectorTask, process.ak4CaloL1FastL2L3ResidualCorrectorTask, process.ak4CaloL1L2L3CorrectorTask, process.ak4CaloL1L2L3ResidualCorrectorTask, process.ak4CaloL2L3CorrectorTask, process.ak4CaloL2L3L6CorrectorTask, process.ak4CaloL2L3ResidualCorrectorTask, process.ak4JPTL1FastL2L3CorrectorTask, process.ak4JPTL1FastL2L3ResidualCorrectorTask, process.ak4JPTL1L2L3CorrectorTask, process.ak4JPTL1L2L3ResidualCorrectorTask, process.ak4JPTL2L3CorrectorTask, process.ak4JPTL2L3ResidualCorrectorTask, process.ak4L1JPTFastjetCorrectorTask, process.ak4L1JPTOffsetCorrectorTask, process.ak4PFCHSL1FastL2L3CorrectorTask, process.ak4PFCHSL1FastL2L3ResidualCorrectorTask, process.ak4PFCHSL1L2L3CorrectorTask, process.ak4PFCHSL1L2L3ResidualCorrectorTask, process.ak4PFCHSL2L3CorrectorTask, process.ak4PFCHSL2L3ResidualCorrectorTask, process.ak4PFL1FastL2L3CorrectorTask, process.ak4PFL1FastL2L3L6CorrectorTask, process.ak4PFL1FastL2L3ResidualCorrectorTask, process.ak4PFL1L2L3CorrectorTask, process.ak4PFL1L2L3ResidualCorrectorTask, process.ak4PFL2L3CorrectorTask, process.ak4PFL2L3L6CorrectorTask, process.ak4PFL2L3ResidualCorrectorTask, process.ak4PFPuppiL1FastL2L3CorrectorTask, process.ak4PFPuppiL1FastL2L3ResidualCorrectorTask, process.ak4PFPuppiL1L2L3CorrectorTask, process.ak4PFPuppiL1L2L3ResidualCorrectorTask, process.ak4PFPuppiL2L3CorrectorTask, process.ak4PFPuppiL2L3ResidualCorrectorTask, process.ak4TrackL2L3CorrectorTask)


process.patAlgosToolsTask = cms.Task(process.basicJetsForMet, process.cleanedPatJets, process.genMetExtractor, process.jetCorrectorsTask, process.jetSelectorForMet, process.metrawCalo, process.patCHSMet, process.patCaloMet, process.patJetCorrFactorsReapplyJEC, process.patJetsReapplyJEC, process.patPFMetT1, process.patPFMetT1ElectronEnDown, process.patPFMetT1ElectronEnUp, process.patPFMetT1JetEnDown, process.patPFMetT1JetEnUp, process.patPFMetT1JetResDown, process.patPFMetT1JetResUp, process.patPFMetT1MuonEnDown, process.patPFMetT1MuonEnUp, process.patPFMetT1PhotonEnDown, process.patPFMetT1PhotonEnUp, process.patPFMetT1Smear, process.patPFMetT1SmearElectronEnDown, process.patPFMetT1SmearElectronEnUp, process.patPFMetT1SmearJetEnDown, process.patPFMetT1SmearJetEnUp, process.patPFMetT1SmearJetResDown, process.patPFMetT1SmearJetResUp, process.patPFMetT1SmearMuonEnDown, process.patPFMetT1SmearMuonEnUp, process.patPFMetT1SmearPhotonEnDown, process.patPFMetT1SmearPhotonEnUp, process.patPFMetT1SmearTauEnDown, process.patPFMetT1SmearTauEnUp, process.patPFMetT1SmearUnclusteredEnDown, process.patPFMetT1SmearUnclusteredEnUp, process.patPFMetT1TauEnDown, process.patPFMetT1TauEnUp, process.patPFMetT1Txy, process.patPFMetT1UnclusteredEnDown, process.patPFMetT1UnclusteredEnUp, process.patPFMetT2SmearCorrTask, process.patPFMetTxy, process.patPFMetTxyCorrTask, process.patTrkMet, process.pfCHS, process.pfCandsForUnclusteredUnc, process.pfCandsNoJets, process.pfCandsNoJetsNoEle, process.pfCandsNoJetsNoEleNoMu, process.pfCandsNoJetsNoEleNoMuNoTau, process.pfMet, process.pfMetCHS, process.pfMetTrk, process.pfTrk, process.producePatPFMETCorrectionsTask, process.shiftedPatElectronEnDown, process.shiftedPatElectronEnUp, process.shiftedPatJetEnDown, process.shiftedPatJetEnUp, process.shiftedPatJetResDown, process.shiftedPatJetResUp, process.shiftedPatMETCorrElectronEnDown, process.shiftedPatMETCorrElectronEnUp, process.shiftedPatMETCorrJetEnDown, process.shiftedPatMETCorrJetEnUp, process.shiftedPatMETCorrJetResDown, process.shiftedPatMETCorrJetResUp, process.shiftedPatMETCorrMuonEnDown, process.shiftedPatMETCorrMuonEnUp, process.shiftedPatMETCorrPhotonEnDown, process.shiftedPatMETCorrPhotonEnUp, process.shiftedPatMETCorrSmearedJetResDown, process.shiftedPatMETCorrSmearedJetResUp, process.shiftedPatMETCorrTauEnDown, process.shiftedPatMETCorrTauEnUp, process.shiftedPatMETCorrUnclusteredEnDown, process.shiftedPatMETCorrUnclusteredEnUp, process.shiftedPatMuonEnDown, process.shiftedPatMuonEnUp, process.shiftedPatPhotonEnDown, process.shiftedPatPhotonEnUp, process.shiftedPatSmearedJetResDown, process.shiftedPatSmearedJetResUp, process.shiftedPatTauEnDown, process.shiftedPatTauEnUp, process.shiftedPatUnclusteredEnDown, process.shiftedPatUnclusteredEnUp, process.slimmedMETs)


process.ak4JPTL2L3ResidualCorrectorChain = cms.Sequence(process.ak4JPTL2L3ResidualCorrectorTask)


process.type0PFMEtCorrectionPFCandToVertexAssociationForValidation = cms.Sequence(process.type0PFMEtCorrectionPFCandToVertexAssociationTask)


process.ak4PFL1FastL2L3L6CorrectorChain = cms.Sequence(process.ak4PFL1FastL2L3L6CorrectorTask)


process.miniAODrhoSequence = cms.Sequence(process.fixedGridRhoAllMaker+process.fixedGridRhoFastJetAllMaker+process.fixedGridRhoFastJetAllCaloMaker+process.fixedGridRhoFastJetCentralCaloMaker+process.fixedGridRhoFastJetCentralChargedPileUpMaker+process.fixedGridRhoFastJetCentralNeutralMaker+process.fixedGridRhoFastJetAllCentralMaker)


process.ak4L1JPTOffsetCorrectorChain = cms.Sequence(process.ak4L1JPTOffsetCorrectorTask)


process.ak4TrackL2L3CorrectorChain = cms.Sequence(process.ak4TrackL2L3CorrectorTask)


process.patPFMetT1T2CorrSequence = cms.Sequence(cms.Task(process.patPFMetT1T2Corr))


process.ak4PFPuppiL2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFPuppiL2L3ResidualCorrectorTask)


process.ak4L1JPTFastjetCorrectorChain = cms.Sequence(process.ak4L1JPTFastjetCorrectorTask)


process.hltMakerSequence = cms.Sequence(process.hltMaker)


process.ak4PFL2L3L6CorrectorChain = cms.Sequence(process.ak4PFL2L3L6CorrectorTask)


process.ak4PFL1L2L3CorrectorChain = cms.Sequence(process.ak4PFL1L2L3CorrectorTask)


process.ak4PFCHSL1FastL2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFCHSL1FastL2L3ResidualCorrectorTask)


process.type0PFMEtCorrectionPFCandToVertexAssociation = cms.Sequence(process.type0PFMEtCorrectionPFCandToVertexAssociationTask)


process.ak4PFCHSL1FastL2L3CorrectorChain = cms.Sequence(process.ak4PFCHSL1FastL2L3CorrectorTask)


process.ak4CaloL1L2L3CorrectorChain = cms.Sequence(process.ak4CaloL1L2L3CorrectorTask)


process.ak4PFCHSL1L2L3CorrectorChain = cms.Sequence(process.ak4PFCHSL1L2L3CorrectorTask)


process.ak4CaloL1FastL2L3ResidualCorrectorChain = cms.Sequence(process.ak4CaloL1FastL2L3ResidualCorrectorTask)


process.ak4PFL1FastL2L3CorrectorChain = cms.Sequence(process.ak4PFL1FastL2L3CorrectorTask)


process.ak4CaloL2L3L6CorrectorChain = cms.Sequence(process.ak4CaloL2L3L6CorrectorTask)


process.patMetModuleSequence = cms.Sequence(process.pfMet+process.genMetExtractor+process.patJetCorrFactorsReapplyJEC+process.patJetsReapplyJEC+process.basicJetsForMet+process.jetSelectorForMet+process.cleanedPatJets+process.metrawCalo+process.pfCHS+process.pfMetCHS+process.patCHSMet+process.pfTrk+process.pfMetTrk+process.patTrkMet+process.patPFMet)


process.patPFMetT2SmearCorrSequence = cms.Sequence(process.patPFMetT2SmearCorrTask)


process.ak4JPTL2L3CorrectorChain = cms.Sequence(process.ak4JPTL2L3CorrectorTask)


process.ak4CaloL2L3CorrectorChain = cms.Sequence(process.ak4CaloL2L3CorrectorTask)


process.ak4JPTL1L2L3ResidualCorrectorChain = cms.Sequence(process.ak4JPTL1L2L3ResidualCorrectorTask)


process.patPFMetTxyCorrSequence = cms.Sequence(process.patPFMetTxyCorrTask)


process.producePatPFMETCorrectionsUnc = cms.Sequence(process.producePatPFMETCorrectionsUncTask)


process.ak4PFPuppiL1L2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFPuppiL1L2L3ResidualCorrectorTask)


process.patPFMetT1SmearpatShiftedModuleSequence = cms.Sequence(process.patPFMetT1SmearJetResDown+process.patPFMetT1SmearJetResUp+process.patPFMetT1SmearMuonEnUp+process.patPFMetT1SmearMuonEnDown+process.patPFMetT1SmearJetEnUp+process.patPFMetT1SmearJetEnDown+process.patPFMetT1SmearTauEnUp+process.patPFMetT1SmearTauEnDown+process.patPFMetT1SmearPhotonEnUp+process.patPFMetT1SmearPhotonEnDown+process.patPFMetT1SmearElectronEnDown+process.patPFMetT1SmearElectronEnUp+process.patPFMetT1SmearUnclusteredEnUp+process.patPFMetT1SmearUnclusteredEnDown)


process.ak4PFPuppiL1FastL2L3CorrectorChain = cms.Sequence(process.ak4PFPuppiL1FastL2L3CorrectorTask)


process.correctionTermsPfMetType1Type2 = cms.Sequence(process.correctionTermsPfMetType1Type2Task)


process.ak4PFL1FastL2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFL1FastL2L3ResidualCorrectorTask)


process.producePatPFMETCorrections = cms.Sequence(process.producePatPFMETCorrectionsTask)


process.ak4JPTL1L2L3CorrectorChain = cms.Sequence(process.ak4JPTL1L2L3CorrectorTask)


process.ak4PFPuppiL1FastL2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFPuppiL1FastL2L3ResidualCorrectorTask)


process.type0PFMEtCorrection = cms.Sequence(process.type0PFMEtCorrectionTask)


process.ak4JPTL1FastL2L3ResidualCorrectorChain = cms.Sequence(process.ak4JPTL1FastL2L3ResidualCorrectorTask)


process.ak4PFCHSL2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFCHSL2L3ResidualCorrectorTask)


process.egmGsfElectronIDSequence = cms.Sequence(process.electronMVAValueMapProducer+process.egmGsfElectronIDs)


process.patPFMetT1SmearpatMetUncertaintySequence = cms.Sequence(process.shiftedPatSmearedJetResDown+process.shiftedPatMETCorrSmearedJetResDown+process.shiftedPatSmearedJetResUp+process.shiftedPatMETCorrSmearedJetResUp)


process.ak4CaloL1FastL2L3L6CorrectorChain = cms.Sequence(process.ak4CaloL1FastL2L3L6CorrectorTask)


process.type0PFMEtCorrectionPFCandToVertexAssociationForValidationMiniAOD = cms.Sequence(process.type0PFMEtCorrectionPFCandToVertexAssociationTask)


process.ak4PFL2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFL2L3ResidualCorrectorTask)


process.ak4PFCHSL2L3CorrectorChain = cms.Sequence(process.ak4PFCHSL2L3CorrectorTask)


process.patPFMetT0CorrSequence = cms.Sequence(process.patPFMetT0CorrTask)


process.patPFMetSmearCorrSequence = cms.Sequence(process.patPFMetSmearCorrTask)


process.ak4PFPuppiL1L2L3CorrectorChain = cms.Sequence(process.ak4PFPuppiL1L2L3CorrectorTask)


process.ak4CaloL2L3ResidualCorrectorChain = cms.Sequence(process.ak4CaloL2L3ResidualCorrectorTask)


process.ak4PFCHSL1L2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFCHSL1L2L3ResidualCorrectorTask)


process.ak4PFL2L3CorrectorChain = cms.Sequence(process.ak4PFL2L3CorrectorTask)


process.ak4CaloL1L2L3ResidualCorrectorChain = cms.Sequence(process.ak4CaloL1L2L3ResidualCorrectorTask)


process.ak4CaloL1FastL2L3CorrectorChain = cms.Sequence(process.ak4CaloL1FastL2L3CorrectorTask)


process.ak4PFPuppiL2L3CorrectorChain = cms.Sequence(process.ak4PFPuppiL2L3CorrectorTask)


process.patPFMetT2CorrSequence = cms.Sequence(process.patPFMetT2CorrTask)


process.ak4JPTL1FastL2L3CorrectorChain = cms.Sequence(process.ak4JPTL1FastL2L3CorrectorTask)


process.ak4PFL1L2L3ResidualCorrectorChain = cms.Sequence(process.ak4PFL1L2L3ResidualCorrectorTask)


process.patMetCorrectionSequence = cms.Sequence(process.patPFMetT1T2CorrSequence+process.patPFMetT1+process.patPFMetTxyCorrSequence+process.patPFMetT1Txy+process.patPFMetTxy+process.patPFMetSmearCorrSequence+process.patPFMetT1Smear)


process.patMetUncertaintySequence = cms.Sequence(process.ak4PFCHSL1FastL2L3CorrectorChain+process.ak4PFCHSL1FastL2L3ResidualCorrectorChain+(process.shiftedPatJetResDown+process.shiftedPatMETCorrJetResDown+process.shiftedPatJetResUp+process.shiftedPatMETCorrJetResUp+process.pfCandsNoJets+process.pfCandsNoJetsNoEle+process.pfCandsNoJetsNoEleNoMu+process.pfCandsNoJetsNoEleNoMuNoTau+process.pfCandsForUnclusteredUnc+process.shiftedPatMuonEnDown+process.shiftedPatMETCorrMuonEnDown+process.shiftedPatMuonEnUp+process.shiftedPatMETCorrMuonEnUp+process.shiftedPatJetEnDown+process.shiftedPatMETCorrJetEnDown+process.shiftedPatJetEnUp+process.shiftedPatMETCorrJetEnUp+process.shiftedPatTauEnDown+process.shiftedPatMETCorrTauEnDown+process.shiftedPatTauEnUp+process.shiftedPatMETCorrTauEnUp+process.shiftedPatPhotonEnDown+process.shiftedPatMETCorrPhotonEnDown+process.shiftedPatPhotonEnUp+process.shiftedPatMETCorrPhotonEnUp+process.shiftedPatElectronEnDown+process.shiftedPatMETCorrElectronEnDown+process.shiftedPatElectronEnUp+process.shiftedPatMETCorrElectronEnUp+process.shiftedPatUnclusteredEnDown+process.shiftedPatMETCorrUnclusteredEnDown+process.shiftedPatUnclusteredEnUp+process.shiftedPatMETCorrUnclusteredEnUp)+process.patPFMetT1SmearpatMetUncertaintySequence)


process.patShiftedModuleSequence = cms.Sequence(process.patPFMetT1JetResUp+process.patPFMetT1JetResDown+process.patPFMetT1MuonEnUp+process.patPFMetT1MuonEnDown+process.patPFMetT1JetEnUp+process.patPFMetT1JetEnDown+process.patPFMetT1TauEnDown+process.patPFMetT1TauEnUp+process.patPFMetT1PhotonEnUp+process.patPFMetT1PhotonEnDown+process.patPFMetT1ElectronEnUp+process.patPFMetT1ElectronEnDown+process.patPFMetT1UnclusteredEnUp+process.patPFMetT1UnclusteredEnDown+process.patPFMetT1SmearpatShiftedModuleSequence)


process.fullPatMetSequence = cms.Sequence(process.patMetModuleSequence+process.patMetCorrectionSequence+process.patMetUncertaintySequence+process.patShiftedModuleSequence+process.patCaloMet+process.slimmedMETs)


process.p = cms.Path(process.metFilterMaker+process.egmGsfElectronIDSequence+process.vertexMaker+process.secondaryVertexMaker+process.eventMaker+process.pfCandidateMaker+process.isoTrackMaker+process.isoForEle+process.isoForMu+process.electronMaker+process.muonMaker+process.pfJetMaker+process.pfJetPUPPIMaker+process.pfmetMaker+process.pfmetpuppiMaker+process.hltMakerSequence+process.miniAODrhoSequence+process.pftauMaker+process.photonMaker+process.genMaker+process.genJetMaker+process.candToGenAssMaker+process.pdfinfoMaker+process.puSummaryInfoMaker+process.hypDilepMaker)


process.outpath = cms.EndPath(process.out)


