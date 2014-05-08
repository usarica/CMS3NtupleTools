import FWCore.ParameterSet.Config as cms

slimcms2 = cms.untracked.vstring()

### -> pass 0

#PFCandidateMaker
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfcandsposAtEcalp4_CMS2*'))       # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfcandsisMuIso_CMS2*'))           # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfcandsecalE_CMS2*'))             # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfcandshcalE_CMS2*'))             # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfcandsrawEcalE_CMS2*'))          # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfcandsrawHcalE_CMS2*'))          # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfcandspS1E_CMS2*'))              # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfcandspS2E_CMS2*'))              # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfcandsdeltaP_CMS2*'))            # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfcandsmvaepi_CMS2*'))            # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfcandsmvaemu_CMS2*'))            # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfcandsmvapimu_CMS2*'))           # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfcandsmvanothinggamma_CMS2*'))   # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfcandsmvanothingnh_CMS2*'))      # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfcandsflag_CMS2*'))              # The miniAOD doesn't fill these anyway

#TrackMaker
slimcms2.extend(cms.untracked.vstring('drop *_trackMaker_trksvertexp4_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_trackMaker_trksouterp4_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_trackMaker_trksinnerposition_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_trackMaker_trksouterposition_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_trackMaker_trkslayer1layer_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_trackMaker_trkslayer1sizerphi_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_trackMaker_trkslayer1sizerz_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_trackMaker_trkslayer1charge_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_trackMaker_trkslayer1det_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_trackMaker_trksd0corrPhi_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_trackMaker_trksnlayers3D_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_trackMaker_trksnlayersLost_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_trackMaker_trkslostpixelhits_CMS2*'))

#L1 Maker
#slimcms2.extend(cms.untracked.vstring('drop *_l1Maker_*_CMS2*'))

#GenMaker
#slimcms2.extend(cms.untracked.vstring('drop *_genMaker_genps*_CMS2*'))

#CandToGenAssMaker ->drop all for jet and trk
#slimcms2.extend(cms.untracked.vstring('drop *_*_jetsmcdr_CMS2*'))                    # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_jetsmcidx_CMS2*'))                   # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_jetsmcemEnergy_CMS2*'))              # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_jetsmchadEnergy_CMS2*'))             # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_jetsmcinvEnergy_CMS2*'))             # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_jetsmcotherEnergy_CMS2*'))           # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_jetsmcp4_CMS2*'))                    # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_jetsmcgpdr_CMS2*'))                  # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_jetsmcgpidx_CMS2*'))                 # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_jetsmcgpp4_CMS2*'))                  # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_jetsmcid_CMS2*'))                    # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_jetsmcmotherid_CMS2*'))              # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_jetsmcmotherp4_CMS2*'))              # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_jetsmc3dr_CMS2*'))                   # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_jetsmc3idx_CMS2*'))                  # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_jetsmc3id_CMS2*'))                   # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_trkmcid_CMS2*'))                     # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_trkmcmotherid_CMS2*'))               # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_trkmcidx_CMS2*'))                    # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_trkmcp4_CMS2*'))                     # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_trkmcdr_CMS2*'))                     # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_trkmc3id_CMS2*'))                    # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_trkmc3motherid_CMS2*'))              # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_trkmc3idx_CMS2*'))                   # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_trkmc3motheridx_CMS2*'))             # The miniAOD doesn't fill these anyway
#slimcms2.extend(cms.untracked.vstring('drop *_*_trkmc3dr_CMS2*'))                    # The miniAOD doesn't fill these anyway

#CandToGenAssMaker ->keep mcp4, mcid, motherid for els, mus, photons, pfjets
#slimcms2.extend(cms.untracked.vstring('drop *_*_elsmcidx_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_elsmcmotherp4_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_elsmcdr_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_elsmc3id_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_elsmc3motherid_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_elsmc3idx_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_elsmc3motheridx_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_elsmc3dr_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_photonsmcidx_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_photonsmcmotherp4_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_photonsmcdr_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_photonsmc3id_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_photonsmc3motherid_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_photonsmc3idx_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_photonsmc3motheridx_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_photonsmc3dr_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_musmcidx_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_musmcmotherp4_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_musmcdr_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_musmc3id_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_musmc3motherid_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_musmc3idx_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_musmc3motheridx_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_musmc3dr_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfjetsmcidx_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfjetsmcgpdr_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfjetsmcgpidx_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfjetsmcgpp4_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfjetsmcemEnergy_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfjetsmchadEnergy_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfjetsmcinvEnergy_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfjetsmcotherEnergy_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfjetsmcmotherp4_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfjetsmc3dr_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfjetsmc3idx_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_*_pfjetsmc3id_CMS2*'))

### -> pass 1

#Jet Collections
slimcms2.extend(cms.untracked.vstring('drop *_jptMaker_*_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_bTagJPTJetMaker_*_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_jetMaker_*_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_bTagMaker_*_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_trkJetMaker_*_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_bTagTrkMaker_*_CMS2*'))

#bTagPFJetMaker
slimcms2.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetssoftElectronByIP3dBJetTag_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetssoftElectronByPtBJetTag_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetssoftMuonByIP3dBJetTag_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetssoftMuonByPtBJetTag_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetssoftMuonBJetTag_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetssoftElectronTag_CMS2*'))

#Vertex Collections
slimcms2.extend(cms.untracked.vstring('drop *_davertexMaker_*_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_vertexMakerWithBS_*_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_secondaryVertexMaker_*_CMS2*'))

#GSFTrackMaker
slimcms2.extend(cms.untracked.vstring('drop *_*_gsftrksvertexp4_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_gsftrksouterp4_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_gsftrksinnerposition_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_gsftrksouterposition_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_gsftrkslayer1layer_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_gsftrkslayer1sizerphi_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_gsftrkslayer1sizerz_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_gsftrkslayer1charge_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_gsftrkslayer1det_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_gsftrksd0corrPhi_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_gsftrksnlayers3D_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_gsftrksnlayersLost_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_gsftrkslostpixelhits_CMS2*'))

#Other Track-related Makers
slimcms2.extend(cms.untracked.vstring('drop *_trkToVtxAssMaker_*_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_beamHaloMaker_*_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_myTrkJetMetMaker_*_CMS2*'))

#RecoConversionMaker
slimcms2.extend(cms.untracked.vstring('drop *_*_convsrefitPairMomp4_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_convsvtxpos_CMS2*'))

#HypMakers
slimcms2.extend(cms.untracked.vstring('drop *_hypDilepVertexMaker_*_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_hypTrilepMaker_*_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_hypQuadlepMaker_*_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_hypDilepMaker_hypotherjetsp4_CMS2*'))

#HypDilepMaker
slimcms2.extend(cms.untracked.vstring('drop *_*_hypnjets_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypnojets_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypltvalidHits_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypltlostHits_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypltd0_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypltz0_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypltd0corr_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypltz0corr_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypltchi2_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypltndof_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypltd0Err_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypltz0Err_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypltptErr_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypltetaErr_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypltphiErr_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hyplttrkp4_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypllvalidHits_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hyplllostHits_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hyplld0_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypllz0_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hyplld0corr_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypllz0corr_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypllchi2_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypllndof_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hyplld0Err_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypllz0Err_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypllptErr_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hyplletaErr_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypllphiErr_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hyplltrkp4_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypltdPhiunCorrMet_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hyplldPhiunCorrMet_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypltdPhimuCorrMet_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hyplldPhimuCorrMet_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypltdPhitcMet_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hyplldPhitcMet_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypltdPhimetMuonJESCorr_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hyplldPhimetMuonJESCorr_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypdPhinJetunCorrMet_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypdPhinJetmuCorrMet_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypdPhinJettcMet_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypdPhinJetmetMuonJESCorr_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypmt2tcMet_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypmt2muCorrMet_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypmt2metMuonJESCorr_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypjetsidx_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypotherjetsidx_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypjetsp4_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_*_hypotherjetsp4_CMS2*'))

#Other Stuff
slimcms2.extend(cms.untracked.vstring('drop *_electronMaker_elsconvsposp4_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_pfElectronMaker_pfelsposAtEcalp4_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_puSummaryInfoMaker_puInfosumpthighpt_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_puSummaryInfoMaker_puInfozpositions_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_puSummaryInfoMaker_puInfosumptlowpt_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_puSummaryInfoMaker_puInfontrkshighpt_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_puSummaryInfoMaker_puInfontrkslowpt_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_jetToElAssMaker_jetsclosestElectronDR_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_jetToMuAssMaker_jetsclosestMuonDR_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_muonMaker_musgfitouterPosp4_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_muonMaker_musfitfirsthitp4_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_muonMaker_musfitdefaultp4_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_muonMaker_musfitpickyp4_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_muonMaker_musecalposp4_CMS2*'))
slimcms2.extend(cms.untracked.vstring('drop *_scMaker_scsvtxp4_CMS2*'))

### -> pass 2

#slimcms2.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetssimpleSecondaryVertexHighPurBJetTag_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetssimpleSecondaryVertexHighEffBJetTag_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetscombinedSecondaryVertexMVABJetTag_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetscombinedSecondaryVertexBJetTag_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetssimpleSecondaryVertexBJetTag_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetstrackCountingHighPurBJetTag_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetstrackCountingHighEffBJetTag_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetsjetBProbabilityBJetTag_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetsjetProbabilityBJetTag_CMS2*'))

#PFJet Maker
#slimcms2.extend(cms.untracked.vstring('drop *_pfJetMaker_pfjetschargedMultiplicity_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_pfJetMaker_pfjetsneutralMultiplicity_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_pfJetMaker_pfjetsneutralHadronE_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_pfJetMaker_pfjetschargedHadronE_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_pfJetMaker_pfjetsmuonMultiplicity_CMS2*'))

#For testing purposes only!
#slimcms2.extend(cms.untracked.vstring('drop *_pfCandidateMaker_*_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_trackMaker_*_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_hltMaker_*_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_candToGenAssMaker_*_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_pfJetMaker_*_CMS2*'))
#slimcms2.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_*_CMS2*'))
