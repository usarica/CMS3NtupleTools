import FWCore.ParameterSet.Config as cms

slimcms3 = cms.untracked.vstring()

#TrackMaker
slimcms3.extend(cms.untracked.vstring('drop *_trackMaker_trksvertexp4_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_trackMaker_trksouterp4_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_trackMaker_trksinnerposition_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_trackMaker_trksouterposition_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_trackMaker_trkslayer1layer_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_trackMaker_trkslayer1sizerphi_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_trackMaker_trkslayer1sizerz_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_trackMaker_trkslayer1charge_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_trackMaker_trkslayer1det_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_trackMaker_trksd0corrPhi_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_trackMaker_trksnlayers3D_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_trackMaker_trksnlayersLost_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_trackMaker_trkslostpixelhits_CMS3*'))

#Jet Collections
slimcms3.extend(cms.untracked.vstring('drop *_jptMaker_*_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_bTagJPTJetMaker_*_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_jetMaker_*_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_bTagMaker_*_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_trkJetMaker_*_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_bTagTrkMaker_*_CMS3*'))

#bTagPFJetMaker
slimcms3.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetssoftElectronByIP3dBJetTag_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetssoftElectronByPtBJetTag_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetssoftMuonByIP3dBJetTag_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetssoftMuonByPtBJetTag_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetssoftMuonBJetTag_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_bTagPFJetMaker_pfjetssoftElectronTag_CMS3*'))

#Vertex Collections
slimcms3.extend(cms.untracked.vstring('drop *_davertexMaker_*_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_vertexMakerWithBS_*_CMS3*'))

#GSFTrackMaker
slimcms3.extend(cms.untracked.vstring('drop *_*_gsftrksvertexp4_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_gsftrksouterp4_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_gsftrksinnerposition_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_gsftrksouterposition_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_gsftrkslayer1layer_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_gsftrkslayer1sizerphi_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_gsftrkslayer1sizerz_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_gsftrkslayer1charge_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_gsftrkslayer1det_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_gsftrksd0corrPhi_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_gsftrksnlayers3D_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_gsftrksnlayersLost_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_gsftrkslostpixelhits_CMS3*'))

#Other Track-related Makers
slimcms3.extend(cms.untracked.vstring('drop *_trkToVtxAssMaker_*_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_myTrkJetMetMaker_*_CMS3*'))

#RecoConversionMaker
slimcms3.extend(cms.untracked.vstring('drop *_*_convsrefitPairMomp4_CMS3*')) #this one actually does something
slimcms3.extend(cms.untracked.vstring('drop *_*_convsvtxpos_CMS3*')) #this one too

#HypMakers
slimcms3.extend(cms.untracked.vstring('drop *_hypDilepVertexMaker_*_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_hypTrilepMaker_*_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_hypQuadlepMaker_*_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_hypDilepMaker_hypotherjetsp4_CMS3*'))

#HypDilepMaker
slimcms3.extend(cms.untracked.vstring('drop *_*_hypnjets_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypnojets_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypltvalidHits_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypltlostHits_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypltd0_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypltz0_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypltd0corr_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypltz0corr_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypltchi2_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypltndof_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypltd0Err_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypltz0Err_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypltptErr_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypltetaErr_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypltphiErr_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hyplttrkp4_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypllvalidHits_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hyplllostHits_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hyplld0_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypllz0_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hyplld0corr_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypllz0corr_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypllchi2_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypllndof_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hyplld0Err_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypllz0Err_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypllptErr_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hyplletaErr_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypllphiErr_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hyplltrkp4_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypltdPhiunCorrMet_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hyplldPhiunCorrMet_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypltdPhimuCorrMet_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hyplldPhimuCorrMet_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypltdPhitcMet_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hyplldPhitcMet_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypltdPhimetMuonJESCorr_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hyplldPhimetMuonJESCorr_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypdPhinJetunCorrMet_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypdPhinJetmuCorrMet_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypdPhinJettcMet_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypdPhinJetmetMuonJESCorr_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypmt2tcMet_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypmt2muCorrMet_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypmt2metMuonJESCorr_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypjetsidx_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypotherjetsidx_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypjetsp4_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_*_hypotherjetsp4_CMS3*'))

#Other Stuff
slimcms3.extend(cms.untracked.vstring('drop *_electronMaker_elsconvsposp4_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_pfElectronMaker_pfelsposAtEcalp4_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_puSummaryInfoMaker_puInfosumpthighpt_CMS3*')) #this one
slimcms3.extend(cms.untracked.vstring('drop *_puSummaryInfoMaker_puInfozpositions_CMS3*')) #and this
slimcms3.extend(cms.untracked.vstring('drop *_puSummaryInfoMaker_puInfosumptlowpt_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_puSummaryInfoMaker_puInfontrkshighpt_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_puSummaryInfoMaker_puInfontrkslowpt_CMS3*')) #through here
slimcms3.extend(cms.untracked.vstring('drop *_jetToElAssMaker_jetsclosestElectronDR_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_jetToMuAssMaker_jetsclosestMuonDR_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_muonMaker_musgfitouterPosp4_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_muonMaker_musfitfirsthitp4_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_muonMaker_musfitdefaultp4_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_muonMaker_musfitpickyp4_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_muonMaker_musecalposp4_CMS3*'))
slimcms3.extend(cms.untracked.vstring('drop *_scMaker_scsvtxp4_CMS3*'))
