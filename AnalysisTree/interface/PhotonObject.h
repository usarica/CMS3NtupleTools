#ifndef PHOTONOBJECT_H
#define PHOTONOBJECT_H

#include <string>

#include "ParticleObject.h"


#define PHOTON_COMMON_VARIABLES \
PHOTON_VARIABLE(float, etaSC, 0) \
PHOTON_VARIABLE(float, scale_smear_corr, 1) \
PHOTON_VARIABLE(float, scale_smear_corr_scale_totalUp, 1) \
PHOTON_VARIABLE(float, scale_smear_corr_scale_totalDn, 1) \
PHOTON_VARIABLE(float, scale_smear_corr_smear_totalUp, 1) \
PHOTON_VARIABLE(float, scale_smear_corr_smear_totalDn, 1) \
PHOTON_VARIABLE(cms3_egamma_fid_type_mask_t, fid_mask, 0) \
PHOTON_VARIABLE(bool, hasPixelSeed, false) \
PHOTON_VARIABLE(bool, passElectronVeto, false) \
PHOTON_VARIABLE(bool, id_MVA_Fall17V2_pass_wp90, false) \
PHOTON_VARIABLE(bool, id_MVA_Fall17V2_pass_wp80, false) \
PHOTON_VARIABLE(cms3_photon_cutbasedbits_t, id_cutBased_Fall17V2_Loose_Bits, 0) \
PHOTON_VARIABLE(cms3_photon_cutbasedbits_t, id_cutBased_Fall17V2_Medium_Bits, 0) \
PHOTON_VARIABLE(cms3_photon_cutbasedbits_t, id_cutBased_Fall17V2_Tight_Bits, 0) \
PHOTON_VARIABLE(cms3_photon_cutbasedbits_hgg_t, id_cutBased_HGG_Bits, 0) \
PHOTON_VARIABLE(cms3_photon_cutbasedbits_egPFPhoton_t, id_egamma_pfPhoton_Bits, 0) \
PHOTON_VARIABLE(bool, hcal_is_valid, 0) \
PHOTON_VARIABLE(float, hOverE, 0) \
PHOTON_VARIABLE(float, hOverEtowBC, 0) \
PHOTON_VARIABLE(float, r9, 0) \
PHOTON_VARIABLE(float, sigmaIEtaIEta, 0) \
PHOTON_VARIABLE(float, sigmaIPhiIPhi, 0) \
PHOTON_VARIABLE(float, full5x5_r9, 0) \
PHOTON_VARIABLE(float, full5x5_sigmaIEtaIEta, 0) \
PHOTON_VARIABLE(float, full5x5_sigmaIPhiIPhi, 0) \
PHOTON_VARIABLE(float, MIPTotalEnergy, 0) \
PHOTON_VARIABLE(float, E4overE1, 0) \
PHOTON_VARIABLE(float, seedTime, 0) \
PHOTON_VARIABLE(float, pfChargedHadronIso_EAcorr, 0) \
PHOTON_VARIABLE(float, pfNeutralHadronIso_EAcorr, 0) \
PHOTON_VARIABLE(float, pfEMIso_EAcorr, 0) \
PHOTON_VARIABLE(float, pfIso_comb, 0) \
PHOTON_VARIABLE(float, pfWorstChargedHadronIso_allVtxs, 0) \
PHOTON_VARIABLE(float, pfWorstChargedHadronIso_pt0p1_minDR0p02_allVtxs, 0) \
PHOTON_VARIABLE(float, pfWorstChargedHadronIso_pt0p1_minDR0p02_dxy_dz_firstPV_allVtxs, 0) \
PHOTON_VARIABLE(float, trkIso03_hollow, 0) \
PHOTON_VARIABLE(float, trkIso03_solid, 0) \
PHOTON_VARIABLE(float, ecalRecHitIso03, 0) \
PHOTON_VARIABLE(float, hcalTowerIso03, 0) \
PHOTON_VARIABLE(cms3_listSize_t, n_associated_pfcands, 0) \
PHOTON_VARIABLE(cms3_listSize_t, n_associated_pfphotons, 0) \
PHOTON_VARIABLE(float, associated_pfcands_sum_px, 0) \
PHOTON_VARIABLE(float, associated_pfcands_sum_py, 0) \
PHOTON_VARIABLE(float, min_dR_photon_pfphoton_associated, 0) \
PHOTON_VARIABLE(float, closestPFPhoton_associated_px, 0) \
PHOTON_VARIABLE(float, closestPFPhoton_associated_py, 0)

#define PHOTON_GENINFO_VARIABLES \
PHOTON_VARIABLE(bool, is_genMatched, false) \
PHOTON_VARIABLE(bool, is_genMatched_prompt, false)

#define PHOTON_MVAID_EXTRA_VARIABLES \
PHOTON_VARIABLE(float, id_MVA_Fall17V2_Val, 0) \
PHOTON_VARIABLE(cms3_photon_mvacat_t, id_MVA_Fall17V2_Cat, 0)

#define PHOTON_RECO_VARIABLES \
PHOTON_COMMON_VARIABLES \
/*PHOTON_MVAID_EXTRA_VARIABLES*/

#define PHOTON_VARIABLES \
PHOTON_RECO_VARIABLES \
PHOTON_GENINFO_VARIABLES

class PhotonVariables{
public:
#define PHOTON_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  PHOTON_VARIABLES;
#undef PHOTON_VARIABLE

  PhotonVariables();
  PhotonVariables(PhotonVariables const& other);
  PhotonVariables& operator=(const PhotonVariables& other);

  void swap(PhotonVariables& other);

};

class PhotonObject : public ParticleObject{
public:
  PhotonVariables extras;
  float currentSystScale;

  PhotonObject();
  PhotonObject(LorentzVector_t const& mom_);
  PhotonObject(const PhotonObject& other);
  PhotonObject& operator=(const PhotonObject& other);
  ~PhotonObject();

  void swap(PhotonObject& other);

  void makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const&);

  float const& etaSC() const{ return extras.etaSC; }

  bool isEBEEGap() const;
  bool isEB() const;
  bool isEE() const;
  bool isAnyGap() const;
  bool isGap() const{ return this->isAnyGap(); }

  ParticleObject::LorentzVector_t::Scalar uncorrected_pt() const;
  ParticleObject::LorentzVector_t uncorrected_p4() const;

};

#endif
