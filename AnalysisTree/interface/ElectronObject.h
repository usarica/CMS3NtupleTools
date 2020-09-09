#ifndef ELECTRONOBJECT_H
#define ELECTRONOBJECT_H

#include <string>

#include "ParticleObject.h"


#define ELECTRONS_HAVE_FALL17V1_CUTBASED 0

#define ELECTRON_COMMON_IDISO_MOMENTUMSCALE_VARIABLES \
ELECTRON_VARIABLE(float, etaSC, 0) \
ELECTRON_VARIABLE(float, ecalEnergy, 0) \
ELECTRON_VARIABLE(float, thetaSC_pos, 0) \
ELECTRON_VARIABLE(float, scale_smear_corr, 1) \
ELECTRON_VARIABLE(float, scale_smear_corr_scale_totalUp, 1) \
ELECTRON_VARIABLE(float, scale_smear_corr_scale_totalDn, 1) \
ELECTRON_VARIABLE(float, scale_smear_corr_smear_totalUp, 1) \
ELECTRON_VARIABLE(float, scale_smear_corr_smear_totalDn, 1) \
ELECTRON_VARIABLE(bool, conv_vtx_flag, false) \
ELECTRON_VARIABLE(cms3_electron_missinghits_t, n_missing_inner_hits, 0) \
ELECTRON_VARIABLE(cms3_electron_missinghits_t, n_all_missing_inner_hits, 0) \
ELECTRON_VARIABLE(cms3_electron_charge_consistency_bits_t, charge_consistency_bits, 0) \
ELECTRON_VARIABLE(cms3_egamma_fid_type_mask_t, fid_mask, 0) \
ELECTRON_VARIABLE(cms3_egamma_fid_type_mask_t, type_mask, 0) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_Iso_pass_wpLoose, false) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_Iso_pass_wp90, false) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_Iso_pass_wp80, false) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_Iso_pass_wpHZZ, false) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_NoIso_pass_wpLoose, false) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_NoIso_pass_wp90, false) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_NoIso_pass_wp80, false) \
ELECTRON_VARIABLE(bool, id_MVA_HZZRun2Legacy_Iso_pass_wpHZZ, false) \
ELECTRON_VARIABLE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V2_Veto_Bits, 0) \
ELECTRON_VARIABLE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V2_Loose_Bits, 0) \
ELECTRON_VARIABLE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V2_Medium_Bits, 0) \
ELECTRON_VARIABLE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V2_Tight_Bits, 0) \
ELECTRON_VARIABLE(cms3_electron_cutbasedbits_egPFElectron_t, id_egamma_pfElectron_Bits, 0) \
ELECTRON_VARIABLE(cms3_electron_cutbasedbits_triggeremulation_t, id_cutBased_triggerEmulationV1_Bits, 0) \
ELECTRON_VARIABLE(cms3_electron_cutbasedbits_triggeremulation_t, id_cutBased_triggerEmulationV2_Bits, 0) \
ELECTRON_VARIABLE(float, full5x5_sigmaIEtaIEta, 0) \
ELECTRON_VARIABLE(float, full5x5_sigmaIPhiIPhi, 0) \
ELECTRON_VARIABLE(float, full5x5_r9, 0) \
ELECTRON_VARIABLE(float, r9, 0) \
ELECTRON_VARIABLE(float, E4overE1, 0) \
ELECTRON_VARIABLE(float, seedTime, 0) \
ELECTRON_VARIABLE(float, dxy_firstPV, 0) \
ELECTRON_VARIABLE(float, dz_firstPV, 0) \
ELECTRON_VARIABLE(float, SIP3D, 0) \
ELECTRON_VARIABLE(float, pfIso03_sum_charged_nofsr, 0) \
ELECTRON_VARIABLE(float, pfIso03_sum_neutral_nofsr, 0) \
ELECTRON_VARIABLE(float, pfIso03_comb_nofsr, 0) \
ELECTRON_VARIABLE(float, pfIso04_sum_charged_nofsr, 0) \
ELECTRON_VARIABLE(float, pfIso04_sum_neutral_nofsr, 0) \
ELECTRON_VARIABLE(float, pfIso04_comb_nofsr, 0) \
ELECTRON_VARIABLE(float, miniIso_sum_charged_nofsr, 0) \
ELECTRON_VARIABLE(float, miniIso_sum_neutral_nofsr, 0) \
ELECTRON_VARIABLE(float, miniIso_comb_nofsr, 0) \
ELECTRON_VARIABLE(float, miniIso_comb_nofsr_uncorrected, 0) \
ELECTRON_VARIABLE(cms3_listSize_t, n_associated_pfcands, 0) \
ELECTRON_VARIABLE(cms3_listSize_t, n_associated_pfelectrons, 0) \
ELECTRON_VARIABLE(float, associated_pfcands_sum_px, 0) \
ELECTRON_VARIABLE(float, associated_pfcands_sum_py, 0) \
ELECTRON_VARIABLE(float, min_dR_electron_pfelectron_associated, 0)


#define ELECTRON_MVAID_EXTRA_VARIABLES \
ELECTRON_VARIABLE(float, id_MVA_Fall17V2_Iso_Val, 0) \
ELECTRON_VARIABLE(cms3_electron_mvacat_t, id_MVA_Fall17V2_Iso_Cat, 0) \
ELECTRON_VARIABLE(float, id_MVA_Fall17V2_NoIso_Val, 0) \
ELECTRON_VARIABLE(cms3_electron_mvacat_t, id_MVA_Fall17V2_NoIso_Cat, 0) \
ELECTRON_VARIABLE(float, id_MVA_HZZRun2Legacy_Iso_Val, 0) \
ELECTRON_VARIABLE(cms3_electron_mvacat_t, id_MVA_HZZRun2Legacy_Iso_Cat, 0)

#define ELECTRON_FALL17V1_CUTBASED_VARIABLES \
ELECTRON_VARIABLE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V1_Veto_Bits, 0) \
ELECTRON_VARIABLE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V1_Loose_Bits, 0) \
ELECTRON_VARIABLE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V1_Medium_Bits, 0) \
ELECTRON_VARIABLE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V1_Tight_Bits, 0)

#define ELECTRON_GENINFO_VARIABLES \
ELECTRON_VARIABLE(bool, is_genMatched, false) \
ELECTRON_VARIABLE(bool, is_genMatched_prompt, false)

#if ELECTRONS_HAVE_FALL17V1_CUTBASED == 1
#define ELECTRON_COMMON_VARIABLES \
ELECTRON_COMMON_IDISO_MOMENTUMSCALE_VARIABLES \
ELECTRON_FALL17V1_CUTBASED_VARIABLES
#else
#define ELECTRON_COMMON_VARIABLES \
ELECTRON_COMMON_IDISO_MOMENTUMSCALE_VARIABLES
#endif

#define ELECTRON_RECO_VARIABLES \
ELECTRON_COMMON_VARIABLES \
ELECTRON_MVAID_EXTRA_VARIABLES
#define ELECTRON_VARIABLES \
ELECTRON_RECO_VARIABLES \
ELECTRON_GENINFO_VARIABLES


class ElectronVariables{
public:
#define ELECTRON_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  ELECTRON_VARIABLES;
#undef ELECTRON_VARIABLE

  ElectronVariables();
  ElectronVariables(ElectronVariables const& other);
  ElectronVariables& operator=(const ElectronVariables& other);

  void swap(ElectronVariables& other);

};

class ElectronObject : public ParticleObject{
public:
  ElectronVariables extras;
  float currentSystScale;

  ElectronObject();
  ElectronObject(cms3_id_t id_);
  ElectronObject(cms3_id_t id_, LorentzVector_t const& mom_);
  ElectronObject(const ElectronObject& other);
  ElectronObject& operator=(const ElectronObject& other);
  ~ElectronObject();

  void swap(ElectronObject& other);

  void makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const&);
  void applyFSRIsoCorr(ParticleObject::LorentzVector_t::Scalar const& dR_fsr, ParticleObject::LorentzVector_t::Scalar const& pt_fsr);

  float const& etaSC() const{ return extras.etaSC; }

  bool isEBEEGap() const;
  bool isEB() const;
  bool isEE() const;
  bool isAnyGap() const;
  bool isGap() const{ return this->isAnyGap(); }

  ParticleObject::LorentzVector_t::Scalar uncorrected_pt() const;
  ParticleObject::LorentzVector_t uncorrected_p4() const;

  bool hasFSR() const;

};

#endif
