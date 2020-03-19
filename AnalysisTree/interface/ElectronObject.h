#ifndef ELECTRONOBJECT_H
#define ELECTRONOBJECT_H

#include "ParticleObject.h"

#define ELECTRONS_HAVE_FALL17V1_CUTBASED 0

#define ELECTRON_COMMON_VARIABLES \
ELECTRON_VARIABLE(float, etaSC, 0) \
ELECTRON_VARIABLE(float, scale_smear_corr, 1) \
ELECTRON_VARIABLE(float, scale_smear_corr_scale_totalUp, 1) \
ELECTRON_VARIABLE(float, scale_smear_corr_scale_totalDn, 1) \
ELECTRON_VARIABLE(float, scale_smear_corr_smear_totalUp, 1) \
ELECTRON_VARIABLE(float, scale_smear_corr_smear_totalDn, 1) \
ELECTRON_VARIABLE(float, id_MVA_Fall17V2_Iso_Val, 0) \
ELECTRON_VARIABLE(cms3_electron_mvacat_t, id_MVA_Fall17V2_Iso_Cat, 0) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_Iso_pass_wpLoose, false) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_Iso_pass_wp90, false) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_Iso_pass_wp80, false) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_Iso_pass_wpHZZ, false) \
ELECTRON_VARIABLE(float, id_MVA_Fall17V2_NoIso_Val, 0) \
ELECTRON_VARIABLE(cms3_electron_mvacat_t, id_MVA_Fall17V2_NoIso_Cat, 0) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_NoIso_pass_wpLoose, false) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_NoIso_pass_wp90, false) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_NoIso_pass_wp80, false) \
ELECTRON_VARIABLE(float, id_MVA_HZZRun2Legacy_Iso_Val, 0) \
ELECTRON_VARIABLE(cms3_electron_mvacat_t, id_MVA_HZZRun2Legacy_Iso_Cat, 0) \
ELECTRON_VARIABLE(bool, id_MVA_HZZRun2Legacy_Iso_pass_wpHZZ, false) \
ELECTRON_VARIABLE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V2_Veto_Bits, 0) \
ELECTRON_VARIABLE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V2_Loose_Bits, 0) \
ELECTRON_VARIABLE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V2_Medium_Bits, 0) \
ELECTRON_VARIABLE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V2_Tight_Bits, 0) \
ELECTRON_VARIABLE(float, pfIso03_sum_charged_nofsr, 0) \
ELECTRON_VARIABLE(float, pfIso03_sum_neutral_nofsr, 0) \
ELECTRON_VARIABLE(float, pfIso03_comb_nofsr, 0) \
ELECTRON_VARIABLE(float, pfIso04_sum_charged_nofsr, 0) \
ELECTRON_VARIABLE(float, pfIso04_sum_neutral_nofsr, 0) \
ELECTRON_VARIABLE(float, pfIso04_comb_nofsr, 0) \
ELECTRON_VARIABLE(float, miniIso_sum_charged_nofsr, 0) \
ELECTRON_VARIABLE(float, miniIso_sum_neutral_nofsr, 0) \
ELECTRON_VARIABLE(float, miniIso_comb_nofsr, 0) \
ELECTRON_VARIABLE(float, miniIso_comb_nofsr_uncorrected, 0)
#define ELECTRON_FALL17V1_CUTBASED_VARIABLES \
ELECTRON_VARIABLE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V1_Veto_Bits, 0) \
ELECTRON_VARIABLE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V1_Loose_Bits, 0) \
ELECTRON_VARIABLE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V1_Medium_Bits, 0) \
ELECTRON_VARIABLE(cms3_electron_cutbasedbits_t, id_cutBased_Fall17V1_Tight_Bits, 0)

#if ELECTRONS_HAVE_FALL17V1_CUTBASED == 1
#define ELECTRON_VARIABLES \
ELECTRON_COMMON_VARIABLES \
ELECTRON_FALL17V1_CUTBASED_VARIABLES
#else
#define ELECTRON_VARIABLES \
ELECTRON_COMMON_VARIABLES
#endif

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

  float const& etaSC() const{ return extras.etaSC; }

  ParticleObject::LorentzVector_t::Scalar uncorrected_pt() const;

};

#endif
