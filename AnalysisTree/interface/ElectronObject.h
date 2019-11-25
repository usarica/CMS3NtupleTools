#ifndef ELECTRONOBJECT_H
#define ELECTRONOBJECT_H

#include "ParticleObject.h"


#define ELECTRON_VARIABLES \
ELECTRON_VARIABLE(float, scale_smear_corr, 1) \
ELECTRON_VARIABLE(float, scale_smear_corr_scale_totalUp, 1) \
ELECTRON_VARIABLE(float, scale_smear_corr_scale_totalDn, 1) \
ELECTRON_VARIABLE(float, scale_smear_corr_smear_totalUp, 1) \
ELECTRON_VARIABLE(float, scale_smear_corr_smear_totalDn, 1) \
ELECTRON_VARIABLE(float, id_MVA_Fall17V2_Iso_Val, 0) \
ELECTRON_VARIABLE(unsigned int, id_MVA_Fall17V2_Iso_Cat, 0) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_Iso_pass_wpLoose, false) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_Iso_pass_wp90, false) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_Iso_pass_wp80, false) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_Iso_pass_wpHZZ, false) \
ELECTRON_VARIABLE(float, id_MVA_Fall17V2_NoIso_Val, 0) \
ELECTRON_VARIABLE(unsigned int, id_MVA_Fall17V2_NoIso_Cat, 0) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_NoIso_pass_wpLoose, false) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_NoIso_pass_wp90, false) \
ELECTRON_VARIABLE(bool, id_MVA_Fall17V2_NoIso_pass_wp80, false) \
ELECTRON_VARIABLE(float, id_MVA_HZZRun2Legacy_Iso_Val, 0) \
ELECTRON_VARIABLE(unsigned int, id_MVA_HZZRun2Legacy_Iso_Cat, 0) \
ELECTRON_VARIABLE(bool, id_MVA_HZZRun2Legacy_Iso_pass_wpHZZ, false) \
ELECTRON_VARIABLE(unsigned int, id_cutBased_Fall17V2_Veto_Bits, 0) \
ELECTRON_VARIABLE(unsigned int, id_cutBased_Fall17V2_Loose_Bits, 0) \
ELECTRON_VARIABLE(unsigned int, id_cutBased_Fall17V2_Medium_Bits, 0) \
ELECTRON_VARIABLE(unsigned int, id_cutBased_Fall17V2_Tight_Bits, 0) \
ELECTRON_VARIABLE(unsigned int, id_cutBased_Fall17V1_Veto_Bits, 0) \
ELECTRON_VARIABLE(unsigned int, id_cutBased_Fall17V1_Loose_Bits, 0) \
ELECTRON_VARIABLE(unsigned int, id_cutBased_Fall17V1_Medium_Bits, 0) \
ELECTRON_VARIABLE(unsigned int, id_cutBased_Fall17V1_Tight_Bits, 0) \
ELECTRON_VARIABLE(float, pfIso03_comb_nofsr, 0) \
ELECTRON_VARIABLE(float, pfIso04_comb_nofsr, 0) \
ELECTRON_VARIABLE(float, miniIso_comb_nofsr, 0) \
ELECTRON_VARIABLE(float, miniIso_comb_nofsr_uncorrected, 0)


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

  ElectronObject();
  ElectronObject(int id_);
  ElectronObject(int id_, LorentzVector_t const& mom_);
  ElectronObject(const ElectronObject& other);
  ElectronObject& operator=(const ElectronObject& other);
  ~ElectronObject();

  void swap(ElectronObject& other);

  void makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const&);

};

#endif
