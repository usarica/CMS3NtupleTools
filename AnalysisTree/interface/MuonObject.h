#ifndef MUONOBJECT_H
#define MUONOBJECT_H

#include "ParticleObject.h"


#define MUON_VARIABLES \
MUON_VARIABLE(unsigned int, POG_selector_bits, 0) \
MUON_VARIABLE(float, pfIso03_comb_nofsr, 0) \
MUON_VARIABLE(float, pfIso04_comb_nofsr, 0) \
MUON_VARIABLE(float, miniIso_comb_nofsr, 0) \
MUON_VARIABLE(float, miniIso_comb_nofsr_uncorrected, 0) \
MUON_VARIABLE(int, time_comb_ndof, 0) \
MUON_VARIABLE(float, time_comb_IPInOut, 0) \
MUON_VARIABLE(float, time_comb_IPOutIn, 0) \
MUON_VARIABLE(float, time_comb_IPInOutError, 0) \
MUON_VARIABLE(float, time_comb_IPOutInError, 0) \
MUON_VARIABLE(int, time_rpc_ndof, 0) \
MUON_VARIABLE(float, time_rpc_IPInOut, 0) \
MUON_VARIABLE(float, time_rpc_IPOutIn, 0) \
MUON_VARIABLE(float, time_rpc_IPInOutError, 0) \
MUON_VARIABLE(float, time_rpc_IPOutInError, 0) \
MUON_VARIABLE(float, scale_smear_pt_corr, 0) \
MUON_VARIABLE(float, scale_smear_pt_corr_scale_totalUp, 0) \
MUON_VARIABLE(float, scale_smear_pt_corr_scale_totalDn, 0) \
MUON_VARIABLE(float, scale_smear_pt_corr_smear_totalUp, 0) \
MUON_VARIABLE(float, scale_smear_pt_corr_smear_totalDn, 0)


class MuonVariables{
public:
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  MUON_VARIABLES;
#undef MUON_VARIABLE

  MuonVariables();
  MuonVariables(MuonVariables const& other);
  MuonVariables& operator=(const MuonVariables& other);

  void swap(MuonVariables& other);

};

class MuonObject : public ParticleObject{
public:
  MuonVariables extras;
  float currentSystScale;

  MuonObject();
  MuonObject(int id_);
  MuonObject(int id_, LorentzVector_t const& mom_);
  MuonObject(const MuonObject& other);
  MuonObject& operator=(const MuonObject& other);
  ~MuonObject();

  void swap(MuonObject& other);

  void makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const&);

};

#endif
