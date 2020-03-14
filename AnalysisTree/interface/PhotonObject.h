#ifndef PHOTONOBJECT_H
#define PHOTONOBJECT_H

#include "ParticleObject.h"


#define PHOTON_VARIABLES \
PHOTON_VARIABLE(float, etaSC, 0) \
PHOTON_VARIABLE(float, scale_smear_corr, 1) \
PHOTON_VARIABLE(float, scale_smear_corr_scale_totalUp, 1) \
PHOTON_VARIABLE(float, scale_smear_corr_scale_totalDn, 1) \
PHOTON_VARIABLE(float, scale_smear_corr_smear_totalUp, 1) \
PHOTON_VARIABLE(float, scale_smear_corr_smear_totalDn, 1) \
PHOTON_VARIABLE(float, id_MVA_Fall17V2_Val, 0) \
PHOTON_VARIABLE(cms3_photon_mvacat_t, id_MVA_Fall17V2_Cat, 0) \
PHOTON_VARIABLE(bool, id_MVA_Fall17V2_pass_wp90, false) \
PHOTON_VARIABLE(bool, id_MVA_Fall17V2_pass_wp80, false) \
PHOTON_VARIABLE(cms3_photon_cutbasedbits_t, id_cutBased_Fall17V2_Loose_Bits, 0) \
PHOTON_VARIABLE(cms3_photon_cutbasedbits_t, id_cutBased_Fall17V2_Medium_Bits, 0) \
PHOTON_VARIABLE(cms3_photon_cutbasedbits_t, id_cutBased_Fall17V2_Tight_Bits, 0) \
PHOTON_VARIABLE(bool, hasPixelSeed, false) \
PHOTON_VARIABLE(bool, passElectronVeto, false) \
PHOTON_VARIABLE(float, pfIso_comb, 0) \
PHOTON_VARIABLE(float, pfChargedHadronIso_EAcorr, 0)


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

};

#endif
