#ifndef AK8JETOBJECT_H
#define AK8JETOBJECT_H

#include "ParticleObject.h"


#define AK8JET_RECO_VARIABLES \
AK8JET_VARIABLE(float, JECNominal, 0)

#define AK8JET_GENINFO_VARIABLES \
AK8JET_VARIABLE(float, relJECUnc, 0) \
AK8JET_VARIABLE(float, JERNominal, 1) \
AK8JET_VARIABLE(float, JERDn, 1) \
AK8JET_VARIABLE(float, JERUp, 1) \
AK8JET_VARIABLE(bool, is_genMatched, false) \
AK8JET_VARIABLE(bool, is_genMatched_fullCone, false) \
AK8JET_VARIABLE(cms3_jet_genflavor_t, partonFlavour, 0) \
AK8JET_VARIABLE(cms3_jet_genflavor_t, hadronFlavour, 0)

#define AK8JET_VARIABLES \
AK8JET_RECO_VARIABLES \
AK8JET_GENINFO_VARIABLES


class AK8JetVariables{
public:
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  AK8JET_VARIABLES;
#undef AK8JET_VARIABLE

  AK8JetVariables();
  AK8JetVariables(AK8JetVariables const& other);
  AK8JetVariables& operator=(const AK8JetVariables& other);

  void swap(AK8JetVariables& other);

};

class AK8JetObject : public ParticleObject{
protected:
  LorentzVector_t mom_original;

public:
  constexpr static float ConeRadiusConstant = 0.8;

  AK8JetVariables extras;
  float currentSystScale;

  AK8JetObject();
  AK8JetObject(LorentzVector_t const& mom_);
  AK8JetObject(const AK8JetObject& other);
  AK8JetObject& operator=(const AK8JetObject& other);
  ~AK8JetObject();

  void swap(AK8JetObject& other);

  void makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const& syst);

};

#endif
