#ifndef AK4JETOBJECT_H
#define AK4JETOBJECT_H

#include "ParticleObject.h"
#include "BTagCalibrationStandalone.h"


#define AK4JET_CORE_VARIABLES \
AK4JET_VARIABLE(bool, is_genMatched, false) \
AK4JET_VARIABLE(bool, is_genMatched_fullCone, false) \
AK4JET_VARIABLE(bool, pass_looseId, false) \
AK4JET_VARIABLE(bool, pass_tightId, false) \
AK4JET_VARIABLE(bool, pass_leptonVetoId, false) \
AK4JET_VARIABLE(cms3_jet_pujetid_t, pileupJetId, false) \
AK4JET_VARIABLE(cms3_jet_pujetid_t, pileupJetId_default, false) \
/*AK4JET_VARIABLE(size_t, n_pfcands, 0)*/ \
/*AK4JET_VARIABLE(size_t, n_mucands, 0)*/ \
/*AK4JET_VARIABLE(float, area, 0)*/ \
/*AK4JET_VARIABLE(float, pt_resolution, 0)*/ \
/*AK4JET_VARIABLE(float, ptDistribution, 0)*/ \
/*AK4JET_VARIABLE(float, totalMultiplicity, 0)*/ \
/*AK4JET_VARIABLE(float, axis1, 0)*/ \
/*AK4JET_VARIABLE(float, axis2, 0)*/ \
AK4JET_VARIABLE(float, JECNominal, 1) \
AK4JET_VARIABLE(float, JECUp, 1) \
AK4JET_VARIABLE(float, JECDn, 1) \
AK4JET_VARIABLE(float, JERNominal, 1) \
AK4JET_VARIABLE(float, JERUp, 1) \
AK4JET_VARIABLE(float, JERDn, 1) \
AK4JET_VARIABLE(float, JECL1Nominal, 1) \
AK4JET_VARIABLE(float, mucands_sump4_px, 1) \
AK4JET_VARIABLE(float, mucands_sump4_py, 1) \
AK4JET_VARIABLE(float, NEMF, 0) \
AK4JET_VARIABLE(float, CEMF, 0) \
AK4JET_VARIABLE(cms3_jet_genflavor_t, partonFlavour, 0) \
AK4JET_VARIABLE(cms3_jet_genflavor_t, hadronFlavour, 0)

#define AK4JET_BTAGGING_VARIABLES \
AK4JET_VARIABLE(float, deepFlavourprobb, -1) \
AK4JET_VARIABLE(float, deepFlavourprobbb, -1) \
AK4JET_VARIABLE(float, deepFlavourprobc, -1) \
AK4JET_VARIABLE(float, deepFlavourprobg, -1) \
AK4JET_VARIABLE(float, deepFlavourproblepb, -1) \
AK4JET_VARIABLE(float, deepFlavourprobuds, -1) \
AK4JET_VARIABLE(float, deepCSVprobb, -1) \
AK4JET_VARIABLE(float, deepCSVprobbb, -1) \
AK4JET_VARIABLE(float, deepCSVprobc, -1) \
AK4JET_VARIABLE(float, deepCSVprobudsg, -1)

#define AK4JET_VARIABLES \
AK4JET_CORE_VARIABLES \
AK4JET_BTAGGING_VARIABLES


class AK4JetVariables{
public:
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  AK4JET_VARIABLES;
#undef AK4JET_VARIABLE

  AK4JetVariables();
  AK4JetVariables(AK4JetVariables const& other);
  AK4JetVariables& operator=(const AK4JetVariables& other);

  void swap(AK4JetVariables& other);

};

class AK4JetObject : public ParticleObject{
public:
  constexpr static float ConeRadiusConstant = 0.4;

  AK4JetVariables extras;
  float currentJEC_full;
  float currentJEC_L1only;
  float currentJER;
  float currentSystScale;

  AK4JetObject();
  AK4JetObject(LorentzVector_t const& mom_);
  AK4JetObject(const AK4JetObject& other);
  AK4JetObject& operator=(const AK4JetObject& other);
  ~AK4JetObject();

  void swap(AK4JetObject& other);

  void makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const&);

  BTagEntry::JetFlavor getBTagJetFlavor() const;
  float getBtagValue() const;

  float const& getJECValue(bool useL1only=false) const{ return (!useL1only? currentJEC_full : currentJEC_L1only); }
  float const& getJERValue() const{ return currentJER; }

  LorentzVector_t uncorrected_p4() const{ return this->p4()*(1.f/currentSystScale); }
  LorentzVector_t p4_nomu() const;

};

#endif
