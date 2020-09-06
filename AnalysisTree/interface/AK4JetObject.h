#ifndef AK4JETOBJECT_H
#define AK4JETOBJECT_H

#include <string>

#include "ParticleObject.h"
#include "BTagCalibrationStandalone.h"


#define AK4JET_CORE_VARIABLES \
AK4JET_VARIABLE(bool, pass_looseId, false) \
AK4JET_VARIABLE(bool, pass_tightId, false) \
AK4JET_VARIABLE(bool, pass_leptonVetoId, false) \
AK4JET_VARIABLE(cms3_jet_pujetid_t, pileupJetId, 0) \
AK4JET_VARIABLE(cms3_jet_pujetid_t, pileupJetId_default, 0) \
AK4JET_VARIABLE(cms3_metsafety_t, isMETJERCSafe_Bits, 0) \
AK4JET_VARIABLE(cms3_metsafety_t, isMETJERCSafe_p4Preserved_Bits, 0) \
AK4JET_VARIABLE(float, JECNominal, 1) \
AK4JET_VARIABLE(float, JECL1Nominal, 1) \
AK4JET_VARIABLE(float, mucands_sump4_px, 1) \
AK4JET_VARIABLE(float, mucands_sump4_py, 1) \
AK4JET_VARIABLE(float, NEMF, 0) \
AK4JET_VARIABLE(float, CEMF, 0)

#define AK4JET_GENINFO_VARIABLES \
AK4JET_VARIABLE(float, relJECUnc, 0) \
AK4JET_VARIABLE(float, relJECUnc_nomus, 0) \
AK4JET_VARIABLE(float, relJECUnc_nomus_JERNominal, 0) \
AK4JET_VARIABLE(float, JERNominal, 1) \
AK4JET_VARIABLE(float, JERDn, 1) \
AK4JET_VARIABLE(float, JERUp, 1) \
AK4JET_VARIABLE(bool, is_genMatched, false) \
AK4JET_VARIABLE(bool, is_genMatched_fullCone, false) \
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

#define AK4JET_RECO_VARIABLES \
AK4JET_CORE_VARIABLES \
AK4JET_BTAGGING_VARIABLES
#define AK4JET_VARIABLES \
AK4JET_RECO_VARIABLES \
AK4JET_GENINFO_VARIABLES


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
protected:
  LorentzVector_t mom_original;

public:
  constexpr static float ConeRadiusConstant = 0.4;

  static const std::string colName;

  AK4JetVariables extras;
  SystematicsHelpers::SystematicVariationTypes currentSyst;
  float currentJEC;
  float currentJER;
  float currentSystScale;

  AK4JetObject();
  AK4JetObject(LorentzVector_t const& mom_);
  AK4JetObject(const AK4JetObject& other);
  AK4JetObject& operator=(const AK4JetObject& other);
  ~AK4JetObject();

  void swap(AK4JetObject& other);

  void makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const& syst);

  BTagEntry::JetFlavor getBTagJetFlavor() const;
  float getBtagValue() const;

  float const& getJECValue() const{ return currentJEC; }
  float const& getJERValue() const{ return currentJER; }

  LorentzVector_t uncorrected_p4() const{ return this->mom_original; }
  LorentzVector_t p4_mucands() const{ return LorentzVector_t(this->extras.mucands_sump4_px, this->extras.mucands_sump4_py, 0, 0); }
  // Unfortunately, there could be multiple versions of this function. This one is the most straightforward version.
  LorentzVector_t p4_nomus_basic() const{ return this->p4() - this->p4_mucands(); }
  // And here is why:
  bool isMETSafe(SystematicsHelpers::SystematicVariationTypes const& syst, bool useP4Preserved, bool applyJER) const;
  bool isMETSafe(bool useP4Preserved, bool applyJER) const{ return isMETSafe(currentSyst, useP4Preserved, applyJER); }
  // Returns whether a contribution can/should be acquired. Notice p4_metShift=-(corrected - uncorrected)
  bool getT1METShift(SystematicsHelpers::SystematicVariationTypes const& syst, bool useP4Preserved, bool applyJER, ParticleObject::LorentzVector_t& p4_metShift) const;
  bool getT1METShift(bool useP4Preserved, bool applyJER, ParticleObject::LorentzVector_t& p4_metShift) const{ return getT1METShift(currentSyst, useP4Preserved, applyJER, p4_metShift); }

  static LorentzVector_t compute_METShift(
    bool preserve_corrected_jet_p4,
    ParticleObject::LorentzVector_t const& p4_jet_uncorrected, ParticleObject::LorentzVector_t const& p4_mucands,
    float const& JEC_L1L2L3, float const& JEC_L1, float const& JERval,
    char const& iJECshift, float const& nativeRelJECUnc
  );

};

#endif
