#ifndef METOBJECT_H
#define METOBJECT_H

#include "SystematicVariations.h"
#include "ParticleObject.h"


#define MET_RECORDED_VARIABLES \
MET_VARIABLE(float, met_Nominal, 0) \
MET_VARIABLE(float, metPhi_Nominal, 0) \
MET_VARIABLE(float, met_JECUp, 0) \
MET_VARIABLE(float, metPhi_JECUp, 0) \
MET_VARIABLE(float, met_JECDn, 0) \
MET_VARIABLE(float, metPhi_JECDn, 0) \
MET_VARIABLE(float, metSignificance, 0)

#define MET_JERSHIFT_VARIABLES \
MET_VARIABLE(float, metShift_px_JERNominal, 0) \
MET_VARIABLE(float, metShift_py_JERNominal, 0) \
MET_VARIABLE(float, metShift_px_JERUp, 0) \
MET_VARIABLE(float, metShift_py_JERUp, 0) \
MET_VARIABLE(float, metShift_px_JERDn, 0) \
MET_VARIABLE(float, metShift_py_JERDn, 0)

#define MET_EXTRA_PT_VARIABLES \
MET_VARIABLE(float, met_original, 0) \
MET_VARIABLE(float, met_JERUp, 0) \
MET_VARIABLE(float, met_JERDn, 0) \
MET_VARIABLE(float, met_PUUp, 0) \
MET_VARIABLE(float, met_PUDn, 0) \
MET_VARIABLE(float, met_METUp, 0) \
MET_VARIABLE(float, met_METDn, 0) \

#define MET_EXTRA_PHI_VARIABLES \
MET_VARIABLE(float, metPhi_original, 0) \
MET_VARIABLE(float, metPhi_JERUp, 0) \
MET_VARIABLE(float, metPhi_JERDn, 0) \
MET_VARIABLE(float, metPhi_PUUp, 0) \
MET_VARIABLE(float, metPhi_PUDn, 0) \
MET_VARIABLE(float, metPhi_METUp, 0) \
MET_VARIABLE(float, metPhi_METDn, 0)


#define MET_EXTRA_VARIABLES \
MET_EXTRA_PT_VARIABLES \
MET_EXTRA_PHI_VARIABLES

#define MET_VARIABLES \
MET_RECORDED_VARIABLES \
MET_JERSHIFT_VARIABLES \
MET_EXTRA_VARIABLES


class METVariables{
public:
#define MET_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  MET_VARIABLES;
#undef MET_VARIABLE

  METVariables();
  METVariables(METVariables const& other);
  METVariables& operator=(const METVariables& other);

  void swap(METVariables& other);

};

class METObject{
public:
  METVariables extras;

protected:
  SystematicsHelpers::SystematicVariationTypes currentSyst;
  ParticleObject::LorentzVector_t currentXYshift;
  ParticleObject::LorentzVector_t currentJERShift;
  ParticleObject::LorentzVector_t particleMomentumCorrections;

  void setJERShifts();

public:
  METObject();
  METObject(const METObject& other);
  METObject& operator=(const METObject& other);
  ~METObject();

  void swap(METObject& other);

  void setSystematic(SystematicsHelpers::SystematicVariationTypes const&);

  void setXYShift(ParticleObject::LorentzVector_t const& shift){ currentXYshift=shift; }
  void setXYShift(float const& shift_x, float const& shift_y){ this->setXYShift(ParticleObject::LorentzVector_t(shift_x, shift_y, 0., 0.)); }

  void setParticleShifts(ParticleObject::LorentzVector_t const& shift){ particleMomentumCorrections=shift; }
  void setParticleShifts(float const& shift_x, float const& shift_y){ this->setParticleShifts(ParticleObject::LorentzVector_t(shift_x, shift_y, 0., 0.)); }

  void getPtPhi(float& pt, float& phi, bool addXYShifts, bool addJERShifts, bool addParticleShifts) const;
  ParticleObject::LorentzVector_t::Scalar met(bool addXYShifts, bool addJERShifts, bool addParticleShifts) const;
  ParticleObject::LorentzVector_t::Scalar phi(bool addXYShifts, bool addJERShifts, bool addParticleShifts) const;
  ParticleObject::LorentzVector_t::Scalar pt(bool addXYShifts, bool addJERShifts, bool addParticleShifts) const{ return met(addXYShifts, addJERShifts, addParticleShifts); }
  ParticleObject::LorentzVector_t::Scalar px(bool addXYShifts, bool addJERShifts, bool addParticleShifts, float phi_rot=0) const;
  ParticleObject::LorentzVector_t::Scalar py(bool addXYShifts, bool addJERShifts, bool addParticleShifts, float phi_rot=0) const;
  ParticleObject::LorentzVector_t p4(bool addXYShifts, bool addJERShifts, bool addParticleShifts, float phi_rot=0) const;

};

#endif
