#ifndef METOBJECT_H
#define METOBJECT_H

#include "SystematicVariations.h"
#include "ParticleObject.h"


#define MET_RECORDED_CORE_VARIABLES \
MET_VARIABLE(float, met_Nominal, 0) \
MET_VARIABLE(float, metPhi_Nominal, 0) \
MET_VARIABLE(float, met_Raw, 0) \
MET_VARIABLE(float, metPhi_Raw, 0)
#define MET_RECORDED_GENINFO_VARIABLES \
MET_VARIABLE(float, met_JECDn, 0) \
MET_VARIABLE(float, metPhi_JECDn, 0) \
MET_VARIABLE(float, met_JECUp, 0) \
MET_VARIABLE(float, metPhi_JECUp, 0)
#define MET_RECORDED_VARIABLES \
MET_RECORDED_CORE_VARIABLES \
MET_RECORDED_GENINFO_VARIABLES

// JER shifts are not recorded for PUPPI MET because they do not apply.
#define MET_SHIFT_CORE_VARIABLES \
MET_VARIABLE(float, metShift_px_JECNominal, 0) \
MET_VARIABLE(float, metShift_py_JECNominal, 0) \
MET_VARIABLE(float, metShift_px_JECNominal_JERNominal, 0) \
MET_VARIABLE(float, metShift_py_JECNominal_JERNominal, 0) \
MET_VARIABLE(float, metShift_p4Preserved_px_JECNominal, 0) \
MET_VARIABLE(float, metShift_p4Preserved_py_JECNominal, 0) \
MET_VARIABLE(float, metShift_p4Preserved_px_JECNominal_JERNominal, 0) \
MET_VARIABLE(float, metShift_p4Preserved_py_JECNominal_JERNominal, 0)
#define MET_SHIFT_GENINFO_VARIABLES \
MET_VARIABLE(float, metShift_px_JECDn, 0) \
MET_VARIABLE(float, metShift_py_JECDn, 0) \
MET_VARIABLE(float, metShift_px_JECUp, 0) \
MET_VARIABLE(float, metShift_py_JECUp, 0) \
MET_VARIABLE(float, metShift_px_JECNominal_JERDn, 0) \
MET_VARIABLE(float, metShift_py_JECNominal_JERDn, 0) \
MET_VARIABLE(float, metShift_px_JECNominal_JERUp, 0) \
MET_VARIABLE(float, metShift_py_JECNominal_JERUp, 0) \
MET_VARIABLE(float, metShift_px_JECDn_JERNominal, 0) \
MET_VARIABLE(float, metShift_py_JECDn_JERNominal, 0) \
MET_VARIABLE(float, metShift_px_JECUp_JERNominal, 0) \
MET_VARIABLE(float, metShift_py_JECUp_JERNominal, 0) \
MET_VARIABLE(float, metShift_p4Preserved_px_JECDn, 0) \
MET_VARIABLE(float, metShift_p4Preserved_py_JECDn, 0) \
MET_VARIABLE(float, metShift_p4Preserved_px_JECUp, 0) \
MET_VARIABLE(float, metShift_p4Preserved_py_JECUp, 0) \
MET_VARIABLE(float, metShift_p4Preserved_px_JECNominal_JERDn, 0) \
MET_VARIABLE(float, metShift_p4Preserved_py_JECNominal_JERDn, 0) \
MET_VARIABLE(float, metShift_p4Preserved_px_JECNominal_JERUp, 0) \
MET_VARIABLE(float, metShift_p4Preserved_py_JECNominal_JERUp, 0) \
MET_VARIABLE(float, metShift_p4Preserved_px_JECDn_JERNominal, 0) \
MET_VARIABLE(float, metShift_p4Preserved_py_JECDn_JERNominal, 0) \
MET_VARIABLE(float, metShift_p4Preserved_px_JECUp_JERNominal, 0) \
MET_VARIABLE(float, metShift_p4Preserved_py_JECUp_JERNominal, 0)

#define MET_SHIFT_VARIABLES \
MET_SHIFT_CORE_VARIABLES \
MET_SHIFT_GENINFO_VARIABLES

#define MET_CORE_VARIABLES \
MET_RECORDED_CORE_VARIABLES \
MET_SHIFT_CORE_VARIABLES

#define MET_GENINFO_VARIABLES \
MET_RECORDED_GENINFO_VARIABLES \
MET_SHIFT_GENINFO_VARIABLES

#define MET_VARIABLES \
MET_RECORDED_VARIABLES \
MET_SHIFT_VARIABLES


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
  ParticleObject::LorentzVector_t currentMETShift_noJER;
  ParticleObject::LorentzVector_t currentMETShift;
  ParticleObject::LorentzVector_t currentMETShift_p4Preserved_noJER;
  ParticleObject::LorentzVector_t currentMETShift_p4Preserved;
  ParticleObject::LorentzVector_t currentXYshift;
  ParticleObject::LorentzVector_t particleMomentumCorrections;

  std::vector<ParticleObject::LorentzVector_t> currentMETCorrections;
  std::vector<ParticleObject::LorentzVector_t> currentJetOverlapCorrections;

  void setMETShifts();

public:
  METObject();
  METObject(const METObject& other);
  METObject& operator=(const METObject& other);
  ~METObject();

  void swap(METObject& other);

  void setSystematic(SystematicsHelpers::SystematicVariationTypes const&);
  SystematicsHelpers::SystematicVariationTypes const& getCurrentSystematic() const{ return currentSyst; }

  void setXYShift(ParticleObject::LorentzVector_t const& shift){ currentXYshift=shift; }
  void setXYShift(float const& shift_x, float const& shift_y){ this->setXYShift(ParticleObject::LorentzVector_t(shift_x, shift_y, 0., 0.)); }

  void setParticleShifts(ParticleObject::LorentzVector_t const& shift){ particleMomentumCorrections=shift; }
  void setParticleShifts(float const& shift_x, float const& shift_y){ this->setParticleShifts(ParticleObject::LorentzVector_t(shift_x, shift_y, 0., 0.)); }

  void setMETCorrection(ParticleObject::LorentzVector_t const& corr, bool hasXYShifts, bool hasJERShifts, bool addParticleShifts, bool preserveP4);

  void setJetOverlapCorrection(ParticleObject::LorentzVector_t const& corr, bool hasJERShifts, bool preserveP4);

  void getPtPhi(float& pt, float& phi, bool addXYShifts, bool addJERShifts, bool addParticleShifts, bool preserveP4) const;
  ParticleObject::LorentzVector_t::Scalar met(bool addXYShifts, bool addJERShifts, bool addParticleShifts, bool preserveP4) const;
  ParticleObject::LorentzVector_t::Scalar phi(bool addXYShifts, bool addJERShifts, bool addParticleShifts, bool preserveP4) const;
  ParticleObject::LorentzVector_t::Scalar pt(bool addXYShifts, bool addJERShifts, bool addParticleShifts, bool preserveP4) const{ return met(addXYShifts, addJERShifts, addParticleShifts, preserveP4); }
  ParticleObject::LorentzVector_t::Scalar px(bool addXYShifts, bool addJERShifts, bool addParticleShifts, bool preserveP4, float phi_rot=0) const;
  ParticleObject::LorentzVector_t::Scalar py(bool addXYShifts, bool addJERShifts, bool addParticleShifts, bool preserveP4, float phi_rot=0) const;
  ParticleObject::LorentzVector_t p4(bool addXYShifts, bool addJERShifts, bool addParticleShifts, bool preserveP4, float phi_rot=0) const;

};

#endif
