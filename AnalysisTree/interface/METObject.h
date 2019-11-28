#ifndef METOBJECT_H
#define METOBJECT_H

#include "SystematicVariations.h"


#define MET_RECORDED_VARIABLES \
MET_VARIABLE(float, met, 0) \
MET_VARIABLE(float, phi, 0) \
MET_VARIABLE(float, met_JECup, 0) \
MET_VARIABLE(float, phi_JECup, 0) \
MET_VARIABLE(float, met_JECdn, 0) \
MET_VARIABLE(float, phi_JECdn, 0) \
MET_VARIABLE(float, metSignificance, 0)

#define MET_EXTRA_PT_VARIABLES \
MET_VARIABLE(float, met_original, 0) \
MET_VARIABLE(float, met_METup, 0) \
MET_VARIABLE(float, met_METdn, 0)

#define MET_EXTRA_PHI_VARIABLES \
MET_VARIABLE(float, phi_original, 0) \
MET_VARIABLE(float, phi_METup, 0) \
MET_VARIABLE(float, phi_METdn, 0)

#define MET_EXTRA_VARIABLES \
MET_EXTRA_PT_VARIABLES \
MET_EXTRA_PHI_VARIABLES

#define MET_VARIABLES \
MET_RECORDED_VARIABLES \
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

public:
  METObject();
  METObject(const METObject& other);
  METObject& operator=(const METObject& other);
  ~METObject();

  void swap(METObject& other);

  void setSystematic(SystematicsHelpers::SystematicVariationTypes const&);

  void getPtPhi(float const*& pt, float const*& phi) const;
  float const& met() const;
  float const& phi() const;
  float const& pt() const{ return met(); }
  float px(float phi_rot=0) const;
  float py(float phi_rot=0) const;

};

#endif
