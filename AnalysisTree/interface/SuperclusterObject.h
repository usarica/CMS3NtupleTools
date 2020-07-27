#ifndef SUPERCLUSTEROBJECT_H
#define SUPERCLUSTEROBJECT_H

#include "ParticleObject.h"


#define SUPERCLUSTER_COMMON_VARIABLES \
SUPERCLUSTER_VARIABLE(float, energy, 0) \
SUPERCLUSTER_VARIABLE(float, sinTheta_SC_pos, 0) \
SUPERCLUSTER_VARIABLE(cms3_egamma_fid_type_mask_t, fid_mask, 0)

#define SUPERCLUSTER_VARIABLES \
SUPERCLUSTER_COMMON_VARIABLES \

class SuperclusterVariables{
public:
#define SUPERCLUSTER_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  SUPERCLUSTER_VARIABLES;
#undef SUPERCLUSTER_VARIABLE

  SuperclusterVariables();
  SuperclusterVariables(SuperclusterVariables const& other);
  SuperclusterVariables& operator=(const SuperclusterVariables& other);

  void swap(SuperclusterVariables& other);

};

class SuperclusterObject : public ParticleObject{
public:
  SuperclusterVariables extras;

  SuperclusterObject();
  SuperclusterObject(LorentzVector_t const& mom_);
  SuperclusterObject(const SuperclusterObject& other);
  SuperclusterObject& operator=(const SuperclusterObject& other);
  ~SuperclusterObject();

  void swap(SuperclusterObject& other);

  float etaSC() const{ return this->eta(); }

  bool isEBEEGap() const;
  bool isEB() const;
  bool isEE() const;
  bool isAnyGap() const;
  bool isGap() const{ return this->isAnyGap(); }

};

#endif
