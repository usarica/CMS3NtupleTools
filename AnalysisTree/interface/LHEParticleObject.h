#ifndef LHEPARTICLEOBJECT_H
#define LHEPARTICLEOBJECT_H

#include "ParticleObject.h"


#define LHEPARTICLE_ID_MOMENTUM_VARIABLES \
LHEPARTICLE_VARIABLE(cms3_id_t, id, 0) \
LHEPARTICLE_VARIABLE(cms3_genstatus_t, status, 0) \
LHEPARTICLE_VARIABLE(float, px, 0) \
LHEPARTICLE_VARIABLE(float, py, 0) \
LHEPARTICLE_VARIABLE(float, pz, 0) \
LHEPARTICLE_VARIABLE(float, E, 0)

#define LHEPARTICLE_MOTHER_VARIABLES \
LHEPARTICLE_VARIABLE(int, mother0_index, 0) \
LHEPARTICLE_VARIABLE(int, mother1_index, 0)

#define LHEPARTICLE_VARIABLES \
LHEPARTICLE_ID_MOMENTUM_VARIABLES \
LHEPARTICLE_MOTHER_VARIABLES


class LHEParticleObject : public ParticleObject{
protected:
  cms3_genstatus_t st;

public:
  LHEParticleObject();
  LHEParticleObject(cms3_id_t id_, cms3_genstatus_t st_);
  LHEParticleObject(cms3_id_t id_, cms3_genstatus_t st_, LorentzVector_t const& mom_);
  LHEParticleObject(const LHEParticleObject& other);
  LHEParticleObject& operator=(const LHEParticleObject& other);
  ~LHEParticleObject();

  void swap(LHEParticleObject& other);

  cms3_genstatus_t& status(){ return this->st; }
  cms3_genstatus_t const& status() const{ return this->st; }

};

#endif
