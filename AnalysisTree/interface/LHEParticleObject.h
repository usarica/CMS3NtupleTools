#ifndef LHEPARTICLEOBJECT_H
#define LHEPARTICLEOBJECT_H

#include "ParticleObject.h"


#define LHEPARTICLE_ID_MOMENTUM_VARIABLES \
LHEPARTICLE_VARIABLE(int, id, 0) \
LHEPARTICLE_VARIABLE(float, px, 0) \
LHEPARTICLE_VARIABLE(float, py, 0) \
LHEPARTICLE_VARIABLE(float, pz, 0) \
LHEPARTICLE_VARIABLE(float, E, 0)

#define LHEPARTICLE_EXTRA_VARIABLES \
LHEPARTICLE_VARIABLE(int, status, 0) \
LHEPARTICLE_VARIABLE(int, mother0_index, 0) \
LHEPARTICLE_VARIABLE(int, mother1_index, 0)

#define LHEPARTICLE_VARIABLES \
LHEPARTICLE_ID_MOMENTUM_VARIABLES \
LHEPARTICLE_EXTRA_VARIABLES


class LHEParticleObject : public ParticleObject{
protected:
  int st;

public:
  LHEParticleObject();
  LHEParticleObject(int id_, int st_);
  LHEParticleObject(int id_, int st_, LorentzVector_t const& mom_);
  LHEParticleObject(const LHEParticleObject& other);
  LHEParticleObject& operator=(const LHEParticleObject& other);
  ~LHEParticleObject();

  void swap(LHEParticleObject& other);

  int& status(){ return this->st; }
  int const& status() const{ return this->st; }

};

#endif
