#ifndef GENPARTICLEOBJECT_H
#define GENPARTICLEOBJECT_H

#include "ParticleObject.h"


#define GENPARTICLE_ID_MOMENTUM_VARIABLES \
GENPARTICLE_VARIABLE(cms3_id_t, id, 0) \
GENPARTICLE_VARIABLE(cms3_genstatus_t, status, 0) \
GENPARTICLE_VARIABLE(float, pt, 0) \
GENPARTICLE_VARIABLE(float, eta, 0) \
GENPARTICLE_VARIABLE(float, phi, 0) \
GENPARTICLE_VARIABLE(float, mass, 0)

#define GENPARTICLE_MOTHER_VARIABLES \
GENPARTICLE_VARIABLE(int, mom0_index, 0) \
GENPARTICLE_VARIABLE(int, mom1_index, 0)

#define GENPARTICLE_EXTRA_VARIABLES \
GENPARTICLE_VARIABLE(bool, is_packed, false) \
GENPARTICLE_VARIABLE(bool, isPromptFinalState, false) \
GENPARTICLE_VARIABLE(bool, isPromptDecayed, false) \
GENPARTICLE_VARIABLE(bool, isDirectPromptTauDecayProductFinalState, false) \
GENPARTICLE_VARIABLE(bool, isHardProcess, false) \
GENPARTICLE_VARIABLE(bool, fromHardProcessFinalState, false) \
GENPARTICLE_VARIABLE(bool, fromHardProcessDecayed, false) \
GENPARTICLE_VARIABLE(bool, isDirectHardProcessTauDecayProductFinalState, false) \
GENPARTICLE_VARIABLE(bool, fromHardProcessBeforeFSR, false) \
GENPARTICLE_VARIABLE(bool, isLastCopy, false) \
GENPARTICLE_VARIABLE(bool, isLastCopyBeforeFSR, false)

#define GENPARTICLE_VARIABLES \
GENPARTICLE_ID_MOMENTUM_VARIABLES \
GENPARTICLE_MOTHER_VARIABLES \
GENPARTICLE_EXTRA_VARIABLES


class GenParticleVariables{
public:
#define GENPARTICLE_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME;
  GENPARTICLE_EXTRA_VARIABLES;
#undef GENPARTICLE_VARIABLE

  GenParticleVariables();
  GenParticleVariables(GenParticleVariables const& other);
  GenParticleVariables& operator=(const GenParticleVariables& other);

  void swap(GenParticleVariables& other);

};

class GenParticleObject : public ParticleObject{
protected:
  cms3_genstatus_t st;

public:
  GenParticleVariables extras;

  GenParticleObject();
  GenParticleObject(cms3_id_t id_, cms3_genstatus_t st_);
  GenParticleObject(cms3_id_t id_, cms3_genstatus_t st_, LorentzVector_t const& mom_);
  GenParticleObject(const GenParticleObject& other);
  GenParticleObject& operator=(const GenParticleObject& other);
  ~GenParticleObject();

  void swap(GenParticleObject& other);

  cms3_genstatus_t& status(){ return this->st; }
  cms3_genstatus_t const& status() const{ return this->st; }

};

#endif
