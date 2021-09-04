#ifndef PARTICLEOBJECT_H
#define PARTICLEOBJECT_H

#include "IvyParticle.h"
#include <CMS3/Dictionaries/interface/CommonTypedefs.h>
#include "SystematicVariations.h"


class ParticleObject : public IvyParticle{
public:
  ParticleObject();
  ParticleObject(cms3_id_t id_);
  ParticleObject(cms3_id_t id_, LorentzVector_t const& mom_);
  ParticleObject(const ParticleObject& other);
  virtual ~ParticleObject(){}

  // Swap and assignment operators are not virtual; they bring more complication than necessary, so they are implemented in the derived classes.
  void swap(ParticleObject& other);

  virtual void makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const&){}

  static bool checkDeepDaughtership_NoPFCandidates(IvyParticle const* part1, IvyParticle const* part2);

};

#endif
