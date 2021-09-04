#include <algorithm>
#include <utility>
#include "ParticleObject.h"
#include "PFCandidateObject.h"
#include "HelperFunctions.h"
#include "PDGHelpers.h"
#include "TMath.h"


using namespace PDGHelpers;


ParticleObject::ParticleObject() :
  IvyParticle()
{}
ParticleObject::ParticleObject(cms3_id_t id_) :
  IvyParticle(id_)
{}
ParticleObject::ParticleObject(cms3_id_t id_, LorentzVector_t const& momentum_) :
  IvyParticle(id_, momentum_)
{}
ParticleObject::ParticleObject(const ParticleObject& other) :
  IvyParticle(other)
{}

void ParticleObject::swap(ParticleObject& other){ IvyParticle::swap(other); }

bool ParticleObject::checkDeepDaughtership_NoPFCandidates(IvyParticle const* part1, IvyParticle const* part2){
  if (!part1 || !part2) return false;
  if (part1 == part2) return true;
  std::vector<IvyParticle*> const& daughters1 = part1->getDaughters();
  std::vector<IvyParticle*> daughters1_filtered; daughters1_filtered.reserve(daughters1.size());
  for (auto const& dau:daughters1){ if (dynamic_cast<PFCandidateObject*>(dau)==nullptr) daughters1_filtered.push_back(dau); }
  std::vector<IvyParticle*> const& daughters2 = part2->getDaughters();
  std::vector<IvyParticle*> daughters2_filtered; daughters2_filtered.reserve(daughters2.size());
  for (auto const& dau:daughters2){ if (dynamic_cast<PFCandidateObject*>(dau)==nullptr) daughters2_filtered.push_back(dau); }
  return HelperFunctions::hasCommonElements(daughters1_filtered, daughters2_filtered);
}
