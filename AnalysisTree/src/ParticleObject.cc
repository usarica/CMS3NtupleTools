#include <algorithm>
#include <utility>
#include "ParticleObject.h"
#include "HelperFunctions.h"
#include "PDGHelpers.h"
#include "TMath.h"


using namespace PDGHelpers;


ParticleObject::ParticleObject() :
  id(-9000),
  uniqueIdentifier(0),
  selectionBits(0),
  momentum(0, 0, 0, 0)
{}
ParticleObject::ParticleObject(cms3_id_t id_) :
  id(id_),
  uniqueIdentifier(0),
  selectionBits(0),
  momentum(0, 0, 0, 0)
{}
ParticleObject::ParticleObject(cms3_id_t id_, LorentzVector_t const& momentum_) :
  id(id_),
  uniqueIdentifier(0),
  selectionBits(0),
  momentum(momentum_)
{}
ParticleObject::ParticleObject(const ParticleObject& other) :
  id(other.id),
  uniqueIdentifier(other.uniqueIdentifier),
  selectionBits(other.selectionBits),
  momentum(other.momentum),
  mothers(other.mothers),
  daughters(other.daughters)
{}

void ParticleObject::swap(ParticleObject& other){
  std::swap(id, other.id);
  std::swap(uniqueIdentifier, other.uniqueIdentifier);
  std::swap(selectionBits, other.selectionBits);
  std::swap(momentum, other.momentum);
  std::swap(mothers, other.mothers);
  std::swap(daughters, other.daughters);
}

void ParticleObject::setSelectionBit(unsigned int ibit, bool val){ HelperFunctions::set_bit(this->selectionBits, ibit, val); }
bool ParticleObject::testSelectionBit(unsigned int ibit) const{ return HelperFunctions::test_bit(this->selectionBits, ibit); }

float ParticleObject::charge() const{
  float cpos=0;
  if (isAWBoson(id) || abs(id)==37 || abs(id)==2212 || abs(id)==211 || abs(id)==321 || abs(id)==411 || abs(id)==521) cpos = 1.;
  else if (isALepton(id)) cpos = -1.;
  else if (isUpTypeQuark(id)) cpos = 2./3.;
  else if (isDownTypeQuark(id)) cpos = -1./3.;
  if (id<0) cpos *= -1.;
  return cpos;
}
ParticleObject::LorentzVector_t::Scalar ParticleObject::deltaR(LorentzVector_t::Scalar eta_, LorentzVector_t::Scalar phi_) const{
  ParticleObject::LorentzVector_t::Scalar res;
  HelperFunctions::deltaR(eta(), phi(), eta_, phi_, res);
  return res;
}
ParticleObject::LorentzVector_t::Scalar ParticleObject::deltaEta(LorentzVector_t::Scalar eta_) const{
  ParticleObject::LorentzVector_t::Scalar res;
  HelperFunctions::deltaEta(eta(), eta_, res);
  return res;
}
ParticleObject::LorentzVector_t::Scalar ParticleObject::deltaPhi(LorentzVector_t::Scalar phi_) const{
  ParticleObject::LorentzVector_t::Scalar res;
  HelperFunctions::deltaPhi(phi(), phi_, res);
  return res;
}

bool ParticleObject::checkParticleExists(ParticleObject* myParticle, std::vector<ParticleObject*> const& particleArray){ return HelperFunctions::checkListVariable(particleArray, myParticle); }
bool ParticleObject::checkDeepDaughtership(ParticleObject const* part1, ParticleObject const* part2){
  if (!part1 || !part2) return false;
  if (part1 == part2) return true;
  std::vector<ParticleObject*> const& daughters1 = part1->getDaughters();
  std::vector<ParticleObject*> const& daughters2 = part2->getDaughters();
  return HelperFunctions::hasCommonElements(daughters1, daughters2);
}

void ParticleObject::addMother(ParticleObject* myParticle){ if (!checkParticleExists(myParticle, mothers)) mothers.push_back(myParticle); }
void ParticleObject::addDaughter(ParticleObject* myParticle){ if (!checkParticleExists(myParticle, daughters)) daughters.push_back(myParticle); }

bool ParticleObject::hasMother(ParticleObject* part) const{ return checkParticleExists(part, mothers); }
bool ParticleObject::hasDaughter(ParticleObject* part) const{ return checkParticleExists(part, daughters); }

void ParticleObject::getDeepDaughters(std::vector<ParticleObject*>& deepdaus, bool addSelf){
  if (this->getNDaughters()==0){
    if (addSelf && !HelperFunctions::checkListVariable(deepdaus, (ParticleObject*) this)) deepdaus.push_back(this);
  }
  else{
    for (auto const& dau:daughters) dau->getDeepDaughters(deepdaus, true);
  }
}
void ParticleObject::getDeepDaughters(std::vector<ParticleObject const*>& deepdaus, bool addSelf) const{
  if (this->getNDaughters()==0){
    if (addSelf && !HelperFunctions::checkListVariable(deepdaus, (ParticleObject const*) this)) deepdaus.push_back(this);
  }
  else{
    for (auto const& dau:daughters) dau->getDeepDaughters(deepdaus, true);
  }
}
