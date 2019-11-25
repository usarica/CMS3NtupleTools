#include <algorithm>
#include <utility>
#include "ParticleObject.h"
#include "HelperFunctions.h"
#include "PDGHelpers.h"
#include "TMath.h"


using namespace PDGHelpers;


ParticleObject::ParticleObject() :
  id(-9000),
  selectionBits(0),
  momentum(0, 0, 0, 0)
{}
ParticleObject::ParticleObject(int id_) :
  id(id_),
  selectionBits(0),
  momentum(0, 0, 0, 0)
{}
ParticleObject::ParticleObject(int id_, LorentzVector_t const& momentum_) :
  id(id_),
  selectionBits(0),
  momentum(momentum_)
{}
ParticleObject::ParticleObject(const ParticleObject& other) :
  id(other.id),
  selectionBits(other.selectionBits),
  momentum(other.momentum)
{}

void ParticleObject::setSelectionBit(unsigned int ibit){ HelperFunctions::set_bit(this->selectionBits, ibit); }
bool ParticleObject::testSelection(unsigned int ibit) const{ return HelperFunctions::test_bit(this->selectionBits, ibit); }

float ParticleObject::charge()const{
  float cpos=0;
  if (isAWBoson(id) || abs(id)==37 || abs(id)==2212 || abs(id)==211 || abs(id)==321 || abs(id)==411 || abs(id)==521) cpos = 1.;
  else if (isALepton(id)) cpos = -1.;
  else if (isUpTypeQuark(id)) cpos = 2./3.;
  else if (isDownTypeQuark(id)) cpos = -1./3.;
  if (id<0) cpos *= -1.;
  return cpos;
}
float ParticleObject::deltaPhi(float phi_) const{
  float dPhi = phi_-phi();
  if (dPhi>TMath::Pi()) dPhi = 2.*TMath::Pi() - dPhi;
  else if (dPhi<-TMath::Pi()) dPhi = 2.*TMath::Pi() + dPhi;
  return dPhi;
}
