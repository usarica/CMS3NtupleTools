#ifndef PARTICLEOBJECT_H
#define PARTICLEOBJECT_H

#include <DataFormats/Math/interface/deltaR.h>
#include <Math/GenVector/DisplacementVector2D.h>
#include <DataFormats/Math/interface/Vector3D.h>

#include "CMSLorentzVector.h"
#include "SystematicVariations.h"
#include "TLorentzVector.h"


class ParticleObject{
public:
  typedef CMSLorentzVector_d LorentzVector_t;
  typedef math::XYZVectorD Vector3D_t;
  typedef ROOT::Math::DisplacementVector2D< ROOT::Math::Cartesian2D<double> > Vector2D_t;

  int id;
  unsigned long long selectionBits;
  LorentzVector_t momentum;

  ParticleObject();
  ParticleObject(int id_);
  ParticleObject(int id_, LorentzVector_t const& mom_);
  ParticleObject(const ParticleObject& other);
  virtual ~ParticleObject(){}

  // Swap and assignment operators are not virtual; they bring more complication than necessary, so they are implemented in the derived classes.

  void setPdgId(int const& id_){ id=id_; }
  void setP4(LorentzVector_t const& momentum_){ momentum=momentum_; }

  void resetSelectionBits(){ selectionBits=0; }
  void setSelectionBit(unsigned int ibit);
  bool testSelection(unsigned int ibit) const;

  int const& pdgId() const{ return id; }
  int& pdgId(){ return id; }

  LorentzVector_t const& p4() const{ return momentum; }
  LorentzVector_t& p4(){ return momentum; }
  virtual void makeFinalMomentum(SystematicsHelpers::SystematicVariationTypes const&) = 0;

  unsigned long long const& getSelectionBits() const{ return selectionBits; }
  unsigned long long& getSelectionBits(){ return selectionBits; }

  float charge() const;
  float m() const{ return momentum.M(); }
  float x() const{ return momentum.X(); }
  float y() const{ return momentum.Y(); }
  float z() const{ return momentum.Z(); }
  float t() const{ return momentum.T(); }
  float energy() const{ return this->t(); }
  float p() const{ return momentum.P(); }
  float pt() const{ return momentum.Pt(); }
  float eta() const{ return momentum.Eta(); }
  float phi() const{ return momentum.Phi(); }
  float rapidity() const{ return momentum.Rapidity(); }
  float dot(const TLorentzVector& v) const{ return (momentum.T()*v.T()-(momentum.X()*v.X()+momentum.Y()*v.Y()+momentum.Z()*v.Z())); }
  float dot(const LorentzVector_t& v) const{ return (momentum.T()*v.T()-(momentum.X()*v.X()+momentum.Y()*v.Y()+momentum.Z()*v.Z())); }
  float dot(const ParticleObject& part) const{ return dot(part.momentum); }
  float dot(const ParticleObject* part) const{ if (part!=0) return dot(*part); else return 0; }
  float deltaR(const TLorentzVector& v) const{ TLorentzVector tmp(momentum.X(), momentum.Y(), momentum.Z(), momentum.T()); return tmp.DeltaR(v); }
  float deltaR(const LorentzVector_t& v) const{ return reco::deltaR(momentum, v); }
  float deltaR(const ParticleObject& part) const{ return deltaR(part.momentum); }
  float deltaR(const ParticleObject* part) const{ if (part) return deltaR(*part); else return -1; }
  float deltaPhi(float phi_) const;

  Vector3D_t vect() const{ return Vector3D_t(momentum.X(), momentum.Y(), momentum.Z()); }
  Vector2D_t vect_trans() const{ return Vector2D_t(momentum.X(), momentum.Y()); }

};

#endif
