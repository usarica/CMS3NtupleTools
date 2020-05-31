#include <cassert>
#include <cmath>
#include "HLTObjectProperties.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


HLTObjectProperties::HLTObjectProperties() :
  type(HLTObjectProperties::nTriggerObjectTypes)
{
  resetCuts();
}
HLTObjectProperties::HLTObjectProperties(HLTObjectProperties::TriggerObjectType const& type_) :
  type(type_)
{
  resetCuts();
}
HLTObjectProperties::HLTObjectProperties(HLTObjectProperties::TriggerObjectType const& type_, std::vector< std::pair<HLTObjectProperties::CutType, float> > const& cuts_) :
  type(type_)
{
  resetCuts();

  for (auto const& type_val_pair:cuts_){
    float const& cut_val = type_val_pair.second;
    switch (type_val_pair.first){
    case kPt:
      setPtCut(cut_val);
      break;
    case kPtHigh:
      setPtHighCut(cut_val);
      break;
    case kEta:
      setEtaCut(cut_val);
      break;
    case kMass:
      setMassCut(cut_val);
      break;
    default:
      MELAerr << "HLTObjectProperties::HLTObjectProperties: Type " << type_val_pair.first << " is not defined." << endl;
      assert(0);
      break;
    }
  }
}
HLTObjectProperties::HLTObjectProperties(HLTObjectProperties const& other) :
  type(other.type),

  pt_cut(other.pt_cut),
  pt_high_cut(other.pt_high_cut),
  eta_cut(other.eta_cut),
  mass_cut(other.mass_cut),

  has_pt_cut(other.has_pt_cut),
  has_pt_high_cut(other.has_pt_high_cut),
  has_eta_cut(other.has_eta_cut),
  has_mass_cut(other.has_mass_cut)
{}

void HLTObjectProperties::resetCuts(){
  pt_cut = -1;
  pt_high_cut = -1;
  eta_cut = 10;
  mass_cut = -1;

  has_pt_cut = false;
  has_pt_high_cut = false;
  has_eta_cut = false;
  has_mass_cut = false;
}

bool HLTObjectProperties::testCuts(ParticleObject::LorentzVector_t const& p4, HLTObjectProperties::TriggerObjectType const& type_) const{
  if (type != type_ || type == HLTObjectProperties::nTriggerObjectTypes || type_ == HLTObjectProperties::nTriggerObjectTypes) return true;

  bool res = true;

  res &= (!has_pt_cut || p4.Pt()>=pt_cut);
  res &= (!has_pt_high_cut || p4.Pt()<pt_high_cut);
  res &= (!has_eta_cut || std::abs(p4.Eta())<eta_cut);
  res &= (!has_mass_cut || p4.M()>=mass_cut);

  return res;
}

bool HLTObjectProperties::operator == (HLTObjectProperties const& other) const{
  return (
    type == other.type
    &&
    pt_cut == other.pt_cut
    &&
    pt_high_cut == other.pt_high_cut
    &&
    eta_cut == other.eta_cut
    &&
    mass_cut == other.mass_cut
    );
}
bool HLTObjectProperties::operator > (HLTObjectProperties const& other) const{
  // Do not require the upper pt cut in these comparisons
  return (
    type == other.type
    && (
      pt_cut>other.pt_cut
      ||
      (pt_cut==other.pt_cut && eta_cut<other.eta_cut)
      ||
      (pt_cut==other.pt_cut && eta_cut==other.eta_cut && mass_cut>other.mass_cut)
      )
    );
}
