#include <cassert>
#include "HelperFunctions.h"
#include "HLTTriggerPathProperties.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


HLTTriggerPathProperties::HLTTriggerPathProperties(){}
HLTTriggerPathProperties::HLTTriggerPathProperties(std::string const& name_) :
  name(name_)
{
  setupName();
  setup();
}
HLTTriggerPathProperties::HLTTriggerPathProperties(std::string const& name_, std::vector<HLTObjectProperties> const& triggerObjectProperties_) :
  name(name_)
{
  setupName();
  for (auto const& props:triggerObjectProperties_) addObjectProperties(props);
  setup();
}
HLTTriggerPathProperties::HLTTriggerPathProperties(HLTTriggerPathProperties const& other) :
  name(other.name),
  triggerObjectProperties(other.triggerObjectProperties)
{}

void HLTTriggerPathProperties::setupName(){
  size_t pos = name.find("*");
  if (pos!=std::string::npos){
    if (pos != name.length()-1){
      MELAerr << "HLTTriggerPathProperties::HLTTriggerPathProperties: Trigger name " << name << " can only contain * as the last character. Please fix the name passed!" << endl;
      assert(0);
    }
    name = name.substr(0, pos);
    pos = pos-1;
    if (name.at(pos)!='v' || name.at(pos-1)!='_'){
      MELAerr << "HLTTriggerPathProperties::HLTTriggerPathProperties: Trigger name " << name << " has to end with _v. Please fix the name passed!" << endl;
      assert(0);
    }
  }
}

void HLTTriggerPathProperties::setup(){
  for (auto& it:triggerObjectProperties){
    auto& props = it.second;
    HLTObjectProperties::sortByMoreRestrictive(props);
  }
}

void HLTTriggerPathProperties::addObjectProperties(HLTObjectProperties const& props){
  auto const& type=props.getType();
  if (triggerObjectProperties.find(type)==triggerObjectProperties.end()) triggerObjectProperties[type]=std::vector<HLTObjectProperties>();
  triggerObjectProperties[type].push_back(props);
}

bool HLTTriggerPathProperties::isSameTrigger(std::string const& name_) const{
  return (name_.find(name)!=std::string::npos);
}

void HLTTriggerPathProperties::resetCuts(){
  for (auto& it:triggerObjectProperties){
    for (auto& props:it.second) props.resetCuts();
  }
}

bool HLTTriggerPathProperties::testCuts(
  std::vector<MuonObject const*> const& muons,
  std::vector<ElectronObject const*> const& electrons,
  std::vector<PhotonObject const*> const& photons,
  std::vector<AK4JetObject const*> const& ak4jets,
  std::vector<AK8JetObject const*> const& ak8jets,
  ParticleObject::LorentzVector_t const& pfmet_p4,
  ParticleObject::LorentzVector_t const& pfmet_nomu_p4,
  ParticleObject::LorentzVector_t const& ht_p4,
  ParticleObject::LorentzVector_t const& ht_nomu_p4
) const{
  bool res = true;
  for (auto const& it:triggerObjectProperties){
    auto const& prop_type = it.first;
    auto const& props = it.second;
    switch (prop_type){
    case HLTObjectProperties::kMuon:
      res = res && testCutSet(prop_type, props, muons);
      break;
    case HLTObjectProperties::kElectron:
      res = res && testCutSet(prop_type, props, electrons);
      break;
    case HLTObjectProperties::kPhoton:
      res = res && testCutSet(prop_type, props, photons);
      break;
    case HLTObjectProperties::kAK4Jet:
      res = res && testCutSet(prop_type, props, ak4jets);
      break;
    case HLTObjectProperties::kAK8Jet:
      res = res && testCutSet(prop_type, props, ak8jets);
      break;
    case HLTObjectProperties::kHT:
      for (auto const& propcuts:props) res = res && propcuts.testCuts(ht_p4, HLTObjectProperties::kHT);
      break;
    case HLTObjectProperties::kHT_NoMu:
      for (auto const& propcuts:props) res = res && propcuts.testCuts(ht_nomu_p4, HLTObjectProperties::kHT_NoMu);
      break;
    case HLTObjectProperties::kMET:
      for (auto const& propcuts:props) res = res && propcuts.testCuts(pfmet_p4, HLTObjectProperties::kMET);
      break;
    case HLTObjectProperties::kMET_NoMu:
      for (auto const& propcuts:props) res = res && propcuts.testCuts(pfmet_nomu_p4, HLTObjectProperties::kMET_NoMu);
      break;
    default:
      break;
    }
  }
  return res;
}
