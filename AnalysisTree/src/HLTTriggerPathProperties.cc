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
  name(name_),
  triggerObjectProperties(triggerObjectProperties_)
{
  setupName();
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
  HLTObjectProperties::sortByMoreRestrictive(triggerObjectProperties);
}

void HLTTriggerPathProperties::addObjectProperties(HLTObjectProperties const& props){
  triggerObjectProperties.push_back(props);
}

bool HLTTriggerPathProperties::isSameTrigger(std::string const& name_) const{
  return (name_.find(name)!=std::string::npos);
}

bool HLTTriggerPathProperties::testCuts(ParticleObject::LorentzVector_t const& p4, HLTObjectProperties::TriggerObjectType const& type_) const{
  bool res=true;
  for (auto const& propcuts:triggerObjectProperties){
    if (propcuts.hasSameType(type_)) res &= propcuts.testCuts(p4, type_);
  }
  return res;
}
