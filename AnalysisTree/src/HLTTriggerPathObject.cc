#include <algorithm>
#include <utility>
#include <cmath>
#include "HLTTriggerPathObject.h"


HLTTriggerPathObject::HLTTriggerPathObject(){
#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=DEFVAL;
  HLTTRIGGERPATH_VARIABLES;
#undef HLTTRIGGERPATH_VARIABLE
}
HLTTriggerPathObject::HLTTriggerPathObject(HLTTriggerPathObject const& other){
#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) this->NAME=other.NAME;
  HLTTRIGGERPATH_VARIABLES;
#undef HLTTRIGGERPATH_VARIABLE
}
void HLTTriggerPathObject::swap(HLTTriggerPathObject& other){
#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) std::swap(this->NAME, other.NAME);
  HLTTRIGGERPATH_VARIABLES;
#undef HLTTRIGGERPATH_VARIABLE
}
HLTTriggerPathObject& HLTTriggerPathObject::operator=(const HLTTriggerPathObject& other){
  HLTTriggerPathObject tmp(other);
  swap(tmp);
  return *this;
}
