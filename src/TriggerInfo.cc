#include "CMS3/NtupleMaker/interface/TriggerInfo.h"


TriggerInfo::TriggerInfo() :
  name(""),
  passTrigger(false),
  HLTprescale(1),
  L1prescale(1)
{}

TriggerInfo::TriggerInfo(std::string const& name_, bool const& passTrigger_, int const& HLTprescale_, int const& L1prescale_) :
  name(name_),
  passTrigger(passTrigger_),
  HLTprescale(HLTprescale_),
  L1prescale(L1prescale_)
{}

TriggerInfo::TriggerInfo(TriggerInfo const& other) :
  name(other.name),
  passTrigger(other.passTrigger),
  HLTprescale(other.HLTprescale),
  L1prescale(other.L1prescale)
{}
