#ifndef TRIGGERINFO_H
#define TRIGGERINFO_H

#include <string>


struct TriggerInfo{
  std::string name;
  bool passTrigger;
  int HLTprescale;
  int L1prescale;

  TriggerInfo();
  TriggerInfo(std::string const&, bool const&, int const&, int const&);
  TriggerInfo(TriggerInfo const&);

};


#endif
