#ifndef TRIGGERINFO_H
#define TRIGGERINFO_H

#include <string>


struct TriggerInfo{
  std::string name; // Trigger name, including version
  size_t index; // Index at the reference collection
  bool passTrigger; // Pass or fail
  int HLTprescale; // Usually 1
  int L1prescale; // Usually 1

  TriggerInfo();
  TriggerInfo(std::string const&, size_t const&, bool const&, int const&, int const&);
  TriggerInfo(TriggerInfo const&);

};


#endif
