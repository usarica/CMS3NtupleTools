#ifndef TRIGGEROBJECTINFO_H
#define TRIGGEROBJECTINFO_H

#define TRIGGEROBJECTINFO_INDEX_BY_ORIGINAL 0

#include <string>
#include <CMSDataTools/AnalysisTree/interface/CMSLorentzVector.h>


struct TriggerObjectInfo{
  std::vector<unsigned int> triggerCollectionIndices;
  std::vector<bool> passAllTriggerFiltersList;

  size_t triggerObjectCollectionIndex;
  std::vector<int> types;
  CMSLorentzVector_d p4;

  TriggerObjectInfo();
  TriggerObjectInfo(size_t const&, std::vector<int> const&, CMSLorentzVector_d const&);
  TriggerObjectInfo(TriggerObjectInfo const&);

  void addTriggerCollectionIndex(unsigned int const&, bool const&);

  int bestType() const;

};


#endif
