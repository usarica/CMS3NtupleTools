#ifndef TRIGGEROBJECTINFO_H
#define TRIGGEROBJECTINFO_H

#define TRIGGEROBJECTINFO_INDEX_BY_ORIGINAL 0

#include <string>
#include <CMSDataTools/AnalysisTree/interface/CMSLorentzVector.h>
#include <CMS3/Dictionaries/interface/CommonTypedefs.h>


struct TriggerObjectInfo{
  std::vector<cms3_triggerIndex_t> triggerCollectionIndices;
  std::vector<bool> passAllTriggerFiltersList;

  size_t triggerObjectCollectionIndex;
  std::vector<int> types;
  CMSLorentzVector_d p4;

  TriggerObjectInfo();
  TriggerObjectInfo(size_t const&, std::vector<int> const&, CMSLorentzVector_d const&);
  TriggerObjectInfo(TriggerObjectInfo const&);

  void addTriggerCollectionIndex(cms3_triggerIndex_t const&, bool const&);

  int bestType() const;

};


#endif
