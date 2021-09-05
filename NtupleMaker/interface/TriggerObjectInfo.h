#ifndef TRIGGEROBJECTINFO_H
#define TRIGGEROBJECTINFO_H

#define TRIGGEROBJECTINFO_INDEX_BY_ORIGINAL 0

#include <string>
#include <IvyFramework/IvyDataTools/interface/CMSLorentzVector.h>
#include <CMS3/Dictionaries/interface/CommonTypedefs.h>


struct TriggerObjectInfo{
  std::vector<cms3_triggerIndex_t> triggerCollectionIndices;
  std::vector<bool> passAllTriggerFiltersList;

  size_t triggerObjectCollectionIndex;
  std::vector<cms3_triggertype_t> types;
  CMSLorentzVector_d p4;

  TriggerObjectInfo();
  TriggerObjectInfo(size_t const&, std::vector<cms3_triggertype_t> const&, CMSLorentzVector_d const&);
  TriggerObjectInfo(TriggerObjectInfo const&);

  void addTriggerCollectionIndex(cms3_triggerIndex_t const&, bool const&);

  cms3_triggertype_t bestType() const;

};


#endif
