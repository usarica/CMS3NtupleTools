#include <cassert>
#include <algorithm>
#include <utility>
#include <DataFormats/HLTReco/interface/TriggerTypeDefs.h>
#include <CMS3/NtupleMaker/interface/TriggerObjectInfo.h>


TriggerObjectInfo::TriggerObjectInfo() :
  triggerObjectCollectionIndex(0),
  p4(0, 0, 0, 0)
{}

TriggerObjectInfo::TriggerObjectInfo(size_t const& triggerObjectCollectionIndex_, std::vector<cms3_triggertype_t> const& types_, CMSLorentzVector_d const& p4_) :
  triggerObjectCollectionIndex(triggerObjectCollectionIndex_),
  types(types_),
  p4(p4_)
{}

TriggerObjectInfo::TriggerObjectInfo(TriggerObjectInfo const& other) :
  triggerCollectionIndices(other.triggerCollectionIndices),
  triggerObjectCollectionIndex(other.triggerObjectCollectionIndex),
  types(other.types),
  p4(other.p4)
{}

void TriggerObjectInfo::addTriggerCollectionIndex(cms3_triggerIndex_t const& index, bool const& passedAllFilters){
  triggerCollectionIndices.push_back(index);
  passAllTriggerFiltersList.push_back(passedAllFilters);
}

cms3_triggertype_t TriggerObjectInfo::bestType() const{
  cms3_triggertype_t const i_muon = static_cast<cms3_triggertype_t>(trigger::TriggerMuon);
  cms3_triggertype_t const i_electron = static_cast<cms3_triggertype_t>(trigger::TriggerElectron);
  cms3_triggertype_t const i_photon = static_cast<cms3_triggertype_t>(trigger::TriggerPhoton);
  cms3_triggertype_t const i_cluster = static_cast<cms3_triggertype_t>(trigger::TriggerCluster);
  cms3_triggertype_t const i_bjet = static_cast<cms3_triggertype_t>(trigger::TriggerBJet);

  if (std::find(types.cbegin(), types.cend(), i_muon)!=types.cend()) return i_muon;
  else if (std::find(types.cbegin(), types.cend(), i_electron)!=types.cend()) return i_electron;
  else if (std::find(types.cbegin(), types.cend(), i_photon)!=types.cend()) return i_photon;
  else if (std::find(types.cbegin(), types.cend(), i_cluster)!=types.cend()) return i_cluster;
  else if (std::find(types.cbegin(), types.cend(), i_bjet)!=types.cend()) return i_bjet;
  else return types.front();
}
