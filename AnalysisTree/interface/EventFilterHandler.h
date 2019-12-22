#ifndef EVENTFILTERHANDLER_H
#define EVENTFILTERHANDLER_H

#include <vector>
#include <unordered_map>
#include "IvyBase.h"
#include "ElectronObject.h"
#include "PhotonObject.h"
#include "AK4JetObject.h"
#include "AK8JetObject.h"
#include "HLTTriggerPathObject.h"
#include "SystematicVariations.h"


class EventFilterHandler : public IvyBase{
public:
  static const std::string colName_HLTpaths;
  static const std::string colName_metfilters;

protected:
  std::vector<HLTTriggerPathObject*> product_HLTpaths;
  std::unordered_map<std::string, bool> product_metfilters;

  void clear();

  bool constructHLTPaths();
  bool constructMETFilters();

public:
  // Constructors
  EventFilterHandler();

  // Destructors
  ~EventFilterHandler(){ clear(); }

  bool constructFilters();

  bool hasMatchingTriggerPath(std::vector<std::string> const& hltpaths_) const;
  float getTriggerWeight(std::vector<std::string> const& hltpaths_) const;
  bool passMETFilters() const;
  // Special event filters for various specific issues
  bool test2018HEMFilter(
    std::vector<ElectronObject*> const* electrons,
    std::vector<PhotonObject*> const* photons,
    std::vector<AK4JetObject*> const* ak4jets,
    std::vector<AK8JetObject*> const* ak8jets
  ) const;

  std::vector<HLTTriggerPathObject*> const& getHLTPaths() const{ return this->product_HLTpaths; }
  std::unordered_map<std::string, bool> const& getMETFilters() const{ return this->product_metfilters; }

  void bookBranches(BaseTree* intree);
  static std::vector<std::string> acquireMETFilterFlags(BaseTree* intree);

};


#endif
