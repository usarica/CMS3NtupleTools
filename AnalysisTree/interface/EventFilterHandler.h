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

protected:
  std::vector<HLTTriggerPathObject*> product_HLTpaths;

  void clear();

  bool constructHLTPaths();

public:
  // Constructors
  EventFilterHandler();

  // Destructors
  ~EventFilterHandler(){ clear(); }

  bool constructFilters();

  float getTriggerWeight(std::vector<std::string> const& hltpaths_) const;
  // Special event filters for various specific issues
  /*
  // Requires an implementation of Samples.h, so skip for now
  bool test2018HEMFilter(
    std::vector<ElectronObject*> const* electrons,
    std::vector<PhotonObject*> const* photons,
    std::vector<AK4JetObject*> const* ak4jets,
    std::vector<AK8JetObject*> const* ak8jets
  ) const;
  */

  std::vector<HLTTriggerPathObject*> const& getHLTPaths() const{ return this->product_HLTpaths; }

  void bookBranches(BaseTree* intree);
  //static std::vector<TString> acquireMETFilterFlags(BaseTree* intree);
  //static std::unordered_map<TString, std::vector<TString>> acquireHLTPaths(BaseTree* intree);

};


#endif
