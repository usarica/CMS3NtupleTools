#ifndef EVENTFILTERHANDLER_H
#define EVENTFILTERHANDLER_H

#include <vector>
#include <unordered_map>
#include "IvyBase.h"
#include "SimEventHandler.h"
#include "MuonObject.h"
#include "ElectronObject.h"
#include "PhotonObject.h"
#include "AK4JetObject.h"
#include "AK8JetObject.h"
#include "METObject.h"
#include "HLTTriggerPathObject.h"
#include "TriggerObject.h"
#include "TriggerHelpersCore.h"
#include "SystematicVariations.h"


class EventFilterHandler : public IvyBase{
public:
  enum METFilterCutType{
    kMETFilters_Standard = 0,
    kMETFilters_Tight,
    nMETFilterCutTypes
  };

  static const std::string colName_HLTpaths;
  static const std::string colName_triggerobjects;
  static const std::string colName_metfilters;

protected:
  bool trackDataEvents;
  bool checkUniqueDataEvent;
  bool checkHLTPathRunRanges;
  bool trackTriggerObjects;
  bool checkTriggerObjectsForHLTPaths;

  bool product_passCommonSkim;
  bool product_uniqueEvent;

  std::vector<HLTTriggerPathObject*> product_HLTpaths;
  std::vector<TriggerObject*> product_triggerobjects;
  std::unordered_map<std::string, bool> product_metfilters;
  std::unordered_map<unsigned int, std::unordered_map<unsigned int, std::vector<unsigned long long>> > era_dataeventblock_map;

  void clear();

  bool constructCommonSkim();
  bool constructHLTPaths(SimEventHandler const* simEventHandler);
  bool constructTriggerObjects();
  bool constructMETFilters();
  bool accumulateRunLumiEventBlock();

public:
  // Constructors
  EventFilterHandler();

  // Destructors
  ~EventFilterHandler(){ clear(); }

  bool constructFilters(SimEventHandler const* simEventHandler);

  bool hasMatchingTriggerPath(std::vector<std::string> const& hltpaths_) const;
  float getTriggerWeight(std::vector<std::string> const& hltpaths_) const;
  float getTriggerWeight(
    std::vector< std::unordered_map< TriggerHelpers::TriggerType, std::vector<HLTTriggerPathProperties> >::const_iterator > const& hltpathprops_,
    std::vector<MuonObject*> const* muons,
    std::vector<ElectronObject*> const* electrons,
    std::vector<PhotonObject*> const* photons,
    std::vector<AK4JetObject*> const* ak4jets,
    std::vector<AK8JetObject*> const* ak8jets,
    METObject const* pfmet,
    HLTTriggerPathObject const** firstPassingHLTPath = nullptr
  ) const;
  float getTriggerWeight(
    std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > const& hltpathprops_,
    std::vector<MuonObject*> const* muons,
    std::vector<ElectronObject*> const* electrons,
    std::vector<PhotonObject*> const* photons,
    std::vector<AK4JetObject*> const* ak4jets,
    std::vector<AK8JetObject*> const* ak8jets,
    METObject const* pfmet,
    HLTTriggerPathObject const** firstPassingHLTPath = nullptr
  ) const;
  bool passMETFilters(EventFilterHandler::METFilterCutType const& cuttype) const;
  // Special event filters for various specific issues
  bool test2018HEMFilter(
    SimEventHandler const* simEventHandler,
    std::vector<ElectronObject*> const* electrons,
    std::vector<PhotonObject*> const* photons,
    std::vector<AK4JetObject*> const* ak4jets,
    std::vector<AK8JetObject*> const* ak8jets
  ) const;
  bool const& passCommonSkim() const{ return product_passCommonSkim; }
  // For data trees. MC is always true
  bool const& isUniqueDataEvent() const{ return product_uniqueEvent; }

  void setTrackDataEvents(bool flag){ this->trackDataEvents=flag; }
  void setCheckUniqueDataEvent(bool flag){ this->checkUniqueDataEvent=flag; }
  void setCheckHLTPathRunRanges(bool flag){ this->checkHLTPathRunRanges=flag; }
  void setTrackTriggerObjects(bool flag){ this->trackTriggerObjects=flag; }
  void setCheckTriggerObjectsForHLTPaths(bool flag){ this->checkTriggerObjectsForHLTPaths=flag; if (flag) this->trackTriggerObjects=flag; }

  std::vector<HLTTriggerPathObject*> const& getHLTPaths() const{ return this->product_HLTpaths; }
  std::vector<TriggerObject*> const& getTriggerObjects() const{ return this->product_triggerobjects; }
  std::unordered_map<std::string, bool> const& getMETFilters() const{ return this->product_metfilters; }

  void bookBranches(BaseTree* intree);
  static std::vector<std::string> acquireMETFilterFlags(BaseTree* intree, EventFilterHandler::METFilterCutType const& cuttype);

};


#endif
