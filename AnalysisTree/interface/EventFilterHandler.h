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
  static const std::string colName_HLTpaths;
  static const std::string colName_triggerobjects;
  static const std::string colName_metfilters;
  static const std::string colName_vertices;

protected:
  bool trackDataEvents;
  bool checkUniqueDataEvent;
  bool trackTriggerObjects;
  bool checkTriggerObjectsForHLTPaths;

  bool product_passCommonSkim;
  bool product_hasGoodVertex;
  bool product_uniqueEvent;

  std::vector<HLTTriggerPathObject*> product_HLTpaths;
  std::vector<TriggerObject*> product_triggerobjects;
  std::unordered_map<std::string, bool> product_metfilters;
  std::unordered_map<unsigned int, std::unordered_map<unsigned int, std::vector<unsigned long long>> > era_dataeventblock_map;

  void clear();

  bool constructCommonSkim();
  bool constructHLTPaths();
  bool constructTriggerObjects();
  bool constructMETFilters();
  bool constructVertexFilter();
  bool accumulateRunLumiEventBlock();

public:
  // Constructors
  EventFilterHandler();

  // Destructors
  ~EventFilterHandler(){ clear(); }

  bool constructFilters();

  bool hasMatchingTriggerPath(std::vector<std::string> const& hltpaths_) const;
  float getTriggerWeight(std::vector<std::string> const& hltpaths_) const;
  float getTriggerWeight(
    std::vector< std::unordered_map< TriggerHelpers::TriggerType, std::vector<HLTTriggerPathProperties> >::const_iterator > const& hltpathprops_,
    std::vector<MuonObject*> const* muons,
    std::vector<ElectronObject*> const* electrons,
    std::vector<PhotonObject*> const* photons,
    std::vector<AK4JetObject*> const* ak4jets,
    std::vector<AK8JetObject*> const* ak8jets,
    METObject const* pfmet
  ) const;
  bool passMETFilters() const;
  // Special event filters for various specific issues
  bool test2018HEMFilter(
    SimEventHandler const* simEventHandler,
    std::vector<ElectronObject*> const* electrons,
    std::vector<PhotonObject*> const* photons,
    std::vector<AK4JetObject*> const* ak4jets,
    std::vector<AK8JetObject*> const* ak8jets
  ) const;
  bool const& passCommonSkim() const{ return product_passCommonSkim; }
  bool const& hasGoodVertex() const{ return product_hasGoodVertex; }
  // For data trees. MC is always true
  bool const& isUniqueDataEvent() const{ return product_uniqueEvent; }

  void setTrackDataEvents(bool flag){ this->trackDataEvents=flag; }
  void setCheckUniqueDataEvent(bool flag){ this->checkUniqueDataEvent=flag; }
  void setTrackTriggerObjects(bool flag){ this->trackTriggerObjects=flag; }
  void setCheckTriggerObjectsForHLTPaths(bool flag){ this->checkTriggerObjectsForHLTPaths=flag; if (flag) this->trackTriggerObjects=flag; }

  std::vector<HLTTriggerPathObject*> const& getHLTPaths() const{ return this->product_HLTpaths; }
  std::vector<TriggerObject*> const& getTriggerObjects() const{ return this->product_triggerobjects; }
  std::unordered_map<std::string, bool> const& getMETFilters() const{ return this->product_metfilters; }

  void bookBranches(BaseTree* intree);
  static std::vector<std::string> acquireMETFilterFlags(BaseTree* intree);

};


#endif
