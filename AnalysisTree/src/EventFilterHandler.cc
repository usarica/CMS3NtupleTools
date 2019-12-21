#include <cassert>
#include "TRandom3.h"

#include "EventFilterHandler.h"
#include "SamplesCore.h"
#include "HelperFunctionsCore.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


#define VECTOR_ITERATOR_HANDLER_DIRECTIVES_HLTTRIGGERPATHS \
HLTTRIGGERPATH_VARIABLES


const std::string EventFilterHandler::colName_HLTpaths = "triggers";
const std::string EventFilterHandler::colName_metfilters = "metfilter";

EventFilterHandler::EventFilterHandler() :
  IvyBase()
{
#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(EventFilterHandler::colName_HLTpaths + "_" + #NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_HLTTRIGGERPATHS;
#undef HLTTRIGGERPATH_VARIABLE
}

void EventFilterHandler::clear(){
  for (auto*& prod:product_HLTpaths) delete prod;
  product_HLTpaths.clear();
}

bool EventFilterHandler::constructFilters(){
  clear();
  if (!currentTree) return false;

  bool res = this->constructHLTPaths() && this->constructMETFilters();

  return res;
}

bool EventFilterHandler::hasMatchingTriggerPath(std::vector<std::string> const& hltpaths_) const{
  bool res = false;
  for (auto str:hltpaths_){
    HelperFunctions::replaceString(str, "*", "");
    for (auto const* prod:product_HLTpaths){ if (prod->name.find(str)!=std::string::npos){ res = true; break; } }
  }
  return res;
}
float EventFilterHandler::getTriggerWeight(std::vector<std::string> const& hltpaths_) const{
  if (hltpaths_.empty()) return 1;
  float failRate = 1;
  bool foundAtLeastOneTrigger = false;
  for (auto str:hltpaths_){
    HelperFunctions::replaceString(str, "*", "");
    for (auto const* prod:product_HLTpaths){
      if (prod->name.find(str)!=std::string::npos && prod->passTrigger){
        float wgt = 1.f;
        if (prod->L1prescale>0) wgt *= static_cast<float>(prod->L1prescale);
        if (prod->HLTprescale>0) wgt *= static_cast<float>(prod->HLTprescale);
        if (wgt == 1.f) return wgt; // If the event passes an unprescaled trigger, its weight is 1.
        else if (wgt == 0.f) continue;
        foundAtLeastOneTrigger = true;
        failRate *= 1.f-1.f/wgt;
      }
    }
  }
  return (foundAtLeastOneTrigger ? 1.f/(1.f-failRate) : 0.f);
}

bool EventFilterHandler::passMETFilters() const{
  bool res = true;
  for (auto it:product_metfilters) res &= it.second;
  return res;
}

bool EventFilterHandler::constructHLTPaths(){
#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_HLTpaths_##NAME, itEnd_HLTpaths_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_HLTTRIGGERPATHS;
#undef HLTTRIGGERPATH_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(EventFilterHandler::colName_HLTpaths + "_" + #NAME, &itBegin_HLTpaths_##NAME, &itEnd_HLTpaths_##NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_HLTTRIGGERPATHS;
#undef HLTTRIGGERPATH_VARIABLE
  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "EventFilterHandler::constructHLTPaths: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "EventFilterHandler::constructHLTPaths: All variables are set up!" << endl;

  size_t n_HLTpaths = (itEnd_HLTpaths_name - itBegin_HLTpaths_name);
  product_HLTpaths.reserve(n_HLTpaths);
#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) auto it_HLTpaths_##NAME = itBegin_HLTpaths_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_HLTTRIGGERPATHS;
#undef HLTTRIGGERPATH_VARIABLE
  while (it_HLTpaths_name != itEnd_HLTpaths_name){
    product_HLTpaths.push_back(new HLTTriggerPathObject());
    HLTTriggerPathObject*& obj = product_HLTpaths.back();

#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) obj->NAME = *it_HLTpaths_##NAME;
    VECTOR_ITERATOR_HANDLER_DIRECTIVES_HLTTRIGGERPATHS;
#undef HLTTRIGGERPATH_VARIABLE

#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) it_HLTpaths_##NAME++;
    VECTOR_ITERATOR_HANDLER_DIRECTIVES_HLTTRIGGERPATHS;
#undef HLTTRIGGERPATH_VARIABLE
  }

  return true;
}

bool EventFilterHandler::constructMETFilters(){
  auto strmetfilters = EventFilterHandler::acquireMETFilterFlags(currentTree);
  product_metfilters.clear();

  bool allVariablesPresent = true;
  for (auto const& strmetfilter:strmetfilters){
    product_metfilters[strmetfilter] = false;
    allVariablesPresent &= this->getConsumedValue<bool>(strmetfilter, product_metfilters.find(strmetfilter)->second);
  }
  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "EventFilterHandler::constructMETFilters: Not all variables are consumed properly!" << endl;
    assert(0);
  }
  if (this->verbosity>=TVar::DEBUG) MELAout << "EventFilterHandler::constructMETFilters: All variables are set up!" << endl;

  return allVariablesPresent;
}

std::vector<std::string> EventFilterHandler::acquireMETFilterFlags(BaseTree* intree){
  std::vector<std::string> res;

  switch (SampleHelpers::theDataYear){
  case 2016:
  {
    res = std::vector<std::string>{
      "goodVertices",
      "HBHENoiseFilter",
      "HBHENoiseIsoFilter",
      "EcalDeadCellTriggerPrimitiveFilter"
    };
    if (!SampleHelpers::checkSampleIsFastSim(intree->sampleIdentifier)) res.push_back("globalSuperTightHalo2016Filter"); // For data or non-FS MC
    if (SampleHelpers::checkSampleIsData(intree->sampleIdentifier)) res.push_back("eeBadScFilter"); // Only for data

    if (!SampleHelpers::checkSampleIs80X(intree->sampleIdentifier)){ // These MET filters are available in CMSSW_VERSION>=94X
      res.push_back("BadPFMuonFilter");
      //res.push_back("BadChargedCandidateFilter"); // Disabled per https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#2016_data
    }
    // Else need "Bad PF Muon Filter" and "Bad Charged Hadron Filter" to be calculated on the fly for data, MC and FastSim, see https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#How_to_run_the_Bad_Charged_Hadro
    break;
  }
  case 2017:
  {
    res = std::vector<std::string>{
      "goodVertices",
      "HBHENoiseFilter",
      "HBHENoiseIsoFilter",
      "EcalDeadCellTriggerPrimitiveFilter",
      "BadPFMuonFilter",
      //"BadChargedCandidateFilter", // FIXME: To be updated following https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#2017_data
      "ecalBadCalibFilterUpdated"
    };
    if (!SampleHelpers::checkSampleIsFastSim(intree->sampleIdentifier)) res.push_back("globalSuperTightHalo2016Filter"); // For data or non-FS MC
    if (SampleHelpers::checkSampleIsData(intree->sampleIdentifier)) res.push_back("eeBadScFilter"); // Only for data
    break;
  }
  case 2018:
  {
    res = std::vector<std::string>{
      "goodVertices",
      "HBHENoiseFilter",
      "HBHENoiseIsoFilter",
      "EcalDeadCellTriggerPrimitiveFilter",
      "BadPFMuonFilter",
      //"BadChargedCandidateFilter", // FIXME: To be updated following https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#2018_data
      "ecalBadCalibFilterUpdated"
    };
    if (!SampleHelpers::checkSampleIsFastSim(intree->sampleIdentifier)) res.push_back("globalSuperTightHalo2016Filter"); // For data or non-FS MC
    if (SampleHelpers::checkSampleIsData(intree->sampleIdentifier)) res.push_back("eeBadScFilter"); // Only for data
    break;
  }
  default:
    MELAerr << "EventFilterHandler::acquireMETFilterFlags: Data year " << SampleHelpers::theDataYear << " is not implemented!" << endl;
    assert(0);
  }

  for (auto& strmetfilter:res) strmetfilter = EventFilterHandler::colName_metfilters + "_" + strmetfilter;
  return res;
}

void EventFilterHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

  // Book HLT paths
#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(EventFilterHandler::colName_HLTpaths + "_" + #NAME, nullptr);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_HLTTRIGGERPATHS;
#undef HLTTRIGGERPATH_VARIABLE

  // Book MET filters
  auto strmetfilters = EventFilterHandler::acquireMETFilterFlags(tree);
  for (auto const& strmetfilter:strmetfilters){
    tree->bookBranch<bool>(strmetfilter, true); // Default should be true to avoid non-existing branches
    this->addConsumed<bool>(strmetfilter);
    this->defineConsumedSloppy(strmetfilter); // Define as sloppy so that different samples from different years/versions can be processed.
  }
}


#undef VECTOR_ITERATOR_HANDLER_DIRECTIVES_HLTTRIGGERPATHS
