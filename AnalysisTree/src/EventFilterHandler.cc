#include <cassert>
#include "TRandom3.h"

#include "ElectronSelectionHelpers.h"
#include "PhotonSelectionHelpers.h"
#include "AK4JetSelectionHelpers.h"
#include "AK8JetSelectionHelpers.h"
#include "EventFilterHandler.h"
#include "HelperFunctionsCore.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


#define VECTOR_ITERATOR_HANDLER_DIRECTIVES_HLTTRIGGERPATHS \
HLTTRIGGERPATH_VARIABLES


const std::string EventFilterHandler::colName_HLTpaths = "triggers";

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

  bool res = this->constructHLTPaths();

  return res;
}

float EventFilterHandler::testHLTPaths(std::vector<std::string> const& hltpaths_) const{
  float failRate = 1;
  for (auto str:hltpaths_){
    HelperFunctions::replaceString(str, "*", "");
    for (auto const* prod:product_HLTpaths){
      if (prod->name.find(str)!=std::string::npos && prod->passTrigger){
        float wgt=1;
        if (prod->L1prescale!=-1) wgt *= static_cast<float>(prod->L1prescale);
        if (prod->HLTprescale!=-1) wgt *= static_cast<float>(prod->HLTprescale);
        if (wgt == 1.f) return 1; // If the event passes an unprescaled trigger, its weight is 1.
        failRate *= 1.f-1.f/wgt;
      }
    }
  }
  return 1.f/(1.f-failRate);
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
  for (;it_HLTpaths_name != itEnd_HLTpaths_name; it_HLTpaths_name++){
    product_HLTpaths.push_back(new HLTTriggerPathObject());
    HLTTriggerPathObject*& obj = product_HLTpaths.back();

#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) obj->NAME = *it_HLTpaths_##NAME;
    VECTOR_ITERATOR_HANDLER_DIRECTIVES_HLTTRIGGERPATHS;
#undef HLTTRIGGERPATH_VARIABLE
  }

  return true;
}

void EventFilterHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

  // Book HLT paths
#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(EventFilterHandler::colName_HLTpaths + "_" + #NAME, nullptr);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_HLTTRIGGERPATHS
#undef HLTTRIGGERPATH_VARIABLE
}


#undef VECTOR_ITERATOR_HANDLER_DIRECTIVES_HLTTRIGGERPATHS
