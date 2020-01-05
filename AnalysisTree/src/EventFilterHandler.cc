#include <cassert>
#include "TRandom3.h"

#include "EventFilterHandler.h"
#include "SamplesCore.h"
#include "HelperFunctionsCore.h"
#include "RunLumiEventBlock.h"
#include "AK4JetSelectionHelpers.h"
#include "AK8JetSelectionHelpers.h"
#include "ParticleSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


#define VERTEX_VARIABLES \
VERTEX_VARIABLE(unsigned int, nvtxs_good, 0)


const std::string EventFilterHandler::colName_HLTpaths = "triggers";
const std::string EventFilterHandler::colName_metfilters = "metfilter";
const std::string EventFilterHandler::colName_vertices = "vtxs";

EventFilterHandler::EventFilterHandler() :
  IvyBase(),
  product_passCommonSkim(true),
  product_uniqueEvent(true)
{
  // Common skim
  this->addConsumed<bool>("passCommonSkim");

  // HLT trigger variables
#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(EventFilterHandler::colName_HLTpaths + "_" + #NAME);
  HLTTRIGGERPATH_VARIABLES;
#undef HLTTRIGGERPATH_VARIABLE

  // Vertex variables
#define VERTEX_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<TYPE>(EventFilterHandler::colName_vertices + "_" + #NAME);
  VERTEX_VARIABLES;
#undef VERTEX_VARIABLE
}

void EventFilterHandler::clear(){
  product_passCommonSkim = true;
  product_uniqueEvent = true;
  for (auto*& prod:product_HLTpaths) delete prod;
  product_HLTpaths.clear();
}

bool EventFilterHandler::constructFilters(){
  clear();
  if (!currentTree) return false;

  bool res = this->constructCommonSkim() && this->constructHLTPaths() && this->constructMETFilters() && this->constructVertexFilter() && this->accumulateRunLumiEventBlock();

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

bool EventFilterHandler::test2018HEMFilter(
  SimEventHandler const* simEventHandler,
  std::vector<ElectronObject*> const* electrons,
  std::vector<PhotonObject*> const* photons,
  std::vector<AK4JetObject*> const* ak4jets,
  std::vector<AK8JetObject*> const* ak8jets
) const{
  if (SampleHelpers::theDataYear != 2018) return true;
  if (verbosity>=TVar::DEBUG) MELAerr << "Begin EventFilterHandler::test2018HEMFilter..." << endl;

  // Do not run clear because this is a special filter that does not modify the actual class
  if (!currentTree){
    if (verbosity>=TVar::ERROR) MELAerr << "EventFilterHandler::test2018HEMFilter: Current tree is null!" << endl;
    return false;
  }

  if (!SampleHelpers::checkSampleIsData(currentTree->sampleIdentifier)){
    if (!simEventHandler){
      if (verbosity>=TVar::ERROR) MELAerr << "EventFilterHandler::test2018HEMFilter: MC checks require a SimEventHandler!" << endl;
      assert(0);
      return false;
    }
    if (!simEventHandler->getHasHEM2018Issue()) return true;
  }
  else{
    bool allVariablesPresent = true;
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME = DEFVAL; allVariablesPresent &= this->getConsumedValue<TYPE>(#NAME, NAME);
    RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE
    if (!allVariablesPresent){
      if (this->verbosity>=TVar::ERROR) MELAerr << "EventFilterHandler::test2018HEMFilter: Not all variables of the data case are consumed properly!" << endl;
      assert(0);
    }
    if (!SampleHelpers::isHEM2018Affected(RunNumber)) return true; // All ok runs
  }

  // For affected runs, check object presence.
  static const std::pair<float, float> eta_region(-4.7, -1.4);
  static const std::pair<float, float> phi_region(-1.6, -0.8);
  bool doVeto = false;
  if (!doVeto && electrons){
    for (auto const* part:(*electrons)){
      if (!ParticleSelectionHelpers::isVetoParticle(part)) continue;

      float const eta = part->eta();
      float const phi = part->phi();
      if (eta>=eta_region.first && eta<=eta_region.second && phi>=phi_region.first && phi<=phi_region.second){ doVeto=true; break; }
    }
    if (verbosity>=TVar::DEBUG && doVeto) MELAout << "EventFilterHandler::test2018HEMFilter: Found at least one electron satisfying HEM15/16." << endl;
  }
  if (!doVeto && photons){
    for (auto const* part:(*photons)){
      if (!ParticleSelectionHelpers::isVetoParticle(part)) continue;

      float const eta = part->eta();
      float const phi = part->phi();
      if (eta>=eta_region.first && eta<=eta_region.second && phi>=phi_region.first && phi<=phi_region.second){ doVeto=true; break; }
    }
    if (verbosity>=TVar::DEBUG && doVeto) MELAout << "EventFilterHandler::test2018HEMFilter: Found at least one photon satisfying HEM15/16." << endl;
  }
  // Require a pT>30 GeV cut on jets
  if (!doVeto && ak4jets){
    for (auto const* part:(*ak4jets)){
      if (part->pt()<30.f) continue;
      if (!part->testSelectionBit(AK4JetSelectionHelpers::kLooseId)) continue;

      float const eta = part->eta();
      float const phi = part->phi();
      if (eta>=eta_region.first && eta<=eta_region.second && phi>=phi_region.first && phi<=phi_region.second){ doVeto=true; break; }
    }
    if (verbosity>=TVar::DEBUG && doVeto) MELAout << "EventFilterHandler::test2018HEMFilter: Found at least one AK4 jet satisfying HEM15/16." << endl;
  }
  // Be careful! There is no equivalent of tight ID in ak8 jets, so testing is done on all jets
  if (!doVeto && ak8jets){
    for (auto const* part:(*ak8jets)){
      if (part->pt()<30.f) continue;
      if (!part->testSelectionBit(AK8JetSelectionHelpers::kLooseId)) continue;

      float const eta = part->eta();
      float const phi = part->phi();
      if (eta>=eta_region.first && eta<=eta_region.second && phi>=phi_region.first && phi<=phi_region.second){ doVeto=true; break; }
    }
    if (verbosity>=TVar::DEBUG && doVeto) MELAout << "EventFilterHandler::test2018HEMFilter: Found at least one AK8 jet satisfying HEM15/16." << endl;
  }

  if (verbosity>=TVar::DEBUG) MELAerr << "End EventFilterHandler::test2018HEMFilter successfully." << endl;
  return !doVeto;
}

bool EventFilterHandler::constructCommonSkim(){
  if (!this->getConsumedValue<bool>("passCommonSkim", product_passCommonSkim)){
    if (this->verbosity>=TVar::ERROR) MELAerr << "EventFilterHandler::constructCommonSkim: Not all variables are consumed properly!" << endl;
    assert(0);
  }
  if (this->verbosity>=TVar::DEBUG) MELAout << "EventFilterHandler::constructCommonSkim: All variables are set up!" << endl;

  return true;
}

bool EventFilterHandler::constructHLTPaths(){
#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_HLTpaths_##NAME, itEnd_HLTpaths_##NAME;
  HLTTRIGGERPATH_VARIABLES;
#undef HLTTRIGGERPATH_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(EventFilterHandler::colName_HLTpaths + "_" + #NAME, &itBegin_HLTpaths_##NAME, &itEnd_HLTpaths_##NAME);
  HLTTRIGGERPATH_VARIABLES;
#undef HLTTRIGGERPATH_VARIABLE
  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "EventFilterHandler::constructHLTPaths: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "EventFilterHandler::constructHLTPaths: All variables are set up!" << endl;

  size_t n_HLTpaths = (itEnd_HLTpaths_name - itBegin_HLTpaths_name);
  product_HLTpaths.reserve(n_HLTpaths);
#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) auto it_HLTpaths_##NAME = itBegin_HLTpaths_##NAME;
  HLTTRIGGERPATH_VARIABLES;
#undef HLTTRIGGERPATH_VARIABLE
  while (it_HLTpaths_name != itEnd_HLTpaths_name){
    product_HLTpaths.push_back(new HLTTriggerPathObject());
    HLTTriggerPathObject*& obj = product_HLTpaths.back();

#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) obj->NAME = *it_HLTpaths_##NAME;
    HLTTRIGGERPATH_VARIABLES;
#undef HLTTRIGGERPATH_VARIABLE

#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) it_HLTpaths_##NAME++;
    HLTTRIGGERPATH_VARIABLES;
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

bool EventFilterHandler::constructVertexFilter(){
#define VERTEX_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME=DEFVAL;
  VERTEX_VARIABLES;
#undef VERTEX_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define VERTEX_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedValue<TYPE>(EventFilterHandler::colName_vertices + "_" + #NAME, NAME);
  VERTEX_VARIABLES;
#undef VERTEX_VARIABLE
  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "EventFilterHandler::constructVertexFilter: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "EventFilterHandler::constructVertexFilter: All variables are set up!" << endl;

  product_hasGoodVertex = (nvtxs_good>0);

  return true;
}


bool EventFilterHandler::accumulateRunLumiEventBlock(){
  if (!SampleHelpers::checkSampleIsData(currentTree->sampleIdentifier)){
    product_uniqueEvent = true;
    return true;
  }

  bool allVariablesPresent = true;
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME = DEFVAL; allVariablesPresent &= this->getConsumedValue<TYPE>(#NAME, NAME);
  RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE
  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "EventFilterHandler::accumulateRunLumiEventBlock: Not all variables are consumed properly!" << endl;
    assert(0);
  }
  if (this->verbosity>=TVar::DEBUG) MELAout << "EventFilterHandler::accumulateRunLumiEventBlock: All variables are set up!" << endl;

  auto it_run = era_dataeventblock_map.find(RunNumber);
  if (it_run == era_dataeventblock_map.end()){
    era_dataeventblock_map[RunNumber] = std::unordered_map<unsigned int, std::vector<unsigned long long>>();
    it_run = era_dataeventblock_map.find(RunNumber);
  }
  auto it_lumi = it_run->second.find(LuminosityBlock);
  if (it_lumi == it_run->second.end()){
    it_run->second[LuminosityBlock] = std::vector<unsigned long long>();
    it_lumi = it_run->second.find(LuminosityBlock);
  }
  auto it_event = std::find(it_lumi->second.begin(), it_lumi->second.end(), EventNumber);
  if (it_event == it_lumi->second.end()){
    it_lumi->second.push_back(EventNumber);
    product_uniqueEvent = true;
  }
  else product_uniqueEvent = false;

  return true;
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

  // Common skim
  tree->bookBranch<bool>("passCommonSkim", false);

  // Book HLT paths
#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(EventFilterHandler::colName_HLTpaths + "_" + #NAME, nullptr);
  HLTTRIGGERPATH_VARIABLES;
#undef HLTTRIGGERPATH_VARIABLE

  // Vertex variables
#define VERTEX_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<TYPE>(EventFilterHandler::colName_vertices + "_" + #NAME, DEFVAL);
  VERTEX_VARIABLES;
#undef VERTEX_VARIABLE

  // Book MET filters
  auto strmetfilters = EventFilterHandler::acquireMETFilterFlags(tree);
  for (auto const& strmetfilter:strmetfilters){
    tree->bookBranch<bool>(strmetfilter, true); // Default should be true to avoid non-existing branches
    this->addConsumed<bool>(strmetfilter);
    this->defineConsumedSloppy(strmetfilter); // Define as sloppy so that different samples from different years/versions can be processed.
  }

  // Do these for data trees
  if (SampleHelpers::checkSampleIsData(tree->sampleIdentifier)){
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<TYPE>(#NAME, DEFVAL); this->addConsumed<TYPE>(#NAME); this->defineConsumedSloppy(#NAME);
    RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE
  }
}


#undef VERTEX_VARIABLES
