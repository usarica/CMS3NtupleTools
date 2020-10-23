#include <cassert>
#include "TRandom3.h"

#include <CMS3/Dictionaries/interface/GlobalCollectionNames.h>

#include "EventFilterHandler.h"
#include "SamplesCore.h"
#include "HelperFunctionsCore.h"
#include "RunLumiEventBlock.h"
#include "AK4JetSelectionHelpers.h"
#include "AK8JetSelectionHelpers.h"
#include "ParticleSelectionHelpers.h"
#include "ParticleObjectHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


const std::string EventFilterHandler::colName_HLTpaths = GlobalCollectionNames::colName_triggers;
const std::string EventFilterHandler::colName_triggerobjects = GlobalCollectionNames::colName_triggerObjects;
const std::string EventFilterHandler::colName_metfilters = GlobalCollectionNames::colName_metfilter;

EventFilterHandler::EventFilterHandler() :
  IvyBase(),
  trackDataEvents(true),
  checkUniqueDataEvent(true),
  checkHLTPathRunRanges(true),
  trackTriggerObjects(false),
  checkTriggerObjectsForHLTPaths(false),
  product_passCommonSkim(true),
  product_uniqueEvent(true)
{
  // Common skim
  this->addConsumed<bool>("passCommonSkim");

  // HLT trigger variables
#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(EventFilterHandler::colName_HLTpaths + "_" + #NAME);
  HLTTRIGGERPATH_VARIABLES;
#undef HLTTRIGGERPATH_VARIABLE
}

void EventFilterHandler::clear(){
  this->resetCache();

  product_passCommonSkim = true;
  product_uniqueEvent = true;
  for (auto*& prod:product_HLTpaths) delete prod;
  product_HLTpaths.clear();
  for (auto*& prod:product_triggerobjects) delete prod;
  product_triggerobjects.clear();
}

bool EventFilterHandler::constructFilters(SimEventHandler const* simEventHandler){
  if (this->isAlreadyCached()) return true;

  clear();
  if (!currentTree) return false;

  bool res = (
    this->constructCommonSkim()
    &&
    (!trackTriggerObjects || this->constructTriggerObjects())
    &&
    this->constructHLTPaths(simEventHandler)
    &&
    this->constructMETFilters()
    &&
    this->accumulateRunLumiEventBlock()
    );

  if (res) this->cacheEvent();
  return res;
}

bool EventFilterHandler::hasMatchingTriggerPath(std::vector<std::string> const& hltpaths_) const{
  bool res = false;
  for (auto const& str:hltpaths_){
    for (auto const* prod:product_HLTpaths){ if (prod->name.find(str)!=std::string::npos){ res = true; break; } }
  }
  return res;
}
float EventFilterHandler::getTriggerWeight(std::vector<std::string> const& hltpaths_) const{
  if (hltpaths_.empty()) return 0;
  float failRate = 1;
  bool foundAtLeastOneTrigger = false;
  for (auto const& str:hltpaths_){
    for (auto const* prod:product_HLTpaths){
      if (prod->isValid() && prod->name.find(str)!=std::string::npos && prod->passTrigger){
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
float EventFilterHandler::getTriggerWeight(
  std::vector< std::unordered_map< TriggerHelpers::TriggerType, std::vector<HLTTriggerPathProperties> >::const_iterator > const& hltpathprops_,
  std::vector<MuonObject*> const* muons,
  std::vector<ElectronObject*> const* electrons,
  std::vector<PhotonObject*> const* photons,
  std::vector<AK4JetObject*> const* ak4jets,
  std::vector<AK8JetObject*> const* ak8jets,
  METObject const* pfmet,
  HLTTriggerPathObject const** firstPassingHLTPath,
  std::vector<ParticleObject const*>* outparticles_TOmatched
) const{
  unsigned int isize=0;
  for (auto const& p:hltpathprops_) isize += p->second.size();

  std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > hltpathprops_new; hltpathprops_new.reserve(isize);
  for (auto const& p:hltpathprops_){ for (auto const& pp:p->second) hltpathprops_new.emplace_back(p->first, &pp); }

  return getTriggerWeight(
    hltpathprops_new,
    muons, electrons, photons, ak4jets, ak8jets, pfmet,
    firstPassingHLTPath, outparticles_TOmatched
  );
}
float EventFilterHandler::getTriggerWeight(
  std::vector< std::pair<TriggerHelpers::TriggerType, HLTTriggerPathProperties const*> > const& hltpathprops_,
  std::vector<MuonObject*> const* muons,
  std::vector<ElectronObject*> const* electrons,
  std::vector<PhotonObject*> const* photons,
  std::vector<AK4JetObject*> const* ak4jets,
  std::vector<AK8JetObject*> const* ak8jets,
  METObject const* pfmet,
  HLTTriggerPathObject const** firstPassingHLTPath,
  std::vector<ParticleObject const*>* outparticles_TOmatched
) const{
  if (hltpathprops_.empty()) return 0;

  std::vector<MuonObject const*> muons_trigcheck; if (muons){ muons_trigcheck.reserve(muons->size()); for (auto const& part:(*muons)){ if (ParticleSelectionHelpers::isParticleForTriggerChecking(part)) muons_trigcheck.push_back(part); } }
  std::vector<ElectronObject const*> electrons_trigcheck; if (electrons){ electrons_trigcheck.reserve(electrons->size()); for (auto const& part:(*electrons)){ if (ParticleSelectionHelpers::isParticleForTriggerChecking(part)) electrons_trigcheck.push_back(part); } }
  std::vector<PhotonObject const*> photons_trigcheck; if (photons){ photons_trigcheck.reserve(photons->size()); for (auto const& part:(*photons)){ if (ParticleSelectionHelpers::isParticleForTriggerChecking(part)) photons_trigcheck.push_back(part); } }
  std::vector<AK4JetObject const*> ak4jets_trigcheck; if (ak4jets){ ak4jets_trigcheck.reserve(ak4jets->size()); for (auto const& jet:(*ak4jets)){ if (ParticleSelectionHelpers::isJetForTriggerChecking(jet)) ak4jets_trigcheck.push_back(jet); } }
  std::vector<AK8JetObject const*> ak8jets_trigcheck; if (ak8jets){ ak8jets_trigcheck.reserve(ak8jets->size()); for (auto const& jet:(*ak8jets)){ if (ParticleSelectionHelpers::isJetForTriggerChecking(jet)) ak8jets_trigcheck.push_back(jet); } }

  ParticleObject::LorentzVector_t pfmet_p4, pfmet_nomus_p4;
  if (pfmet){
    pfmet_p4 = pfmet->p4(true, true, true, true);
    pfmet_nomus_p4 = pfmet_p4;
    for (auto const& part:muons_trigcheck) pfmet_nomus_p4 += part->p4();
  }

  float ht_pt=0, ht_nomus_pt=0;
  ParticleObject::LorentzVector_t ht_p4, ht_nomus_p4;
  for (auto const& jet:ak4jets_trigcheck){
    auto jet_p4_nomus = jet->p4_nomus_basic();
    auto const& jet_p4 = jet->p4();

    ht_pt += jet_p4.Pt();
    ht_p4 += jet_p4;

    ht_nomus_pt += jet_p4_nomus.Pt();
    ht_nomus_p4 += jet_p4_nomus;
  }
  ht_p4 = ParticleObject::PolarLorentzVector_t(ht_pt, 0, 0, ht_p4.Pt());
  ht_nomus_p4 = ParticleObject::PolarLorentzVector_t(ht_nomus_pt, 0, 0, ht_nomus_p4.Pt());

  for (auto const& enumType_props_pair:hltpathprops_){
    assert(enumType_props_pair.second != nullptr);
    auto const& hltprop = *(enumType_props_pair.second);
    for (auto const& prod:product_HLTpaths){
      if (prod->isValid() && prod->passTrigger && hltprop.isSameTrigger(prod->name)){
        std::vector<MuonObject const*> muons_trigcheck_TOmatched;
        std::vector<ElectronObject const*> electrons_trigcheck_TOmatched;
        std::vector<PhotonObject const*> photons_trigcheck_TOmatched;

        if (checkTriggerObjectsForHLTPaths){
          HLTTriggerPathProperties::TriggerObjectExceptionType const& TOexception = hltprop.getTOException();
          auto const& passedTriggerObjects = prod->getPassedTriggerObjects();

          bool hasTORecovery = false;
          std::vector<TriggerObject const*> passedTriggerObjectsWithRecovery;
          if (TOexception == HLTTriggerPathProperties::toRecoverObjectsFromFailing){
            // Determine what to recover
            unsigned short n_TOmuons = 0;
            unsigned short n_TOelectrons = 0;
            unsigned short n_TOphotons = 0;
            for (auto const& TOobj:passedTriggerObjects){
              if (TOobj->isTriggerObjectType(trigger::TriggerMuon)) n_TOmuons++;
              else if (TOobj->isTriggerObjectType(trigger::TriggerElectron)) n_TOelectrons++;
              else if (TOobj->isTriggerObjectType(trigger::TriggerPhoton) || TOobj->isTriggerObjectType(trigger::TriggerCluster)){
                n_TOelectrons++;
                n_TOphotons++;
              }
            }
            auto const& props_map = hltprop.getObjectProperties();
            unsigned short n_TOmuons_req = (props_map.find(HLTObjectProperties::kMuon)!=props_map.end() ? props_map.find(HLTObjectProperties::kMuon)->second.size() : 0);
            unsigned short n_TOelectrons_req = (props_map.find(HLTObjectProperties::kElectron)!=props_map.end() ? props_map.find(HLTObjectProperties::kElectron)->second.size() : 0);
            unsigned short n_TOphotons_req = (props_map.find(HLTObjectProperties::kPhoton)!=props_map.end() ? props_map.find(HLTObjectProperties::kPhoton)->second.size() : 0);
            bool needMuonRecovery = (n_TOmuons<n_TOmuons_req);
            bool needElectronRecovery = (n_TOelectrons<n_TOelectrons_req);
            bool needPhotonRecovery = (n_TOphotons<n_TOphotons_req);
            // Add existing passing objects
            for (auto const& TOobj:passedTriggerObjects) passedTriggerObjectsWithRecovery.push_back(TOobj);
            // Add from failing objects
            for (auto const& TOobj:prod->getFailedTriggerObjects()){
              if (
                (needMuonRecovery && TOobj->isTriggerObjectType(trigger::TriggerMuon))
                ||
                (needElectronRecovery && (TOobj->isTriggerObjectType(trigger::TriggerElectron) || TOobj->isTriggerObjectType(trigger::TriggerPhoton) || TOobj->isTriggerObjectType(trigger::TriggerCluster)))
                ||
                (needPhotonRecovery && (TOobj->isTriggerObjectType(trigger::TriggerPhoton) || TOobj->isTriggerObjectType(trigger::TriggerCluster)))
                ) passedTriggerObjectsWithRecovery.push_back(TOobj);
            }
            hasTORecovery = true;
          }

          std::vector<TriggerObject const*> const& passedTriggerObjects_final = (!hasTORecovery ? passedTriggerObjects : passedTriggerObjectsWithRecovery);

          /*
          // Consistency check
          unsigned short n_TOmuons = 0;
          unsigned short n_TOelectrons = 0;
          unsigned short n_TOphotons = 0;
          for (auto const& TOobj:passedTriggerObjects_final){
            if (TOobj->isTriggerObjectType(trigger::TriggerMuon)) n_TOmuons++;
            else if (TOobj->isTriggerObjectType(trigger::TriggerElectron)) n_TOelectrons++;
            else if (TOobj->isTriggerObjectType(trigger::TriggerPhoton) || TOobj->isTriggerObjectType(trigger::TriggerCluster)){
              n_TOelectrons++;
              n_TOphotons++;
            }
          }
          auto const& props_map = hltprop.getObjectProperties();
          unsigned short n_TOmuons_req = (props_map.find(HLTObjectProperties::kMuon)!=props_map.end() ? props_map.find(HLTObjectProperties::kMuon)->second.size() : 0);
          unsigned short n_TOelectrons_req = (props_map.find(HLTObjectProperties::kElectron)!=props_map.end() ? props_map.find(HLTObjectProperties::kElectron)->second.size() : 0);
          unsigned short n_TOphotons_req = (props_map.find(HLTObjectProperties::kPhoton)!=props_map.end() ? props_map.find(HLTObjectProperties::kPhoton)->second.size() : 0);
          if (this->verbosity>=TVar::ERROR){
            bool hasError = false;
            if (n_TOmuons<n_TOmuons_req){
              MELAerr << "EventFilterHandler::getTriggerWeight[" << prod->name << "]: Number of muons " << n_TOmuons << " < number of req. muons " << n_TOmuons_req << endl;
              hasError = true;
            }
            if (n_TOelectrons<n_TOelectrons_req){
              MELAerr << "EventFilterHandler::getTriggerWeight[" << prod->name << "]: Number of electrons " << n_TOelectrons << " < number of req. electrons " << n_TOelectrons_req << endl;
              hasError = true;
            }
            if (n_TOphotons<n_TOphotons_req){
              MELAerr << "EventFilterHandler::getTriggerWeight[" << prod->name << "]: Number of photons " << n_TOphotons << " < number of req. photons " << n_TOphotons_req << endl;
              hasError = true;
            }
            if (hasError){
              MELAout << "\t\t- Passing trigger objects: ";
              for (auto const& TOobj:passedTriggerObjects){
                MELAout << "[" << TOobj->getTriggerObjectType() << "] ( " << TOobj->pt() << ", " << TOobj->eta() << ", " << TOobj->phi() << endl;
              }
              MELAout << "\t\t- Failing trigger objects: ";
              for (auto const& TOobj:prod->getFailedTriggerObjects()){
                MELAout << "[" << TOobj->getTriggerObjectType() << "] ( " << TOobj->pt() << ", " << TOobj->eta() << ", " << TOobj->phi() << endl;
              }
              MELAout << "\t\t- Leptons: ";
              for (auto const& part:muons_trigcheck){
                MELAout << "[" << part->pdgId() << "] ( " << part->pt() << ", " << part->eta() << ", " << part->phi() << endl;
              }
              for (auto const& part:electrons_trigcheck){
                MELAout << "[" << part->pdgId() << "] ( " << part->pt() << ", " << part->eta() << ", " << part->phi() << endl;
              }
              for (auto const& part:photons_trigcheck){
                MELAout << "[" << part->pdgId() << "] ( " << part->pt() << ", " << part->eta() << ", " << part->phi() << endl;
              }
            }
          }
          */

          if (this->verbosity>=TVar::DEBUG){
            MELAout << "EventFilterHandler::getTriggerWeight: Checking " << prod->name << " trigger objects:" << endl;
            MELAout << "\t- Number of passed trigger objects: " << passedTriggerObjects_final.size() << endl;
            MELAout << "\t\t- Trigger object types: ";
            std::vector<trigger::TriggerObjectType> TOtypes;
            for (auto const& TOobj:passedTriggerObjects_final) TOtypes.push_back(TOobj->getTriggerObjectType());
            MELAout << TOtypes << endl;
            MELAout << "\t- Number of muons: " << muons_trigcheck.size() << endl;
            MELAout << "\t- Number of electrons: " << electrons_trigcheck.size() << endl;
            MELAout << "\t- Number of photons: " << photons_trigcheck.size() << endl;
          }

          TriggerObject::getMatchedPhysicsObjects(
            passedTriggerObjects_final, { trigger::TriggerMuon }, 0.2,
            muons_trigcheck, muons_trigcheck_TOmatched
          );
          TriggerObject::getMatchedPhysicsObjects(
            passedTriggerObjects_final, { trigger::TriggerElectron, trigger::TriggerPhoton, trigger::TriggerCluster }, 0.2,
            electrons_trigcheck, electrons_trigcheck_TOmatched
          );
          TriggerObject::getMatchedPhysicsObjects(
            passedTriggerObjects_final, { trigger::TriggerPhoton, trigger::TriggerCluster }, 0.2,
            photons_trigcheck, photons_trigcheck_TOmatched
          );

          if (this->verbosity>=TVar::DEBUG){
            MELAout << "\t- Number of matched muons: " << muons_trigcheck_TOmatched.size() << " / " << muons_trigcheck.size() << endl;
            MELAout << "\t- Number of matched electrons: " << electrons_trigcheck_TOmatched.size() << " / " << electrons_trigcheck.size() << endl;
            MELAout << "\t- Number of matched photons: " << photons_trigcheck_TOmatched.size() << " / " << photons_trigcheck.size() << endl;
          }
        }

        if (
          hltprop.testCuts(
            (!checkTriggerObjectsForHLTPaths ? muons_trigcheck : muons_trigcheck_TOmatched),
            (!checkTriggerObjectsForHLTPaths ? electrons_trigcheck : electrons_trigcheck_TOmatched),
            (!checkTriggerObjectsForHLTPaths ? photons_trigcheck : photons_trigcheck_TOmatched),
            ak4jets_trigcheck,
            ak8jets_trigcheck,
            pfmet_p4,
            pfmet_nomus_p4,
            ht_p4,
            ht_nomus_p4
          )
          ){
          float wgt = 1.f;
          if (prod->L1prescale>0) wgt *= static_cast<float>(prod->L1prescale);
          if (prod->HLTprescale>0) wgt *= static_cast<float>(prod->HLTprescale);
          if (wgt == 0.f) continue;
          else{
            if (firstPassingHLTPath) *firstPassingHLTPath = prod;
            if (outparticles_TOmatched && checkTriggerObjectsForHLTPaths){
              outparticles_TOmatched->clear();
              outparticles_TOmatched->reserve(muons_trigcheck_TOmatched.size() + electrons_trigcheck_TOmatched.size() + photons_trigcheck_TOmatched.size());
              for (auto const& part:muons_trigcheck_TOmatched) outparticles_TOmatched->push_back(part);
              for (auto const& part:electrons_trigcheck_TOmatched) outparticles_TOmatched->push_back(part);
              for (auto const& part:photons_trigcheck_TOmatched) outparticles_TOmatched->push_back(part);
            }
            return wgt; // Take the first trigger that passed.
          }
        }
      }
    }
  }

  return 0;
}

bool EventFilterHandler::passMETFilters(EventFilterHandler::METFilterCutType const& cuttype) const{
  bool res = true;
  auto strmetfilters = EventFilterHandler::acquireMETFilterFlags(currentTree, cuttype);
  for (auto const& it:product_metfilters){ if (HelperFunctions::checkListVariable(strmetfilters, it.first)) res &= it.second; }
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
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) TYPE const* NAME = nullptr; allVariablesPresent &= this->getConsumed(#NAME, NAME);
    RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE
    if (!allVariablesPresent){
      if (this->verbosity>=TVar::ERROR) MELAerr << "EventFilterHandler::test2018HEMFilter: Not all variables of the data case are consumed properly!" << endl;
      assert(0);
    }
    if (!SampleHelpers::isHEM2018Affected(*RunNumber)) return true; // All ok runs
  }

  // For affected runs, check object presence.
  static const std::pair<float, float> eta_region(-3.0, -1.4);
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
      if (!ParticleSelectionHelpers::isJetForHEMVeto(part)) continue;

      float const eta = part->eta();
      float const phi = part->phi();
      if (eta>=eta_region.first && eta<=eta_region.second && phi>=phi_region.first && phi<=phi_region.second){ doVeto=true; break; }
    }
    if (verbosity>=TVar::DEBUG && doVeto) MELAout << "EventFilterHandler::test2018HEMFilter: Found at least one AK4 jet satisfying HEM15/16." << endl;
  }
  // Be careful! There is no equivalent of tight ID in ak8 jets, so testing is done on all jets
  if (!doVeto && ak8jets){
    for (auto const* part:(*ak8jets)){
      if (!ParticleSelectionHelpers::isJetForHEMVeto(part)) continue;

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

bool EventFilterHandler::constructHLTPaths(SimEventHandler const* simEventHandler){
  bool isData = SampleHelpers::checkSampleIsData(currentTree->sampleIdentifier);
  if (!isData && simEventHandler){
    if (!simEventHandler->isAlreadyCached()){
      if (this->verbosity>=TVar::ERROR) MELAerr << "EventFilterHandler::constructHLTPaths: Need to update the SimEventHandler object first!" << endl;
      assert(0);
    }
  }

#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_HLTpaths_##NAME, itEnd_HLTpaths_##NAME;
  HLTTRIGGERPATH_VARIABLES;
#undef HLTTRIGGERPATH_VARIABLE
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) TYPE const* NAME = nullptr;
  RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(EventFilterHandler::colName_HLTpaths + "_" + #NAME, &itBegin_HLTpaths_##NAME, &itEnd_HLTpaths_##NAME);
  HLTTRIGGERPATH_VARIABLES;
#undef HLTTRIGGERPATH_VARIABLE
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumed(#NAME, NAME);
  if (isData){
    RUNLUMIEVENT_VARIABLES;
  }
#undef RUNLUMIEVENT_VARIABLE
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
  {
    cms3_triggerIndex_t itrig=0;
    unsigned int RunNumber_sim=0;
    bool set_RunNumber_sim=false;
    while (it_HLTpaths_name != itEnd_HLTpaths_name){
      product_HLTpaths.push_back(new HLTTriggerPathObject());
      HLTTriggerPathObject*& obj = product_HLTpaths.back();

#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) obj->NAME = *it_HLTpaths_##NAME;
      HLTTRIGGERPATH_VARIABLES;
#undef HLTTRIGGERPATH_VARIABLE

      // Set list index as its unique identifier
      obj->setUniqueIdentifier(itrig);

      // Associate trigger objects
      obj->setTriggerObjects(product_triggerobjects);

      bool isValid = true;
      HLTTriggerPathProperties const* hltprop = nullptr;
      bool hasRunRangeExclusions = checkHLTPathRunRanges && TriggerHelpers::hasRunRangeExclusions(obj->name, &hltprop);
      if (hasRunRangeExclusions){
        assert(hltprop!=nullptr);
        if (!isData && !set_RunNumber_sim){
          if (!simEventHandler){
            if (this->verbosity>=TVar::ERROR) MELAerr << "EventFilterHandler::constructHLTPaths: simEventHandler is needed to determine run range exclusions!" << endl;
            assert(0);
          }
          RunNumber_sim = simEventHandler->getChosenRunNumber();
          set_RunNumber_sim = true;
        }
        isValid = hltprop->testRun((isData ? *RunNumber : RunNumber_sim));
      }
      obj->setValid(isValid);

      itrig++;
#define HLTTRIGGERPATH_VARIABLE(TYPE, NAME, DEFVAL) it_HLTpaths_##NAME++;
      HLTTRIGGERPATH_VARIABLES;
#undef HLTTRIGGERPATH_VARIABLE
    }
  }

  return true;
}

bool EventFilterHandler::constructTriggerObjects(){
#define TRIGGEROBJECT_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME;
  TRIGGEROBJECT_VARIABLES;
#undef TRIGGEROBJECT_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define TRIGGEROBJECT_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(EventFilterHandler::colName_triggerobjects + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  TRIGGEROBJECT_VARIABLES;
#undef TRIGGEROBJECT_VARIABLE
  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "EventFilterHandler::constructTriggerObjects: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "EventFilterHandler::constructTriggerObjects: All variables are set up!" << endl;

  size_t n_products = (itEnd_type - itBegin_type);
  product_triggerobjects.reserve(n_products);
#define TRIGGEROBJECT_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  TRIGGEROBJECT_VARIABLES;
#undef TRIGGEROBJECT_VARIABLE
  {
    size_t ip=0;
    while (it_type != itEnd_type){
      if (this->verbosity>=TVar::DEBUG) MELAout << "EventFilterHandler::constructTriggerObjects: Attempting trigger object " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_pt, *it_eta, *it_phi, *it_mass);
      product_triggerobjects.push_back(new TriggerObject(*it_type, momentum));
      TriggerObject* const& obj = product_triggerobjects.back();

#define TRIGGEROBJECT_VARIABLE(TYPE, NAME, DEFVAL) obj->extras.NAME = *it_##NAME;
      TRIGGEROBJECT_EXTRA_VARIABLES;
#undef TRIGGEROBJECT_VARIABLE

      // Set particle index as its unique identifier
      obj->setUniqueIdentifier(ip);

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define TRIGGEROBJECT_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      TRIGGEROBJECT_VARIABLES;
#undef TRIGGEROBJECT_VARIABLE
    }
  }

  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(product_triggerobjects);

  return true;
}

bool EventFilterHandler::constructMETFilters(){
  auto strmetfilters = EventFilterHandler::acquireMETFilterFlags(currentTree, EventFilterHandler::nMETFilterCutTypes);
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

bool EventFilterHandler::accumulateRunLumiEventBlock(){
  bool isData = SampleHelpers::checkSampleIsData(currentTree->sampleIdentifier);
  if (!isData || !trackDataEvents || !checkUniqueDataEvent){
    product_uniqueEvent = true;
    if (!isData || !trackDataEvents) return true;
  }

  bool allVariablesPresent = true;
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) TYPE const* NAME = nullptr; allVariablesPresent &= this->getConsumed(#NAME, NAME);
  RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE
  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "EventFilterHandler::accumulateRunLumiEventBlock: Not all variables are consumed properly!" << endl;
    assert(0);
  }
  if (this->verbosity>=TVar::DEBUG) MELAout << "EventFilterHandler::accumulateRunLumiEventBlock: All variables are set up!" << endl;

  auto it_run = era_dataeventblock_map.find(*RunNumber);
  if (it_run == era_dataeventblock_map.end()){
    if (this->verbosity>=TVar::DEBUG) MELAout << "EventFilterHandler::accumulateRunLumiEventBlock: Run " << *RunNumber << " is new." << endl;
    era_dataeventblock_map[*RunNumber] = std::unordered_map<unsigned int, std::vector<unsigned long long>>();
    it_run = era_dataeventblock_map.find(*RunNumber);
  }
  else if (this->verbosity>=TVar::DEBUG) MELAout << "EventFilterHandler::accumulateRunLumiEventBlock: Run " << *RunNumber << " is already in the tracked list." << endl;

  auto it_lumi = it_run->second.find(*LuminosityBlock);
  if (it_lumi == it_run->second.end()){
    if (this->verbosity>=TVar::DEBUG) MELAout << "EventFilterHandler::accumulateRunLumiEventBlock: Lumi. section " << *LuminosityBlock << " is new." << endl;
    it_run->second[*LuminosityBlock] = std::vector<unsigned long long>();
    it_lumi = it_run->second.find(*LuminosityBlock);
  }
  else if (this->verbosity>=TVar::DEBUG) MELAout << "EventFilterHandler::accumulateRunLumiEventBlock: Lumi. section " << *LuminosityBlock << " is already in the tracked list." << endl;

  if (checkUniqueDataEvent){
    if (this->verbosity>=TVar::DEBUG) MELAout << "EventFilterHandler::accumulateRunLumiEventBlock: Checking if the event is unique..." << endl;
    auto it_event = std::find(it_lumi->second.begin(), it_lumi->second.end(), *EventNumber);
    if (it_event == it_lumi->second.end()){
      if (this->verbosity>=TVar::DEBUG) MELAout << "EventFilterHandler::accumulateRunLumiEventBlock: Event " << *EventNumber << " is unique." << endl;
      it_lumi->second.push_back(*EventNumber);
      product_uniqueEvent = true;
    }
    else{
      if (this->verbosity>=TVar::DEBUG) MELAout << "EventFilterHandler::accumulateRunLumiEventBlock: Event " << *EventNumber << " is already covered." << endl;
      product_uniqueEvent = false;
    }
  }
  else{
    if (this->verbosity>=TVar::DEBUG) MELAout << "EventFilterHandler::accumulateRunLumiEventBlock: No checking is performed for event uniquenes." << endl;
    it_lumi->second.push_back(*EventNumber);
  }

  return true;
}

std::vector<std::string> EventFilterHandler::acquireMETFilterFlags(BaseTree* intree, EventFilterHandler::METFilterCutType const& cuttype){
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
    if (!SampleHelpers::checkSampleIsFastSim(intree->sampleIdentifier)){ // For data or non-FS MC
      if (cuttype>=EventFilterHandler::kMETFilters_Standard) res.push_back("globalSuperTightHalo2016Filter");
      if (cuttype>=EventFilterHandler::kMETFilters_Tight) res.push_back("globalTightHalo2016Filter");
    }
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
    if (!SampleHelpers::checkSampleIsFastSim(intree->sampleIdentifier)){ // For data or non-FS MC
      if (cuttype>=EventFilterHandler::kMETFilters_Standard) res.push_back("globalSuperTightHalo2016Filter");
      if (cuttype>=EventFilterHandler::kMETFilters_Tight) res.push_back("globalTightHalo2016Filter");
    }
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
    if (!SampleHelpers::checkSampleIsFastSim(intree->sampleIdentifier)){ // For data or non-FS MC
      if (cuttype>=EventFilterHandler::kMETFilters_Standard) res.push_back("globalSuperTightHalo2016Filter");
      if (cuttype>=EventFilterHandler::kMETFilters_Tight) res.push_back("globalTightHalo2016Filter");
    }
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

  // Trigger objects
  if (trackTriggerObjects){
#define TRIGGEROBJECT_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(EventFilterHandler::colName_triggerobjects + "_" + #NAME, nullptr); this->addConsumed<std::vector<TYPE>*>(EventFilterHandler::colName_triggerobjects + "_" + #NAME);
    TRIGGEROBJECT_VARIABLES;
#undef TRIGGEROBJECT_VARIABLE
  }

  // Book MET filters
  auto strmetfilters = EventFilterHandler::acquireMETFilterFlags(tree, EventFilterHandler::nMETFilterCutTypes);
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
