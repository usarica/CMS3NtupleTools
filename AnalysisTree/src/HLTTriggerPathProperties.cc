#include <cassert>
#include "HelperFunctions.h"
#include "SamplesCore.h"
#include "HLTTriggerPathProperties.h"
#include <CMS3/Dictionaries/interface/CMS3StreamHelpers.h>


using namespace std;
using namespace IvyStreamHelpers;


HLTTriggerPathProperties::HLTTriggerPathProperties(){}
HLTTriggerPathProperties::HLTTriggerPathProperties(std::string const& name_) :
  name(name_),
  hasCompositeFilters(false),
  TOexception(nTriggerObjectExceptionTypes)
{
  setupName();
  setup();
}
HLTTriggerPathProperties::HLTTriggerPathProperties(std::string const& name_, std::vector<HLTObjectProperties> const& triggerObjectProperties_) :
  name(name_),
  hasCompositeFilters(false),
  TOexception(nTriggerObjectExceptionTypes)
{
  setupName();
  for (auto const& props:triggerObjectProperties_) addObjectProperties(props);
  setup();
}
HLTTriggerPathProperties::HLTTriggerPathProperties(HLTTriggerPathProperties const& other) :
  name(other.name),
  triggerObjectProperties(other.triggerObjectProperties),
  hasCompositeFilters(other.hasCompositeFilters),
  TOexception(other.TOexception),
  excluded_runRange_list(other.excluded_runRange_list)
{}

void HLTTriggerPathProperties::setupName(){
  size_t pos = name.find("*");
  if (pos!=std::string::npos){
    if (pos != name.length()-1){
      IVYerr << "HLTTriggerPathProperties::HLTTriggerPathProperties: Trigger name " << name << " can only contain * as the last character. Please fix the name passed!" << endl;
      assert(0);
    }
    name = name.substr(0, pos);
    pos = pos-1;
    if (name.at(pos)!='v' || name.at(pos-1)!='_'){
      IVYerr << "HLTTriggerPathProperties::HLTTriggerPathProperties: Trigger name " << name << " has to end with _v. Please fix the name passed!" << endl;
      assert(0);
    }
  }
}

void HLTTriggerPathProperties::setup(){
  for (auto& it:triggerObjectProperties){
    auto& props = it.second;
    HLTObjectProperties::sortByMoreRestrictive(props);
    // The list can be expanded
    if (it.first == HLTObjectProperties::kAK4DiJetSumWithDEtaDPhi){
      hasCompositeFilters = true;
    }
  }
}

void HLTTriggerPathProperties::addObjectProperties(HLTObjectProperties const& props){
  auto const& type=props.getType();
  if (triggerObjectProperties.find(type)==triggerObjectProperties.end()) triggerObjectProperties[type]=std::vector<HLTObjectProperties>();
  triggerObjectProperties[type].push_back(props);
}

bool HLTTriggerPathProperties::isSameTrigger(std::string const& name_) const{
  return (name_.find(name)!=std::string::npos);
}

void HLTTriggerPathProperties::resetCuts(){
  for (auto& it:triggerObjectProperties){
    for (auto& props:it.second) props.resetCuts();
  }
}

bool HLTTriggerPathProperties::testCuts(
  std::vector<MuonObject const*> const& muons,
  std::vector<ElectronObject const*> const& electrons,
  std::vector<PhotonObject const*> const& photons,
  std::vector<AK4JetObject const*> const& ak4jets,
  std::vector<AK8JetObject const*> const& ak8jets,
  ParticleObject::LorentzVector_t const& pfmet_p4,
  ParticleObject::LorentzVector_t const& pfmet_nomu_p4,
  ParticleObject::LorentzVector_t const& ht_p4,
  ParticleObject::LorentzVector_t const& ht_nomu_p4
) const{
  bool res = true;
  std::vector<std::vector<AK4JetObject const*>> perms_ak4jets_passing;
  for (auto const& it:triggerObjectProperties){
    auto const& prop_type = it.first;
    auto const& props = it.second;
    switch (prop_type){
    case HLTObjectProperties::kMuon:
      res = res && testCutSet(prop_type, props, muons);
      break;
    case HLTObjectProperties::kElectron:
      res = res && testCutSet(prop_type, props, electrons);
      break;
    case HLTObjectProperties::kPhoton:
      res = res && testCutSet(prop_type, props, photons);
      break;
    case HLTObjectProperties::kAK4Jet:
      res = res && testCutSet(prop_type, props, ak4jets, (hasCompositeFilters ? &perms_ak4jets_passing : nullptr));
      break;
    case HLTObjectProperties::kAK4DiJetSumWithDEtaDPhi:
      // hasCompositeFilters is already set. A separate check follows below, but it needs the check on HLTObjectProperties::kAK4Jet and the perms_ak4jets_passing vector filled first.
      break;
    case HLTObjectProperties::kAK8Jet:
      res = res && testCutSet(prop_type, props, ak8jets);
      break;
    case HLTObjectProperties::kHT:
      for (auto const& propcuts:props) res = res && propcuts.testCuts(ht_p4, HLTObjectProperties::kHT);
      break;
    case HLTObjectProperties::kHT_NoMu:
      for (auto const& propcuts:props) res = res && propcuts.testCuts(ht_nomu_p4, HLTObjectProperties::kHT_NoMu);
      break;
    case HLTObjectProperties::kMET:
      for (auto const& propcuts:props) res = res && propcuts.testCuts(pfmet_p4, HLTObjectProperties::kMET);
      break;
    case HLTObjectProperties::kMET_NoMu:
      for (auto const& propcuts:props) res = res && propcuts.testCuts(pfmet_nomu_p4, HLTObjectProperties::kMET_NoMu);
      break;
    default:
      break;
    }
  }
  // Composite filters
  if (hasCompositeFilters && res){
    for (auto const& it:triggerObjectProperties){
      auto const& prop_type = it.first;
      auto const& props = it.second;
      switch (prop_type){
      case HLTObjectProperties::kAK4DiJetSumWithDEtaDPhi:
      {
        bool passAtLeastOnePerm = false;
        for (auto const& perm_ak4jets_passing:perms_ak4jets_passing){
          if (perm_ak4jets_passing.size()<2) continue;
          std::vector<std::vector<int>> perms_dijet;
          TNumericUtil::PermutationGenerator(perm_ak4jets_passing.size(), 2, perms_dijet, 0, 1);
          for (auto const& perm_dijet:perms_dijet){
            std::vector<AK4JetObject const*> tmpvec_ak4jets; tmpvec_ak4jets.reserve(2);
            for (auto const& ipos:perm_dijet) tmpvec_ak4jets.push_back(perm_ak4jets_passing.at(ipos));
            assert(tmpvec_ak4jets.size()==2);
            ParticleObject::LorentzVector_t p4_jetsum; p4_jetsum += tmpvec_ak4jets.front()->p4(); p4_jetsum += tmpvec_ak4jets.back()->p4();
            p4_jetsum = ParticleObject::PolarLorentzVector_t(p4_jetsum.Pt(), tmpvec_ak4jets.front()->eta() - tmpvec_ak4jets.back()->eta(), tmpvec_ak4jets.front()->phi() - tmpvec_ak4jets.back()->phi(), p4_jetsum.M());
            ParticleObject const tmppart(0, p4_jetsum);
            passAtLeastOnePerm = testCutSet(prop_type, props, std::vector<ParticleObject const*>{ &tmppart });
          }
        }
        res &= passAtLeastOnePerm;
        break;
      }
      default:
        break;
      }
    }
  }
  return res;
}

bool HLTTriggerPathProperties::testRun(unsigned int const& run) const{
  for (auto const& rr:excluded_runRange_list){
    if (run>=rr.first && run<=rr.second) return false;
  }
  return true;
}
