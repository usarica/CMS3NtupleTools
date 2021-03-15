#include <cassert>

#include <CMS3/Dictionaries/interface/GlobalCollectionNames.h>
#include "ParticleObjectHelpers.h"

#include "SamplesCore.h"
#include "HelperFunctions.h"
#include "ReweightingFunctions.h"
#include "GenInfoHandler.h"
#include "GenParticleSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


const std::string GenInfoHandler::colName_lheparticles = GlobalCollectionNames::colName_lheparticles;
const std::string GenInfoHandler::colName_genparticles = GlobalCollectionNames::colName_genparticles;
const std::string GenInfoHandler::colName_genak4jets = GlobalCollectionNames::colName_genak4jets;
const std::string GenInfoHandler::colName_genak8jets = GlobalCollectionNames::colName_genak8jets;

GenInfoHandler::GenInfoHandler() :
  IvyBase(),

  acquireCoreGenInfo(true),
  acquireLHEMEWeights(true),
  acquireLHEParticles(true),
  acquireGenParticles(true),
  acquireGenAK4Jets(false),
  acquireGenAK8Jets(false),
  doGenJetsCleaning(true),
  doGenJetsVDecayCleaning(false),
  allowLargeGenWeightRemoval(false),

  genWeightException(SampleHelpers::nGenWeightExceptionType),
  abs_genWeight_default_thr(nullptr),

  genInfo(nullptr)
{}

void GenInfoHandler::clear(){
  this->resetCache();

  delete genInfo; genInfo=nullptr;

  for (auto& part:lheparticles) delete part;
  lheparticles.clear();

  for (auto& part:genparticles) delete part;
  genparticles.clear();

  for (auto& part:genak4jets) delete part;
  genak4jets.clear();

  for (auto& part:genak8jets) delete part;
  genak8jets.clear();
}

bool GenInfoHandler::constructGenInfo(SystematicsHelpers::SystematicVariationTypes const& syst){
  if (this->isAlreadyCached()) return true;

  clear();
  if (!currentTree) return false;

  bool res = constructCoreGenInfo(syst) && constructLHEParticles() && constructGenParticles() && constructGenAK4Jets() && constructGenAK8Jets() && applyGenJetCleaning();

  if (res) this->cacheEvent();
  return res;
}

bool GenInfoHandler::constructCoreGenInfo(SystematicsHelpers::SystematicVariationTypes const& syst){
  // Use non-const pointer here because we might have to modify some of the weights
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) TYPE* NAME = nullptr;
  GENINFO_VARIABLES;
#undef GENINFO_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumed(#NAME, NAME);
  if (acquireCoreGenInfo){
    GENINFO_VARIABLES;
  }
#undef GENINFO_VARIABLE

  std::unordered_map<TString, float const*> kfactorlist;
  for (TString const& strkfactor:tree_kfactorlist_map[currentTree]){
    kfactorlist[strkfactor] = nullptr;
    allVariablesPresent &= this->getConsumed(strkfactor, kfactorlist.find(strkfactor)->second);
    if (!(kfactorlist.find(strkfactor)->second)){
      if (this->verbosity>=TVar::ERROR) MELAerr << "GenInfoHandler::constructCoreGenInfo: K factor handle for " << strkfactor << " is null!" << endl;
      assert(0);
    }
  }

  std::unordered_map<TString, float const*> MElist;
  if (acquireLHEMEWeights){
    for (TString const& strme:tree_MElist_map[currentTree]){
      MElist[strme] = nullptr;
      allVariablesPresent &= this->getConsumed(strme, MElist.find(strme)->second);
      if (!(MElist.find(strme)->second)){
        if (this->verbosity>=TVar::ERROR) MELAerr << "GenInfoHandler::constructCoreGenInfo: ME handle for " << strme << " is null!" << endl;
        assert(0);
      }
    }
  }

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "GenInfoHandler::constructCoreGenInfo: Not all variables are consumed properly!" << endl;
    assert(0);
  }
  if (this->verbosity>=TVar::DEBUG) MELAout << "GenInfoHandler::constructCoreGenInfo: All variables are set up!" << endl;

  genInfo = new GenInfoObject();
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) if (NAME) genInfo->extras.NAME = *NAME;
  if (acquireCoreGenInfo){
    GENINFO_VARIABLES;
  }
#undef GENINFO_VARIABLE

  genInfo->setSystematic(syst);

  for (auto it:kfactorlist) genInfo->extras.Kfactors[it.first] = (it.second ? *(it.second) : 1.f);
  for (auto it:MElist) genInfo->extras.LHE_ME_weights[it.first] = (it.second ? *(it.second) : 0.f);

  if (
    genWeightException == SampleHelpers::kLargeDefaultGenWeight
    &&
    (
      (abs_genWeight_default_thr && *abs_genWeight_default_thr>0.f && std::abs(genInfo->extras.genHEPMCweight_default)>*abs_genWeight_default_thr)
      ||
      (LHEweight_scaledOriginalWeight_default/* && *LHEweight_scaledOriginalWeight_default != 0.f*/ && *LHEweight_scaledOriginalWeight_default != *genHEPMCweight_default)
      )
    ){
    if (this->verbosity>=TVar::INFO) MELAout
      << "GenInfoHandler::constructCoreGenInfo: genHEPMCweight_default = " << genInfo->extras.genHEPMCweight_default
      << " (original genHEPMCweight_NNPDF30 = " << genInfo->extras.genHEPMCweight_NNPDF30 << ")"
      << " is invalid! A threshold of " << *abs_genWeight_default_thr << " may be applied." << endl;
    // Attempt to fix the weights
    if (LHEweight_scaledOriginalWeight_default){
      *genHEPMCweight_default = genInfo->extras.genHEPMCweight_default = *LHEweight_scaledOriginalWeight_default;
      *genHEPMCweight_NNPDF30 = genInfo->extras.genHEPMCweight_NNPDF30 = *LHEweight_scaledOriginalWeight_NNPDF30;
    }
    else{ *genHEPMCweight_default = *genHEPMCweight_NNPDF30 = genInfo->extras.genHEPMCweight_default = genInfo->extras.genHEPMCweight_NNPDF30 = 0.f; }
  }

  return true;
}
bool GenInfoHandler::constructLHEParticles(){
  bool doLHEParticles = acquireLHEParticles && tree_lheparticles_present_map[currentTree];
  if (!doLHEParticles) return true;

#define LHEPARTICLE_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_lheparticles_##NAME, itEnd_lheparticles_##NAME;
  LHEPARTICLE_VARIABLES;
#undef LHEPARTICLE_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define LHEPARTICLE_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(GenInfoHandler::colName_lheparticles + "_" + #NAME, &itBegin_lheparticles_##NAME, &itEnd_lheparticles_##NAME);
  LHEPARTICLE_VARIABLES;
#undef LHEPARTICLE_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "GenInfoHandler::constructLHEParticles: Not all variables are consumed properly!" << endl;
    assert(0);
  }
  if (this->verbosity>=TVar::DEBUG) MELAout << "GenInfoHandler::constructLHEParticles: All variables are set up!" << endl;

  size_t nlheparticles = (itEnd_lheparticles_id - itBegin_lheparticles_id);
  lheparticles.reserve(nlheparticles);
  std::vector<std::pair<int, int>> mother_index_pairs; mother_index_pairs.reserve(nlheparticles);
#define LHEPARTICLE_VARIABLE(TYPE, NAME, DEFVAL) auto it_lheparticles_##NAME = itBegin_lheparticles_##NAME;
  LHEPARTICLE_VARIABLES;
#undef LHEPARTICLE_VARIABLE
  {
    size_t ip=0;
    while (it_lheparticles_id != itEnd_lheparticles_id){
      if (this->verbosity>=TVar::DEBUG) MELAout << "GenInfoHandler::constructLHEParticles: Attempting LHE particle " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum(*it_lheparticles_px, *it_lheparticles_py, *it_lheparticles_pz, *it_lheparticles_E);
      lheparticles.push_back(new LHEParticleObject(*it_lheparticles_id, *it_lheparticles_status, momentum));
      mother_index_pairs.emplace_back(*it_lheparticles_mother0_index, *it_lheparticles_mother1_index);

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define LHEPARTICLE_VARIABLE(TYPE, NAME, DEFVAL) it_lheparticles_##NAME++;
      LHEPARTICLE_VARIABLES;
#undef LHEPARTICLE_VARIABLE
    }
  }

  {
    assert(mother_index_pairs.size() == lheparticles.size());

    auto it_lheparticles = lheparticles.begin();
    auto it_mother_index_pairs = mother_index_pairs.begin();
    while (it_lheparticles != lheparticles.end()){
      auto*& part = *it_lheparticles;

      int const& imom = it_mother_index_pairs->first;
      if (imom>=0){
        assert((unsigned int) imom<nlheparticles);
        auto*& mother = lheparticles.at(imom);
        part->addMother(mother);
        mother->addDaughter(part);
      }

      int const& jmom = it_mother_index_pairs->second;
      if (jmom>=0){
        assert((unsigned int) jmom<nlheparticles);
        auto*& mother = lheparticles.at(jmom);
        part->addMother(mother);
        mother->addDaughter(part);
      }

      it_lheparticles++;
      it_mother_index_pairs++;
    }
  }

  return true;
}
bool GenInfoHandler::constructGenParticles(){
  if (!acquireGenParticles) return true;

#define GENPARTICLE_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_genparticles_##NAME, itEnd_genparticles_##NAME;
  GENPARTICLE_VARIABLES;
#undef GENPARTICLE_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define GENPARTICLE_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(GenInfoHandler::colName_genparticles + "_" + #NAME, &itBegin_genparticles_##NAME, &itEnd_genparticles_##NAME);
  GENPARTICLE_VARIABLES;
#undef GENPARTICLE_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "GenInfoHandler::constructGenParticles: Not all variables are consumed properly!" << endl;
    assert(0);
  }
  if (this->verbosity>=TVar::DEBUG) MELAout << "GenInfoHandler::constructGenParticles: All variables are set up!" << endl;

  size_t ngenparticles = (itEnd_genparticles_id - itBegin_genparticles_id);
  genparticles.reserve(ngenparticles);
  std::vector<std::pair<int, int>> mother_index_pairs; mother_index_pairs.reserve(ngenparticles);
#define GENPARTICLE_VARIABLE(TYPE, NAME, DEFVAL) auto it_genparticles_##NAME = itBegin_genparticles_##NAME;
  GENPARTICLE_VARIABLES;
#undef GENPARTICLE_VARIABLE
  {
    size_t ip=0;
    while (it_genparticles_id != itEnd_genparticles_id){
      if (this->verbosity>=TVar::DEBUG) MELAout << "GenInfoHandler::constructGenParticles: Attempting LHE particle " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_genparticles_pt, *it_genparticles_eta, *it_genparticles_phi, *it_genparticles_mass); // Yes you have to do this on a separate line because CMSSW...
      genparticles.push_back(new GenParticleObject(*it_genparticles_id, *it_genparticles_status, momentum));
      mother_index_pairs.emplace_back(*it_genparticles_mom0_index, *it_genparticles_mom1_index);
      GenParticleObject*& obj = genparticles.back();

      // Set extras
#define GENPARTICLE_VARIABLE(TYPE, NAME, DEFVAL) obj->extras.NAME = *it_genparticles_##NAME;
      GENPARTICLE_EXTRA_VARIABLES;
#undef GENPARTICLE_VARIABLE

      // Set selection bits
      GenParticleSelectionHelpers::setSelectionBits(*obj);

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define GENPARTICLE_VARIABLE(TYPE, NAME, DEFVAL) it_genparticles_##NAME++;
      GENPARTICLE_VARIABLES;
#undef GENPARTICLE_VARIABLE
    }
  }

  {
    assert(mother_index_pairs.size() == genparticles.size());

    auto it_genparticles = genparticles.begin();
    auto it_mother_index_pairs = mother_index_pairs.begin();
    while (it_genparticles != genparticles.end()){
      auto*& part = *it_genparticles;

      int const& imom = it_mother_index_pairs->first;
      if (imom>=0){
        assert((unsigned int) imom<ngenparticles);
        auto*& mother = genparticles.at(imom);
        part->addMother(mother);
        mother->addDaughter(part);
      }

      int const& jmom = it_mother_index_pairs->second;
      if (jmom>=0){
        assert((unsigned int) jmom<ngenparticles);
        auto*& mother = genparticles.at(jmom);
        part->addMother(mother);
        mother->addDaughter(part);
      }

      it_genparticles++;
      it_mother_index_pairs++;
    }
  }

  return true;
}

bool GenInfoHandler::constructGenAK4Jets(){
  if (!acquireGenAK4Jets) return true;

#define GENJET_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_genak4jets_##NAME, itEnd_genak4jets_##NAME;
  GENJET_VARIABLES;
#undef GENJET_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define GENJET_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(GenInfoHandler::colName_genak4jets + "_" + #NAME, &itBegin_genak4jets_##NAME, &itEnd_genak4jets_##NAME);
  GENJET_VARIABLES;
#undef GENJET_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "GenInfoHandler::constructGenAK4Jets: Not all variables are consumed properly!" << endl;
    assert(0);
  }
  if (this->verbosity>=TVar::DEBUG) MELAout << "GenInfoHandler::constructGenAK4Jets: All variables are set up!" << endl;

  size_t ngenak4jets = (itEnd_genak4jets_pt - itBegin_genak4jets_pt);
  genak4jets.reserve(ngenak4jets);
#define GENJET_VARIABLE(TYPE, NAME, DEFVAL) auto it_genak4jets_##NAME = itBegin_genak4jets_##NAME;
  GENJET_VARIABLES;
#undef GENJET_VARIABLE
  {
    size_t ip=0;
    while (it_genak4jets_pt != itEnd_genak4jets_pt){
      if (this->verbosity>=TVar::DEBUG) MELAout << "GenInfoHandler::constructGenAK4Jets: Attempting LHE particle " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_genak4jets_pt, *it_genak4jets_eta, *it_genak4jets_phi, *it_genak4jets_mass); // Yes you have to do this on a separate line because CMSSW...
      genak4jets.push_back(new GenJetObject(0, momentum));

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define GENJET_VARIABLE(TYPE, NAME, DEFVAL) it_genak4jets_##NAME++;
      GENJET_VARIABLES;
#undef GENJET_VARIABLE
    }
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(genak4jets);

  return true;
}
bool GenInfoHandler::constructGenAK8Jets(){
  if (!acquireGenAK8Jets) return true;

#define GENJET_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_genak8jets_##NAME, itEnd_genak8jets_##NAME;
  GENJET_VARIABLES;
#undef GENJET_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define GENJET_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(GenInfoHandler::colName_genak8jets + "_" + #NAME, &itBegin_genak8jets_##NAME, &itEnd_genak8jets_##NAME);
  GENJET_VARIABLES;
#undef GENJET_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "GenInfoHandler::constructGenAK8Jets: Not all variables are consumed properly!" << endl;
    assert(0);
  }
  if (this->verbosity>=TVar::DEBUG) MELAout << "GenInfoHandler::constructGenAK8Jets: All variables are set up!" << endl;

  size_t ngenak8jets = (itEnd_genak8jets_pt - itBegin_genak8jets_pt);
  genak8jets.reserve(ngenak8jets);
#define GENJET_VARIABLE(TYPE, NAME, DEFVAL) auto it_genak8jets_##NAME = itBegin_genak8jets_##NAME;
  GENJET_VARIABLES;
#undef GENJET_VARIABLE
  {
    size_t ip=0;
    while (it_genak8jets_pt != itEnd_genak8jets_pt){
      if (this->verbosity>=TVar::DEBUG) MELAout << "GenInfoHandler::constructGenAK8Jets: Attempting LHE particle " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_genak8jets_pt, *it_genak8jets_eta, *it_genak8jets_phi, *it_genak8jets_mass); // Yes you have to do this on a separate line because CMSSW...
      genak8jets.push_back(new GenJetObject(0, momentum));

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define GENJET_VARIABLE(TYPE, NAME, DEFVAL) it_genak8jets_##NAME++;
      GENJET_VARIABLES;
#undef GENJET_VARIABLE
    }
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(genak8jets);

  return true;
}

bool GenInfoHandler::applyGenJetCleaning(){
  if (!acquireGenAK4Jets && !acquireGenAK8Jets) return true;
  if (!doGenJetsCleaning) return true;

  std::vector<ParticleObject*> cleaningparticles; cleaningparticles.reserve(genparticles.size());
  // Accumulate all prompt leptons and photons, or hard-process taus
  // Require the kHardPromptFinalParticle selection bit in order to avoid soft particles,
  // which could bias jet cleaning through combinatorics.
  // Also add hard-process taus to cleaning.
  for (auto const& part:genparticles){
    auto const& extras = part->extras;
    cms3_id_t const& part_id = part->pdgId();
    ParticleObject::LorentzVector_t::Scalar const part_pt = part->pt();
    if (
      part->testSelectionBit(GenParticleSelectionHelpers::kHardPromptFinalVisibleParticle)
      ||
      (extras.isHardProcess && std::abs(part_id)==15 && part_pt>=GenParticleSelectionHelpers::ptThr_hardparticle_lepton)
      ) cleaningparticles.push_back(part);
  }
  // Accumulate all V->qq decays from LHE level
  if (doGenJetsVDecayCleaning){
    for (auto const& part:lheparticles){
      cms3_id_t const& part_id = part->pdgId();
      cms3_genstatus_t const& part_st = part->status();
      bool hasVmothers = false;
      for (auto const& mother:part->getMothers()){
        if (PDGHelpers::isAWBoson(mother->pdgId()) || PDGHelpers::isAZBoson(mother->pdgId())){
          hasVmothers = true;
          break;
        }
      }
      if (part_st==1 && PDGHelpers::isAQuark(part_id) && hasVmothers) cleaningparticles.push_back(part);
    }
  }

  // Clean gen. ak4 jets
  std::vector<GenJetObject*> genak4jets_new; genak4jets_new.reserve(genak4jets.size());
  for (auto& jet:genak4jets){
    bool doSkip=false;
    for (auto const& part:cleaningparticles){
      if (jet->deltaR(part)<0.4){ doSkip=true; break; }
    }
    if (!doSkip) genak4jets_new.push_back(jet);
    else delete jet;
  }
  genak4jets = genak4jets_new;

  // Clean gen. ak8 jets
  std::vector<GenJetObject*> genak8jets_new; genak8jets_new.reserve(genak8jets.size());
  for (auto& jet:genak8jets){
    bool doSkip=false;
    for (auto const& part:cleaningparticles){
      if (jet->deltaR(part)<0.8){ doSkip=true; break; }
    }
    if (!doSkip) genak8jets_new.push_back(jet);
    else delete jet;
  }
  genak8jets = genak8jets_new;

  return true;
}

bool GenInfoHandler::wrapTree(BaseTree* tree){
  if (!tree) return false;

  abs_genWeight_default_thr = nullptr;
  if (SampleHelpers::hasGenWeightException(tree->sampleIdentifier, SampleHelpers::theDataYear, this->genWeightException)) MELAout
    << "GenInfoHandler::wrapTree: Warning! Sample " << tree->sampleIdentifier << " has a gen. weight exception of type " << this->genWeightException << "."
    << endl;

  bool res = IvyBase::wrapTree(tree);

  if (this->genWeightException == SampleHelpers::kLargeDefaultGenWeight){
    auto it_thr = abs_genWeight_default_thr_map.find(currentTree);
    if (it_thr != abs_genWeight_default_thr_map.cend()) abs_genWeight_default_thr = &(it_thr->second);
    else res &= determineWeightThresholds();
  }

  return res;
}

void GenInfoHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

  std::vector<TString> allbranchnames; tree->getValidBranchNamesWithoutAlias(allbranchnames, false);

  if (acquireCoreGenInfo){
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<TYPE>(#NAME, DEFVAL); this->addConsumed<TYPE>(#NAME);
    GENINFO_CORE_VARIABLES;
#undef GENINFO_VARIABLE
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) \
if (HelperFunctions::checkListVariable<TString>(allbranchnames, #NAME)) tree->bookBranch<TYPE>(#NAME, DEFVAL); \
this->addConsumed<TYPE>(#NAME); this->defineConsumedSloppy(#NAME);
    GENINFO_EXTRA_VARIABLES;
#undef GENINFO_VARIABLE
  }

  // K factor and ME reweighting branches are defined as sloppy
  std::vector<TString> kfactorlist;
  std::vector<TString> melist;
  bool has_lheparticles=false;
  for (TString const& bname:allbranchnames){
    if (bname.Contains("KFactor")){
      tree->bookBranch<float>(bname, 1.f);
      this->addConsumed<float>(bname);
      this->defineConsumedSloppy(bname);
      kfactorlist.push_back(bname);
    }
    if (acquireLHEMEWeights && (bname.Contains("p_Gen") || bname.Contains("LHECandMass"))){
      tree->bookBranch<float>(bname, 0.f);
      this->addConsumed<float>(bname);
      this->defineConsumedSloppy(bname);
      melist.push_back(bname);
    }
    else if (acquireLHEParticles && bname.Contains(colName_lheparticles)) has_lheparticles = true;
  }
  tree_kfactorlist_map[tree] = kfactorlist;
  tree_MElist_map[tree] = melist;
  tree_lheparticles_present_map[tree] = has_lheparticles;

  if (has_lheparticles){
#define LHEPARTICLE_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(colName_lheparticles+ "_" + #NAME, nullptr);
    LHEPARTICLE_VARIABLES;
#undef LHEPARTICLE_VARIABLE
  }
#define LHEPARTICLE_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(colName_lheparticles+ "_" + #NAME); this->defineConsumedSloppy(colName_lheparticles+ "_" + #NAME);
  LHEPARTICLE_VARIABLES;
#undef LHEPARTICLE_VARIABLE

  if (acquireGenParticles){
#define GENPARTICLE_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(colName_genparticles+ "_" + #NAME, nullptr);
    GENPARTICLE_VARIABLES;
#undef GENPARTICLE_VARIABLE
  }
#define GENPARTICLE_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(colName_genparticles+ "_" + #NAME); this->defineConsumedSloppy(colName_genparticles+ "_" + #NAME);
  GENPARTICLE_VARIABLES;
#undef GENPARTICLE_VARIABLE

  if (acquireGenAK4Jets){
#define GENJET_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(colName_genak4jets+ "_" + #NAME, nullptr);
    GENJET_VARIABLES;
#undef GENJET_VARIABLE
  }
#define GENJET_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(colName_genak4jets+ "_" + #NAME); this->defineConsumedSloppy(colName_genak4jets+ "_" + #NAME);
  GENJET_VARIABLES;
#undef GENJET_VARIABLE

  if (acquireGenAK8Jets){
#define GENJET_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(colName_genak8jets+ "_" + #NAME, nullptr);
    GENJET_VARIABLES;
#undef GENJET_VARIABLE
  }
#define GENJET_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(colName_genak8jets+ "_" + #NAME); this->defineConsumedSloppy(colName_genak8jets+ "_" + #NAME);
  GENJET_VARIABLES;
#undef GENJET_VARIABLE
}

bool GenInfoHandler::determineWeightThresholds(){
  if (!allowLargeGenWeightRemoval) return true;

  if (!currentTree) return false;
  abs_genWeight_default_thr = nullptr;
  if (!acquireCoreGenInfo){
    if (this->verbosity>=TVar::ERROR) MELAerr << "GenInfoHandler::determineWeightThresholds: In order to determine weight thresholds, you need to set acquireCoreGenInfo=true." << endl;
    assert(0);
    return false;
  }

  float* genHEPMCweight_default = nullptr;

  bool allVariablesPresent = this->getConsumed("genHEPMCweight_default", genHEPMCweight_default);
  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "GenInfoHandler::determineWeightThresholds: Not all variables are consumed properly!" << endl;
    assert(0);
  }
  if (this->verbosity>=TVar::DEBUG) MELAout << "GenInfoHandler::determineWeightThresholds: All variables are set up!" << endl;

  std::vector<float*> tmpvec_wgts{ genHEPMCweight_default };
  abs_genWeight_default_thr_map[currentTree] = ReweightingFunctions::getAbsWeightThresholdByNeff(currentTree, tmpvec_wgts, ReweightingFunctions::getSimpleWeight, 250000., this->verbosity);
  abs_genWeight_default_thr = &(abs_genWeight_default_thr_map[currentTree]);

  return true;
}
