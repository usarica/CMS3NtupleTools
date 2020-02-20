#include <cassert>
#include "GenInfoHandler.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


const std::string GenInfoHandler::colName_lheparticles = "lheparticles";
const std::string GenInfoHandler::colName_genparticles = "genparticles";

GenInfoHandler::GenInfoHandler() :
  IvyBase(),

  acquireCoreGenInfo(true),
  acquireLHEMEWeights(true),
  acquireLHEParticles(true),
  acquireGenParticles(true),

  genInfo(nullptr)
{}

void GenInfoHandler::clear(){
  delete genInfo; genInfo=nullptr;

  for (auto*& part:lheparticles) delete part;
  lheparticles.clear();

  for (auto*& part:genparticles) delete part;
  genparticles.clear();
}

bool GenInfoHandler::constructGenInfo(SystematicsHelpers::SystematicVariationTypes const& syst){
  clear();
  if (!currentTree) return false;

  bool res = constructCoreGenInfo(syst) && constructLHEParticles() && constructGenParticles();
  return res;
}

bool GenInfoHandler::constructCoreGenInfo(SystematicsHelpers::SystematicVariationTypes const& syst){
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) TYPE const* NAME = nullptr;
  GENINFO_VARIABLES;
#undef GENINFO_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
  if (acquireCoreGenInfo){
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumed(#NAME, NAME);
    GENINFO_VARIABLES;
#undef GENINFO_VARIABLE
  }

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
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) genInfo->extras.NAME = (NAME ? *NAME : TYPE(DEFVAL));
  GENINFO_VARIABLES;
#undef GENINFO_VARIABLE
  genInfo->setSystematic(syst);

  for (auto it:kfactorlist) genInfo->extras.Kfactors[it.first] = (it.second ? *(it.second) : 1.f);
  for (auto it:MElist) genInfo->extras.LHE_ME_weights[it.first] = (it.second ? *(it.second) : 0.f);

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

void GenInfoHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

  if (acquireCoreGenInfo){
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<TYPE>(#NAME, DEFVAL); this->addConsumed<TYPE>(#NAME);
    GENINFO_VARIABLES;
#undef GENINFO_VARIABLE
  }

  // K factor and ME reweighting branches are defined as sloppy
  std::vector<TString> allbranchnames; tree->getValidBranchNamesWithoutAlias(allbranchnames, false);
  std::vector<TString> kfactorlist;
  std::vector<TString> melist;
  bool has_lheparticles=false;
  for (TString const& bname : allbranchnames){
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
}
