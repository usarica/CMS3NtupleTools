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


bool GenInfoHandler::constructGenInfo(SystematicsHelpers::SystematicVariationTypes const& syst){
  clear();
  if (!currentTree) return false;

  bool res = constructCoreGenInfo(syst) && constructLHEParticles();
  return res;
}

bool GenInfoHandler::constructCoreGenInfo(SystematicsHelpers::SystematicVariationTypes const& syst){
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME = DEFVAL;
  GENINFO_VARIABLES;
#undef GENINFO_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
  if (acquireCoreGenInfo){
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedValue(#NAME, NAME);
    GENINFO_VARIABLES;
#undef GENINFO_VARIABLE
  }

  std::unordered_map<TString, float> MElist;
  if (acquireLHEMEWeights){
    for (TString const& strme:tree_MElist_map[currentTree]){
      MElist[strme] = 0;
      allVariablesPresent &= this->getConsumedValue(strme, MElist.find(strme)->second);
    }
  }

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "GenInfoHandler::constructCoreGenInfo: Not all variables are consumed properly!" << endl;
    assert(0);
  }
  if (this->verbosity>=TVar::DEBUG) MELAout << "GenInfoHandler::constructCoreGenInfo: All variables are set up!" << endl;

  genInfo = new GenInfoObject();
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) genInfo->extras.NAME = NAME;
  GENINFO_VARIABLES;
#undef GENINFO_VARIABLE
  genInfo->setSystematic(syst);

  for (auto it:MElist) genInfo->extras.LHE_ME_weights[it.first] = it.second;

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

void GenInfoHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

  if (acquireCoreGenInfo){
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<TYPE>(#NAME, DEFVAL); this->addConsumed<TYPE>(#NAME);
    GENINFO_VARIABLES;
#undef GENINFO_VARIABLE
  }

  // ME reweighting branches are defined as sloppy
  std::vector<TString> allbranchnames; tree->getValidBranchNamesWithoutAlias(allbranchnames, false);
  std::vector<TString> melist;
  bool has_lheparticles=false;
  for (TString const& bname : allbranchnames){
    if (acquireLHEMEWeights && (bname.Contains("p_Gen") || bname.Contains("LHECandMass"))){
      tree->bookBranch<float>(bname, 0.f);
      this->addConsumed<float>(bname);
      this->defineConsumedSloppy(bname);
      melist.push_back(bname);
    }
    else if (acquireLHEParticles && bname.Contains(colName_lheparticles)) has_lheparticles = true;
  }
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
}
