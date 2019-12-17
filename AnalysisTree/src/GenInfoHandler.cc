#include <cassert>
#include "GenInfoHandler.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


GenInfoHandler::GenInfoHandler() :
  IvyBase(),
  genInfo(nullptr)
{
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<TYPE>(#NAME);
  GENINFO_VARIABLES
#undef GENINFO_VARIABLE
}


bool GenInfoHandler::constructGenInfo(SystematicsHelpers::SystematicVariationTypes const& syst){
  clear();
  if (!currentTree) return false;

  bool doLHEParticles = tree_lheparticles_present_map[currentTree];

#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME = DEFVAL;
  GENINFO_VARIABLES;
#undef GENINFO_VARIABLE
#define GENINFO_VECTOR_VARIABLE(TYPE, NAME) TYPE* NAME = nullptr;
  GENINFO_VECTOR_VARIABLES;
#undef GENINFO_VECTOR_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedValue(#NAME, NAME);
  GENINFO_VARIABLES;
#undef GENINFO_VARIABLE

  std::unordered_map<TString, float> MElist;
  for (TString const& strme:tree_MElist_map[currentTree]){
    MElist[strme] = 0;
    allVariablesPresent &= this->getConsumedValue(strme, MElist.find(strme)->second);
  }

  if (doLHEParticles){
#define GENINFO_VECTOR_VARIABLE(TYPE, NAME) allVariablesPresent &= this->getConsumedValue(#NAME, NAME);
    GENINFO_VECTOR_VARIABLES;
#undef GENINFO_VECTOR_VARIABLE
  }

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "GenInfoHandler::constructGenInfos: Not all variables are consumed properly!" << endl;
    assert(0);
  }
  if (this->verbosity>=TVar::DEBUG) MELAout << "GenInfoHandler::constructGenInfo: All variables are set up!" << endl;

  genInfo = new GenInfoObject();
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) genInfo->extras.NAME = NAME;
  GENINFO_VARIABLES;
#undef GENINFO_VARIABLE

  for (auto it:MElist) genInfo->extras.LHE_ME_weights[it.first] = it.second;

  if (doLHEParticles){
#define GENINFO_VECTOR_VARIABLE(TYPE, NAME) genInfo->extras.NAME = *NAME;
    GENINFO_VECTOR_VARIABLES;
#undef GENINFO_VECTOR_VARIABLE
  }

  genInfo->setSystematic(syst);

  return true;
}

void GenInfoHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<TYPE>(#NAME, DEFVAL);
  GENINFO_VARIABLES
#undef GENINFO_VARIABLE

  // ME reweighting branches are defined as sloppy
  std::vector<TString> allbranchnames; tree->getValidBranchNamesWithoutAlias(allbranchnames, false);
  std::vector<TString> melist;
  bool has_lheparticles=false;
  for (TString const& bname : allbranchnames){
    if (bname.Contains("p_Gen")){
      tree->bookBranch<float>(bname, 0.f);
      this->addConsumed<float>(bname);
      this->defineConsumedSloppy(bname);
      melist.push_back(bname);
    }
    else if (bname.Contains("lheparticles")) has_lheparticles = true;
  }
  tree_MElist_map[tree] = melist;
  tree_lheparticles_present_map[tree] = has_lheparticles;

#define GENINFO_VECTOR_VARIABLE(TYPE, NAME) \
if (has_lheparticles){ tree->bookBranch<TYPE*>(#NAME, nullptr); } \
this->addConsumed<TYPE*>(#NAME); \
this->defineConsumedSloppy(#NAME);

  GENINFO_VECTOR_VARIABLES;
#undef GENINFO_VECTOR_VARIABLE
}
