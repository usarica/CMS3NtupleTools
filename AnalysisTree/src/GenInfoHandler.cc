#include <cassert>
#include "GenInfoHandler.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


GenInfoHandler::GenInfoHandler() : IvyBase()
{
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<TYPE>(#NAME);
  GENINFO_VARIABLES;
#undef GENINFO_VARIABLE
}


bool GenInfoHandler::constructGenInfo(SystematicsHelpers::SystematicVariationTypes const& syst){
  clear();
  if (!currentTree) return false;

#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) TYPE NAME = DEFVAL;
  GENINFO_VARIABLES;
#undef GENINFO_VARIABLE

    // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedValue<TYPE>(#NAME, NAME);
  GENINFO_VARIABLES;
#undef GENINFO_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "GenInfoHandler::constructGenInfos: Not all variables are consumed properly!" << endl;
    assert(0);
  }
  if (this->verbosity>=TVar::DEBUG) MELAout << "GenInfoHandler::constructGenInfo: All variables are set up!" << endl;

  genInfo = new GenInfoObject();
#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) genInfo->extras.NAME = NAME;
  GENINFO_VARIABLES;
#undef GENINFO_VARIABLE
  genInfo->setSystematic(syst);

  return true;
}

void GenInfoHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

#define GENINFO_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<TYPE>(#NAME, DEFVAL);
  GENINFO_VARIABLES
#undef GENINFO_VARIABLE
}
