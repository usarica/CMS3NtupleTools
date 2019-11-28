#include <cassert>
#include "ParticleObjectHelpers.h"
#include "MuonHandler.h"
//#include "MuonSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


#define VECTOR_ITERATOR_HANDLER_DIRECTIVES \
MUON_VARIABLE(float, pt, 0) \
MUON_VARIABLE(float, eta, 0) \
MUON_VARIABLE(float, phi, 0) \
MUON_VARIABLE(float, mass, 0) \
MUON_VARIABLE(int, charge, 0) \
MUON_VARIABLES


const std::string MuonHandler::colName = "muons";

MuonHandler::MuonHandler() : IvyBase()
{
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(MuonHandler::colName + "_" + #NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef MUON_VARIABLE
}


bool MuonHandler::constructMuons(SystematicsHelpers::SystematicVariationTypes const& syst){
  clear();
  if (!currentTree) return false;

#define MUON_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef MUON_VARIABLE

    // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(MuonHandler::colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef MUON_VARIABLE

  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "MuonHandler::constructMuons: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "MuonHandler::constructMuons: All variables are set up!" << endl;

  if (itBegin_charge == itEnd_charge) return true; // Construction is successful, it is just that no muons exist.

  size_t nProducts = (itEnd_charge - itBegin_charge);
  productList.reserve(nProducts);
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef MUON_VARIABLE
  {
    size_t ip=0;
    while (it_charge != itEnd_charge){
      if (this->verbosity>=TVar::DEBUG) MELAout << "MuonHandler::constructMuons: Attempting muon " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_pt, *it_eta, *it_phi, *it_mass); // Yes you have to do this on a separate line because CMSSW...
      productList.push_back(new MuonObject(-13*(*it_charge>0 ? 1 : -1), momentum));
      MuonObject*& obj = productList.back();

      // Set extras
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) obj->extras.NAME = *it_##NAME;
      MUON_VARIABLES;
#undef MUON_VARIABLE

      // Replace momentum
      obj->makeFinalMomentum(syst);

      // Set the selection bits
      //MuonSelectionHelpers::setSelectionBits(*obj);

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      VECTOR_ITERATOR_HANDLER_DIRECTIVES;
#undef MUON_VARIABLE
    }
  }
    // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(productList);

  return true;
}

void MuonHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

#define MUON_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(MuonHandler::colName + "_" + #NAME, nullptr);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES
#undef MUON_VARIABLE
}


#undef VECTOR_ITERATOR_HANDLER_DIRECTIVES
