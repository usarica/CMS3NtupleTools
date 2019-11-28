#include <cassert>
#include "ParticleObjectHelpers.h"
#include "JetMETHandler.h"
//#include "MuonSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


#define VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK4JETS \
AK4JET_VARIABLE(float, pt, 0) \
AK4JET_VARIABLE(float, eta, 0) \
AK4JET_VARIABLE(float, phi, 0) \
AK4JET_VARIABLE(float, mass, 0) \
AK4JET_VARIABLES
#define VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK8JETS \
AK8JET_VARIABLE(float, pt, 0) \
AK8JET_VARIABLE(float, eta, 0) \
AK8JET_VARIABLE(float, phi, 0) \
AK8JET_VARIABLE(float, mass, 0) \
AK8JET_VARIABLES

const std::string JetMETHandler::colName_ak4jets = "ak4jets";
const std::string JetMETHandler::colName_ak8jets = "ak8jets";


JetMETHandler::JetMETHandler() : IvyBase()
{
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(JetMETHandler::colName_ak4jets + "_" + #NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK4JETS;
#undef AK4JET_VARIABLE
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(JetMETHandler::colName_ak8jets + "_" + #NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK8JETS;
#undef AK8JET_VARIABLE
}

void JetMETHandler::clear(){
  for (auto*& prod:ak4jets) delete prod;
  ak4jets.clear();
  for (auto*& prod:ak8jets) delete prod;
  ak8jets.clear();
}

bool JetMETHandler::constructProducts(SystematicsHelpers::SystematicVariationTypes const& syst){
  clear();
  if (!currentTree) return false;

#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_ak4jets_##NAME, itEnd_ak4jets_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK4JETS;
#undef AK4JET_VARIABLE
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_ak8jets_##NAME, itEnd_ak8jets_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK8JETS;
#undef AK8JET_VARIABLE

    // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(JetMETHandler::colName_ak4jets + "_" + #NAME, &itBegin_ak4jets_##NAME, &itEnd_ak4jets_##NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK4JETS;
#undef AK4JET_VARIABLE
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(JetMETHandler::colName_ak8jets + "_" + #NAME, &itBegin_ak8jets_##NAME, &itEnd_ak8jets_##NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK8JETS;
#undef AK8JET_VARIABLE

  if (!allVariablesPresent && this->verbosity>=TVar::ERROR){
    MELAerr << "JetMETHandler::constructProducts: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructProducts: All variables are set up!" << endl;

  /************/
  /* ak4 jets */
  /************/
  size_t nak4jets = (itEnd_ak4jets_pt - itBegin_ak4jets_pt);
  ak4jets.reserve(nak4jets);
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) auto it_ak4jets_##NAME = itBegin_ak4jets_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK4JETS;
#undef AK4JET_VARIABLE
  {
    size_t ip=0;
    while (it_ak4jets_pt != itEnd_ak4jets_pt){
      if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructProducts: Attempting ak4 jet " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_ak4jets_pt, *it_ak4jets_eta, *it_ak4jets_phi, *it_ak4jets_mass); // Yes you have to do this on a separate line because CMSSW...
      ak4jets.push_back(new AK4JetObject(momentum));
      AK4JetObject*& obj = ak4jets.back();

      // Set extras
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) obj->extras.NAME = *it_ak4jets_##NAME;
      AK4JET_VARIABLES;
#undef AK4JET_VARIABLE

      // Replace momentum
      obj->makeFinalMomentum(syst);

      // Set the selection bits
      //AK4JetSelectionHelpers::setSelectionBits(*obj);

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) it_ak4jets_##NAME++;
      VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK4JETS;
#undef AK4JET_VARIABLE
    }
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(ak4jets);

  /************/
  /* ak8 jets */
  /************/
  size_t nak8jets = (itEnd_ak8jets_pt - itBegin_ak8jets_pt);
  ak8jets.reserve(nak8jets);
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) auto it_ak8jets_##NAME = itBegin_ak8jets_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK8JETS;
#undef AK8JET_VARIABLE
  {
    size_t ip=0;
    while (it_ak8jets_pt != itEnd_ak8jets_pt){
      if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructProducts: Attempting ak8 jet " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_ak8jets_pt, *it_ak8jets_eta, *it_ak8jets_phi, *it_ak8jets_mass); // Yes you have to do this on a separate line because CMSSW...
      ak8jets.push_back(new AK8JetObject(momentum));
      AK8JetObject*& obj = ak8jets.back();

      // Set extras
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) obj->extras.NAME = *it_ak8jets_##NAME;
      AK8JET_VARIABLES;
#undef AK8JET_VARIABLE

      // Replace momentum
      obj->makeFinalMomentum(syst);

      // Set the selection bits
      //AK8JetSelectionHelpers::setSelectionBits(*obj);

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) it_ak8jets_##NAME++;
      VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK8JETS;
#undef AK8JET_VARIABLE
    }
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(ak8jets);

  return true;
}

void JetMETHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(JetMETHandler::colName_ak4jets + "_" + #NAME, nullptr);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK4JETS
#undef AK4JET_VARIABLE
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(JetMETHandler::colName_ak8jets + "_" + #NAME, nullptr);
    VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK8JETS
#undef AK8JET_VARIABLE
}


#undef VECTOR_ITERATOR_HANDLER_DIRECTIVES
