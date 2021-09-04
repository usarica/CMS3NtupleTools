#include <cassert>

#include <CMS3/Dictionaries/interface/GlobalCollectionNames.h>

#include "ParticleObjectHelpers.h"
#include "MuonHandler.h"
#include "MuonSelectionHelpers.h"
#include <CMS3/Dictionaries/interface/CMS3StreamHelpers.h>


using namespace std;
using namespace IvyStreamHelpers;


#define MUON_MOMENTUM_VARIABLES \
MUON_VARIABLE(float, pt, 0) \
MUON_VARIABLE(float, eta, 0) \
MUON_VARIABLE(float, phi, 0) \
MUON_VARIABLE(float, mass, 0) \
MUON_VARIABLE(cms3_charge_t, charge, 0)


const std::string MuonHandler::colName = GlobalCollectionNames::colName_muons;

MuonHandler::MuonHandler() :
  IvyBase(),
  has_precomputed_timing(false),
  has_genmatching(false),
  hasOverlapMaps(false),
  overlapMap_muons_ak4jets(nullptr),
  overlapMap_muons_ak8jets(nullptr)
{
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(MuonHandler::colName + "_" + #NAME);
  MUON_MOMENTUM_VARIABLES;
  MUON_IDISO_VARIABLES;
  MUON_MOMENTUMSCALE_VARIABLES;
#undef MUON_VARIABLE
}


bool MuonHandler::constructMuons(SystematicsHelpers::SystematicVariationTypes const& syst, std::vector<PFCandidateObject*> const* pfcandidates){
  if (this->isAlreadyCached()) return true;

  clear();
  if (!currentTree) return false;

  bool res = (constructMuonObjects(syst) && associatePFCandidates(pfcandidates) && linkOverlapElements());

  if (res) this->cacheEvent();
  return res;
}

bool MuonHandler::associatePFCandidates(std::vector<PFCandidateObject*> const* pfcandidates) const{
  if (!pfcandidates) return true;

  for (auto const& pfcand:(*pfcandidates)){
    auto const& associated_particle_indices = pfcand->extras.matched_muon_index_list;
    for (auto const& part:productList){
      if (HelperFunctions::checkListVariable(associated_particle_indices, part->getUniqueIdentifier())){
        part->addDaughter(pfcand);
        pfcand->addMother(part);
      }
    }
  }

  return true;
}

bool MuonHandler::linkOverlapElements() const{
  if (!hasOverlapMaps) return true;

  overlapMap_muons_ak4jets->constructOverlapMaps();
  for (auto const& ome:overlapMap_muons_ak4jets->getProducts()) ome->linkFirstElement(productList);

  overlapMap_muons_ak8jets->constructOverlapMaps();
  for (auto const& ome:overlapMap_muons_ak8jets->getProducts()) ome->linkFirstElement(productList);

  return true;
}

bool MuonHandler::constructMuonObjects(SystematicsHelpers::SystematicVariationTypes const& syst){
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME;
  MUON_MOMENTUM_VARIABLES;
  MUON_VARIABLES;
#undef MUON_VARIABLE

    // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(MuonHandler::colName + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  MUON_MOMENTUM_VARIABLES;
  MUON_IDISO_VARIABLES;
  MUON_MOMENTUMSCALE_VARIABLES;
  if (this->has_precomputed_timing){
    MUON_PRETESTED_VARIABLES;
  }
  else{
    MUON_FULLTIMING_VARIABLES;
  }
  if (this->has_genmatching){
    MUON_GENINFO_VARIABLES;
  }
#undef MUON_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=MiscUtils::ERROR) IVYerr << "MuonHandler::constructMuonObjects: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=MiscUtils::DEBUG) IVYout << "MuonHandler::constructMuonObjects: All variables are set up!" << endl;

  if (itBegin_charge == itEnd_charge) return true; // Construction is successful, it is just that no muons exist.

  size_t nProducts = (itEnd_charge - itBegin_charge);
  productList.reserve(nProducts);
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  MUON_MOMENTUM_VARIABLES;
  MUON_VARIABLES;
#undef MUON_VARIABLE
  {
    size_t ip=0;
    while (it_charge != itEnd_charge){
      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "MuonHandler::constructMuonObjects: Attempting muon " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_pt, *it_eta, *it_phi, *it_mass); // Yes you have to do this on a separate line because CMSSW...
      productList.push_back(new MuonObject(-13*(*it_charge>0 ? 1 : -1), momentum));
      MuonObject*& obj = productList.back();

      // Set extras
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) obj->extras.NAME = *it_##NAME;
      MUON_IDISO_VARIABLES;
      MUON_MOMENTUMSCALE_VARIABLES;
      if (this->has_precomputed_timing){
        MUON_PRETESTED_VARIABLES;
      }
      else{
        MUON_FULLTIMING_VARIABLES;
      }
      if (this->has_genmatching){
        MUON_GENINFO_VARIABLES;
      }
#undef MUON_VARIABLE

      // Set particle index as its unique identifier
      obj->setUniqueIdentifier(ip);

      // Replace momentum
      obj->makeFinalMomentum(syst);

      // Set the selection bits
      MuonSelectionHelpers::setSelectionBits(*obj);

      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "\t- Success!" << endl;

      ip++;
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      MUON_MOMENTUM_VARIABLES;
      MUON_IDISO_VARIABLES;
      MUON_MOMENTUMSCALE_VARIABLES;
      if (this->has_precomputed_timing){
        MUON_PRETESTED_VARIABLES;
      }
      else{
        MUON_FULLTIMING_VARIABLES;
      }
      if (this->has_genmatching){
        MUON_GENINFO_VARIABLES;
      }
#undef MUON_VARIABLE
    }
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(productList);

  return true;
}

void MuonHandler::checkOptionalInfo(BaseTree* tree, bool& flag_precomputed_timing, bool& flag_genmatching){
  flag_precomputed_timing = flag_genmatching = true;

  std::vector<TString> bnames;
  tree->getValidBranchNamesWithoutAlias(bnames, false);

#define MUON_VARIABLE(TYPE, NAME, DEFVAL) flag_precomputed_timing &= (std::find(bnames.cbegin(), bnames.cend(), MuonHandler::colName + "_" + #NAME)!=bnames.cend());
  MUON_PRETESTED_VARIABLES;
#undef MUON_VARIABLE
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) flag_genmatching &= (std::find(bnames.cbegin(), bnames.cend(), MuonHandler::colName + "_" + #NAME)!=bnames.cend());
  MUON_GENINFO_VARIABLES;
#undef MUON_VARIABLE
}

bool MuonHandler::wrapTree(BaseTree* tree){
  if (!tree) return false;

  MuonHandler::checkOptionalInfo(tree, this->has_precomputed_timing, this->has_genmatching);

  return IvyBase::wrapTree(tree);
}


void MuonHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

  MuonHandler::checkOptionalInfo(tree, this->has_precomputed_timing, this->has_genmatching);
#define MUON_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(MuonHandler::colName + "_" + #NAME); this->defineConsumedSloppy(MuonHandler::colName + "_" + #NAME);
  if (this->has_precomputed_timing){
    MUON_PRETESTED_VARIABLES;
  }
  else{
    MUON_FULLTIMING_VARIABLES;
  }
  if (this->has_genmatching){
    MUON_GENINFO_VARIABLES;
  }
#undef MUON_VARIABLE

#define MUON_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(MuonHandler::colName + "_" + #NAME, nullptr);
  MUON_MOMENTUM_VARIABLES;
  MUON_IDISO_VARIABLES;
  MUON_MOMENTUMSCALE_VARIABLES;
  if (this->has_precomputed_timing){
    MUON_PRETESTED_VARIABLES;
  }
  else{
    MUON_FULLTIMING_VARIABLES;
  }
  if (this->has_genmatching){
    MUON_GENINFO_VARIABLES;
  }
#undef MUON_VARIABLE
}

void MuonHandler::registerOverlapMaps(
  OverlapMapHandler<MuonObject, AK4JetObject>& overlapMap_muons_ak4jets_,
  OverlapMapHandler<MuonObject, AK8JetObject>& overlapMap_muons_ak8jets_
){
  overlapMap_muons_ak4jets = &overlapMap_muons_ak4jets_;
  overlapMap_muons_ak8jets = &overlapMap_muons_ak8jets_;
  hasOverlapMaps = true;
}


#undef MUON_MOMENTUM_VARIABLES
