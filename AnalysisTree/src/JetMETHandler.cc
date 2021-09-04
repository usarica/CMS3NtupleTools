#include <cassert>

#include <CMS3/Dictionaries/interface/GlobalCollectionNames.h>
#include <CMS3/Dictionaries/interface/JetMETEnums.h>

#include "RunLumiEventBlock.h"
#include "ParticleObjectHelpers.h"
#include "SamplesCore.h"
#include "FSRObject.h"
#include "JetMETHandler.h"
#include "MuonSelectionHelpers.h"
#include "ElectronSelectionHelpers.h"
#include "PhotonSelectionHelpers.h"
#include "AK4JetSelectionHelpers.h"
#include "AK8JetSelectionHelpers.h"
#include "ParticleSelectionHelpers.h"
#include <CMS3/Dictionaries/interface/CMS3StreamHelpers.h>


using namespace std;
using namespace IvyStreamHelpers;


#define VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK4JETS \
AK4JET_VARIABLE(float, pt, 0) \
AK4JET_VARIABLE(float, eta, 0) \
AK4JET_VARIABLE(float, phi, 0) \
AK4JET_VARIABLE(float, mass, 0) \
AK4JET_RECO_VARIABLES
#define VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK8JETS \
AK8JET_VARIABLE(float, pt, 0) \
AK8JET_VARIABLE(float, eta, 0) \
AK8JET_VARIABLE(float, phi, 0) \
AK8JET_VARIABLE(float, mass, 0) \
AK8JET_RECO_VARIABLES
#define JETMET_METXY_VERTEX_VARIABLES \
JETMET_METXY_VERTEX_VARIABLE(cms3_vtxs_nvtxs_t, nvtxs, 0)


const std::string JetMETHandler::colName_ak4jets = GlobalCollectionNames::colName_ak4jets;
const std::string JetMETHandler::colName_ak8jets = GlobalCollectionNames::colName_ak8jets;
const std::string JetMETHandler::colName_pfmet = GlobalCollectionNames::colName_pfmet;
const std::string JetMETHandler::colName_pfpuppimet = GlobalCollectionNames::colName_puppimet;
const std::string JetMETHandler::colName_vertices = GlobalCollectionNames::colName_vtxs;

JetMETHandler::JetMETHandler() :
  IvyBase(),

  flag_hasRevertMETFixVariables(false),
  hasOverlapMaps(false),
  overlapMap_muons_ak4jets(nullptr),
  overlapMap_muons_ak8jets(nullptr),
  overlapMap_electrons_ak4jets(nullptr),
  overlapMap_electrons_ak8jets(nullptr),
  overlapMap_photons_ak4jets(nullptr),
  overlapMap_photons_ak8jets(nullptr),

  pfmet(nullptr),
  pfpuppimet(nullptr),

  pfmet_XYcorr_xCoeffA(0),
  pfmet_XYcorr_xCoeffB(0),
  pfmet_XYcorr_yCoeffA(0),
  pfmet_XYcorr_yCoeffB(0)
{
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(JetMETHandler::colName_ak4jets + "_" + #NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK4JETS;
#undef AK4JET_VARIABLE
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(JetMETHandler::colName_ak8jets + "_" + #NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK8JETS;
#undef AK8JET_VARIABLE
#define MET_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<TYPE>(JetMETHandler::colName_pfmet + "_" + #NAME);
  MET_CORE_VARIABLES;
#undef MET_VARIABLE
#define MET_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<TYPE>(JetMETHandler::colName_pfpuppimet + "_" + #NAME);
  MET_RECORDED_CORE_VARIABLES;
#undef MET_VARIABLE

  // Vertex variables
#define JETMET_METXY_VERTEX_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<TYPE>(JetMETHandler::colName_vertices + "_" + #NAME);
  JETMET_METXY_VERTEX_VARIABLES;
#undef JETMET_METXY_VERTEX_VARIABLE
}

void JetMETHandler::clear(){
  this->resetCache();

  for (auto*& prod:ak4jets) delete prod;
  for (auto*& prod:ak4jets_masked) delete prod;
  ak4jets.clear();
  ak4jets_masked.clear();
  for (auto*& prod:ak8jets) delete prod;
  for (auto*& prod:ak8jets_masked) delete prod;
  ak8jets.clear();
  ak8jets_masked.clear();
  delete pfmet; pfmet=nullptr;
  delete pfpuppimet; pfpuppimet=nullptr;
}

JetMETHandler::METFixMode JetMETHandler::getMETFixMode(SimEventHandler const* simEventHandler) const{
  JetMETHandler::METFixMode res = kMETFix_NoPatch;

  int const& year = SampleHelpers::getDataYear();
  if (year!=2017 || !flag_hasRevertMETFixVariables) return res;

  TString effDataPeriod;
  if (!SampleHelpers::checkSampleIsData(currentTree->sampleIdentifier)){
    if (!simEventHandler){
      if (verbosity>=MiscUtils::ERROR) IVYerr << "JetMETHandler::getMETFixMode: MC checks require a SimEventHandler!" << endl;
      assert(0);
      return res;
    }
    effDataPeriod = simEventHandler->getChosenDataPeriod();
  }
  else{
    bool allVariablesPresent = true;
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) TYPE const* NAME = nullptr; allVariablesPresent &= this->getConsumed(#NAME, NAME);
    RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE
    if (!allVariablesPresent){
      if (this->verbosity>=MiscUtils::ERROR) IVYerr << "JetMETHandler::getMETFixMode: Not all variables of the data case are consumed properly!" << endl;
      assert(0);
    }
    effDataPeriod = SampleHelpers::getDataPeriodFromRunNumber(*RunNumber);
  }

  if (effDataPeriod == "2017E" || effDataPeriod == "2017F") res = kMETFix_RecoverGoodJets;
  else res = kMETFix_RevertMETFix;

  return res;
}

bool JetMETHandler::constructJetMET(
  SimEventHandler const* simEventHandler,
  SystematicsHelpers::SystematicVariationTypes const& syst,
  std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons, std::vector<PhotonObject*> const* photons,
  std::vector<PFCandidateObject*> const* pfcandidates
){
  if (this->isAlreadyCached()) return true;

  clear();
  if (!currentTree) return false;

  METFixMode const mode_metfix = this->getMETFixMode(simEventHandler);
  if (verbosity>=MiscUtils::DEBUG) IVYout << "JetMETHandler::constructJetMET: MET fix mode: " << mode_metfix << endl;

  bool res = (
    constructAK4Jets(syst, mode_metfix) && constructAK8Jets(syst) && constructMET(syst, mode_metfix)
    &&
    associatePFCandidates(pfcandidates) && linkOverlapElements()
    &&
    applyJetCleaning(hasOverlapMaps && (pfcandidates!=nullptr), mode_metfix, muons, electrons, photons)
    &&
    assignMETXYShifts(syst) && applyMETParticleShifts(hasOverlapMaps && (pfcandidates!=nullptr), muons, electrons, photons)
    );

  if (res) this->cacheEvent();
  return res;
}

bool JetMETHandler::constructAK4Jets(SystematicsHelpers::SystematicVariationTypes const& syst, JetMETHandler::METFixMode const& mode_metfix){
  bool const isData = SampleHelpers::checkSampleIsData(currentTree->sampleIdentifier);
  bool const doRevertMETFix = (mode_metfix == kMETFix_RevertMETFix);

#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK4JETS;
  AK4JET_GENINFO_VARIABLES;
  AK4JET_REVERTMETFIX_VARIABLES;
#undef AK4JET_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(JetMETHandler::colName_ak4jets + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK4JETS;
  if (!isData){
    AK4JET_GENINFO_VARIABLES;
  }
  if (doRevertMETFix){
    AK4JET_REVERTMETFIX_VARIABLES;
  }
#undef AK4JET_VARIABLE
  if (!allVariablesPresent){
    if (this->verbosity>=MiscUtils::ERROR) IVYerr << "JetMETHandler::constructAK4Jets: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=MiscUtils::DEBUG) IVYout << "JetMETHandler::constructAK4Jets: All variables are set up!" << endl;

  /************/
  /* ak4 jets */
  /************/
  size_t n_objects = (itEnd_pt - itBegin_pt);
  ak4jets.reserve(n_objects);
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK4JETS;
  AK4JET_GENINFO_VARIABLES;
  AK4JET_REVERTMETFIX_VARIABLES;
#undef AK4JET_VARIABLE
  {
    size_t ip=0;
    while (it_pt != itEnd_pt){
      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "JetMETHandler::constructAK4Jets: Attempting ak4 jet " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_pt, *it_eta, *it_phi, *it_mass); // Yes you have to do this on a separate line because CMSSW...
      ak4jets.push_back(new AK4JetObject(momentum));
      AK4JetObject*& obj = ak4jets.back();

      // Set extras
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) obj->extras.NAME = *it_##NAME;
      AK4JET_RECO_VARIABLES;
      if (!isData){
        AK4JET_GENINFO_VARIABLES;
      }
#undef AK4JET_VARIABLE

      // Override MET flags applying an OR with RevertMETFix versions before setting selection bits
      if (doRevertMETFix){
        if (verbosity>=MiscUtils::DEBUG){
          IVYout
            << "\t- isMETJERCSafe_Bits = " << HelperFunctions::displayBitsAsString(obj->extras.isMETJERCSafe_Bits, JetMETEnums::nMETShiftTypes)
            << " |= " << HelperFunctions::displayBitsAsString(*it_isMETJERCSafe_RevertMETFix_Bits, JetMETEnums::nMETShiftTypes)
            << " = " << HelperFunctions::displayBitsAsString(obj->extras.isMETJERCSafe_Bits | *it_isMETJERCSafe_RevertMETFix_Bits, JetMETEnums::nMETShiftTypes)
            << endl;
          IVYout
            << "\t- isMETJERCSafe_p4Preserved_Bits = " << HelperFunctions::displayBitsAsString(obj->extras.isMETJERCSafe_p4Preserved_Bits, JetMETEnums::nMETShiftTypes)
            << " |= " << HelperFunctions::displayBitsAsString(*it_isMETJERCSafe_p4Preserved_RevertMETFix_Bits, JetMETEnums::nMETShiftTypes)
            << " = " << HelperFunctions::displayBitsAsString(obj->extras.isMETJERCSafe_p4Preserved_Bits | *it_isMETJERCSafe_p4Preserved_RevertMETFix_Bits, JetMETEnums::nMETShiftTypes)
            << endl;
        }
        obj->extras.isMETJERCSafe_Bits = obj->extras.isMETJERCSafe_Bits | *it_isMETJERCSafe_RevertMETFix_Bits;
        obj->extras.isMETJERCSafe_p4Preserved_Bits = obj->extras.isMETJERCSafe_p4Preserved_Bits | *it_isMETJERCSafe_p4Preserved_RevertMETFix_Bits;
      }

      // Replace momentum
      obj->makeFinalMomentum(syst);

      // Set the selection bits
      AK4JetSelectionHelpers::setSelectionBits(*obj);

      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "\t- Success!" << endl;

      ip++;
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK4JETS;
      if (!isData){
        AK4JET_GENINFO_VARIABLES;
      }
      if (doRevertMETFix){
        AK4JET_REVERTMETFIX_VARIABLES;
      }
#undef AK4JET_VARIABLE
    }
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(ak4jets);

  return true;
}
bool JetMETHandler::constructAK8Jets(SystematicsHelpers::SystematicVariationTypes const& syst){
  bool const isData = SampleHelpers::checkSampleIsData(currentTree->sampleIdentifier);

#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_##NAME, itEnd_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK8JETS;
  AK8JET_GENINFO_VARIABLES;
#undef AK8JET_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(JetMETHandler::colName_ak8jets + "_" + #NAME, &itBegin_##NAME, &itEnd_##NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK8JETS;
  if (!isData){
    AK8JET_GENINFO_VARIABLES;
  }
#undef AK8JET_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=MiscUtils::ERROR) IVYerr << "JetMETHandler::constructAK8Jets: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=MiscUtils::DEBUG) IVYout << "JetMETHandler::constructAK8Jets: All variables are set up!" << endl;

  /************/
  /* ak8 jets */
  /************/
  size_t n_objects = (itEnd_pt - itBegin_pt);
  ak8jets.reserve(n_objects);
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) auto it_##NAME = itBegin_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK8JETS;
  AK8JET_GENINFO_VARIABLES;
#undef AK8JET_VARIABLE
  {
    size_t ip=0;
    while (it_pt != itEnd_pt){
      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "JetMETHandler::constructAK8Jets: Attempting ak8 jet " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_pt, *it_eta, *it_phi, *it_mass); // Yes you have to do this on a separate line because CMSSW...
      ak8jets.push_back(new AK8JetObject(momentum));
      AK8JetObject*& obj = ak8jets.back();

      // Set extras
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) obj->extras.NAME = *it_##NAME;
      AK8JET_RECO_VARIABLES;
      if (!isData){
        AK8JET_GENINFO_VARIABLES;
      }
#undef AK8JET_VARIABLE

      // Replace momentum
      obj->makeFinalMomentum(syst);

      // Set the selection bits
      AK8JetSelectionHelpers::setSelectionBits(*obj);

      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "\t- Success!" << endl;

      ip++;
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) it_##NAME++;
      VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK8JETS;
      if (!isData){
        AK8JET_GENINFO_VARIABLES;
      }
#undef AK8JET_VARIABLE
    }
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(ak8jets);

  return true;
}
bool JetMETHandler::associatePFCandidates(std::vector<PFCandidateObject*> const* pfcandidates) const{
  if (!pfcandidates) return true;

  for (auto const& part:(*pfcandidates)){
    auto const& associated_ak4jet_indices = part->extras.matched_ak4jet_index_list;
    for (auto const& jet:ak4jets){
      if (HelperFunctions::checkListVariable(associated_ak4jet_indices, jet->getUniqueIdentifier())){
        jet->addDaughter(part);
        part->addMother(jet);
      }
    }
    auto const& associated_ak8jet_indices = part->extras.matched_ak8jet_index_list;
    for (auto const& jet:ak8jets){
      if (HelperFunctions::checkListVariable(associated_ak8jet_indices, jet->getUniqueIdentifier())){
        jet->addDaughter(part);
        part->addMother(jet);
      }
    }
  }

  return true;
}

bool JetMETHandler::linkOverlapElements() const{
  if (!hasOverlapMaps) return true;

  overlapMap_muons_ak4jets->constructOverlapMaps();
  overlapMap_electrons_ak4jets->constructOverlapMaps();
  overlapMap_photons_ak4jets->constructOverlapMaps();

  overlapMap_muons_ak8jets->constructOverlapMaps();
  overlapMap_electrons_ak8jets->constructOverlapMaps();
  overlapMap_photons_ak8jets->constructOverlapMaps();

  for (auto const& ome:overlapMap_muons_ak4jets->getProducts()) ome->linkSecondElement(ak4jets);
  for (auto const& ome:overlapMap_electrons_ak4jets->getProducts()) ome->linkSecondElement(ak4jets);
  for (auto const& ome:overlapMap_photons_ak4jets->getProducts()) ome->linkSecondElement(ak4jets);

  for (auto const& ome:overlapMap_muons_ak8jets->getProducts()) ome->linkSecondElement(ak8jets);
  for (auto const& ome:overlapMap_electrons_ak8jets->getProducts()) ome->linkSecondElement(ak8jets);
  for (auto const& ome:overlapMap_photons_ak8jets->getProducts()) ome->linkSecondElement(ak8jets);

  return true;
}

bool JetMETHandler::applyJetCleaning(bool usePFCandidates, JetMETHandler::METFixMode const& mode_metfix, std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons, std::vector<PhotonObject*> const* photons){
  // Baseline correction to be computed before cleaning
  ParticleObject::LorentzVector_t pfmet_METCorrection_baseline[4];
  // Only compute recovery when there is a MET fix applied but we need to add back good jets.
  if (mode_metfix == kMETFix_RecoverGoodJets){
    if (this->verbosity>=MiscUtils::DEBUG) IVYout << "JetMETHandler::applyJetCleaning: Recovering good jets from the MET fix:" << endl;

    for (auto const& jet:ak4jets){
      if (!jet->testSelectionBit(AK4JetSelectionHelpers::kMETFixRecoverySuitable)) continue;
      // The uncorrected p4 is to be subtracted straightaway.
      ParticleObject::LorentzVector_t const p4_uncorrected = jet->uncorrected_p4();
      bool const isMETJERCSafe_Base = jet->testSelectionBit(AK4JetSelectionHelpers::kMETJERCSafe_Base);
      if (this->verbosity>=MiscUtils::DEBUG) IVYout
        << "=> Jet with corrected, uncorrected pT = " << jet->pt() << ", " << p4_uncorrected.Pt()
        << " is suitable. JERC safety = " << isMETJERCSafe_Base
        << ", JEC_L1L2L3 = " << jet->extras.JECNominal
        << ", JEC_L1 = " << jet->extras.JECL1Nominal
        << ", relJECUnc = " << jet->extras.relJECUnc
        << ", relJECUnc_nomus = " << jet->extras.relJECUnc_nomus
        << ", relJECUnc_nomus_JERNominal = " << jet->extras.relJECUnc_nomus_JERNominal
        << ", JERNominal = " << jet->extras.JERNominal
        << ", isMETJERCSafe_Bits = " << HelperFunctions::displayBitsAsString(jet->extras.isMETJERCSafe_Bits, JetMETEnums::nMETShiftTypes)
        << ", isMETJERCSafe_p4Preserved_Bits = " << HelperFunctions::displayBitsAsString(jet->extras.isMETJERCSafe_p4Preserved_Bits, JetMETEnums::nMETShiftTypes)
        << endl;
      for (unsigned char ijer=0; ijer<2; ijer++){
        for (unsigned char ip4=0; ip4<2; ip4++){
          ParticleObject::LorentzVector_t p4_corr = -p4_uncorrected;
          if (isMETJERCSafe_Base){
            // p4_diff already has the (-) sign
            ParticleObject::LorentzVector_t p4_diff;
            jet->getT1METShift(ip4, ijer, p4_diff, true);
            p4_corr += p4_diff;
          }
          pfmet_METCorrection_baseline[2*ip4 + ijer] += p4_corr;
        }
      }
    }
    if (this->verbosity>=MiscUtils::DEBUG){
      IVYout << "=> Adding the following baseline vectors to PF MET:" << endl;
      for (unsigned char ijer=0; ijer<2; ijer++){
        for (unsigned char ip4=0; ip4<2; ip4++){
          IVYout
            << "\t- " << (ijer==0 ? "[No JER]" : "[JER]") << (ip4==0 ? "[p4 default]" : "[p4-preserved]")
            << " (pT, phi) = (" << pfmet_METCorrection_baseline[2*ip4 + ijer].Pt() << ", " << pfmet_METCorrection_baseline[2*ip4 + ijer].Phi() << ")"
            << endl;
        }
      }
    }
  }

  std::vector<AK4JetObject*> ak4jets_new; ak4jets_new.reserve(ak4jets.size()); ak4jets_masked.reserve(ak4jets.size());
  std::vector<AK8JetObject*> ak8jets_new; ak8jets_new.reserve(ak8jets.size()); ak8jets_masked.reserve(ak8jets.size());

  std::vector<MuonObject*> muons_jetcleaning;
  std::vector<ElectronObject*> electrons_jetcleaning;
  std::vector<PhotonObject*> photons_jetcleaning;
  if (muons){
    for (auto const& part:(*muons)){
      if (!ParticleSelectionHelpers::isParticleForJetCleaning(part)) continue;
      muons_jetcleaning.push_back(part);
    }
  }
  if (electrons){
    for (auto const& part:(*electrons)){
      if (!ParticleSelectionHelpers::isParticleForJetCleaning(part)) continue;
      electrons_jetcleaning.push_back(part);
    }
  }
  if (photons){
    for (auto const& part:(*photons)){
      if (!ParticleSelectionHelpers::isParticleForJetCleaning(part)) continue;
      photons_jetcleaning.push_back(part);
    }
  }

  if (!usePFCandidates){
    constexpr bool undoT1METCorrFromCleaned = true; // Gives slightly better performance, also more consistent with particle momentum corrections
    // In this scenario, a simple delta-R cleaning is done.
    for (auto*& jet:ak4jets){
      bool doSkip=false;
      if (!doSkip){
        for (auto const& part:muons_jetcleaning){
          if (!ParticleSelectionHelpers::isParticleForJetCleaning(part)) continue;
          if (jet->deltaR(part)<jet->ConeRadiusConstant){ doSkip=true; break; }
        }
      }
      if (!doSkip){
        for (auto const& part:electrons_jetcleaning){
          if (!ParticleSelectionHelpers::isParticleForJetCleaning(part)) continue;
          if (jet->deltaR(part)<jet->ConeRadiusConstant){ doSkip=true; break; }
        }
      }
      if (!doSkip){
        for (auto const& part:photons_jetcleaning){
          if (!ParticleSelectionHelpers::isParticleForJetCleaning(part)) continue;
          if (jet->deltaR(part)<jet->ConeRadiusConstant){ doSkip=true; break; }
        }
      }

      if (!doSkip) ak4jets_new.push_back(jet);
      else ak4jets_masked.push_back(jet);
    }
    ak4jets = ak4jets_new;

    // Correct T1 PF MET for the cleaned jets
    if (undoT1METCorrFromCleaned){
      ParticleObject::LorentzVector_t sump4_METContribution[4];
      for (auto const& jet:ak4jets_masked){
        bool const overrideBits = (mode_metfix == kMETFix_RecoverGoodJets && jet->testSelectionBit(AK4JetSelectionHelpers::kMETFixRecoverySuitable));
        bool const isMETJERCSafe_Base = jet->testSelectionBit(AK4JetSelectionHelpers::kMETJERCSafe_Base);
        for (unsigned char ijer=0; ijer<2; ijer++){
          for (unsigned char ip4=0; ip4<2; ip4++){
            if (!overrideBits || isMETJERCSafe_Base){
              ParticleObject::LorentzVector_t p4_METContribution;
              jet->getT1METShift(ip4, ijer, p4_METContribution, overrideBits);
              sump4_METContribution[2*ip4 + ijer] += p4_METContribution;
            }
          }
        }
      }
      for (unsigned char ijer=0; ijer<2; ijer++){
        for (unsigned char ip4=0; ip4<2; ip4++){
          ParticleObject::LorentzVector_t p4_corr = -sump4_METContribution[2*ip4 + ijer] + pfmet_METCorrection_baseline[2*ip4 + ijer];
          if (this->verbosity>=MiscUtils::DEBUG) IVYout
            << "JetMETHandler::applyJetCleaning: Total MET overlap correction "
            << (ijer==0 ? "without" : "with") << " JER, "
            << (ip4==0 ? "without" : "with") << " p4 preservation = " << p4_corr
            << endl;
          pfmet->setJetOverlapCorrection(p4_corr, ijer, ip4);
          //// Set the same correction for PUPPI
          //pfpuppimet->setJetOverlapCorrection(p4_corr, ijer, ip4);
        }
      }
    }

    for (auto*& jet:ak8jets){
      bool doSkip=false;
      if (!doSkip){
        for (auto const& part:muons_jetcleaning){
          if (!ParticleSelectionHelpers::isParticleForJetCleaning(part)) continue;
          if (jet->deltaR(part)<jet->ConeRadiusConstant){ doSkip=true; break; }
        }
      }
      if (!doSkip){
        for (auto const& part:electrons_jetcleaning){
          if (!ParticleSelectionHelpers::isParticleForJetCleaning(part)) continue;
          if (jet->deltaR(part)<jet->ConeRadiusConstant){ doSkip=true; break; }
        }
      }
      if (!doSkip){
        for (auto const& part:photons_jetcleaning){
          if (!ParticleSelectionHelpers::isParticleForJetCleaning(part)) continue;
          if (jet->deltaR(part)<jet->ConeRadiusConstant){ doSkip=true; break; }
        }
      }
      if (!doSkip) ak8jets_new.push_back(jet);
      else ak8jets_masked.push_back(jet);
    }
    ak8jets = ak8jets_new;
  }
  else{
    // In this scenario, overlaps are checked explicitly.
    // No jets are skipped. They are modified instead.
    // Modifications are propagated to MET!
    constexpr bool applyConeVetoToStripping = true;
    ParticleObject::LorentzVector_t sump4_METContribution_old[4];
    ParticleObject::LorentzVector_t sump4_METContribution_new[4];

    // ak4 jets
    for (auto*& jet:ak4jets){
      std::vector<PFCandidateObject*> common_pfcands;
      ParticleObject::LorentzVector_t sump4_overlaps;
      ParticleObject::LorentzVector_t sump4_overlaps_mucands;
      bool hasCorrections = false;
      AK4JetObject* oldjet = nullptr;

      {
        for (auto const& part:muons_jetcleaning){
          if (!ParticleSelectionHelpers::isParticleForJetCleaning(part)) continue;
          if (applyConeVetoToStripping && jet->deltaR(part)>=jet->ConeRadiusConstant) continue;
          MuonObject* thePart = nullptr;
          FSRObject* theFSR = nullptr;
          if (part->hasFSR()){
            for (auto const& dau:part->getDaughters()){
              MuonObject* thePartDau = dynamic_cast<MuonObject*>(dau);
              FSRObject* theFSRDau = dynamic_cast<FSRObject*>(dau);
              if (!thePart && thePartDau) thePart = thePartDau;
              if (!theFSR && theFSRDau) theFSR = theFSRDau;
            }
          }
          else thePart = part;
          auto overlapElement = overlapMap_muons_ak4jets->getMatchingOverlapMap(thePart, jet);
          // If an overlap element is found, it means the particle overlaps with the jet.
          if (overlapElement){
            hasCorrections = true;
            if (!oldjet) oldjet = new AK4JetObject(*jet);
            oldjet->addDaughter(thePart);
            ParticleObject::LorentzVector_t p4_overlap = overlapElement->p4_common();
            ParticleObject::LorentzVector_t p4_overlap_mucands = overlapElement->p4_commonMuCands_goodMET();
            std::vector<IvyParticle*> daughters_part;
            thePart->getDeepDaughters(daughters_part, false);
            auto const& daughters_jet = jet->getDaughters();
            for (auto const& daughter_part:daughters_part){
              PFCandidateObject* pfcand = dynamic_cast<PFCandidateObject*>(daughter_part);
              if (pfcand){
                if (HelperFunctions::checkListVariable(daughters_jet, daughter_part)){
                  if (!HelperFunctions::checkListVariable(common_pfcands, pfcand)) common_pfcands.push_back(pfcand);
                  p4_overlap -= pfcand->p4();
                  if (MuonSelectionHelpers::testGoodMETPFMuon(*pfcand)) p4_overlap_mucands -= pfcand->p4();
                }
              }
            }
            sump4_overlaps += p4_overlap;
            sump4_overlaps_mucands += p4_overlap_mucands;
            if (this->verbosity>=MiscUtils::DEBUG) IVYout
              << "JetMETHandler::applyJetCleaning: ak4 jet " << jet->getUniqueIdentifier() << " has overlap with muon " << thePart->getUniqueIdentifier() << ":"
              << "\n\t- Jet p4 = " << jet->p4()
              << "\n\t- Particle p4 = " << thePart->p4()
              << "\n\t- Overlap p4 = " << p4_overlap
              << "\n\t- Overlap muons p4 = " << p4_overlap_mucands
              << endl;
          }
          if (theFSR){
            if (HelperFunctions::checkListVariable(theFSR->extras.fsrMatch_ak4jet_index_list, jet->getUniqueIdentifier())){
              hasCorrections = true;
              if (!oldjet) oldjet = new AK4JetObject(*jet);
              oldjet->addDaughter(theFSR);
              sump4_overlaps += theFSR->p4();
              if (this->verbosity>=MiscUtils::DEBUG) IVYout
                << "JetMETHandler::applyJetCleaning: ak4 jet " << jet->getUniqueIdentifier() << " has overlap with muon FSR " << theFSR->getUniqueIdentifier() << ":"
                << "\n\t- Jet p4 = " << jet->p4()
                << "\n\t- FSR p4 = " << theFSR->p4()
                << endl;
            }
          }
        }
      }
      {
        for (auto const& part:electrons_jetcleaning){
          if (!ParticleSelectionHelpers::isParticleForJetCleaning(part)) continue;
          if (applyConeVetoToStripping && jet->deltaR(part)>=jet->ConeRadiusConstant) continue;
          ElectronObject* thePart = nullptr;
          FSRObject* theFSR = nullptr;
          if (part->hasFSR()){
            for (auto const& dau:part->getDaughters()){
              ElectronObject* thePartDau = dynamic_cast<ElectronObject*>(dau);
              FSRObject* theFSRDau = dynamic_cast<FSRObject*>(dau);
              if (!thePart && thePartDau) thePart = thePartDau;
              if (!theFSR && theFSRDau) theFSR = theFSRDau;
            }
          }
          else thePart = part;
          auto overlapElement = overlapMap_electrons_ak4jets->getMatchingOverlapMap(thePart, jet);
          // If an overlap element is found, it means the particle overlaps with the jet.
          if (overlapElement){
            hasCorrections = true;
            if (!oldjet) oldjet = new AK4JetObject(*jet);
            oldjet->addDaughter(thePart);
            ParticleObject::LorentzVector_t p4_overlap = overlapElement->p4_common();
            ParticleObject::LorentzVector_t p4_overlap_mucands = overlapElement->p4_commonMuCands_goodMET();
            std::vector<IvyParticle*> daughters_part;
            thePart->getDeepDaughters(daughters_part, false);
            auto const& daughters_jet = jet->getDaughters();
            for (auto const& daughter_part:daughters_part){
              PFCandidateObject* pfcand = dynamic_cast<PFCandidateObject*>(daughter_part);
              if (pfcand){
                if (HelperFunctions::checkListVariable(daughters_jet, daughter_part)){
                  if (!HelperFunctions::checkListVariable(common_pfcands, pfcand)) common_pfcands.push_back(pfcand);
                  p4_overlap -= pfcand->p4();
                  if (MuonSelectionHelpers::testGoodMETPFMuon(*pfcand)) p4_overlap_mucands -= pfcand->p4();
                }
              }
            }
            sump4_overlaps += p4_overlap;
            sump4_overlaps_mucands += p4_overlap_mucands;
            if (this->verbosity>=MiscUtils::DEBUG) IVYout
              << "JetMETHandler::applyJetCleaning: ak4 jet " << jet->getUniqueIdentifier() << " has overlap with electron " << thePart->getUniqueIdentifier() << ":"
              << "\n\t- Jet p4 = " << jet->p4()
              << "\n\t- Particle p4 = " << thePart->p4()
              << "\n\t- Overlap p4 = " << p4_overlap
              << "\n\t- Overlap muons p4 = " << p4_overlap_mucands
              << endl;
          }
          if (theFSR){
            if (HelperFunctions::checkListVariable(theFSR->extras.fsrMatch_ak4jet_index_list, jet->getUniqueIdentifier())){
              hasCorrections = true;
              if (!oldjet) oldjet = new AK4JetObject(*jet);
              oldjet->addDaughter(theFSR);
              sump4_overlaps += theFSR->p4();
              if (this->verbosity>=MiscUtils::DEBUG) IVYout
                << "JetMETHandler::applyJetCleaning: ak4 jet " << jet->getUniqueIdentifier() << " has overlap with electron FSR " << theFSR->getUniqueIdentifier() << ":"
                << "\n\t- Jet p4 = " << jet->p4()
                << "\n\t- FSR p4 = " << theFSR->p4()
                << endl;
            }
          }
        }
      }
      {
        for (auto const& part:photons_jetcleaning){
          if (!ParticleSelectionHelpers::isParticleForJetCleaning(part)) continue;
          if (applyConeVetoToStripping && jet->deltaR(part)>=jet->ConeRadiusConstant) continue;
          auto overlapElement = overlapMap_photons_ak4jets->getMatchingOverlapMap(part, jet);
          // If an overlap element is found, it means the particle overlaps with the jet.
          if (overlapElement){
            hasCorrections = true;
            if (!oldjet) oldjet = new AK4JetObject(*jet);
            oldjet->addDaughter(part);
            ParticleObject::LorentzVector_t p4_overlap = overlapElement->p4_common();
            ParticleObject::LorentzVector_t p4_overlap_mucands = overlapElement->p4_commonMuCands_goodMET();
            std::vector<IvyParticle*> daughters_part;
            part->getDeepDaughters(daughters_part, false);
            auto const& daughters_jet = jet->getDaughters();
            for (auto const& daughter_part:daughters_part){
              PFCandidateObject* pfcand = dynamic_cast<PFCandidateObject*>(daughter_part);
              if (pfcand){
                if (HelperFunctions::checkListVariable(daughters_jet, daughter_part)){
                  if (!HelperFunctions::checkListVariable(common_pfcands, pfcand)) common_pfcands.push_back(pfcand);
                  p4_overlap -= pfcand->p4();
                  if (MuonSelectionHelpers::testGoodMETPFMuon(*pfcand)) p4_overlap_mucands -= pfcand->p4();
                }
              }
            }
            sump4_overlaps += p4_overlap;
            sump4_overlaps_mucands += p4_overlap_mucands;
            if (this->verbosity>=MiscUtils::DEBUG) IVYout
              << "JetMETHandler::applyJetCleaning: ak4 jet " << jet->getUniqueIdentifier() << " has overlap with photon " << part->getUniqueIdentifier() << ":"
              << "\n\t- Jet p4 = " << jet->p4()
              << "\n\t- Particle p4 = " << part->p4()
              << "\n\t- Overlap p4 = " << p4_overlap
              << "\n\t- Overlap muons p4 = " << p4_overlap_mucands
              << endl;
          }
        }
      }

      if (hasCorrections){
        bool const overrideBits = (mode_metfix == kMETFix_RecoverGoodJets && jet->testSelectionBit(AK4JetSelectionHelpers::kMETFixRecoverySuitable));
        bool const isMETJERCSafe_Base = jet->testSelectionBit(AK4JetSelectionHelpers::kMETJERCSafe_Base);

        ParticleObject::LorentzVector_t p4_METContribution_old[4];
        ParticleObject::LorentzVector_t p4_METContribution_new[4];
        for (unsigned char ijer=0; ijer<2; ijer++){
          for (unsigned char ip4=0; ip4<2; ip4++){
            if (!overrideBits || isMETJERCSafe_Base) jet->getT1METShift(ip4, ijer, p4_METContribution_old[2*ip4 + ijer], overrideBits);
          }
        }

        ParticleObject::LorentzVector_t p4_jet_uncorrected_old = jet->uncorrected_p4();
        ParticleObject::LorentzVector_t p4_jet_uncorrected_new = p4_jet_uncorrected_old - sump4_overlaps;
        ParticleObject::LorentzVector_t p4_jet_uncorrected_mucands_old = jet->p4_mucands();
        ParticleObject::LorentzVector_t p4_jet_uncorrected_mucands_new = p4_jet_uncorrected_mucands_old - sump4_overlaps_mucands;
        for (auto const& pfcand:common_pfcands){
          p4_jet_uncorrected_new -= pfcand->uncorrected_p4();
          if (MuonSelectionHelpers::testGoodMETPFMuon(*pfcand)) p4_jet_uncorrected_mucands_new -= pfcand->uncorrected_p4();
        }
        if (this->verbosity>=MiscUtils::DEBUG){
          IVYout << "\t- Old uncorrected p4 = " << p4_jet_uncorrected_old << endl;
          IVYout << "\t- Old uncorrected mu candidates p4 = " << p4_jet_uncorrected_mucands_old << endl;
          IVYout << "\t- New uncorrected p4 = " << p4_jet_uncorrected_new << endl;
          IVYout << "\t- New uncorrected mu candidates p4 = " << p4_jet_uncorrected_mucands_new << endl;
        }

        assert(oldjet!=nullptr);
        oldjet->addDaughter(jet); jet->addMother(oldjet);

        jet->reset_uncorrected_p4(p4_jet_uncorrected_new);
        jet->reset_p4_mucands(p4_jet_uncorrected_mucands_new);
        jet->makeFinalMomentum(jet->getCurrentSyst());
        AK4JetSelectionHelpers::setSelectionBits(*jet, false, true);

        for (unsigned char ijer=0; ijer<2; ijer++){
          for (unsigned char ip4=0; ip4<2; ip4++){
            if (!overrideBits || isMETJERCSafe_Base) jet->getT1METShift(ip4, ijer, p4_METContribution_new[2*ip4 + ijer], overrideBits);
          }
        }

        // Add the old and new MET contributions
        for (unsigned char ijer=0; ijer<2; ijer++){
          for (unsigned char ip4=0; ip4<2; ip4++){
            sump4_METContribution_old[2*ip4 + ijer] += p4_METContribution_old[2*ip4 + ijer];
            sump4_METContribution_new[2*ip4 + ijer] += p4_METContribution_new[2*ip4 + ijer];

            if (this->verbosity>=MiscUtils::DEBUG){
              IVYout << "\t- Old MET contribution " << (ijer==0 ? "without" : "with") << " JER, " << (ip4==0 ? "without" : "with") << " p4 preservation = " << p4_METContribution_old[2*ip4 + ijer] << endl;
              IVYout << "\t- New MET contribution " << (ijer==0 ? "without" : "with") << " JER, " << (ip4==0 ? "without" : "with") << " p4 preservation = " << p4_METContribution_new[2*ip4 + ijer] << endl;
            }
          }
        }
      }

      // Never skip jets in this mode of operation
      ak4jets_new.push_back(jet);
      if (oldjet) ak4jets_masked.push_back(oldjet);
    }
    ak4jets = ak4jets_new;
    ParticleObjectHelpers::sortByGreaterPt(ak4jets);

    // Propagate MET corrections from overlap removal
    for (unsigned char ijer=0; ijer<2; ijer++){
      for (unsigned char ip4=0; ip4<2; ip4++){
        ParticleObject::LorentzVector_t p4_corr = sump4_METContribution_new[2*ip4 + ijer] - sump4_METContribution_old[2*ip4 + ijer] + pfmet_METCorrection_baseline[2*ip4 + ijer];
        if (this->verbosity>=MiscUtils::DEBUG) IVYout
          << "Total MET overlap correction "
          << (ijer==0 ? "without" : "with") << " JER, "
          << (ip4==0 ? "without" : "with") << " p4 preservation = " << p4_corr
          << endl;
        pfmet->setJetOverlapCorrection(p4_corr, ijer, ip4);
        //// Set the same correction for PUPPI
        //pfpuppimet->setJetOverlapCorrection(p4_corr, ijer, ip4);
      }
    }

    // ak8 jets
    for (auto*& jet:ak8jets){
      std::vector<PFCandidateObject*> common_pfcands;
      ParticleObject::LorentzVector_t sump4_overlaps;
      bool hasCorrections = false;
      AK8JetObject* oldjet = nullptr;

      {
        for (auto const& part:muons_jetcleaning){
          if (!ParticleSelectionHelpers::isParticleForJetCleaning(part)) continue;
          if (applyConeVetoToStripping && jet->deltaR(part)>=jet->ConeRadiusConstant) continue;
          MuonObject* thePart = nullptr;
          FSRObject* theFSR = nullptr;
          if (part->hasFSR()){
            for (auto const& dau:part->getDaughters()){
              MuonObject* thePartDau = dynamic_cast<MuonObject*>(dau);
              FSRObject* theFSRDau = dynamic_cast<FSRObject*>(dau);
              if (!thePart && thePartDau) thePart = thePartDau;
              if (!theFSR && theFSRDau) theFSR = theFSRDau;
            }
          }
          else thePart = part;
          auto overlapElement = overlapMap_muons_ak8jets->getMatchingOverlapMap(thePart, jet);
          // If an overlap element is found, it means the particle overlaps with the jet.
          if (overlapElement){
            hasCorrections = true;
            if (!oldjet) oldjet = new AK8JetObject(*jet);
            oldjet->addDaughter(thePart);
            ParticleObject::LorentzVector_t p4_overlap = overlapElement->p4_common();
            std::vector<IvyParticle*> daughters_part;
            thePart->getDeepDaughters(daughters_part, false);
            auto const& daughters_jet = jet->getDaughters();
            for (auto const& daughter_part:daughters_part){
              PFCandidateObject* pfcand = dynamic_cast<PFCandidateObject*>(daughter_part);
              if (pfcand){
                if (HelperFunctions::checkListVariable(daughters_jet, daughter_part)){
                  if (!HelperFunctions::checkListVariable(common_pfcands, pfcand)) common_pfcands.push_back(pfcand);
                  p4_overlap -= pfcand->p4();
                }
              }
            }
            sump4_overlaps += p4_overlap;
          }
          if (theFSR){
            if (HelperFunctions::checkListVariable(theFSR->extras.fsrMatch_ak8jet_index_list, jet->getUniqueIdentifier())){
              hasCorrections = true;
              if (!oldjet) oldjet = new AK8JetObject(*jet);
              oldjet->addDaughter(theFSR);
              sump4_overlaps += theFSR->p4();
            }
          }
        }
      }
      {
        for (auto const& part:electrons_jetcleaning){
          if (!ParticleSelectionHelpers::isParticleForJetCleaning(part)) continue;
          if (applyConeVetoToStripping && jet->deltaR(part)>=jet->ConeRadiusConstant) continue;
          ElectronObject* thePart = nullptr;
          FSRObject* theFSR = nullptr;
          if (part->hasFSR()){
            for (auto const& dau:part->getDaughters()){
              ElectronObject* thePartDau = dynamic_cast<ElectronObject*>(dau);
              FSRObject* theFSRDau = dynamic_cast<FSRObject*>(dau);
              if (!thePart && thePartDau) thePart = thePartDau;
              if (!theFSR && theFSRDau) theFSR = theFSRDau;
            }
          }
          else thePart = part;
          auto overlapElement = overlapMap_electrons_ak8jets->getMatchingOverlapMap(thePart, jet);
          // If an overlap element is found, it means the particle overlaps with the jet.
          if (overlapElement){
            hasCorrections = true;
            if (!oldjet) oldjet = new AK8JetObject(*jet);
            oldjet->addDaughter(thePart);
            ParticleObject::LorentzVector_t p4_overlap = overlapElement->p4_common();
            std::vector<IvyParticle*> daughters_part;
            thePart->getDeepDaughters(daughters_part, false);
            auto const& daughters_jet = jet->getDaughters();
            for (auto const& daughter_part:daughters_part){
              PFCandidateObject* pfcand = dynamic_cast<PFCandidateObject*>(daughter_part);
              if (pfcand){
                if (HelperFunctions::checkListVariable(daughters_jet, daughter_part)){
                  if (!HelperFunctions::checkListVariable(common_pfcands, pfcand)) common_pfcands.push_back(pfcand);
                  p4_overlap -= pfcand->p4();
                }
              }
            }
            sump4_overlaps += p4_overlap;
          }
          if (theFSR){
            if (HelperFunctions::checkListVariable(theFSR->extras.fsrMatch_ak8jet_index_list, jet->getUniqueIdentifier())){
              hasCorrections = true;
              if (!oldjet) oldjet = new AK8JetObject(*jet);
              oldjet->addDaughter(theFSR);
              sump4_overlaps += theFSR->p4();
            }
          }
        }
      }
      {
        for (auto const& part:photons_jetcleaning){
          if (!ParticleSelectionHelpers::isParticleForJetCleaning(part)) continue;
          if (applyConeVetoToStripping && jet->deltaR(part)>=jet->ConeRadiusConstant) continue;
          auto overlapElement = overlapMap_photons_ak8jets->getMatchingOverlapMap(part, jet);
          // If an overlap element is found, it means the particle overlaps with the jet.
          if (overlapElement){
            hasCorrections = true;
            if (!oldjet) oldjet = new AK8JetObject(*jet);
            oldjet->addDaughter(part);
            ParticleObject::LorentzVector_t p4_overlap = overlapElement->p4_common();
            std::vector<IvyParticle*> daughters_part;
            part->getDeepDaughters(daughters_part, false);
            auto const& daughters_jet = jet->getDaughters();
            for (auto const& daughter_part:daughters_part){
              PFCandidateObject* pfcand = dynamic_cast<PFCandidateObject*>(daughter_part);
              if (pfcand){
                if (HelperFunctions::checkListVariable(daughters_jet, daughter_part)){
                  if (!HelperFunctions::checkListVariable(common_pfcands, pfcand)) common_pfcands.push_back(pfcand);
                  p4_overlap -= pfcand->p4();
                }
              }
            }
            sump4_overlaps += p4_overlap;
          }
        }
      }

      if (hasCorrections){
        ParticleObject::LorentzVector_t p4_jet_uncorrected_old = jet->uncorrected_p4();
        ParticleObject::LorentzVector_t p4_jet_uncorrected_new = p4_jet_uncorrected_old - sump4_overlaps;
        for (auto const& pfcand:common_pfcands) p4_jet_uncorrected_new -= pfcand->uncorrected_p4();

        assert(oldjet!=nullptr);
        oldjet->addDaughter(jet); jet->addMother(oldjet);

        jet->reset_uncorrected_p4(p4_jet_uncorrected_new);
        jet->makeFinalMomentum(jet->getCurrentSyst());
        AK8JetSelectionHelpers::setSelectionBits(*jet, false, true);
      }

      // Never skip jets in this mode of operation
      ak8jets_new.push_back(jet);
      if (oldjet) ak8jets_masked.push_back(oldjet);
    }
    ak8jets = ak8jets_new;
    ParticleObjectHelpers::sortByGreaterPt(ak8jets);
  }

  return true;
}

bool JetMETHandler::constructMET(SystematicsHelpers::SystematicVariationTypes const& syst, JetMETHandler::METFixMode const& mode_metfix){
  bool const isData = SampleHelpers::checkSampleIsData(currentTree->sampleIdentifier);
  bool const doRevertMETFix = (mode_metfix == kMETFix_RevertMETFix);

#define MET_VARIABLE(TYPE, NAME, DEFVAL) TYPE const* pfmet_##NAME = nullptr;
  MET_VARIABLES;
  MET_RECORDED_REVERTMETFIX_VARIABLES;
  MET_SHIFT_RECORDED_REVERTMETFIX_CORE_VARIABLES;
  MET_SHIFT_RECORDED_REVERTMETFIX_GENINFO_VARIABLES;
#undef MET_VARIABLE
#define MET_VARIABLE(TYPE, NAME, DEFVAL) TYPE const* pfpuppimet_##NAME = nullptr;
  MET_RECORDED_VARIABLES;
#undef MET_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define MET_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumed(JetMETHandler::colName_pfmet + "_" + #NAME, pfmet_##NAME);
  MET_CORE_VARIABLES;
  if (!isData){
    MET_GENINFO_VARIABLES;
  }
  if (doRevertMETFix){
    MET_RECORDED_REVERTMETFIX_VARIABLES;
    MET_SHIFT_RECORDED_REVERTMETFIX_CORE_VARIABLES;
    if (!isData){
      MET_SHIFT_RECORDED_REVERTMETFIX_GENINFO_VARIABLES;
    }
  }
#undef MET_VARIABLE
#define MET_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumed(JetMETHandler::colName_pfpuppimet + "_" + #NAME, pfpuppimet_##NAME);
  MET_RECORDED_CORE_VARIABLES;
  if (!isData){
    MET_RECORDED_GENINFO_VARIABLES;
  }
#undef MET_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=MiscUtils::ERROR) IVYerr << "JetMETHandler::constructMET: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=MiscUtils::DEBUG) IVYout << "JetMETHandler::constructMET: All variables are set up!" << endl;

  /**********/
  /* PF MET */
  /**********/
  pfmet = new METObject();
#define MET_VARIABLE(TYPE, NAME, DEFVAL) pfmet->extras.NAME = *pfmet_##NAME;
  MET_CORE_VARIABLES;
  if (!isData){
    MET_GENINFO_VARIABLES;
  }
#undef MET_VARIABLE

  // Apply patch to revert MET fixes
  if (doRevertMETFix){
    ParticleObject::LorentzVector_t p4_Nominal = METObject::constructP4FromPtPhi(*pfmet_met_Nominal, *pfmet_metPhi_Nominal);
    ParticleObject::LorentzVector_t p4_Raw = METObject::constructP4FromPtPhi(*pfmet_met_Raw, *pfmet_metPhi_Raw);
    ParticleObject::LorentzVector_t p4_Raw_Default = METObject::constructP4FromPtPhi(*pfmet_met_Raw_Default, *pfmet_metPhi_Raw_Default);
    ParticleObject::LorentzVector_t p4_extrashift_Nominal = METObject::constructP4FromPxPy(*pfmet_metShift_RevertMETFix_px_JECNominal, *pfmet_metShift_RevertMETFix_py_JECNominal);
    // There is no need to subtract anything from p4_extrashift_Nominal because it is already a product of an if-else construct (else part).

    // New nominal should have the difference in raw MET and the additional shift due to accepted jets added.
    ParticleObject::LorentzVector_t p4_Nominal_Reverted = p4_Nominal + p4_Raw_Default - p4_Raw + p4_extrashift_Nominal;

    if (verbosity>=MiscUtils::DEBUG){
      IVYout << "JetMETHandler::constructMET: Applying a reversion to MET fix:" << endl;
      IVYout
        << "\t- PF MET raw MET-fixed -> MET default (pT, phi) = (" << *pfmet_met_Raw << ", " << *pfmet_metPhi_Raw
        << ") -> (" << *pfmet_met_Raw_Default << ", " << *pfmet_metPhi_Raw_Default << ")"
        << endl;
      IVYout
        << "\t- PF MET extra Nominal JERC shift to revert MET fix (pT, phi) = (" << p4_extrashift_Nominal.Pt() << ", " << p4_extrashift_Nominal.Phi() << ")"
        << endl;
      IVYout
        << "\t- PF MET nominal MET-fixed -> MET default (pT, phi) = (" << *pfmet_met_Nominal << ", " << *pfmet_metPhi_Nominal
        << ") -> (" << p4_Nominal_Reverted.Pt() << ", " << p4_Nominal_Reverted.Phi() << ")"
        << endl;
    }

    // Begin to override variables
    pfmet->extras.met_Nominal = p4_Nominal_Reverted.Pt(); pfmet->extras.metPhi_Nominal = p4_Nominal_Reverted.Phi();
    pfmet->extras.met_Raw = *pfmet_met_Raw_Default; pfmet->extras.metPhi_Raw = *pfmet_metPhi_Raw_Default;

#define REVERTMETFIX_COMMON_COMMANDS \
REVERTMETFIX_COMMAND(JECNominal)
#define REVERTMETFIX_GENINFO_COMMANDS \
REVERTMETFIX_COMMAND(JECDn) \
REVERTMETFIX_COMMAND(JECUp) \
REVERTMETFIX_COMMAND(JECNominal_JERNominal) \
REVERTMETFIX_COMMAND(JECNominal_JERDn) \
REVERTMETFIX_COMMAND(JECNominal_JERUp) \
REVERTMETFIX_COMMAND(JECDn_JERNominal) \
REVERTMETFIX_COMMAND(JECUp_JERNominal)

#define REVERTMETFIX_COMMAND(SYST) \
IVYout << "\t- metShift p4 default before correction [" << #SYST << "] = (" << pfmet->extras.metShift_px_##SYST << ", " << pfmet->extras.metShift_py_##SYST << ")" << endl; \
IVYout << "\t- metShift p4-preserved before correction [" << #SYST << "] = (" << pfmet->extras.metShift_p4Preserved_px_##SYST << ", " << pfmet->extras.metShift_p4Preserved_py_##SYST << ")" << endl;
    if (verbosity>=MiscUtils::DEBUG){
      REVERTMETFIX_COMMON_COMMANDS;
      if (!isData){
        REVERTMETFIX_GENINFO_COMMANDS;
      }
    }
#undef REVERTMETFIX_COMMAND

    // Special treatment for the JECNominal case:
    // extras::metShift_p*_JECNominal are used for bookkeeping. They are not actually added because extras::met*_Nominal already contain them.
    pfmet->extras.metShift_px_JECNominal += *pfmet_metShift_RevertMETFix_px_JECNominal;
    pfmet->extras.metShift_py_JECNominal += *pfmet_metShift_RevertMETFix_py_JECNominal;
    // All other variables need to subtract the MET reversion shift already applied on extras::met*_Nominal.
    pfmet->extras.metShift_p4Preserved_px_JECNominal += *pfmet_metShift_p4Preserved_RevertMETFix_px_JECNominal - *pfmet_metShift_RevertMETFix_px_JECNominal;
    pfmet->extras.metShift_p4Preserved_py_JECNominal += *pfmet_metShift_p4Preserved_RevertMETFix_py_JECNominal - *pfmet_metShift_RevertMETFix_py_JECNominal;
    // COMMON metShift variables (JECNominal) are corrected above, so correct only GENINFO versions as above.
#define REVERTMETFIX_COMMAND(SYST) \
pfmet->extras.metShift_px_##SYST += *pfmet_metShift_RevertMETFix_px_##SYST - *pfmet_metShift_RevertMETFix_px_JECNominal; \
pfmet->extras.metShift_py_##SYST += *pfmet_metShift_RevertMETFix_py_##SYST - *pfmet_metShift_RevertMETFix_py_JECNominal; \
pfmet->extras.metShift_p4Preserved_px_##SYST += *pfmet_metShift_p4Preserved_RevertMETFix_px_##SYST - *pfmet_metShift_RevertMETFix_px_##SYST; \
pfmet->extras.metShift_p4Preserved_py_##SYST += *pfmet_metShift_p4Preserved_RevertMETFix_py_##SYST - *pfmet_metShift_RevertMETFix_py_##SYST;
    if (!isData){
      REVERTMETFIX_GENINFO_COMMANDS;
    }
#undef REVERTMETFIX_COMMAND

#define REVERTMETFIX_COMMAND(SYST) \
IVYout << "\t- metShift p4 default after correction [" << #SYST << "] = (" << pfmet->extras.metShift_px_##SYST << ", " << pfmet->extras.metShift_py_##SYST << ")" << endl; \
IVYout << "\t- metShift p4-preserved after correction [" << #SYST << "] = (" << pfmet->extras.metShift_p4Preserved_px_##SYST << ", " << pfmet->extras.metShift_p4Preserved_py_##SYST << ")" << endl;
    if (verbosity>=MiscUtils::DEBUG){
      REVERTMETFIX_COMMON_COMMANDS;
      if (!isData){
        REVERTMETFIX_GENINFO_COMMANDS;
      }
    }
#undef REVERTMETFIX_COMMAND

#undef REVERTMETFIX_GENINFO_COMMANDS
#undef REVERTMETFIX_COMMON_COMMANDS
  }

  // Apply fix to JECNominal_JERNominal variables
  if (isData){
    pfmet->extras.metShift_p4Preserved_px_JECNominal_JERNominal = pfmet->extras.metShift_p4Preserved_px_JECNominal;
    pfmet->extras.metShift_p4Preserved_py_JECNominal_JERNominal = pfmet->extras.metShift_p4Preserved_py_JECNominal;
  }

  pfmet->setSystematic(syst);

  /*************/
  /* PUPPI MET */
  /*************/
  pfpuppimet = new METObject();
#define MET_VARIABLE(TYPE, NAME, DEFVAL) pfpuppimet->extras.NAME = *pfpuppimet_##NAME;
  MET_RECORDED_CORE_VARIABLES;
  if (!isData){
    MET_RECORDED_GENINFO_VARIABLES;
  }
#undef MET_VARIABLE

  pfpuppimet->setSystematic(syst);

  return true;
}
bool JetMETHandler::assignMETXYShifts(SystematicsHelpers::SystematicVariationTypes const& syst){
  // 200314: Only PF MET has XY shifts
#define JETMET_METXY_VERTEX_VARIABLE(TYPE, NAME, DEFVAL) TYPE const* NAME = nullptr;
  JETMET_METXY_VERTEX_VARIABLES;
#undef JETMET_METXY_VERTEX_VARIABLE

  bool allVariablesPresent = true;
#define JETMET_METXY_VERTEX_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumed(JetMETHandler::colName_vertices + "_" + #NAME, NAME);
  JETMET_METXY_VERTEX_VARIABLES;
#undef JETMET_METXY_VERTEX_VARIABLE
  if (!allVariablesPresent){
    if (this->verbosity>=MiscUtils::ERROR) IVYerr << "JetMETHandler::assignMETXYShifts: Not all variables are consumed properly!" << endl;
    assert(0);
    return false;
  }

  float const npv = std::min(*nvtxs, (cms3_vtxs_nvtxs_t) 100); // Effective number of primary vertices
  float METxcorr = -(pfmet_XYcorr_xCoeffA*npv + pfmet_XYcorr_xCoeffB);
  float METycorr = -(pfmet_XYcorr_yCoeffA*npv + pfmet_XYcorr_yCoeffB);
  pfmet->setXYShift(METxcorr, METycorr);
  if (this->verbosity>=MiscUtils::DEBUG) IVYerr << "JetMETHandler::assignMETXYShifts: Applying an additional XY shift of ( " << METxcorr << ", " << METycorr << " )" << endl;

  return true;
}
bool JetMETHandler::applyMETParticleShifts(bool usePFCandidates, std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons, std::vector<PhotonObject*> const* photons){
  ParticleObject::LorentzVector_t pfmet_particleShift(0, 0, 0, 0);
  if (muons){
    for (auto const* part:(*muons)){
      if (!ParticleSelectionHelpers::isGoodMETParticle(part)) continue;
      ParticleObject::LorentzVector_t p4_uncorrected = part->uncorrected_p4();
      pfmet_particleShift += -(part->p4() - p4_uncorrected);
      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "JetMETHandler::applyMETParticleShifts: MET shift due to muon " << part->getUniqueIdentifier() << " p4 = " << -(part->p4() - p4_uncorrected) << endl;
    }
  }
  if (electrons){
    for (auto const* part:(*electrons)){
      if (!ParticleSelectionHelpers::isGoodMETParticle(part)) continue;
      // PF MET uses PF candidates directly.
      /*
      ParticleObject::LorentzVector_t p4_uncorrected = part->uncorrected_p4();
      if (!(part->testSelection(ElectronSelectionHelpers::kPFElectronId) && part->testSelection(ElectronSelectionHelpers::kPFMETSafe))){
        p4_uncorrected = ParticleObject::LorentzVector_t(part->extras.associated_pfcands_sum_px, part->extras.associated_pfcands_sum_py, 0, 0);
      }
      */
      ParticleObject::LorentzVector_t p4_uncorrected = ParticleObject::LorentzVector_t(part->extras.associated_pfcands_sum_px, part->extras.associated_pfcands_sum_py, 0, 0);
      pfmet_particleShift += -(part->p4() - p4_uncorrected);
      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "JetMETHandler::applyMETParticleShifts: MET shift due to electron " << part->getUniqueIdentifier() << " p4 = " << -(part->p4() - p4_uncorrected) << endl;
    }
  }
  if (photons){
    for (auto const* part:(*photons)){
      if (!ParticleSelectionHelpers::isGoodMETParticle(part)) continue;
      // PF MET uses PF candidates directly.
      /*
      ParticleObject::LorentzVector_t p4_uncorrected = part->uncorrected_p4();
      if (!(part->testSelection(PhotonSelectionHelpers::kPFPhotonId) && part->testSelection(PhotonSelectionHelpers::kPFMETSafe))){
        p4_uncorrected = ParticleObject::LorentzVector_t(part->extras.associated_pfcands_sum_px, part->extras.associated_pfcands_sum_py, 0, 0);
      }
      */
      ParticleObject::LorentzVector_t p4_uncorrected = ParticleObject::LorentzVector_t(part->extras.associated_pfcands_sum_px, part->extras.associated_pfcands_sum_py, 0, 0);
      pfmet_particleShift += -(part->p4() - p4_uncorrected);
      if (this->verbosity>=MiscUtils::DEBUG) IVYout << "JetMETHandler::applyMETParticleShifts: MET shift due to photon " << part->getUniqueIdentifier() << " p4 = " << -(part->p4() - p4_uncorrected) << endl;
    }
  }

  pfmet->setParticleShifts(pfmet_particleShift);
  pfpuppimet->setParticleShifts(pfmet_particleShift); // Particle shifts are the same for PF and PUPPI MET

  return true;
}

bool JetMETHandler::checkRevertMETFixVariables(BaseTree* tree) const{
  bool const isData = SampleHelpers::checkSampleIsData(tree->sampleIdentifier);

  // Only samples which may need MET fix reversion have these variables
  bool res = true;
  std::vector<TString> bnames;
  tree->getValidBranchNamesWithoutAlias(bnames, false);
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) res &= (std::find(bnames.cbegin(), bnames.cend(), JetMETHandler::colName_ak4jets + "_" + #NAME)!=bnames.cend());
  AK4JET_REVERTMETFIX_VARIABLES;
#undef AK4JET_VARIABLE
#define MET_VARIABLE(TYPE, NAME, DEFVAL) res &= (std::find(bnames.cbegin(), bnames.cend(), JetMETHandler::colName_pfmet + "_" + #NAME)!=bnames.cend());
  MET_RECORDED_REVERTMETFIX_VARIABLES;
  MET_SHIFT_RECORDED_REVERTMETFIX_CORE_VARIABLES;
  if (!isData){
    MET_SHIFT_RECORDED_REVERTMETFIX_GENINFO_VARIABLES;
  }
#undef MET_VARIABLE
  return res;
}

bool JetMETHandler::wrapTree(BaseTree* tree){
  if (!tree) return false;

  // 200314: The following is taken from https://lathomas.web.cern.ch/lathomas/METStuff/XYCorrections/XYMETCorrection.h
  // The formula is corr = -(A*npv + B).
  TString theDP = SampleHelpers::getDataPeriod();
  bool const isData = SampleHelpers::checkSampleIsData(tree->sampleIdentifier, &theDP);
  if (isData){
    if (theDP == "2016B"){
      pfmet_XYcorr_xCoeffA = -0.0478335; pfmet_XYcorr_xCoeffB = -0.108032;
      pfmet_XYcorr_yCoeffA = 0.125148; pfmet_XYcorr_yCoeffB = 0.355672;
    }
    else if (theDP == "2016C"){
      pfmet_XYcorr_xCoeffA = -0.0916985; pfmet_XYcorr_xCoeffB = 0.393247;
      pfmet_XYcorr_yCoeffA = 0.151445; pfmet_XYcorr_yCoeffB = 0.114491;
    }
    else if (theDP == "2016D"){
      pfmet_XYcorr_xCoeffA = -0.0581169; pfmet_XYcorr_xCoeffB = 0.567316;
      pfmet_XYcorr_yCoeffA = 0.147549; pfmet_XYcorr_yCoeffB = 0.403088;
    }
    else if (theDP == "2016E"){
      pfmet_XYcorr_xCoeffA = -0.065622; pfmet_XYcorr_xCoeffB = 0.536856;
      pfmet_XYcorr_yCoeffA = 0.188532; pfmet_XYcorr_yCoeffB = 0.495346;
    }
    else if (theDP == "2016F"){
      pfmet_XYcorr_xCoeffA = -0.0313322; pfmet_XYcorr_xCoeffB = 0.39866;
      pfmet_XYcorr_yCoeffA = 0.16081; pfmet_XYcorr_yCoeffB = 0.960177;
    }
    else if (theDP == "2016G"){
      pfmet_XYcorr_xCoeffA = 0.040803; pfmet_XYcorr_xCoeffB = -0.290384;
      pfmet_XYcorr_yCoeffA = 0.0961935; pfmet_XYcorr_yCoeffB = 0.666096;
    }
    else if (theDP == "2016H"){
      pfmet_XYcorr_xCoeffA = 0.0330868; pfmet_XYcorr_xCoeffB = -0.209534;
      pfmet_XYcorr_yCoeffA = 0.141513; pfmet_XYcorr_yCoeffB = 0.816732;
    }
    /*
    else if (theDP == "2017B"){
    pfmet_XYcorr_xCoeffA = -0.259456; pfmet_XYcorr_xCoeffB = 1.95372;
    pfmet_XYcorr_yCoeffA = 0.353928; pfmet_XYcorr_yCoeffB = -2.46685;
    }
    else if (theDP == "2017C"){
    pfmet_XYcorr_xCoeffA = -0.232763; pfmet_XYcorr_xCoeffB = 1.08318;
    pfmet_XYcorr_yCoeffA = 0.257719; pfmet_XYcorr_yCoeffB = -1.1745;
    }
    else if (theDP == "2017D"){
    pfmet_XYcorr_xCoeffA = -0.238067; pfmet_XYcorr_xCoeffB = 1.80541;
    pfmet_XYcorr_yCoeffA = 0.235989; pfmet_XYcorr_yCoeffB = -1.44354;
    }
    else if (theDP == "2017E"){
    pfmet_XYcorr_xCoeffA = -0.212352; pfmet_XYcorr_xCoeffB = 1.851;
    pfmet_XYcorr_yCoeffA = 0.157759; pfmet_XYcorr_yCoeffB = -0.478139;
    }
    else if (theDP == "2017F"){
    pfmet_XYcorr_xCoeffA = -0.232733; pfmet_XYcorr_xCoeffB = 2.24134;
    pfmet_XYcorr_yCoeffA = 0.213341; pfmet_XYcorr_yCoeffB = 0.684588;
    }
    */
    else if (theDP == "2017B"){
      pfmet_XYcorr_xCoeffA = -0.19563; pfmet_XYcorr_xCoeffB = 1.51859;
      pfmet_XYcorr_yCoeffA = 0.306987; pfmet_XYcorr_yCoeffB = -1.84713;
    }
    else if (theDP == "2017C"){
      pfmet_XYcorr_xCoeffA = -0.161661; pfmet_XYcorr_xCoeffB = 0.589933;
      pfmet_XYcorr_yCoeffA = 0.233569; pfmet_XYcorr_yCoeffB = -0.995546;
    }
    else if (theDP == "2017D"){
      pfmet_XYcorr_xCoeffA = -0.180911; pfmet_XYcorr_xCoeffB = 1.23553;
      pfmet_XYcorr_yCoeffA = 0.240155; pfmet_XYcorr_yCoeffB = -1.27449;
    }
    else if (theDP == "2017E"){
      pfmet_XYcorr_xCoeffA = -0.149494; pfmet_XYcorr_xCoeffB = 0.901305;
      pfmet_XYcorr_yCoeffA = 0.178212; pfmet_XYcorr_yCoeffB = -0.535537;
    }
    else if (theDP == "2017F"){
      pfmet_XYcorr_xCoeffA = -0.165154; pfmet_XYcorr_xCoeffB = 1.02018;
      pfmet_XYcorr_yCoeffA = 0.253794; pfmet_XYcorr_yCoeffB = 0.75776;
    }
    else if (theDP == "2018A"){
      pfmet_XYcorr_xCoeffA = 0.362865; pfmet_XYcorr_xCoeffB = -1.94505;
      pfmet_XYcorr_yCoeffA = 0.0709085; pfmet_XYcorr_yCoeffB = -0.307365;
    }
    else if (theDP == "2018B"){
      pfmet_XYcorr_xCoeffA = 0.492083; pfmet_XYcorr_xCoeffB = -2.93552;
      pfmet_XYcorr_yCoeffA = 0.17874; pfmet_XYcorr_yCoeffB = -0.786844;
    }
    else if (theDP == "2018C"){
      pfmet_XYcorr_xCoeffA = 0.521349; pfmet_XYcorr_xCoeffB = -1.44544;
      pfmet_XYcorr_yCoeffA = 0.118956; pfmet_XYcorr_yCoeffB = -1.96434;
    }
    else if (theDP == "2018D"){
      pfmet_XYcorr_xCoeffA = 0.531151; pfmet_XYcorr_xCoeffB = -1.37568;
      pfmet_XYcorr_yCoeffA = 0.0884639; pfmet_XYcorr_yCoeffB = -1.57089;
    }
    else{
      IVYerr << "JetMETHandler::wrapTree: Data period " << theDP << " is undefined for the data MET corrections." << endl;
      return false;
    }
  }
  else{
    switch (SampleHelpers::getDataYear()){
    case 2016:
      pfmet_XYcorr_xCoeffA = -0.195191; pfmet_XYcorr_xCoeffB = -0.170948;
      pfmet_XYcorr_yCoeffA = -0.0311891; pfmet_XYcorr_yCoeffB = 0.787627;
      break;
    case 2017:
      pfmet_XYcorr_xCoeffA = -0.217714; pfmet_XYcorr_xCoeffB = 0.493361;
      pfmet_XYcorr_yCoeffA = 0.177058; pfmet_XYcorr_yCoeffB = -0.336648;
      break;
    case 2018:
      pfmet_XYcorr_xCoeffA = 0.296713; pfmet_XYcorr_xCoeffB = -0.141506;
      pfmet_XYcorr_yCoeffA = 0.115685; pfmet_XYcorr_yCoeffB = 0.0128193;
      break;
    default:
      IVYerr << "JetMETHandler::wrapTree: Year " << SampleHelpers::getDataYear() << " is undefined for the MC MET corrections." << endl;
      return false;
      break;
    }
  }

  // Only samples which may need MET fix reversion have these variables
  flag_hasRevertMETFixVariables = checkRevertMETFixVariables(tree);

  return IvyBase::wrapTree(tree);
}

void JetMETHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

  bool const isData = SampleHelpers::checkSampleIsData(tree->sampleIdentifier);
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) \
this->addConsumed<std::vector<TYPE>*>(JetMETHandler::colName_ak4jets + "_" + #NAME); \
this->defineConsumedSloppy(JetMETHandler::colName_ak4jets + "_" + #NAME); \
tree->bookBranch<std::vector<TYPE>*>(JetMETHandler::colName_ak4jets + "_" + #NAME, nullptr);
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) \
this->addConsumed<std::vector<TYPE>*>(JetMETHandler::colName_ak8jets + "_" + #NAME); \
this->defineConsumedSloppy(JetMETHandler::colName_ak8jets + "_" + #NAME); \
tree->bookBranch<std::vector<TYPE>*>(JetMETHandler::colName_ak8jets + "_" + #NAME, nullptr);
  if (!isData){
    AK4JET_GENINFO_VARIABLES;
    AK8JET_GENINFO_VARIABLES;
  }
#undef AK8JET_VARIABLE
#undef AK4JET_VARIABLE

#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(JetMETHandler::colName_ak4jets + "_" + #NAME, nullptr);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK4JETS;
#undef AK4JET_VARIABLE
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<std::vector<TYPE>*>(JetMETHandler::colName_ak8jets + "_" + #NAME, nullptr);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK8JETS;
#undef AK8JET_VARIABLE

#define MET_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<TYPE>(JetMETHandler::colName_pfmet + "_" + #NAME, DEFVAL);
  MET_CORE_VARIABLES;
#undef MET_VARIABLE
#define MET_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<TYPE>(JetMETHandler::colName_pfpuppimet + "_" + #NAME, DEFVAL);
  MET_RECORDED_CORE_VARIABLES;
#undef MET_VARIABLE
  if (!isData){
#define MET_VARIABLE(TYPE, NAME, DEFVAL) \
this->addConsumed<TYPE>(JetMETHandler::colName_pfmet + "_" + #NAME); \
this->defineConsumedSloppy(JetMETHandler::colName_pfmet + "_" + #NAME); \
tree->bookBranch<TYPE>(JetMETHandler::colName_pfmet + "_" + #NAME, DEFVAL);
    MET_GENINFO_VARIABLES;
#undef MET_VARIABLE
#define MET_VARIABLE(TYPE, NAME, DEFVAL) \
this->addConsumed<TYPE>(JetMETHandler::colName_pfpuppimet + "_" + #NAME); \
this->defineConsumedSloppy(JetMETHandler::colName_pfpuppimet + "_" + #NAME); \
tree->bookBranch<TYPE>(JetMETHandler::colName_pfpuppimet + "_" + #NAME, DEFVAL);
    MET_RECORDED_GENINFO_VARIABLES;
#undef MET_VARIABLE
  }

  // Vertex variables
#define JETMET_METXY_VERTEX_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<TYPE>(JetMETHandler::colName_vertices + "_" + #NAME, DEFVAL);
  JETMET_METXY_VERTEX_VARIABLES;
#undef JETMET_METXY_VERTEX_VARIABLE

  // Only samples which may need MET fix reversion have these variables
  flag_hasRevertMETFixVariables = checkRevertMETFixVariables(tree);
  if (flag_hasRevertMETFixVariables){
    IVYout << "JetMETHandler::bookBranch: Booking MET fix reversion variables..." << endl;
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) \
this->addConsumed<std::vector<TYPE>*>(JetMETHandler::colName_ak4jets + "_" + #NAME); \
this->defineConsumedSloppy(JetMETHandler::colName_ak4jets + "_" + #NAME); \
tree->bookBranch<std::vector<TYPE>*>(JetMETHandler::colName_ak4jets + "_" + #NAME, nullptr);
    AK4JET_REVERTMETFIX_VARIABLES;
#undef AK4JET_VARIABLE
#define MET_VARIABLE(TYPE, NAME, DEFVAL) \
this->addConsumed<TYPE>(JetMETHandler::colName_pfmet + "_" + #NAME); \
this->defineConsumedSloppy(JetMETHandler::colName_pfmet + "_" + #NAME); \
tree->bookBranch<TYPE>(JetMETHandler::colName_pfmet + "_" + #NAME, DEFVAL);
    MET_RECORDED_REVERTMETFIX_VARIABLES;
    MET_SHIFT_RECORDED_REVERTMETFIX_CORE_VARIABLES;
    if (!isData){
      MET_SHIFT_RECORDED_REVERTMETFIX_GENINFO_VARIABLES;
    }
#undef MET_VARIABLE

    if (isData){
#define RUNLUMIEVENT_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<TYPE>(#NAME, DEFVAL); this->addConsumed<TYPE>(#NAME); this->defineConsumedSloppy(#NAME);
      RUNLUMIEVENT_VARIABLES;
#undef RUNLUMIEVENT_VARIABLE
    }
  }
}

void JetMETHandler::registerOverlapMaps(
  OverlapMapHandler<MuonObject, AK4JetObject>& overlapMap_muons_ak4jets_,
  OverlapMapHandler<MuonObject, AK8JetObject>& overlapMap_muons_ak8jets_,
  OverlapMapHandler<ElectronObject, AK4JetObject>& overlapMap_electrons_ak4jets_,
  OverlapMapHandler<ElectronObject, AK8JetObject>& overlapMap_electrons_ak8jets_,
  OverlapMapHandler<PhotonObject, AK4JetObject>& overlapMap_photons_ak4jets_,
  OverlapMapHandler<PhotonObject, AK8JetObject>& overlapMap_photons_ak8jets_
){
  overlapMap_muons_ak4jets = &overlapMap_muons_ak4jets_;
  overlapMap_muons_ak8jets = &overlapMap_muons_ak8jets_;
  overlapMap_electrons_ak4jets = &overlapMap_electrons_ak4jets_;
  overlapMap_electrons_ak8jets = &overlapMap_electrons_ak8jets_;
  overlapMap_photons_ak4jets = &overlapMap_photons_ak4jets_;
  overlapMap_photons_ak8jets = &overlapMap_photons_ak8jets_;
  hasOverlapMaps = true;
}


#undef VECTOR_ITERATOR_HANDLER_DIRECTIVES
