#include <cassert>
#include "RunLumiEventBlock.h"
#include "ParticleObjectHelpers.h"
#include "SamplesCore.h"
#include "JetMETHandler.h"
#include "AK4JetSelectionHelpers.h"
#include "AK8JetSelectionHelpers.h"
#include "ParticleSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace MELAStreamHelpers;


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
JETMET_METXY_VERTEX_VARIABLE(unsigned int, nvtxs, 0)


const std::string JetMETHandler::colName_ak4jets = "ak4jets";
const std::string JetMETHandler::colName_ak8jets = "ak8jets";
const std::string JetMETHandler::colName_pfmet = "pfmet";
const std::string JetMETHandler::colName_pfpuppimet = "puppimet";
const std::string JetMETHandler::colName_vertices = "vtxs";

JetMETHandler::JetMETHandler() :
  IvyBase(),
  pfmet(nullptr),
  pfpuppimet(nullptr),

  pfmet_XYcorr_xCoeffA(0),
  pfmet_XYcorr_xCoeffB(0),
  pfmet_XYcorr_yCoeffA(0),
  pfmet_XYcorr_yCoeffB(0),

  has_genmatching(false)
{
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(JetMETHandler::colName_ak4jets + "_" + #NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK4JETS;
#undef AK4JET_VARIABLE
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(JetMETHandler::colName_ak8jets + "_" + #NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK8JETS;
#undef AK8JET_VARIABLE
#define MET_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<TYPE>(JetMETHandler::colName_pfmet + "_" + #NAME);
  MET_RECORDED_VARIABLES;
  MET_JERSHIFT_VARIABLES;
#undef MET_VARIABLE
#define MET_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<TYPE>(JetMETHandler::colName_pfpuppimet + "_" + #NAME);
  MET_RECORDED_VARIABLES;
#undef MET_VARIABLE

  // Vertex variables
#define JETMET_METXY_VERTEX_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<TYPE>(JetMETHandler::colName_vertices + "_" + #NAME);
  JETMET_METXY_VERTEX_VARIABLES;
#undef JETMET_METXY_VERTEX_VARIABLE
}

void JetMETHandler::clear(){
  for (auto*& prod:ak4jets) delete prod;
  ak4jets.clear();
  for (auto*& prod:ak8jets) delete prod;
  ak8jets.clear();
  delete pfmet; pfmet=nullptr;
  delete pfpuppimet; pfpuppimet=nullptr;
}

bool JetMETHandler::constructJetMET(SystematicsHelpers::SystematicVariationTypes const& syst, std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons, std::vector<PhotonObject*> const* photons){
  clear();
  if (!currentTree) return false;

  bool res = (
    constructAK4Jets(syst) && constructAK8Jets(syst)
    &&
    constructMET(syst) && assignMETXYShifts(syst) && applyMETParticleShifts(muons, electrons, photons)
    &&
    applyJetCleaning(muons, electrons, photons)
    );

  return res;
}

bool JetMETHandler::constructAK4Jets(SystematicsHelpers::SystematicVariationTypes const& syst){
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_ak4jets_##NAME, itEnd_ak4jets_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK4JETS;
  AK4JET_GENINFO_VARIABLES;
#undef AK4JET_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(JetMETHandler::colName_ak4jets + "_" + #NAME, &itBegin_ak4jets_##NAME, &itEnd_ak4jets_##NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK4JETS;
  if (this->has_genmatching){
    AK4JET_GENINFO_VARIABLES;
  }
#undef AK4JET_VARIABLE
  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "JetMETHandler::constructAK4Jets: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructAK4Jets: All variables are set up!" << endl;

  /************/
  /* ak4 jets */
  /************/
  size_t nak4jets = (itEnd_ak4jets_pt - itBegin_ak4jets_pt);
  ak4jets.reserve(nak4jets);
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) auto it_ak4jets_##NAME = itBegin_ak4jets_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK4JETS;
  AK4JET_GENINFO_VARIABLES;
#undef AK4JET_VARIABLE
  {
    size_t ip=0;
    while (it_ak4jets_pt != itEnd_ak4jets_pt){
      if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructAK4Jets: Attempting ak4 jet " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_ak4jets_pt, *it_ak4jets_eta, *it_ak4jets_phi, *it_ak4jets_mass); // Yes you have to do this on a separate line because CMSSW...
      ak4jets.push_back(new AK4JetObject(momentum));
      AK4JetObject*& obj = ak4jets.back();

      // Set extras
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) obj->extras.NAME = *it_ak4jets_##NAME;
      AK4JET_RECO_VARIABLES;
      if (this->has_genmatching){
        AK4JET_GENINFO_VARIABLES;
      }
#undef AK4JET_VARIABLE

      // Replace momentum
      obj->makeFinalMomentum(syst);

      // Set the selection bits
      AK4JetSelectionHelpers::setSelectionBits(*obj);

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) it_ak4jets_##NAME++;
      VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK4JETS;
      if (this->has_genmatching){
        AK4JET_GENINFO_VARIABLES;
      }
#undef AK4JET_VARIABLE
    }
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(ak4jets);

  return true;
}
bool JetMETHandler::constructAK8Jets(SystematicsHelpers::SystematicVariationTypes const& syst){
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) std::vector<TYPE>::const_iterator itBegin_ak8jets_##NAME, itEnd_ak8jets_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK8JETS;
  AK8JET_GENINFO_VARIABLES;
#undef AK8JET_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumedCIterators<std::vector<TYPE>>(JetMETHandler::colName_ak8jets + "_" + #NAME, &itBegin_ak8jets_##NAME, &itEnd_ak8jets_##NAME);
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK8JETS;
  if (this->has_genmatching){
    AK8JET_GENINFO_VARIABLES;
  }
#undef AK8JET_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "JetMETHandler::constructAK8Jets: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructAK8Jets: All variables are set up!" << endl;

  /************/
  /* ak8 jets */
  /************/
  size_t nak8jets = (itEnd_ak8jets_pt - itBegin_ak8jets_pt);
  ak8jets.reserve(nak8jets);
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) auto it_ak8jets_##NAME = itBegin_ak8jets_##NAME;
  VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK8JETS;
  AK8JET_GENINFO_VARIABLES;
#undef AK8JET_VARIABLE
  {
    size_t ip=0;
    while (it_ak8jets_pt != itEnd_ak8jets_pt){
      if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructAK8Jets: Attempting ak8 jet " << ip << "..." << endl;

      ParticleObject::LorentzVector_t momentum;
      momentum = ParticleObject::PolarLorentzVector_t(*it_ak8jets_pt, *it_ak8jets_eta, *it_ak8jets_phi, *it_ak8jets_mass); // Yes you have to do this on a separate line because CMSSW...
      ak8jets.push_back(new AK8JetObject(momentum));
      AK8JetObject*& obj = ak8jets.back();

      // Set extras
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) obj->extras.NAME = *it_ak8jets_##NAME;
      AK8JET_RECO_VARIABLES;
      if (this->has_genmatching){
        AK8JET_GENINFO_VARIABLES;
      }
#undef AK8JET_VARIABLE

      // Replace momentum
      obj->makeFinalMomentum(syst);

      // Set the selection bits
      AK8JetSelectionHelpers::setSelectionBits(*obj);

      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Success!" << endl;

      ip++;
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) it_ak8jets_##NAME++;
      VECTOR_ITERATOR_HANDLER_DIRECTIVES_AK8JETS;
      if (this->has_genmatching){
        AK8JET_GENINFO_VARIABLES;
      }
#undef AK8JET_VARIABLE
    }
  }
  // Sort particles
  ParticleObjectHelpers::sortByGreaterPt(ak8jets);

  return true;
}
bool JetMETHandler::constructMET(SystematicsHelpers::SystematicVariationTypes const& syst){
#define MET_VARIABLE(TYPE, NAME, DEFVAL) TYPE const* pfmet_##NAME = nullptr;
  MET_RECORDED_VARIABLES;
  MET_JERSHIFT_VARIABLES;
#undef MET_VARIABLE
#define MET_VARIABLE(TYPE, NAME, DEFVAL) TYPE const* pfpuppimet_##NAME = nullptr;
  MET_RECORDED_VARIABLES;
#undef MET_VARIABLE

  // Beyond this point starts checks and selection
  bool allVariablesPresent = true;
#define MET_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumed(JetMETHandler::colName_pfmet + "_" + #NAME, pfmet_##NAME);
  MET_RECORDED_VARIABLES;
  MET_JERSHIFT_VARIABLES;
#undef MET_VARIABLE
#define MET_VARIABLE(TYPE, NAME, DEFVAL) allVariablesPresent &= this->getConsumed(JetMETHandler::colName_pfpuppimet + "_" + #NAME, pfpuppimet_##NAME);
  MET_RECORDED_VARIABLES;
#undef MET_VARIABLE

  if (!allVariablesPresent){
    if (this->verbosity>=TVar::ERROR) MELAerr << "JetMETHandler::constructMET: Not all variables are consumed properly!" << endl;
    assert(0);
  }

  if (this->verbosity>=TVar::DEBUG) MELAout << "JetMETHandler::constructMET: All variables are set up!" << endl;

  /**********/
  /* PF MET */
  /**********/
  pfmet = new METObject();
#define MET_VARIABLE(TYPE, NAME, DEFVAL) pfmet->extras.NAME = *pfmet_##NAME;
  MET_RECORDED_VARIABLES;
  MET_JERSHIFT_VARIABLES;
#undef MET_VARIABLE
#define MET_VARIABLE(TYPE, NAME, DEFVAL) pfmet->extras.NAME = pfmet->extras.met_Nominal;
  MET_EXTRA_PT_VARIABLES;
#undef MET_VARIABLE
#define MET_VARIABLE(TYPE, NAME, DEFVAL) pfmet->extras.NAME = pfmet->extras.metPhi_Nominal;
  MET_EXTRA_PHI_VARIABLES;
#undef MET_VARIABLE
  pfmet->setSystematic(syst);

  /*************/
  /* PUPPI MET */
  /*************/
  pfpuppimet = new METObject();
#define MET_VARIABLE(TYPE, NAME, DEFVAL) pfpuppimet->extras.NAME = *pfpuppimet_##NAME;
  MET_RECORDED_VARIABLES;
#undef MET_VARIABLE
#define MET_VARIABLE(TYPE, NAME, DEFVAL) pfpuppimet->extras.NAME = pfpuppimet->extras.met_Nominal;
  MET_EXTRA_PT_VARIABLES;
#undef MET_VARIABLE
#define MET_VARIABLE(TYPE, NAME, DEFVAL) pfpuppimet->extras.NAME = pfpuppimet->extras.metPhi_Nominal;
  MET_EXTRA_PHI_VARIABLES;
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
    if (this->verbosity>=TVar::ERROR) MELAerr << "JetMETHandler::assignMETXYShifts: Not all variables are consumed properly!" << endl;
    assert(0);
    return false;
  }

  float const npv = std::min(*nvtxs, (unsigned int) 100); // Effective number of primary vertices
  float METxcorr = -(pfmet_XYcorr_xCoeffA*npv + pfmet_XYcorr_xCoeffB);
  float METycorr = -(pfmet_XYcorr_yCoeffA*npv + pfmet_XYcorr_yCoeffB);
  pfmet->setXYShift(METxcorr, METycorr);

  return true;
}

bool JetMETHandler::applyMETParticleShifts(std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons, std::vector<PhotonObject*> const* photons){
  ParticleObject::LorentzVector_t pfmet_particleShift(0, 0, 0, 0);
  if (muons){
    for (auto const* part:*(muons)){
      if (!ParticleSelectionHelpers::isLooseParticle(part)) continue;
      ParticleObject::LorentzVector_t p4_uncorrected; p4_uncorrected = ParticleObject::PolarLorentzVector_t(part->uncorrected_pt(), part->eta(), part->phi(), part->mass());
      pfmet_particleShift += -(part->p4() - p4_uncorrected);
    }
  }
  if (electrons){
    for (auto const* part:*(electrons)){
      if (!ParticleSelectionHelpers::isLooseParticle(part)) continue;
      ParticleObject::LorentzVector_t::Scalar diffCorrScale = 1.f - 1.f/part->currentSystScale;
      pfmet_particleShift += -part->p4()*diffCorrScale;
    }
  }
  if (photons){
    for (auto const* part:*(photons)){
      if (!ParticleSelectionHelpers::isLooseParticle(part)) continue;
      ParticleObject::LorentzVector_t::Scalar diffCorrScale = 1.f - 1.f/part->currentSystScale;
      pfmet_particleShift += -part->p4()*diffCorrScale;
    }
  }

  pfmet->setParticleShifts(pfmet_particleShift);
  pfpuppimet->setParticleShifts(pfmet_particleShift); // Particle shofts are the same for PF and PUPPI MET

  return true;
}

bool JetMETHandler::applyJetCleaning(std::vector<MuonObject*> const* muons, std::vector<ElectronObject*> const* electrons, std::vector<PhotonObject*> const* photons){
  std::vector<AK4JetObject*> ak4jets_new; ak4jets_new.reserve(ak4jets.size());
  for (auto*& jet:ak4jets){
    bool doSkip=false;
    if (muons){
      for (auto const* part:*(muons)){
        if (!ParticleSelectionHelpers::isLooseParticle(part)) continue;
        if (reco::deltaR(jet->p4(), part->p4())<jet->ConeRadiusConstant){ doSkip=true; break; }
      }
    }
    if (electrons){
      for (auto const* part:*(electrons)){
        if (!ParticleSelectionHelpers::isLooseParticle(part)) continue;
        if (reco::deltaR(jet->p4(), part->p4())<jet->ConeRadiusConstant){ doSkip=true; break; }
      }
    }
    if (photons){
      for (auto const* part:*(photons)){
        if (!ParticleSelectionHelpers::isLooseParticle(part)) continue;
        if (reco::deltaR(jet->p4(), part->p4())<jet->ConeRadiusConstant){ doSkip=true; break; }
      }
    }
    if (!doSkip) ak4jets_new.push_back(jet);
    else delete jet;
  }
  ak4jets = ak4jets_new;

  std::vector<AK8JetObject*> ak8jets_new; ak8jets_new.reserve(ak8jets.size());
  for (auto*& jet:ak8jets){
    bool doSkip=false;
    if (muons){
      for (auto const* part:*(muons)){
        if (!ParticleSelectionHelpers::isLooseParticle(part)) continue;
        if (reco::deltaR(jet->p4(), part->p4())<jet->ConeRadiusConstant){ doSkip=true; break; }
      }
    }
    if (electrons){
      for (auto const* part:*(electrons)){
        if (!ParticleSelectionHelpers::isLooseParticle(part)) continue;
        if (reco::deltaR(jet->p4(), part->p4())<jet->ConeRadiusConstant){ doSkip=true; break; }
      }
    }
    if (photons){
      for (auto const* part:*(photons)){
        if (!ParticleSelectionHelpers::isLooseParticle(part)) continue;
        if (reco::deltaR(jet->p4(), part->p4())<jet->ConeRadiusConstant){ doSkip=true; break; }
      }
    }
    if (!doSkip) ak8jets_new.push_back(jet);
    else delete jet;
  }
  ak8jets = ak8jets_new;

  return true;
}

void JetMETHandler::checkOptionalInfo(BaseTree* tree, bool& flag_genmatching){
  flag_genmatching = true;

  std::vector<TString> bnames;
  tree->getValidBranchNamesWithoutAlias(bnames, false);

#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) flag_genmatching &= (std::find(bnames.cbegin(), bnames.cend(), JetMETHandler::colName_ak4jets + "_" + #NAME)!=bnames.cend());
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) flag_genmatching &= (std::find(bnames.cbegin(), bnames.cend(), JetMETHandler::colName_ak8jets + "_" + #NAME)!=bnames.cend());
  AK4JET_GENINFO_VARIABLES;
  AK8JET_GENINFO_VARIABLES;
#undef AK8JET_VARIABLE
#undef AK4JET_VARIABLE
}

bool JetMETHandler::wrapTree(BaseTree* tree){
  if (!tree) return false;

  // 200314: The following is taken from https://lathomas.web.cern.ch/lathomas/METStuff/XYCorrections/XYMETCorrection.h
  // The formula is corr = -(A*npv + B).
  TString theDP = SampleHelpers::theDataPeriod;
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
      MELAerr << "JetMETHandler::wrapTree: Data period " << theDP << " is undefined for the data MET corrections." << endl;
      return false;
    }
  }
  else{
    switch (SampleHelpers::theDataYear){
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
      MELAerr << "JetMETHandler::wrapTree: Year " << SampleHelpers::theDataYear << " is undefined for the MC MET corrections." << endl;
      return false;
      break;
    }
  }

  JetMETHandler::checkOptionalInfo(tree, this->has_genmatching);

  return IvyBase::wrapTree(tree);
}

void JetMETHandler::bookBranches(BaseTree* tree){
  if (!tree) return;

  JetMETHandler::checkOptionalInfo(tree, this->has_genmatching);
#define AK4JET_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(JetMETHandler::colName_ak4jets + "_" + #NAME); this->defineConsumedSloppy(#NAME); tree->bookBranch<std::vector<TYPE>*>(JetMETHandler::colName_ak4jets + "_" + #NAME, nullptr);
#define AK8JET_VARIABLE(TYPE, NAME, DEFVAL) this->addConsumed<std::vector<TYPE>*>(JetMETHandler::colName_ak8jets + "_" + #NAME); this->defineConsumedSloppy(#NAME); tree->bookBranch<std::vector<TYPE>*>(JetMETHandler::colName_ak8jets + "_" + #NAME, nullptr);
  if (this->has_genmatching){
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
  MET_RECORDED_VARIABLES;
  MET_JERSHIFT_VARIABLES;
#undef MET_VARIABLE
#define MET_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<TYPE>(JetMETHandler::colName_pfpuppimet + "_" + #NAME, DEFVAL);
  MET_RECORDED_VARIABLES;
#undef MET_VARIABLE

  // Vertex variables
#define JETMET_METXY_VERTEX_VARIABLE(TYPE, NAME, DEFVAL) tree->bookBranch<TYPE>(JetMETHandler::colName_vertices + "_" + #NAME, DEFVAL);
  JETMET_METXY_VERTEX_VARIABLES;
#undef JETMET_METXY_VERTEX_VARIABLE
}


#undef VECTOR_ITERATOR_HANDLER_DIRECTIVES
