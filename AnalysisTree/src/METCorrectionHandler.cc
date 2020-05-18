#include <cassert>
#include "HostHelpersCore.h"
#include "SampleHelpersCore.h"
#include "SamplesCore.h"
#include "METCorrectionHandler.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace SampleHelpers;
using namespace MELAStreamHelpers;


METCorrectionParameters::METCorrectionParameters(){
  sigmas_nominal.reserve(4);
  sigmas_dn.reserve(4);
  sigmas_up.reserve(4);
  fracs.reserve(3);
}
METCorrectionParameters::METCorrectionParameters(METCorrectionParameters const& other) :
  sigmas_nominal(other.sigmas_nominal),
  sigmas_dn(other.sigmas_dn),
  sigmas_up(other.sigmas_up),
  fracs(other.fracs)
{}
METCorrectionParameters& METCorrectionParameters::operator=(METCorrectionParameters const& other){
  METCorrectionParameters tmp(other);
  swap(tmp);
  return *this;
}
void METCorrectionParameters::swap(METCorrectionParameters& other){
  std::swap(sigmas_nominal, other.sigmas_nominal);
  std::swap(sigmas_dn, other.sigmas_dn);
  std::swap(sigmas_up, other.sigmas_up);
  std::swap(fracs, other.fracs);
}

void METCorrectionParameters::translateFractions(){
  // Convert fractions from recursive to non-recursive
  for (unsigned int i=0; i<fracs.size(); i++){
    float& lastFrac = fracs.at(i);
    float prevFracProd = 1;
    for (unsigned int j=0; j<i; j++) prevFracProd *= (1. - fracs.at(j));
    lastFrac *= prevFracProd;
  }
  // Accumulate
  for (unsigned int i=1; i<fracs.size(); i++) fracs.at(i) += fracs.at(i-1);
}


METCorrectionHandler::METCorrectionHandler() : ScaleFactorHandlerBase()
{
  setup();
}

METCorrectionHandler::~METCorrectionHandler(){ this->reset(); }

bool METCorrectionHandler::setup(){
  bool res = true;
  this->reset();

  TString strinputcore = ANALYSISTREEPKGDATAPATH + Form("ScaleFactors/METResolution/%i/", theDataYear);
  HostHelpers::ExpandEnvironmentVariables(strinputcore);

  auto const allowedDataPeriods = SampleHelpers::getValidDataPeriods();
  static const std::vector<SystematicsHelpers::SystematicVariationTypes> allowedSysts{
    SystematicsHelpers::sNominal,
    SystematicsHelpers::eJECDn, SystematicsHelpers::eJECUp,
    SystematicsHelpers::eJERDn, SystematicsHelpers::eJERUp,
    SystematicsHelpers::ePUDn, SystematicsHelpers::ePUUp
  };
  for (auto const& period:allowedDataPeriods){
    pfmet_XY_JER_map[period]=std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>();
    pfmet_JER_map[period]=std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>();
    pfmet_XY_map[period]=std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>();
    pfmet_map[period]=std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>();
    puppimet_map[period]=std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>();
    for (auto const& syst:allowedSysts){
      pfmet_XY_JER_map[period][syst]=METCorrectionParameters();
      pfmet_JER_map[period][syst]=METCorrectionParameters();
      pfmet_XY_map[period][syst]=METCorrectionParameters();
      pfmet_map[period][syst]=METCorrectionParameters();
      puppimet_map[period][syst]=METCorrectionParameters();
    }
    readFile(strinputcore + "fitparameters_pfmet_JEC_XY_JER_PartMomShifts_" + period + ".txt", pfmet_XY_JER_map[period]);
    readFile(strinputcore + "fitparameters_pfmet_JEC_JER_PartMomShifts_" + period + ".txt", pfmet_JER_map[period]);
    readFile(strinputcore + "fitparameters_pfmet_JEC_XY_PartMomShifts_" + period + ".txt", pfmet_XY_map[period]);
    readFile(strinputcore + "fitparameters_pfmet_JEC_PartMomShifts_" + period + ".txt", pfmet_map[period]);
    readFile(strinputcore + "fitparameters_puppimet_JEC_PartMomShifts_" + period + ".txt", puppimet_map[period]);
  }

  return res;
}
void METCorrectionHandler::readFile(TString const& strinput, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>& pars){
  if (!HostHelpers::FileReadable(strinput.Data())){
    MELAerr << "METCorrectionHandler::setup: File " << strinput << " is not readable." << endl;
    assert(0);
  }

  SystematicsHelpers::SystematicVariationTypes currentSyst = SystematicsHelpers::nSystematicVariations; // Data

  ifstream fin;
  fin.open(strinput.Data(), std::ios_base::in);
  if (fin.good()){
    string str_in;
    while (!fin.eof()){
      getline(fin, str_in);
      if (str_in.find("Final fitted parameters")!=std::string::npos){
        if (str_in.find("data")!=std::string::npos) currentSyst = SystematicsHelpers::nSystematicVariations;
        else if (str_in.find("MC")!=std::string::npos && str_in.find("Nominal")!=std::string::npos) currentSyst = SystematicsHelpers::sNominal;
        else if (str_in.find("MC")!=std::string::npos && str_in.find("JECDn ")!=std::string::npos) currentSyst = SystematicsHelpers::eJECDn;
        else if (str_in.find("MC")!=std::string::npos && str_in.find("JECUp ")!=std::string::npos) currentSyst = SystematicsHelpers::eJECUp;
        else if (str_in.find("MC")!=std::string::npos && str_in.find("JERDn ")!=std::string::npos) currentSyst = SystematicsHelpers::eJERDn;
        else if (str_in.find("MC")!=std::string::npos && str_in.find("JERUp ")!=std::string::npos) currentSyst = SystematicsHelpers::eJERUp;
        else if (str_in.find("MC")!=std::string::npos && str_in.find("PUDn ")!=std::string::npos) currentSyst = SystematicsHelpers::ePUDn;
        else if (str_in.find("MC")!=std::string::npos && str_in.find("PUUp ")!=std::string::npos) currentSyst = SystematicsHelpers::ePUUp;
        else{
          MELAerr << "METCorrectionHandler::readFile: Cannot determine the systematic for line " << str_in << endl;
          assert(0);
        }
      }
      else if (str_in.find("Sigma")!=std::string::npos || str_in.find("Frac")!=std::string::npos){
        std::string dummy;
        float val, errlow, errhigh;
        stringstream ss(str_in);
        ss >> dummy >> dummy >> dummy >> val >> errlow >> errhigh;
        if (str_in.find("Sigma")!=std::string::npos){
          pars[currentSyst].sigmas_nominal.push_back(val);
          pars[currentSyst].sigmas_dn.push_back(val - std::abs(errlow));
          pars[currentSyst].sigmas_up.push_back(val + std::abs(errhigh));
        }
        else pars[currentSyst].fracs.push_back(val);
      }
    }
  }
  fin.close();

  for (auto it:pars) it.second.translateFractions();
}

void METCorrectionHandler::reset(){
  pfmet_XY_JER_map.clear();
  pfmet_JER_map.clear();
  pfmet_XY_map.clear();
  pfmet_map.clear();
  puppimet_map.clear();
}

void METCorrectionHandler::applyCorrections(
  TString const& effDataPeriod, SystematicsHelpers::SystematicVariationTypes const& syst,
  float const& genMET, float const& genMETPhi,
  METObject* obj, bool isPFMET
) const{
  static const std::vector<SystematicsHelpers::SystematicVariationTypes> allowedSysts{
    SystematicsHelpers::sNominal,
    SystematicsHelpers::eJECDn, SystematicsHelpers::eJECUp,
    SystematicsHelpers::eJERDn, SystematicsHelpers::eJERUp,
    SystematicsHelpers::ePUDn, SystematicsHelpers::ePUUp
  };
  SystematicsHelpers::SystematicVariationTypes effSyst = SystematicsHelpers::sNominal;
  if (HelperFunctions::checkListVariable(allowedSysts, syst)) effSyst = syst;

  TRandom3 rand; rand.SetSeed(std::abs(static_cast<int>(sin(genMETPhi)*100000.)));
  float const frac_x = rand.Uniform();
  ParticleObject::LorentzVector_t const genmet_p4(genMET*std::cos(genMETPhi), genMET*std::sin(genMETPhi), 0, 0);

  for (unsigned short iXY=0; iXY<2; iXY++){
    for (unsigned short iJER=0; iJER<2; iJER++){
      std::unordered_map<
        TString,
        std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>
      > const* theMap = nullptr;
      if (!isPFMET) theMap = &puppimet_map;
      else{
        if (iXY==0 && iJER==0) theMap = &pfmet_map;
        else if (iXY==0 && iJER==1) theMap = &pfmet_JER_map;
        else if (iXY==1 && iJER==0) theMap = &pfmet_XY_map;
        else theMap = &pfmet_XY_JER_map;
      }

      auto it = theMap->find(effDataPeriod);
      if (it==theMap->cend()){
        MELAerr << "METCorrectionHandler::applyCorrections: The data period " << effDataPeriod << " for (XY, JER) = (" << iXY << ", " << iJER << ") cannot be found!" << endl;
        assert(0);
      }
      METCorrectionParameters const& pars_data = it->second.find(SystematicsHelpers::nSystematicVariations)->second;
      METCorrectionParameters const& pars_MC = it->second.find(effSyst)->second;

      unsigned int iGaussian = pars_data.sigmas_nominal.size()-1;
      for (unsigned int ifrac=0; ifrac<pars_data.fracs.size(); ifrac++){
        if (frac_x<pars_data.fracs.at(ifrac)){
          iGaussian = ifrac;
          break;
        }
      }

      float const& sigma_nominal_data = pars_data.sigmas_nominal.at(iGaussian);
      float const& sigma_dn_data = pars_data.sigmas_dn.at(iGaussian);
      float const& sigma_up_data = pars_data.sigmas_up.at(iGaussian);

      float const& sigma_nominal_MC = pars_MC.sigmas_nominal.at(iGaussian);
      float const& sigma_dn_MC = pars_MC.sigmas_dn.at(iGaussian);
      float const& sigma_up_MC = pars_MC.sigmas_up.at(iGaussian);

      float SF = 1;
      if (syst == SystematicsHelpers::eMETDn) SF = sigma_dn_data / sigma_up_MC;
      else if (syst == SystematicsHelpers::eMETUp) SF = sigma_up_data / sigma_dn_MC;
      else SF = sigma_nominal_data / sigma_nominal_MC;

      for (unsigned short iPMS=0; iPMS<2; iPMS++){
        ParticleObject::LorentzVector_t met_p4_diff = obj->p4((bool) iXY, (bool) iJER, (bool) iPMS) - genmet_p4;
        met_p4_diff *= SF;
        obj->setMETCorrection(met_p4_diff, (bool) iXY, (bool) iJER, (bool) iPMS);
      }
    }
  }
}

void METCorrectionHandler::printParameters() const{
  MELAout << "****************************" << endl;
  MELAout << "METCorrectionHandler::printParameters: Parameters of pfmet_XY_JER corrections:" << endl;
  for (auto it:pfmet_XY_JER_map){
    for (auto it2:it.second){
      MELAout << "\t- Period " << it.first << ", systematic " << (it2.first!=SystematicsHelpers::nSystematicVariations ? SystematicsHelpers::getSystName(it2.first) : "data") << ": " << endl;
      MELAout << "\t\t- Fractions: " << it2.second.fracs << endl;
      MELAout << "\t\t- Nominal sigmas: " << it2.second.sigmas_nominal << endl;
      MELAout << "\t\t- MET dn sigmas: " << it2.second.sigmas_dn << endl;
      MELAout << "\t\t- MET up sigmas: " << it2.second.sigmas_up << endl;
    }
  }
  MELAout << "****************************" << endl;
  MELAout << "METCorrectionHandler::printParameters: Parameters of pfmet_JER corrections:" << endl;
  for (auto it:pfmet_JER_map){
    for (auto it2:it.second){
      MELAout << "\t- Period " << it.first << ", systematic " << (it2.first!=SystematicsHelpers::nSystematicVariations ? SystematicsHelpers::getSystName(it2.first) : "data") << ": " << endl;
      MELAout << "\t\t- Fractions: " << it2.second.fracs << endl;
      MELAout << "\t\t- Nominal sigmas: " << it2.second.sigmas_nominal << endl;
      MELAout << "\t\t- MET dn sigmas: " << it2.second.sigmas_dn << endl;
      MELAout << "\t\t- MET up sigmas: " << it2.second.sigmas_up << endl;
    }
  }
  MELAout << "****************************" << endl;
  MELAout << "METCorrectionHandler::printParameters: Parameters of pfmet_XY corrections:" << endl;
  for (auto it:pfmet_XY_map){
    for (auto it2:it.second){
      MELAout << "\t- Period " << it.first << ", systematic " << (it2.first!=SystematicsHelpers::nSystematicVariations ? SystematicsHelpers::getSystName(it2.first) : "data") << ": " << endl;
      MELAout << "\t\t- Fractions: " << it2.second.fracs << endl;
      MELAout << "\t\t- Nominal sigmas: " << it2.second.sigmas_nominal << endl;
      MELAout << "\t\t- MET dn sigmas: " << it2.second.sigmas_dn << endl;
      MELAout << "\t\t- MET up sigmas: " << it2.second.sigmas_up << endl;
    }
  }
  MELAout << "****************************" << endl;
  MELAout << "METCorrectionHandler::printParameters: Parameters of pfmet bare corrections:" << endl;
  for (auto it:pfmet_map){
    for (auto it2:it.second){
      MELAout << "\t- Period " << it.first << ", systematic " << (it2.first!=SystematicsHelpers::nSystematicVariations ? SystematicsHelpers::getSystName(it2.first) : "data") << ": " << endl;
      MELAout << "\t\t- Fractions: " << it2.second.fracs << endl;
      MELAout << "\t\t- Nominal sigmas: " << it2.second.sigmas_nominal << endl;
      MELAout << "\t\t- MET dn sigmas: " << it2.second.sigmas_dn << endl;
      MELAout << "\t\t- MET up sigmas: " << it2.second.sigmas_up << endl;
    }
  }
  MELAout << "****************************" << endl;
  MELAout << "METCorrectionHandler::printParameters: Parameters of puppimet corrections:" << endl;
  for (auto it:puppimet_map){
    for (auto it2:it.second){
      MELAout << "\t- Period " << it.first << ", systematic " << (it2.first!=SystematicsHelpers::nSystematicVariations ? SystematicsHelpers::getSystName(it2.first) : "data") << ": " << endl;
      MELAout << "\t\t- Fractions: " << it2.second.fracs << endl;
      MELAout << "\t\t- Nominal sigmas: " << it2.second.sigmas_nominal << endl;
      MELAout << "\t\t- MET dn sigmas: " << it2.second.sigmas_dn << endl;
      MELAout << "\t\t- MET up sigmas: " << it2.second.sigmas_up << endl;
    }
  }
  MELAout << "****************************" << endl;
}

