#include <cassert>
#include "TRandom3.h"
#include "HostHelpersCore.h"
#include "SampleHelpersCore.h"
#include "SamplesCore.h"
#include "METCorrectionHandler.h"
#include <CMS3/Dictionaries/interface/CMS3StreamHelpers.h>


using namespace std;
using namespace SampleHelpers;
using namespace IvyStreamHelpers;


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
    pfmet_XY_JER_p4Preserved_map[period]=std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>();
    pfmet_JER_p4Preserved_map[period]=std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>();
    pfmet_XY_p4Preserved_map[period]=std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>();
    pfmet_p4Preserved_map[period]=std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>();
    puppimet_p4Preserved_map[period]=std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>();
    pfmet_XY_JER_map[period]=std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>();
    pfmet_JER_map[period]=std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>();
    pfmet_XY_map[period]=std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>();
    pfmet_map[period]=std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>();
    puppimet_map[period]=std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>();
    for (auto const& syst:allowedSysts){
      pfmet_XY_JER_p4Preserved_map[period][syst]=METCorrectionParameters();
      pfmet_JER_p4Preserved_map[period][syst]=METCorrectionParameters();
      pfmet_XY_p4Preserved_map[period][syst]=METCorrectionParameters();
      pfmet_p4Preserved_map[period][syst]=METCorrectionParameters();
      puppimet_p4Preserved_map[period][syst]=METCorrectionParameters();
      pfmet_XY_JER_map[period][syst]=METCorrectionParameters();
      pfmet_JER_map[period][syst]=METCorrectionParameters();
      pfmet_XY_map[period][syst]=METCorrectionParameters();
      pfmet_map[period][syst]=METCorrectionParameters();
      puppimet_map[period][syst]=METCorrectionParameters();
    }
    readFile(strinputcore + "fitparameters_pfmet_JEC_XY_JER_PartMomShifts_p4Preserved_abseta_lt_4p7_" + period + ".txt", pfmet_XY_JER_p4Preserved_map[period]);
    readFile(strinputcore + "fitparameters_pfmet_JEC_JER_PartMomShifts_p4Preserved_abseta_lt_4p7_" + period + ".txt", pfmet_JER_p4Preserved_map[period]);
    readFile(strinputcore + "fitparameters_pfmet_JEC_XY_PartMomShifts_p4Preserved_abseta_lt_4p7_" + period + ".txt", pfmet_XY_p4Preserved_map[period]);
    readFile(strinputcore + "fitparameters_pfmet_JEC_PartMomShifts_p4Preserved_abseta_lt_4p7_" + period + ".txt", pfmet_p4Preserved_map[period]);
    readFile(strinputcore + "fitparameters_puppimet_JEC_PartMomShifts_p4Preserved_abseta_lt_4p7_" + period + ".txt", puppimet_p4Preserved_map[period]);
    readFile(strinputcore + "fitparameters_pfmet_JEC_XY_JER_PartMomShifts_abseta_lt_4p7_" + period + ".txt", pfmet_XY_JER_map[period]);
    readFile(strinputcore + "fitparameters_pfmet_JEC_JER_PartMomShifts_abseta_lt_4p7_" + period + ".txt", pfmet_JER_map[period]);
    readFile(strinputcore + "fitparameters_pfmet_JEC_XY_PartMomShifts_abseta_lt_4p7_" + period + ".txt", pfmet_XY_map[period]);
    readFile(strinputcore + "fitparameters_pfmet_JEC_PartMomShifts_abseta_lt_4p7_" + period + ".txt", pfmet_map[period]);
    readFile(strinputcore + "fitparameters_puppimet_JEC_PartMomShifts_abseta_lt_4p7_" + period + ".txt", puppimet_map[period]);
  }

  // Test all cases
  for (unsigned short isPFMET=0; isPFMET<2; isPFMET++){
    for (unsigned short iXY=0; iXY<2; iXY++){
      for (unsigned short iJER=0; iJER<2; iJER++){
        for (unsigned short iP4Preserve=0; iP4Preserve<2; iP4Preserve++){
          std::unordered_map<
            TString,
            std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>
          > const* theMap = getCorrectionMap(bool(isPFMET), iXY, iJER, iP4Preserve);
          for (auto const& period:allowedDataPeriods){
            auto it = theMap->find(period);

            if (it == theMap->end()){
              IVYerr << "METCorrectionHandler::setup: "
                << (isPFMET==1 ? "pfmet_JEC" : "puppimet_JEC")
                << (iXY==1 ? "_XY" : "")
                << (iJER==1 ? "_JER" : "")
                << "_PartMomShifts"
                << (iP4Preserve==1 ? "_p4Preserved" : "")
                << " did not set up properly for period " << period << "."
                << endl;
              res = false;
              continue;
            }
            for (auto const& syst:allowedSysts){
              auto it_syst_data = it->second.find(SystematicsHelpers::nSystematicVariations);
              auto it_syst_MC = it->second.find(syst);
              if (it_syst_data == it->second.end()){
                IVYerr << "METCorrectionHandler::setup: "
                  << (isPFMET==1 ? "pfmet_JEC" : "puppimet_JEC")
                  << (iXY==1 ? "_XY" : "")
                  << (iJER==1 ? "_JER" : "")
                  << "_PartMomShifts"
                  << (iP4Preserve==1 ? "_p4Preserved" : "")
                  << " did not set up properly for period " << period << ". "
                  << "Data is missing."
                  << endl;
                res = false;
              }
              if (it_syst_MC == it->second.end()){
                IVYerr << "METCorrectionHandler::setup: "
                  << (isPFMET==1 ? "pfmet_JEC" : "puppimet_JEC")
                  << (iXY==1 ? "_XY" : "")
                  << (iJER==1 ? "_JER" : "")
                  << "_PartMomShifts"
                  << (iP4Preserve==1 ? "_p4Preserved" : "")
                  << " did not set up properly for period " << period << ". "
                  << "Systematic " << SystematicsHelpers::getSystName(syst) << " is missing for the MC."
                  << endl;
                res = false;
              }
            }
          }
        }
      }
    }
  }

  assert(res);

  return res;
}
void METCorrectionHandler::readFile(TString const& strinput, std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>& pars){
  if (!HostHelpers::FileReadable(strinput.Data())){
    IVYerr << "METCorrectionHandler::setup: File " << strinput << " is not readable." << endl;
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
          IVYerr << "METCorrectionHandler::readFile: Cannot determine the systematic for line " << str_in << endl;
          assert(0);
        }
        pars[currentSyst].sigmas_nominal.clear();
        pars[currentSyst].sigmas_dn.clear();
        pars[currentSyst].sigmas_up.clear();
        pars[currentSyst].fracs.clear();
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
  pfmet_XY_JER_p4Preserved_map.clear();
  pfmet_JER_p4Preserved_map.clear();
  pfmet_XY_p4Preserved_map.clear();
  pfmet_p4Preserved_map.clear();
  puppimet_p4Preserved_map.clear();
  pfmet_XY_JER_map.clear();
  pfmet_JER_map.clear();
  pfmet_XY_map.clear();
  pfmet_map.clear();
  puppimet_map.clear();
}

std::unordered_map<
  TString,
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>
> const* METCorrectionHandler::getCorrectionMap(bool const& isPFMET, unsigned short const& iXY, unsigned short const& iJER, unsigned short const& iP4Preserve) const{
  std::unordered_map<
    TString,
    std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>
  > const* res = nullptr;
  if (!isPFMET) res = (iP4Preserve==0 ? &puppimet_map : &puppimet_p4Preserved_map);
  else{
    if (iXY==0 && iJER==0) res = (iP4Preserve==0 ? &pfmet_map : &pfmet_p4Preserved_map);
    else if (iXY==0 && iJER==1) res = (iP4Preserve==0 ? &pfmet_JER_map : &pfmet_JER_p4Preserved_map);
    else if (iXY==1 && iJER==0) res = (iP4Preserve==0 ? &pfmet_XY_map : &pfmet_XY_p4Preserved_map);
    else res = (iP4Preserve==0 ? &pfmet_XY_JER_map : &pfmet_XY_JER_p4Preserved_map);
  }
  return res;
}

void METCorrectionHandler::applyCorrections(
  TString const& effDataPeriod,
  float const& genMET, float const& genMETPhi,
  METObject* obj, bool isPFMET,
  double const* inputRndNum
) const{
  static const std::vector<SystematicsHelpers::SystematicVariationTypes> allowedSysts{
    SystematicsHelpers::sNominal,
    SystematicsHelpers::eJECDn, SystematicsHelpers::eJECUp,
    SystematicsHelpers::eJERDn, SystematicsHelpers::eJERUp,
    SystematicsHelpers::ePUDn, SystematicsHelpers::ePUUp
  };
  SystematicsHelpers::SystematicVariationTypes const& syst = obj->getCurrentSystematic();
  SystematicsHelpers::SystematicVariationTypes effSyst = SystematicsHelpers::sNominal;
  if (HelperFunctions::checkListVariable(allowedSysts, syst)) effSyst = syst;

  double frac_x = -1;
  if (inputRndNum) frac_x = *inputRndNum;
  else{
    TRandom3 rand;
    rand.SetSeed(static_cast<unsigned long long>(std::abs(std::sin(genMETPhi)*100000.)));
    frac_x = rand.Uniform();
  }
  ParticleObject::LorentzVector_t const genmet_p4(genMET*std::cos(genMETPhi), genMET*std::sin(genMETPhi), 0, 0);

  for (unsigned short iXY=0; iXY<2; iXY++){
    for (unsigned short iJER=0; iJER<2; iJER++){
      for (unsigned short iP4Preserve=0; iP4Preserve<2; iP4Preserve++){
        std::unordered_map<
          TString,
          std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>
        > const* theMap = getCorrectionMap(isPFMET, iXY, iJER, iP4Preserve);

        auto it = theMap->find(effDataPeriod);
        if (it==theMap->cend()){
          IVYerr << "METCorrectionHandler::applyCorrections: The data period " << effDataPeriod << " for (XY, JER, p4Preserve) = (" << iXY << ", " << iJER << ", " << iP4Preserve << ") cannot be found!" << endl;
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

        float const SF_nominal = sigma_nominal_data / sigma_nominal_MC;
        float SF = SF_nominal;
        if (syst == SystematicsHelpers::eMETDn) SF = SF_nominal - std::sqrt(std::pow(sigma_dn_data / sigma_nominal_MC - SF_nominal, 2) + std::pow(sigma_nominal_data / sigma_up_MC - SF_nominal, 2));
        else if (syst == SystematicsHelpers::eMETUp) SF = SF_nominal + std::sqrt(std::pow(sigma_up_data / sigma_nominal_MC - SF_nominal, 2) + std::pow(sigma_nominal_data / sigma_dn_MC - SF_nominal, 2));

        for (unsigned short iPMS=0; iPMS<2; iPMS++){
          ParticleObject::LorentzVector_t met_p4_diff = obj->p4((bool) iXY, (bool) iJER, (bool) iPMS, (bool) iP4Preserve) - genmet_p4;
          met_p4_diff *= SF-1.f;
          obj->setMETCorrection(met_p4_diff, (bool) iXY, (bool) iJER, (bool) iPMS, (bool) iP4Preserve);
        }
      }
    }
  }
}

void METCorrectionHandler::printParameters(
  std::unordered_map<
  TString,
  std::unordered_map<SystematicsHelpers::SystematicVariationTypes, METCorrectionParameters>
  > const& met_map,
  TString const& mname
) const{
  IVYout << "****************************" << endl;
  IVYout << "METCorrectionHandler::printParameters: Parameters of " << mname << " corrections:" << endl;
  for (auto const& it:met_map){
    for (auto const& it2:it.second){
      IVYout << "\t- Period " << it.first << ", systematic " << (it2.first!=SystematicsHelpers::nSystematicVariations ? SystematicsHelpers::getSystName(it2.first) : "data") << ": " << endl;
      IVYout << "\t\t- Fractions: " << it2.second.fracs << endl;
      IVYout << "\t\t- Nominal sigmas: " << it2.second.sigmas_nominal << endl;
      IVYout << "\t\t- MET dn sigmas: " << it2.second.sigmas_dn << endl;
      IVYout << "\t\t- MET up sigmas: " << it2.second.sigmas_up << endl;
    }
  }
  IVYout << "****************************" << endl;
}


void METCorrectionHandler::printParameters() const{
  printParameters(pfmet_XY_JER_p4Preserved_map, "pfmet_XY_JER_p4Preserved");
  printParameters(pfmet_JER_p4Preserved_map, "pfmet_JER_p4Preserved");
  printParameters(pfmet_XY_p4Preserved_map, "pfmet_XY_p4Preserved");
  printParameters(pfmet_p4Preserved_map, "pfmet_p4Preserved");
  printParameters(puppimet_p4Preserved_map, "puppimet_p4Preserved");
  printParameters(pfmet_XY_JER_map, "pfmet_XY_JER");
  printParameters(pfmet_JER_map, "pfmet_JER");
  printParameters(pfmet_XY_map, "pfmet_XY");
  printParameters(pfmet_map, "pfmet");
  printParameters(puppimet_map, "puppimet");
}

