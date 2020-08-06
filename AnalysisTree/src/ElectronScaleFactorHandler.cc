#include "ElectronScaleFactorHandler.h"
#include "SamplesCore.h"
#include "HelperFunctions.h"
#include "TDirectory.h"
#include "TFile.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace SampleHelpers;
using namespace MELAStreamHelpers;


ElectronScaleFactorHandler::ElectronScaleFactorHandler() : ScaleFactorHandlerBase()
{
  this->setup();
}

ElectronScaleFactorHandler::~ElectronScaleFactorHandler(){ this->reset(); }


bool ElectronScaleFactorHandler::setup(){
  using namespace SystematicsHelpers;

  bool res = true;
  this->reset();

  TDirectory* curdir = gDirectory;

  TString cinput_main = ANALYSISTREEPKGDATAPATH+Form("ScaleFactors/Electrons/%i/", SampleHelpers::theDataYear);
  HostHelpers::ExpandEnvironmentVariables(cinput_main);
  MELAout << "ElectronScaleFactorHandler::setup: Setting up efficiency and SF histograms for year " << SampleHelpers::theDataYear << endl;

  std::vector<SystematicVariationTypes> const allowedSysts={
    sNominal,
    eEleEffStatDn, eEleEffStatUp,
    eEleEffSystDn, eEleEffSystUp,
    eEleEffAltMCDn, eEleEffAltMCUp,
    ePUDn, ePUUp
  };
  constexpr unsigned int n_non_gap_gap = 2;
  for (auto const& syst:allowedSysts){
    std::vector<ExtendedHistogram_2D> tmpvec(n_non_gap_gap, ExtendedHistogram_2D());

    if (syst == sNominal){
      eff_mc_reco_hists = tmpvec;
      eff_mc_id_hists = tmpvec;
      eff_mc_iso_loose_hists = tmpvec;
      eff_mc_iso_tight_hists = tmpvec;
    }

    syst_SF_reco_map[syst] = tmpvec;
    syst_SF_id_map[syst] = tmpvec;
    syst_SF_iso_loose_map[syst] = tmpvec;
    syst_SF_iso_tight_map[syst] = tmpvec;
  }

  {
    TString cinput = cinput_main + "Efficiencies_ee_nongap_gap_tracking.root";
    if (!HostHelpers::FileReadable(cinput.Data())){
      MELAerr << "ElectronScaleFactorHandler::setup: File " << cinput << " is not readable." << endl;
      assert(0);
    }
    TFile* finput = TFile::Open(cinput, "read"); curdir->cd();
    for (unsigned int igap=0; igap<n_non_gap_gap; igap++){
      res &= getHistogram<TH2F, ExtendedHistogram_2D>(syst_SF_id_map[sNominal].at(igap), finput, "EGamma_SF2D");
      res &= getHistogram<TH2F, ExtendedHistogram_2D>(eff_mc_reco_hists.at(igap), finput, "EGamma_EffMC2D");
    }
    ScaleFactorHandlerBase::closeFile(finput); curdir->cd();
  }

  for (unsigned int igap=0; igap<n_non_gap_gap; igap++){
    TString cinput = cinput_main + Form("Efficiencies_ee_%s_id_looseIso_tightIso.root", (igap==0 ? "nongap" : "gap"));
    if (!HostHelpers::FileReadable(cinput.Data())){
      MELAerr << "ElectronScaleFactorHandler::setup: File " << cinput << " is not readable." << endl;
      assert(0);
    }

    TFile* finput = TFile::Open(cinput, "read"); curdir->cd();
    for (auto const& syst:allowedSysts){
      TString systname;
      switch (syst){
      case eEleEffStatDn:
        systname = "StatDn";
        break;
      case eEleEffStatUp:
        systname = "StatUp";
        break;
      case eEleEffSystDn:
        systname = "SystDn";
        break;
      case eEleEffSystUp:
        systname = "SystUp";
        break;
      case eEleEffAltMCDn:
        systname = "AltMCDn";
        break;
      case eEleEffAltMCUp:
        systname = "AltMCUp";
        break;
      case ePUDn:
        systname = "PUDn";
        break;
      case ePUUp:
        systname = "PUUp";
        break;
      default:
        systname = "Nominal";
        break;
      }

      TString str_SF_id = Form("SF_%s_passId", systname.Data());
      TString str_SF_iso_loose = Form("SF_%s_passId_passLooseIso", systname.Data());
      TString str_SF_iso_tight = Form("SF_%s_passId_passTightIso", systname.Data());
      res &= getHistogram<TH2D, ExtendedHistogram_2D>(syst_SF_id_map[syst].at(igap), finput, str_SF_id);
      res &= getHistogram<TH2D, ExtendedHistogram_2D>(syst_SF_iso_loose_map[syst].at(igap), finput, str_SF_iso_loose);
      res &= getHistogram<TH2D, ExtendedHistogram_2D>(syst_SF_iso_tight_map[syst].at(igap), finput, str_SF_iso_tight);

      if (syst == sNominal){
        TString str_eff_mc_id = Form("eff_MC_%s_passId", systname.Data());
        TString str_eff_mc_iso_loose = Form("eff_MC_%s_passId_passLooseIso", systname.Data());
        TString str_eff_mc_iso_tight = Form("eff_MC_%s_passId_passTightIso", systname.Data());
        res &= getHistogram<TH2D, ExtendedHistogram_2D>(eff_mc_id_hists.at(igap), finput, str_eff_mc_id);
        res &= getHistogram<TH2D, ExtendedHistogram_2D>(eff_mc_iso_loose_hists.at(igap), finput, str_eff_mc_iso_loose);
        res &= getHistogram<TH2D, ExtendedHistogram_2D>(eff_mc_iso_tight_hists.at(igap), finput, str_eff_mc_iso_tight);
      }
    }
    ScaleFactorHandlerBase::closeFile(finput); curdir->cd();
  }

  return res;
}
void ElectronScaleFactorHandler::reset(){
  eff_mc_reco_hists.clear();
  eff_mc_id_hists.clear();
  eff_mc_iso_loose_hists.clear();
  eff_mc_iso_tight_hists.clear();

  syst_SF_reco_map.clear();
  syst_SF_id_map.clear();
  syst_SF_iso_loose_map.clear();
  syst_SF_iso_tight_map.clear();
}

void ElectronScaleFactorHandler::evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, float const& pt, float const& etaSC, ExtendedHistogram_2D const& hist, bool etaOnY, bool useAbsEta) const{
  TH2F const* hh = hist.getHistogram();
  if (!hh) return;

  float const abs_eta = std::abs(etaSC);
  float const& eta_to_use = (!useAbsEta ? etaSC : abs_eta);

  int ix, iy;
  int nbinsx = hh->GetNbinsX();
  int nbinsy = hh->GetNbinsY();
  if (!etaOnY){
    ix = hh->GetXaxis()->FindBin(eta_to_use);
    iy = hh->GetYaxis()->FindBin(pt);
  }
  else{
    ix = hh->GetXaxis()->FindBin(pt);
    iy = hh->GetYaxis()->FindBin(eta_to_use);
  }

  bool out_of_bounds = false;
  if (ix==0){ ix=1; out_of_bounds=true; }
  else if (ix==nbinsx+1){ ix=nbinsx; out_of_bounds=true; }
  if (iy==0){ iy=1; out_of_bounds=true; }
  else if (iy==nbinsy+1){ iy=nbinsy; out_of_bounds=true; }

  float bc = hh->GetBinContent(ix, iy);
  float be = hh->GetBinError(ix, iy);
  if (bc!=0.f) be /= bc;
  if (be<0.f) be=0;

  if (out_of_bounds) be *= 1.5;

  theSF *= bc; theSFRelErr = sqrt(pow(theSFRelErr, 2)+pow(be, 2));
}

void ElectronScaleFactorHandler::getIdIsoSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, float const& pt, float const& etaSC, bool const& is_gap, ElectronScaleFactorHandler::EfficiencyType type, float& val, float* effval) const{
  using namespace SystematicsHelpers;

  val = 1;
  if (effval) *effval = 1;

  std::vector<SystematicVariationTypes> const allowedSysts={
    sNominal,
    eEleEffStatDn, eEleEffStatUp,
    eEleEffSystDn, eEleEffSystUp,
    eEleEffAltMCDn, eEleEffAltMCUp,
    ePUDn, ePUUp
  };

  std::vector<SystematicVariationTypes> activeSysts={ sNominal };
  if (syst == eEleEffDn) activeSysts = std::vector<SystematicVariationTypes>{ eEleEffStatDn, eEleEffSystDn, eEleEffAltMCDn };
  else if (syst == eEleEffUp) activeSysts = std::vector<SystematicVariationTypes>{ eEleEffStatUp, eEleEffSystUp, eEleEffAltMCUp };
  else if (HelperFunctions::checkListVariable(allowedSysts, syst)) activeSysts = std::vector<SystematicVariationTypes>{ syst };

  std::vector<ExtendedHistogram_2D const*> hlist_eff_mc; hlist_eff_mc.reserve(kAllEffs);
  std::vector<ExtendedHistogram_2D const*> hlist_SF_nominal; hlist_SF_nominal.reserve(kAllEffs);
  {
    auto it_syst_SF_reco_map = syst_SF_reco_map.find(sNominal);
    auto it_syst_SF_id_map = syst_SF_id_map.find(sNominal);
    auto it_syst_SF_iso_loose_map = syst_SF_iso_loose_map.find(sNominal);
    auto it_syst_SF_iso_tight_map = syst_SF_iso_tight_map.find(sNominal);
    if ((type == kAllEffs || type == kTrackingEff) && it_syst_SF_reco_map!=syst_SF_reco_map.cend()){
      hlist_eff_mc.push_back(&(eff_mc_reco_hists.at((unsigned int) is_gap)));
      hlist_SF_nominal.push_back(&(it_syst_SF_reco_map->second.at((unsigned int) is_gap)));
    }
    else{
      hlist_eff_mc.push_back(nullptr);
      hlist_SF_nominal.push_back(nullptr);
    }
    if ((type == kAllEffs || type == kIdEff) && it_syst_SF_id_map!=syst_SF_id_map.cend()){
      hlist_eff_mc.push_back(&(eff_mc_id_hists.at((unsigned int) is_gap)));
      hlist_SF_nominal.push_back(&(it_syst_SF_id_map->second.at((unsigned int) is_gap)));
    }
    else{
      hlist_eff_mc.push_back(nullptr);
      hlist_SF_nominal.push_back(nullptr);
    }
    if ((type == kAllEffs || type == kLooseIsoEff) && it_syst_SF_iso_loose_map!=syst_SF_iso_loose_map.cend()){
      hlist_eff_mc.push_back(&(eff_mc_iso_loose_hists.at((unsigned int) is_gap)));
      hlist_SF_nominal.push_back(&(it_syst_SF_iso_loose_map->second.at((unsigned int) is_gap)));
    }
    else{
      hlist_eff_mc.push_back(nullptr);
      hlist_SF_nominal.push_back(nullptr);
    }
    if ((type == kAllEffs || type == kTightIsoEff) && it_syst_SF_iso_tight_map!=syst_SF_iso_tight_map.cend()){
      hlist_eff_mc.push_back(&(eff_mc_iso_tight_hists.at((unsigned int) is_gap)));
      hlist_SF_nominal.push_back(&(it_syst_SF_iso_tight_map->second.at((unsigned int) is_gap)));
    }
    else{
      hlist_eff_mc.push_back(nullptr);
      hlist_SF_nominal.push_back(nullptr);
    }
  }
  std::vector<float> SF_nominal_list; SF_nominal_list.reserve(kAllEffs);
  {
    auto it_SF = hlist_SF_nominal.begin();
    auto it_eff_mc = hlist_eff_mc.begin();
    while (it_SF != hlist_SF_nominal.end()){
      float SF_val = 1;
      if (*it_SF){
        float SF_err = 0;
        evalScaleFactorFromHistogram(SF_val, SF_err, pt, etaSC, **it_SF, false, false);
        if (effval){
          float eff_err = 0;
          evalScaleFactorFromHistogram(*effval, eff_err, pt, etaSC, **it_eff_mc, false, false);
        }
      }

      SF_nominal_list.push_back(SF_val);

      it_SF++;
      it_eff_mc++;
    }
  }

  float SF_err_val = 0;
  float SF_nominal_val = 1;
  for (auto const& v:SF_nominal_list) SF_nominal_val *= v;
  if (!(activeSysts.size() == 1 && activeSysts.front() == sNominal)){
    std::vector< std::vector<float> > systs_SF_list; systs_SF_list.reserve(activeSysts.size());
    for (auto const& asyst:activeSysts){
      auto it_syst_SF_reco_map = syst_SF_reco_map.find(sNominal/*asyst*/); // Look for sNominal only (FIXME after generating own reco SFs)
      auto it_syst_SF_id_map = syst_SF_id_map.find(asyst);
      auto it_syst_SF_iso_loose_map = syst_SF_iso_loose_map.find(asyst);
      auto it_syst_SF_iso_tight_map = syst_SF_iso_tight_map.find(asyst);

      std::vector<ExtendedHistogram_2D const*> hlist_SF; hlist_SF.reserve(kAllEffs);
      if ((type == kAllEffs || type == kTrackingEff) && it_syst_SF_reco_map!=syst_SF_reco_map.cend()) hlist_SF.push_back(&(it_syst_SF_reco_map->second.at((unsigned int) is_gap)));
      else hlist_SF.push_back(nullptr);
      if ((type == kAllEffs || type == kIdEff) && it_syst_SF_id_map!=syst_SF_id_map.cend()) hlist_SF.push_back(&(it_syst_SF_id_map->second.at((unsigned int) is_gap)));
      else hlist_SF.push_back(nullptr);
      if ((type == kAllEffs || type == kLooseIsoEff) && it_syst_SF_iso_loose_map!=syst_SF_iso_loose_map.cend()) hlist_SF.push_back(&(it_syst_SF_iso_loose_map->second.at((unsigned int) is_gap)));
      else hlist_SF.push_back(nullptr);
      if ((type == kAllEffs || type == kTightIsoEff) && it_syst_SF_iso_tight_map!=syst_SF_iso_tight_map.cend()) hlist_SF.push_back(&(it_syst_SF_iso_tight_map->second.at((unsigned int) is_gap)));
      else hlist_SF.push_back(nullptr);

      systs_SF_list.push_back(std::vector<float>(hlist_SF.size(), 1));
      std::vector<float>& current_SF_list = systs_SF_list.back();

      auto it_SF = hlist_SF.begin();
      auto it_SF_val = current_SF_list.begin();
      while (it_SF != hlist_SF.end()){
        if (*it_SF){
          float val_err = 0;
          evalScaleFactorFromHistogram(*it_SF_val, val_err, pt, etaSC, **it_SF, false, false);

          // Old-style histograms read errors from histogram bin error
          if ((type == kAllEffs || type == kTrackingEff) && it_SF == hlist_SF.begin()){ // Special treatment for tracking efficiency (FIXME: Need new format at some point)
            if (asyst == eEleEffSystDn) *it_SF_val = std::max(0.f, *it_SF_val - val_err);
            else if (asyst == eEleEffSystUp) *it_SF_val = std::max(0.f, *it_SF_val + val_err);
          }
        }
        else{
          *it_SF_val = SF_nominal_list.at(it_SF - hlist_SF.begin());
        }
        it_SF++;
        it_SF_val++;
      }
    }
    for (auto const& syst_SF_list:systs_SF_list){
      float SF_syst=1;
      for (auto const& v:syst_SF_list) SF_syst *= v;
      if (activeSysts.size() == 1) SF_err_val = SF_syst - SF_nominal_val;
      else{
        SF_err_val = std::sqrt(std::pow(SF_err_val, 2) + std::pow(SF_syst - SF_nominal_val, 2));
        if (syst == eEleEffDn || syst == eEleEffStatDn || syst == eEleEffSystDn || syst == eEleEffAltMCDn) SF_err_val *= -1.;
      }
    }
  }

  val = SF_nominal_val + SF_err_val;
  if (effval) *effval = std::min(1.f, (*effval)*val);
}
void ElectronScaleFactorHandler::getIdIsoSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, ElectronObject const* obj, ElectronScaleFactorHandler::EfficiencyType type, float& val, float* effval) const{
  val = 1;
  if (effval) *effval = 0;

  if (!obj) return;

  getIdIsoSFAndEff(syst, obj->pt(), obj->etaSC(), obj->isAnyGap(), type, val, effval);
}
