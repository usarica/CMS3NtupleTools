#include <cassert>
#include "TriggerScaleFactorHandler.h"
#include "MuonSelectionHelpers.h"
#include "ElectronSelectionHelpers.h"
#include "ParticleSelectionHelpers.h"
#include "SampleHelpersCore.h"
#include "SamplesCore.h"
#include "HelperFunctions.h"
#include "TDirectory.h"
#include "TFile.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace SampleHelpers;
using namespace MELAStreamHelpers;


TriggerScaleFactorHandler::TriggerScaleFactorHandler() : ScaleFactorHandlerBase()
{
  this->setup();
}

TriggerScaleFactorHandler::~TriggerScaleFactorHandler(){ this->reset(); }


bool TriggerScaleFactorHandler::setup(){
  using namespace SystematicsHelpers;

  bool res = true;
  this->reset();

  if (verbosity>=TVar::INFO) MELAout << "TriggerScaleFactorHandler::setup: Setting up efficiency and SF histograms for year " << SampleHelpers::getDataYear() << endl;

  TDirectory* curdir = gDirectory;
  TDirectory* uppermostdir = SampleHelpers::rootTDirectory;

  TString cinput_main = ANALYSISTREEPKGDATAPATH+Form("ScaleFactors/Trigger/%i/", SampleHelpers::getDataYear());
  HostHelpers::ExpandEnvironmentVariables(cinput_main);

  std::vector<SystematicVariationTypes> const allowedSysts={ sNominal, eTriggerEffDn, eTriggerEffUp, ePUDn, ePUUp };
  std::vector<SystematicVariationTypes> const allowedSysts_eff={ sNominal, ePUDn, ePUUp };

  {
    TString cinput = cinput_main + "trigger_efficiencies_leptons.root";
    if (!HostHelpers::FileReadable(cinput.Data())){
      MELAerr << "TriggerScaleFactorHandler::setup: File " << cinput << " is not readable." << endl;
      assert(0);
    }

    TFile* finput = TFile::Open(cinput, "read"); uppermostdir->cd();
    TDirectory* dir_Dilepton_Combined = (TDirectory*) finput->Get("Dilepton_Combined");
    TDirectory* dir_SingleLepton_Combined = (TDirectory*) finput->Get("SingleLepton_Combined");
    uppermostdir->cd();
    for (auto const& syst:allowedSysts){
      TString strSyst = SystematicsHelpers::getSystName(syst).data();
      TString systname;
      switch (syst){
      case eTriggerEffDn:
        systname = "dn";
        strSyst = SystematicsHelpers::getSystName(SystematicsHelpers::sNominal).data();
        break;
      case eTriggerEffUp:
        systname = "up";
        strSyst = SystematicsHelpers::getSystName(SystematicsHelpers::sNominal).data();
        break;
      default:
        systname = "nominal";
        break;
      }

      uppermostdir->cd();
      std::vector<ExtendedHistogram_2D> tmpvec(4, ExtendedHistogram_2D());
      syst_SF_Dilepton_SingleLepton_mumu_map[syst] = tmpvec;
      syst_SF_Dilepton_SingleLepton_ee_map[syst] = tmpvec;
      syst_SF_Dilepton_SingleLepton_mue_map[syst] = tmpvec;
      syst_SF_SingleMuon_map[syst] = tmpvec.front();
      syst_SF_SingleElectron_map[syst] = tmpvec.front();

      if (HelperFunctions::checkListVariable(allowedSysts_eff, syst)){
        syst_eff_mc_Dilepton_SingleLepton_mumu_map[syst] = tmpvec;
        syst_eff_mc_Dilepton_SingleLepton_ee_map[syst] = tmpvec;
        syst_eff_mc_Dilepton_SingleLepton_mue_map[syst] = tmpvec;
        syst_eff_mc_SingleMuon_map[syst] = tmpvec.front();
        syst_eff_mc_SingleElectron_map[syst] = tmpvec.front();
      }

      // Combined dilepton trigger efficiencies
      std::vector<TString> const benames{ "barrel", "endcap" };
      for (unsigned short ibe=0; ibe<2; ibe++){
        for (unsigned short jbe=0; jbe<2; jbe++){
          unsigned short const idx_hist = 2*ibe+jbe;
          TString hname;

          hname = Form("h_Combined_SF_pt25avg_%s_wcuts_mumu_%s_%s_%s", systname.Data(), benames.at(ibe).Data(), benames.at(jbe).Data(), strSyst.Data());
          res &= getHistogram<TH2F, ExtendedHistogram_2D>(syst_SF_Dilepton_SingleLepton_mumu_map[syst].at(idx_hist), dir_Dilepton_Combined, hname);
          uppermostdir->cd();
          hname = Form("h_Combined_SF_pt25avg_%s_wcuts_ee_%s_%s_%s", systname.Data(), benames.at(ibe).Data(), benames.at(jbe).Data(), strSyst.Data());
          res &= getHistogram<TH2F, ExtendedHistogram_2D>(syst_SF_Dilepton_SingleLepton_ee_map[syst].at(idx_hist), dir_Dilepton_Combined, hname);
          uppermostdir->cd();
          hname = Form("h_Combined_SF_pt25avg_%s_wcuts_mue_%s_%s_%s", systname.Data(), benames.at(ibe).Data(), benames.at(jbe).Data(), strSyst.Data());
          res &= getHistogram<TH2F, ExtendedHistogram_2D>(syst_SF_Dilepton_SingleLepton_mue_map[syst].at(idx_hist), dir_Dilepton_Combined, hname);
          uppermostdir->cd();
          if (HelperFunctions::checkListVariable(allowedSysts_eff, syst)){
            hname = Form("h_Combined_eff_%s_wcuts_mumu_%s_%s_MC_%s", systname.Data(), benames.at(ibe).Data(), benames.at(jbe).Data(), strSyst.Data());
            res &= getHistogram<TH2F, ExtendedHistogram_2D>(syst_eff_mc_Dilepton_SingleLepton_mumu_map[syst].at(idx_hist), dir_Dilepton_Combined, hname);
            uppermostdir->cd();
            hname = Form("h_Combined_eff_%s_wcuts_ee_%s_%s_MC_%s", systname.Data(), benames.at(ibe).Data(), benames.at(jbe).Data(), strSyst.Data());
            res &= getHistogram<TH2F, ExtendedHistogram_2D>(syst_eff_mc_Dilepton_SingleLepton_ee_map[syst].at(idx_hist), dir_Dilepton_Combined, hname);
            uppermostdir->cd();
            hname = Form("h_Combined_eff_%s_wcuts_mue_%s_%s_MC_%s", systname.Data(), benames.at(ibe).Data(), benames.at(jbe).Data(), strSyst.Data());
            res &= getHistogram<TH2F, ExtendedHistogram_2D>(syst_eff_mc_Dilepton_SingleLepton_mue_map[syst].at(idx_hist), dir_Dilepton_Combined, hname);
            uppermostdir->cd();
          }
        }
      }

      {
        TString hname;

        hname = Form("h_SingleMuon_SF_%s_%s", systname.Data(), strSyst.Data());
        res &= getHistogram<TH2F, ExtendedHistogram_2D>(syst_SF_SingleMuon_map[syst], dir_SingleLepton_Combined, hname);
        uppermostdir->cd();
        hname = Form("h_SingleElectron_SF_%s_%s", systname.Data(), strSyst.Data());
        res &= getHistogram<TH2F, ExtendedHistogram_2D>(syst_SF_SingleElectron_map[syst], dir_SingleLepton_Combined, hname);
        uppermostdir->cd();
        if (HelperFunctions::checkListVariable(allowedSysts_eff, syst)){
          hname = Form("h_SingleMuon_eff_%s_MC_%s", systname.Data(), strSyst.Data());
          res &= getHistogram<TH2F, ExtendedHistogram_2D>(syst_eff_mc_SingleMuon_map[syst], dir_SingleLepton_Combined, hname);
          uppermostdir->cd();
          hname = Form("h_SingleElectron_eff_%s_MC_%s", systname.Data(), strSyst.Data());
          res &= getHistogram<TH2F, ExtendedHistogram_2D>(syst_eff_mc_SingleElectron_map[syst], dir_SingleLepton_Combined, hname);
          uppermostdir->cd();
        }
      }

    }
    ScaleFactorHandlerBase::closeFile(finput); curdir->cd();
  }

  return res;
}
void TriggerScaleFactorHandler::reset(){
  syst_eff_mc_Dilepton_SingleLepton_mumu_map.clear();
  syst_eff_mc_Dilepton_SingleLepton_ee_map.clear();
  syst_eff_mc_Dilepton_SingleLepton_mue_map.clear();
  syst_eff_mc_SingleMuon_map.clear();
  syst_eff_mc_SingleElectron_map.clear();

  syst_SF_Dilepton_SingleLepton_mumu_map.clear();
  syst_SF_Dilepton_SingleLepton_ee_map.clear();
  syst_SF_Dilepton_SingleLepton_mue_map.clear();
  syst_SF_SingleMuon_map.clear();
  syst_SF_SingleElectron_map.clear();
}

void TriggerScaleFactorHandler::evalScaleFactorFromHistogram_PtEta(float& theSF, float const& pt, float const& eta, ExtendedHistogram_2D const& hist, bool etaOnY, bool useAbsEta) const{
  TH2F const* hh = hist.getHistogram();
  if (!hh) return;

  float const abs_eta = std::abs(eta);
  float const& eta_to_use = (!useAbsEta ? eta : abs_eta);

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

  if (ix==0) ix=1;
  else if (ix==nbinsx+1) ix=nbinsx;
  if (iy==0) iy=1;
  else if (iy==nbinsy+1) iy=nbinsy;

  theSF *= hh->GetBinContent(ix, iy);
}
void TriggerScaleFactorHandler::evalScaleFactorFromHistogram_PtPt(float& theSF, float const& pt1, float const& pt2, ExtendedHistogram_2D const& hist) const{
  TH2F const* hh = hist.getHistogram();
  if (!hh) MELAerr << "Histogram is NULL!" << endl;
  if (!hh) return;

  int ix, iy;
  int nbinsx = hh->GetNbinsX();
  int nbinsy = hh->GetNbinsY();
  ix = hh->GetXaxis()->FindBin(pt1);
  iy = hh->GetYaxis()->FindBin(pt2);

  if (ix==0) ix=1;
  else if (ix==nbinsx+1) ix=nbinsx;
  if (iy==0) iy=1;
  else if (iy==nbinsy+1) iy=nbinsy;

  theSF *= hh->GetBinContent(ix, iy);
}


void TriggerScaleFactorHandler::getCombinedDileptonSFAndEff(
  SystematicsHelpers::SystematicVariationTypes const& syst,
  float pt1, float eta1, cms3_id_t id1,
  float pt2, float eta2, cms3_id_t id2,
  bool passTrigger,
  float& val, float* effval
) const{
  using namespace SystematicsHelpers;

  if (verbosity>=TVar::DEBUG) MELAout
    << "TriggerScaleFactorHandler::getCombinedDileptonSFAndEff: Evaluating " << (effval ? "SFs and efficiencies" : "SFs")
    << " for pT1, pT2 = " << pt1 << ", " << pt2
    << "; eta1, eta2 = " << eta1 << ", " << eta2
    << "; id1, id2 = " << id1 << ", " << id2
    << "; passTrigger ?= " << passTrigger
    << endl;

  val = 1;
  if (effval) *effval = 1;

  std::vector<SystematicVariationTypes> const allowedSysts={ sNominal, eTriggerEffDn, eTriggerEffUp, ePUDn, ePUUp };
  std::vector<SystematicVariationTypes> const allowedSysts_eff={ sNominal, ePUDn, ePUUp };

  SystematicVariationTypes activeSyst_eff_nominal = sNominal;
  if (HelperFunctions::checkListVariable(allowedSysts_eff, syst)) activeSyst_eff_nominal = syst;

  SystematicVariationTypes activeSyst = sNominal;
  if (HelperFunctions::checkListVariable(allowedSysts, syst)) activeSyst = syst;
  if (verbosity>=TVar::DEBUG) MELAout << "\t- Active systematics: " << activeSyst << " / " << activeSyst_eff_nominal << endl;

  // Order objects 1, 2
  cms3_id_t dilepton_id = id1*id2;
  bool is_mue = false;
  bool is_mumu = false;
  if (std::abs(dilepton_id)==143){
    is_mue = true;
    // DF case: Swap only if the order is e-mu
    if (std::abs(id1)==11){
      std::swap(id1, id2);
      std::swap(pt1, pt2);
      std::swap(eta1, eta2);
    }
  }
  else if (pt1<pt2){
    is_mumu = (std::abs(dilepton_id)==169);
    // SF case: Swap only if pt1<pt2
    std::swap(id1, id2);
    std::swap(pt1, pt2);
    std::swap(eta1, eta2);
  }

  // Find B/E index
  bool const isBarrel1 = std::abs(eta1)<(std::abs(id1)==13 ? 1.2 : 1.479);
  bool const isBarrel2 = std::abs(eta2)<(std::abs(id2)==13 ? 1.2 : 1.479);
  unsigned short idx_BE = (isBarrel1 ? 0 : 2) + (isBarrel2 ? 0 : 1);

  // Obtain histograms
  ExtendedHistogram_2D const* h_eff_mc = nullptr;
  ExtendedHistogram_2D const* h_SF = nullptr;
  if (is_mue){
    h_eff_mc = &(syst_eff_mc_Dilepton_SingleLepton_mue_map.find(activeSyst_eff_nominal)->second.at(idx_BE));
    h_SF = &(syst_SF_Dilepton_SingleLepton_mue_map.find(activeSyst)->second.at(idx_BE));
  }
  else if (is_mumu){
    h_eff_mc = &(syst_eff_mc_Dilepton_SingleLepton_mumu_map.find(activeSyst_eff_nominal)->second.at(idx_BE));
    h_SF = &(syst_SF_Dilepton_SingleLepton_mumu_map.find(activeSyst)->second.at(idx_BE));
  }
  else{
    h_eff_mc = &(syst_eff_mc_Dilepton_SingleLepton_ee_map.find(activeSyst_eff_nominal)->second.at(idx_BE));
    h_SF = &(syst_SF_Dilepton_SingleLepton_ee_map.find(activeSyst)->second.at(idx_BE));
  }

  float eff_nominal_unscaled=1, eff_scaled=1;
  evalScaleFactorFromHistogram_PtPt(eff_nominal_unscaled, pt1, pt2, *h_eff_mc);

  float SF_val = 1;
  evalScaleFactorFromHistogram_PtPt(SF_val, pt1, pt2, *h_SF);
  eff_scaled = std::min(1.f, std::max(0.f, eff_nominal_unscaled * SF_val));
  if (passTrigger){
    val = SF_val;
    if (effval) *effval = eff_scaled;
  }
  else{
    eff_scaled = 1.f - eff_scaled;
    eff_nominal_unscaled = 1.f - eff_nominal_unscaled;

    if (eff_nominal_unscaled>0.f){
      val = eff_scaled / eff_nominal_unscaled;
      if (effval) *effval = eff_scaled;
    }
    else{
      val = 0;
      if (effval) *effval = 0;
    }
  }
}
void TriggerScaleFactorHandler::getCombinedSingleLeptonSFAndEff(
  SystematicsHelpers::SystematicVariationTypes const& syst,
  float const& pt, float const& eta, cms3_id_t const& partId,
  bool passTrigger,
  float& val, float* effval
) const{
  using namespace SystematicsHelpers;

  if (verbosity>=TVar::DEBUG) MELAout
    << "TriggerScaleFactorHandler::getCombinedSingleLeptonSFAndEff: Evaluating " << (effval ? "SFs and efficiencies" : "SFs")
    << " for pT, eta, id = " << pt << ", " << eta << ", " << partId
    << ", passTrigger ?= " << passTrigger
    << endl;

  val = 1;
  if (effval) *effval = 1;

  std::vector<SystematicVariationTypes> const allowedSysts={ sNominal, eTriggerEffDn, eTriggerEffUp, ePUDn, ePUUp };
  std::vector<SystematicVariationTypes> const allowedSysts_eff={ sNominal, ePUDn, ePUUp };

  SystematicVariationTypes activeSyst_eff_nominal = sNominal;
  if (HelperFunctions::checkListVariable(allowedSysts_eff, syst)) activeSyst_eff_nominal = syst;

  SystematicVariationTypes activeSyst = sNominal;
  if (HelperFunctions::checkListVariable(allowedSysts, syst)) activeSyst = syst;
  if (verbosity>=TVar::DEBUG) MELAout << "\t- Active systematics: " << activeSyst << " / " << activeSyst_eff_nominal << endl;

  bool is_mu = std::abs(partId)==13;

  // Obtain histograms
  ExtendedHistogram_2D const* h_eff_mc = nullptr;
  ExtendedHistogram_2D const* h_SF = nullptr;
  if (is_mu){
    h_eff_mc = &(syst_eff_mc_SingleMuon_map.find(activeSyst_eff_nominal)->second);
    h_SF = &(syst_SF_SingleMuon_map.find(activeSyst)->second);
  }
  else{
    h_eff_mc = &(syst_eff_mc_SingleElectron_map.find(activeSyst_eff_nominal)->second);
    h_SF = &(syst_SF_SingleElectron_map.find(activeSyst)->second);
  }

  float eff_nominal_unscaled=1, eff_scaled=1;
  evalScaleFactorFromHistogram_PtEta(eff_nominal_unscaled, pt, eta, *h_eff_mc, false, true);

  float SF_val = 1;
  evalScaleFactorFromHistogram_PtEta(SF_val, pt, eta, *h_SF, false, true);
  eff_scaled = std::min(1.f, std::max(0.f, eff_nominal_unscaled * SF_val));
  if (passTrigger){
    val = SF_val;
    if (effval) *effval = eff_scaled;
  }
  else{
    eff_scaled = 1.f - eff_scaled;
    eff_nominal_unscaled = 1.f - eff_nominal_unscaled;

    if (eff_nominal_unscaled>0.f){
      val = eff_scaled / eff_nominal_unscaled;
      if (effval) *effval = eff_scaled;
    }
    else{
      val = 0;
      if (effval) *effval = 0;
    }
  }
}

void TriggerScaleFactorHandler::getCombinedDileptonSFAndEff(
  SystematicsHelpers::SystematicVariationTypes const& syst,
  ParticleObject const* obj1, ParticleObject const* obj2,
  bool passTrigger,
  float& val, float* effval
) const{
  val = 1;
  if (effval) *effval = 1;

  if (!obj1 || !ParticleSelectionHelpers::isParticleForTriggerChecking(obj1)) return;
  if (!obj2 || !ParticleSelectionHelpers::isParticleForTriggerChecking(obj2)) return;

  float pt1 = obj1->pt();
  float pt2 = obj2->pt();
  cms3_id_t id1 = obj1->pdgId();
  cms3_id_t id2 = obj2->pdgId();
  float eta1, eta2;
  if (std::abs(id1)==11) eta1 = dynamic_cast<ElectronObject const*>(obj1)->etaSC();
  else eta1 = obj1->eta();
  if (std::abs(id2)==11) eta2 = dynamic_cast<ElectronObject const*>(obj2)->etaSC();
  else eta2 = obj2->eta();
  getCombinedDileptonSFAndEff(
    syst,
    pt1, eta1, id1,
    pt2, eta2, id2,
    passTrigger,
    val, effval
  );
}
void TriggerScaleFactorHandler::getCombinedSingleLeptonSFAndEff(
  SystematicsHelpers::SystematicVariationTypes const& syst,
  ParticleObject const* obj,
  bool passTrigger,
  float& val, float* effval
) const{
  val = 1;
  if (effval) *effval = 1;

  if (!obj || !ParticleSelectionHelpers::isParticleForTriggerChecking(obj)) return;

  float pt = obj->pt();
  cms3_id_t partId = obj->pdgId();
  float eta;
  if (std::abs(partId)==11) eta = dynamic_cast<ElectronObject const*>(obj)->etaSC();
  else eta = obj->eta();
  getCombinedSingleLeptonSFAndEff(
    syst,
    pt, eta, partId,
    passTrigger,
    val, effval
  );
}
