#include <cassert>
#include "PUJetIdScaleFactorHandler.h"
#include "AK4JetSelectionHelpers.h"
#include "SamplesCore.h"
#include "HelperFunctions.h"
#include "TDirectory.h"
#include "TFile.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace SampleHelpers;
using namespace MELAStreamHelpers;


PUJetIdScaleFactorHandler::PUJetIdScaleFactorHandler() : ScaleFactorHandlerBase()
{
  this->setup();
}

PUJetIdScaleFactorHandler::~PUJetIdScaleFactorHandler(){ this->reset(); }


bool PUJetIdScaleFactorHandler::setup(){
  using namespace SystematicsHelpers;

  bool res = true;
  this->reset();

  TDirectory* curdir = gDirectory;

  TString cinput_main = ANALYSISTREEPKGDATAPATH+Form("ScaleFactors/PUJetId/%i/", SampleHelpers::theDataYear);
  HostHelpers::ExpandEnvironmentVariables(cinput_main);
  MELAout << "PUJetIdScaleFactorHandler::setup: Setting up efficiency and SF histograms for year " << SampleHelpers::getDataYear() << endl;

  std::vector<TString> pujetidwpnames{
    "PUJetId_LooseJets",
    "PUJetId_MediumJets",
    "PUJetId_TightJets"
  };
  std::vector<TString> pujetidwpsfnames{
    "h2_<EFFMISTAG>_sf<YEAR>_L",
    "h2_<EFFMISTAG>_sf<YEAR>_M",
    "h2_<EFFMISTAG>_sf<YEAR>_T"
  };
  std::vector<TString> strmatches{ "mistagged", "matched" };
  std::vector<SystematicVariationTypes> const allowedSysts={
    sNominal,
    eJECDn, eJECUp,
    eJERDn, eJERUp,
    ePUDn, ePUUp
  };

  for (auto const& syst:allowedSysts){
    std::vector<ExtendedHistogram_2D> tmplist(pujetidwpnames.size(), ExtendedHistogram_2D());

    syst_pujetidwp_effs_map_mistagged[syst] = tmplist;
    syst_pujetidwp_effs_map_matched[syst] = tmplist;
    pujetidwp_SFs_map_mistagged = tmplist;
    pujetidwp_SFs_map_matched = tmplist;
  }

  // Get efficiencies
  {
    TString cinput = cinput_main + "Final_PUJetId_Efficiencies_AllMC.root";
    if (!HostHelpers::FileReadable(cinput.Data())){
      MELAerr << "PUJetIdScaleFactorHandler::setup: File " << cinput << " is not readable." << endl;
      assert(0);
    }

    TFile* finput = TFile::Open(cinput, "read"); curdir->cd();
    for (auto const& syst:allowedSysts){
      TString systname = SystematicsHelpers::getSystName(syst).data();
      unsigned short iwp=0;
      for (auto const& pujetidwpname:pujetidwpnames){
        TString hname;
        hname = Form("%s_%s_%s", pujetidwpname.Data(), strmatches.front().Data(), systname.Data());
        res &= getHistogram<TH2F, ExtendedHistogram_2D>(syst_pujetidwp_effs_map_mistagged[syst].at(iwp), finput, hname);
        hname = Form("%s_%s_%s", pujetidwpname.Data(), strmatches.back().Data(), systname.Data());
        res &= getHistogram<TH2F, ExtendedHistogram_2D>(syst_pujetidwp_effs_map_matched[syst].at(iwp), finput, hname);
        iwp++;
      }
    }
    ScaleFactorHandlerBase::closeFile(finput); curdir->cd();
  }

  // Get SFs
  {
    TString cinput = cinput_main + "scalefactorsPUID_81Xtraining.root";
    if (!HostHelpers::FileReadable(cinput.Data())){
      MELAerr << "PUJetIdScaleFactorHandler::setup: File " << cinput << " is not readable." << endl;
      assert(0);
    }

    TFile* finput = TFile::Open(cinput, "read"); curdir->cd();
    {
      unsigned short iwp=0;
      for (auto pujetidwpsfname:pujetidwpsfnames){
        HelperFunctions::replaceString<TString, const char*>(pujetidwpsfname, "<YEAR>", (const char*) Form("%i", SampleHelpers::getDataYear()));
        TString strmistagged = pujetidwpsfname; HelperFunctions::replaceString<TString, const char*>(strmistagged, "<EFFMISTAG>", "mistag");
        TString strmatched = pujetidwpsfname; HelperFunctions::replaceString<TString, const char*>(strmatched, "<EFFMISTAG>", "eff");
        res &= getHistogramWithUncertainy<TH2F, ExtendedHistogram_2D>(pujetidwp_SFs_map_mistagged.at(iwp), finput, strmistagged, strmistagged+"_Systuncty");
        res &= getHistogramWithUncertainy<TH2F, ExtendedHistogram_2D>(pujetidwp_SFs_map_matched.at(iwp), finput, strmatched, strmatched+"_Systuncty");
        iwp++;
      }
    }
    ScaleFactorHandlerBase::closeFile(finput); curdir->cd();
  }

  return res;
}
void PUJetIdScaleFactorHandler::reset(){
  syst_pujetidwp_effs_map_mistagged.clear();
  syst_pujetidwp_effs_map_matched.clear();

  pujetidwp_SFs_map_mistagged.clear();
  pujetidwp_SFs_map_matched.clear();
}

void PUJetIdScaleFactorHandler::evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, float const& pt, float const& eta, ExtendedHistogram_2D const& hist, bool etaOnY, bool useAbsEta) const{
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

void PUJetIdScaleFactorHandler::getSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, float const& pt, float const& eta, bool const& isMatched, bool const& isLoose, bool const& isMedium, bool const& isTight, float& val, float* effval) const{
  using namespace SystematicsHelpers;

  if (verbosity>=TVar::DEBUG) MELAout
    << "PUJetIdScaleFactorHandler::getSFAndEff: Evaluating " << (effval ? "SFs and efficiencies" : "SFs")
    << " for pT=" << pt
    << ", eta=" << eta
    << ", isMatched=" << isMatched
    << ", isLoose=" << isLoose
    << ", isMedium=" << isMedium
    << ", isTight=" << isTight
    << endl;

  val = 1;
  if (effval) *effval = 1;

  std::vector<SystematicVariationTypes> const allowedSysts={
    sNominal,
    eJECDn, eJECUp,
    eJERDn, eJERUp,
    ePUDn, ePUUp
  };

  SystematicVariationTypes activeSyst = sNominal;
  if (HelperFunctions::checkListVariable(allowedSysts, syst)) activeSyst = syst;
  if (verbosity>=TVar::DEBUG) MELAout << "\t- Active systematic: " << activeSyst << endl;

  std::vector<ExtendedHistogram_2D> const& hlist_eff_MC = (isMatched ? syst_pujetidwp_effs_map_matched.find(activeSyst)->second : syst_pujetidwp_effs_map_mistagged.find(activeSyst)->second);
  std::vector<ExtendedHistogram_2D> const& hlist_SF = (isMatched ? pujetidwp_SFs_map_matched : pujetidwp_SFs_map_mistagged);

  std::vector<float> eff_vals_uncorrected(hlist_eff_MC.size(), 1);
  std::vector<float> SF_vals(hlist_SF.size(), 1);
  {
    auto it_eff_hist = hlist_eff_MC.cbegin();
    auto it_SF_hist = hlist_SF.cbegin();
    auto it_eff_val = eff_vals_uncorrected.begin();
    auto it_SF_val = SF_vals.begin();
    while (it_eff_hist!=hlist_eff_MC.cend()){
      float eff_err = 0;
      float SF_err = 0;

      evalScaleFactorFromHistogram(*it_eff_val, eff_err, pt, eta, *it_eff_hist, false, false);
      evalScaleFactorFromHistogram(*it_SF_val, SF_err, pt, eta, *it_SF_hist, false, false);
      if (syst==ePUJetIdEffUp) *it_SF_val += SF_err;
      else if (syst==ePUJetIdEffDn) *it_SF_val -= SF_err;

      it_eff_hist++;
      it_SF_hist++;
      it_eff_val++;
      it_SF_val++;
    }
  }
  std::vector<float> eff_vals_corrected = eff_vals_uncorrected;
  for (unsigned short iwp=0; iwp<hlist_eff_MC.size(); iwp++){
    float SF_eff = SF_vals.at(iwp)/(iwp==0 ? 1.f : SF_vals.at(iwp-1));
    eff_vals_corrected.at(iwp) *= SF_eff;
  }

  float eff_uncorrected = 1;
  float eff_corrected = 1;
  for (unsigned int iwp=0; iwp<hlist_eff_MC.size(); iwp++){
    if (!isLoose && iwp>0) continue;
    else if (!isMedium && iwp>1) continue;

    bool checkFlag = false;
    if (iwp==0) checkFlag = isLoose;
    else if (iwp==1) checkFlag = isMedium;
    else checkFlag = isTight;

    if (checkFlag){
      eff_uncorrected *= eff_vals_uncorrected.at(iwp);
      eff_corrected *= eff_vals_corrected.at(iwp);
    }
    else{
      eff_uncorrected *= 1.f - eff_vals_uncorrected.at(iwp);
      eff_corrected *= 1.f - eff_vals_corrected.at(iwp);
    }
  }
  eff_uncorrected = std::max(0.f, std::min(eff_uncorrected, 1.f));
  eff_corrected = std::max(0.f, std::min(eff_corrected, 1.f));
  if (eff_uncorrected>0.f){
    val = eff_corrected / eff_uncorrected;
    if (effval) *effval = eff_corrected;
  }
  else{
    val = 0;
    if (effval) *effval = 0;
  }
}
void PUJetIdScaleFactorHandler::getSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, AK4JetObject const* obj, float& val, float* effval) const{
  val = 1;
  if (effval) *effval = 1;

  if (!obj) return;
  if (
    !(
      obj->testSelectionBit(AK4JetSelectionHelpers::kTightId)
      &&
      (!AK4JetSelectionHelpers::getApplyTightLeptonVetoIdToJetsFlag() || obj->testSelectionBit(AK4JetSelectionHelpers::kTightLeptonVetoId))
      )
    ||
    !
    (
      obj->pt()<AK4JetSelectionHelpers::ptThr_PUId
      &&
      std::abs(obj->eta())<AK4JetSelectionHelpers::etaThr_PUId
      )
    ) return;

  bool isMatched = obj->extras.is_genMatched_fullCone;
  bool isLoose = obj->testSelectionBit(AK4JetSelectionHelpers::kLoosePUJetId);
  bool isMedium = obj->testSelectionBit(AK4JetSelectionHelpers::kMediumPUJetId);
  bool isTight = obj->testSelectionBit(AK4JetSelectionHelpers::kTightPUJetId);

  getSFAndEff(syst, obj->pt(), obj->eta(), isMatched, isLoose, isMedium, isTight, val, effval);
}
