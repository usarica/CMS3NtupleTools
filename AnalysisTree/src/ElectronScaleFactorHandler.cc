#include "ElectronScaleFactorHandler.h"
#include "SamplesCore.h"
#include "ElectronSelectionHelpers.h"
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
  constexpr unsigned int n_non_gap_gap = 3;
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
      res &= getHistogram<TH2F, ExtendedHistogram_2D>(syst_SF_reco_map[sNominal].at(igap), finput, "EGamma_SF2D");
      res &= getHistogram<TH2F, ExtendedHistogram_2D>(eff_mc_reco_hists.at(igap), finput, "EGamma_EffMC2D");
    }
    ScaleFactorHandlerBase::closeFile(finput); curdir->cd();
  }

  for (unsigned int igap=0; igap<n_non_gap_gap; igap++){
    TString cinput = cinput_main + Form("Efficiencies_ee_%s_id_looseIso_tightIso.root", (igap==0 ? "nongap" : (igap==1 ? "gap" : "nongap_gap")));
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

void ElectronScaleFactorHandler::getIdIsoSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, float const& pt, float const& etaSC, unsigned short const& idx_gap, bool const& passId, bool const& passLooseIso, bool const& passTightIso, float& val, float* effval) const{
  using namespace SystematicsHelpers;

  if (verbosity>=TVar::DEBUG) MELAout
    << "ElectronScaleFactorHandler::getIdIsoSFAndEff: Evaluating " << (effval ? "SFs and efficiencies" : "SFs")
    << " for pT=" << pt
    << ", etaSC=" << etaSC
    << ", idx_gap=" << idx_gap
    << ", passId=" << passId
    << ", passLooseIso=" << passLooseIso
    << ", passTightIso=" << passTightIso
    << endl;

  val = 1;
  if (effval) *effval = 1;

  constexpr unsigned int n_ID_iso_types = 3;
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
  if (verbosity>=TVar::DEBUG) MELAout << "\t- Active systematics: " << activeSysts << endl;

  std::vector<ExtendedHistogram_2D const*> hlist_eff_mc; hlist_eff_mc.reserve(kAllEffs);
  std::vector<ExtendedHistogram_2D const*> hlist_SF_nominal; hlist_SF_nominal.reserve(kAllEffs);
  {
    auto it_syst_SF_reco_map = syst_SF_reco_map.find(sNominal);
    auto it_syst_SF_id_map = syst_SF_id_map.find(sNominal);
    auto it_syst_SF_iso_loose_map = syst_SF_iso_loose_map.find(sNominal);
    auto it_syst_SF_iso_tight_map = syst_SF_iso_tight_map.find(sNominal);
    if (it_syst_SF_reco_map!=syst_SF_reco_map.cend()){
      hlist_eff_mc.push_back(&(eff_mc_reco_hists.at(idx_gap)));
      hlist_SF_nominal.push_back(&(it_syst_SF_reco_map->second.at(idx_gap)));
    }
    else{
      hlist_eff_mc.push_back(nullptr);
      hlist_SF_nominal.push_back(nullptr);
    }
    if (it_syst_SF_id_map!=syst_SF_id_map.cend()){
      hlist_eff_mc.push_back(&eff_mc_id_hists.at(idx_gap));
      hlist_SF_nominal.push_back(&(it_syst_SF_id_map->second.at(idx_gap)));
    }
    else{
      hlist_eff_mc.push_back(nullptr);
      hlist_SF_nominal.push_back(nullptr);
    }
    if (it_syst_SF_iso_loose_map!=syst_SF_iso_loose_map.cend()){
      hlist_eff_mc.push_back(&eff_mc_iso_loose_hists.at(idx_gap));
      hlist_SF_nominal.push_back(&(it_syst_SF_iso_loose_map->second.at(idx_gap)));
    }
    else{
      hlist_eff_mc.push_back(nullptr);
      hlist_SF_nominal.push_back(nullptr);
    }
    if (it_syst_SF_iso_tight_map!=syst_SF_iso_tight_map.cend()){
      hlist_eff_mc.push_back(&eff_mc_iso_tight_hists.at(idx_gap));
      hlist_SF_nominal.push_back(&(it_syst_SF_iso_tight_map->second.at(idx_gap)));
    }
    else{
      hlist_eff_mc.push_back(nullptr);
      hlist_SF_nominal.push_back(nullptr);
    }
  }
  // Obtain nominal efficiencies before ID/iso assessment
  std::vector<float> eff_nominal_unscaled_list(n_ID_iso_types+1, 1);
  std::vector<float> eff_nominal_scaled_list(n_ID_iso_types+1, 1);
  {
    auto it_SF = hlist_SF_nominal.begin();
    auto it_eff_mc = hlist_eff_mc.begin();
    unsigned int isel = 0;
    while (it_SF != hlist_SF_nominal.end()){
      assert(isel<n_ID_iso_types+1);

      float SF_val = 1;
      if (*it_SF){
        float SF_err = 0;
        evalScaleFactorFromHistogram(SF_val, SF_err, pt, etaSC, **it_SF, false, false);
        float eff_err = 0;
        evalScaleFactorFromHistogram(eff_nominal_unscaled_list.at(isel), eff_err, pt, etaSC, **it_eff_mc, false, false);
      }

      eff_nominal_scaled_list.at(isel) = std::max(0.f, std::min(1.f, SF_val * eff_nominal_unscaled_list.at(isel)));

      it_SF++;
      it_eff_mc++;
      isel++;
    }
  }
  // Calculate the actual ID+iso efficiency based on the three flags
  float SF_err_val = 0;
  float SF_nominal_val = 1;
  float eff_nominal_unscaled_val = 1;
  float eff_nominal_scaled_val = 1;
  for (unsigned int isel=0; isel<n_ID_iso_types+1; isel++){
    if (!passId && isel>1) continue;
    if (!passLooseIso && isel>2) continue;
    float tmp_eff_unscaled=1;
    float tmp_eff_scaled=1;
    bool checkFlag = false;
    if (isel==0) checkFlag = true; // Tracking efficiency always true
    else if (isel==1) checkFlag = passId;
    else if (isel==2) checkFlag = passLooseIso;
    else checkFlag = passTightIso;
    if (checkFlag){
      tmp_eff_unscaled = eff_nominal_unscaled_list.at(isel);
      tmp_eff_scaled = eff_nominal_scaled_list.at(isel);
    }
    else{
      tmp_eff_unscaled = 1. - eff_nominal_unscaled_list.at(isel);
      tmp_eff_scaled = 1. - eff_nominal_scaled_list.at(isel);
    }
    if (verbosity>=TVar::DEBUG){
      if (isel==0) MELAout << "\t- Nominal tracking efficiency ";
      else if (isel==1) MELAout << "\t- Nominal ID efficiency ";
      else if (isel==2) MELAout << "\t- Nominal loose iso. efficiency ";
      else MELAout << "\t- Nominal tight iso. efficiency ";
      MELAout << "(pass/fail, unscaled, scaled, SF) = (" << checkFlag << ", " << tmp_eff_unscaled << ", " << tmp_eff_scaled << ", " << tmp_eff_scaled / tmp_eff_unscaled << ")" << endl;
    }

    SF_nominal_val *= tmp_eff_scaled / tmp_eff_unscaled;
    eff_nominal_unscaled_val *= tmp_eff_unscaled;
    eff_nominal_scaled_val *= tmp_eff_scaled;
  }

  if (!(activeSysts.size() == 1 && activeSysts.front() == sNominal)){
    std::vector< std::vector<float> > eff_syst_scaled_lists(activeSysts.size(), std::vector<float>(n_ID_iso_types+1, 1));
    std::vector< std::vector<float> > eff_syst_scaled_complement_lists(activeSysts.size(), std::vector<float>(n_ID_iso_types+1, 1));
    {
      unsigned int ias = 0;
      for (auto const& asyst:activeSysts){
        SystematicVariationTypes asyst_cpl = getSystComplement(asyst);

        std::vector<float>& eff_syst_scaled_list = eff_syst_scaled_lists.at(ias);
        std::vector<float>& eff_syst_scaled_complement_list = eff_syst_scaled_complement_lists.at(ias);

        std::vector<ExtendedHistogram_2D const*> hlist_SF; hlist_SF.reserve(n_ID_iso_types+1);
        std::vector<ExtendedHistogram_2D const*> hlist_SF_cpl; hlist_SF_cpl.reserve(n_ID_iso_types+1);
        {
          auto it_syst_SF_reco_map = syst_SF_reco_map.find(sNominal);
          auto it_syst_SF_id_map = syst_SF_id_map.find(asyst);
          auto it_syst_SF_iso_loose_map = syst_SF_iso_loose_map.find(asyst);
          auto it_syst_SF_iso_tight_map = syst_SF_iso_tight_map.find(asyst);

          if (it_syst_SF_reco_map!=syst_SF_reco_map.cend()){
            hlist_SF.push_back(&(it_syst_SF_reco_map->second.at(idx_gap)));
            if (activeSysts.size() > 1) hlist_SF_cpl.push_back(&(it_syst_SF_reco_map->second.at(idx_gap)));
          }
          else{
            hlist_SF.push_back(nullptr);
            if (activeSysts.size() > 1) hlist_SF_cpl.push_back(nullptr);
          }
          if (it_syst_SF_id_map!=syst_SF_id_map.cend()) hlist_SF.push_back(&(it_syst_SF_id_map->second.at(idx_gap)));
          else hlist_SF.push_back(nullptr);
          if (it_syst_SF_iso_loose_map!=syst_SF_iso_loose_map.cend()) hlist_SF.push_back(&(it_syst_SF_iso_loose_map->second.at(idx_gap)));
          else hlist_SF.push_back(nullptr);
          if (it_syst_SF_iso_tight_map!=syst_SF_iso_tight_map.cend()) hlist_SF.push_back(&(it_syst_SF_iso_tight_map->second.at(idx_gap)));
          else hlist_SF.push_back(nullptr);
        }

        // Get complementary syst
        if (activeSysts.size() > 1){
          auto it_syst_SF_id_map = syst_SF_id_map.find(asyst_cpl);
          auto it_syst_SF_iso_loose_map = syst_SF_iso_loose_map.find(asyst_cpl);
          auto it_syst_SF_iso_tight_map = syst_SF_iso_tight_map.find(asyst_cpl);

          if (it_syst_SF_id_map!=syst_SF_id_map.cend()) hlist_SF_cpl.push_back(&(it_syst_SF_id_map->second.at(idx_gap)));
          else hlist_SF_cpl.push_back(nullptr);
          if (it_syst_SF_iso_loose_map!=syst_SF_iso_loose_map.cend()) hlist_SF_cpl.push_back(&(it_syst_SF_iso_loose_map->second.at(idx_gap)));
          else hlist_SF_cpl.push_back(nullptr);
          if (it_syst_SF_iso_tight_map!=syst_SF_iso_tight_map.cend()) hlist_SF_cpl.push_back(&(it_syst_SF_iso_tight_map->second.at(idx_gap)));
          else hlist_SF_cpl.push_back(nullptr);
        }

        {
          auto it_SF = hlist_SF.begin(); // Input
          auto it_SF_cpl = hlist_SF_cpl.begin(); // Input
          auto it_eff_nominal_unscaled_val = eff_nominal_unscaled_list.begin(); // Input
          auto it_eff_nominal_scaled_val = eff_nominal_scaled_list.begin(); // Input
          auto it_eff_syst_scaled_val = eff_syst_scaled_list.begin(); // Output
          auto it_eff_syst_scaled_cpl_val = eff_syst_scaled_complement_list.begin(); // Output
          unsigned short ihist = 0;
          while (it_SF != hlist_SF.end()){
            if (*it_SF){
              float SF_val = 1;
              float val_err = 0;
              evalScaleFactorFromHistogram(SF_val, val_err, pt, etaSC, **it_SF, false, false);
              if (ihist>0) *it_eff_syst_scaled_val = std::max(0.f, std::min(1.f, SF_val * (*it_eff_nominal_unscaled_val)));
              else *it_eff_syst_scaled_val = std::max(0.f, std::min(1.f, (SF_val + val_err*(1.f*(asyst==eEleEffSystUp)-1.f*(asyst==eEleEffSystDn))) * (*it_eff_nominal_unscaled_val)));

              if (verbosity>=TVar::DEBUG) MELAout
                << "\t\t- Evaluating SF for syst " << asyst << ". SF_syst = " << SF_val << ", eff = " << *it_eff_syst_scaled_val
                << " from histogram " << (*it_SF)->getName() << endl;
            }
            else *it_eff_syst_scaled_val = *it_eff_nominal_scaled_val;

            if (!hlist_SF_cpl.empty()){
              if (*it_SF_cpl){
                float SF_val = 1;
                float val_err = 0;
                evalScaleFactorFromHistogram(SF_val, val_err, pt, etaSC, **it_SF_cpl, false, false);
                if (ihist>0) *it_eff_syst_scaled_cpl_val = std::max(0.f, std::min(1.f, SF_val * (*it_eff_nominal_unscaled_val)));
                else *it_eff_syst_scaled_cpl_val = std::max(0.f, std::min(1.f, (SF_val + val_err*(1.f*(asyst==eEleEffSystUp)-1.f*(asyst==eEleEffSystDn))) * (*it_eff_nominal_unscaled_val)));

                if (verbosity>=TVar::DEBUG) MELAout
                  << "\t\t- Evaluating complementary SF for syst " << asyst << ". SF_syst = " << SF_val << ", eff = " << *it_eff_syst_scaled_cpl_val
                  << " from histogram " << (*it_SF_cpl)->getName() << endl;
              }
              else *it_eff_syst_scaled_cpl_val = *it_eff_nominal_scaled_val;

              it_SF_cpl++;
            }
            else *it_eff_syst_scaled_cpl_val = *it_eff_syst_scaled_val;

            it_SF++;
            it_eff_nominal_unscaled_val++;
            it_eff_nominal_scaled_val++;
            it_eff_syst_scaled_val++;
            it_eff_syst_scaled_cpl_val++;
            ihist++;
          }
        }

        ias++;
      }
    }
    for (unsigned int ias=0; ias<activeSysts.size(); ias++){
      auto const& asyst = activeSysts.at(ias);
      std::vector<float> const& eff_syst_scaled_list = eff_syst_scaled_lists.at(ias);
      std::vector<float> const& eff_syst_scaled_complement_list = eff_syst_scaled_complement_lists.at(ias);

      float eff_syst = 1;
      float eff_syst_cpl = 1;
      for (unsigned int isel=0; isel<n_ID_iso_types+1; isel++){
        if (!passId && isel>1) continue;
        if (!passLooseIso && isel>2) continue;

        float tmp_eff_syst=1;
        float tmp_eff_syst_cpl=1;

        bool checkFlag = false;
        if (isel==0) checkFlag = true;
        else if (isel==1) checkFlag = passId;
        else if (isel==2) checkFlag = passLooseIso;
        else checkFlag = passTightIso;

        tmp_eff_syst = eff_syst_scaled_list.at(isel);
        if (activeSysts.size()>1) tmp_eff_syst_cpl = eff_syst_scaled_complement_list.at(isel);

        if (!checkFlag){
          tmp_eff_syst = 1. - tmp_eff_syst;
          tmp_eff_syst_cpl = 1. - tmp_eff_syst_cpl;
        }

        eff_syst *= tmp_eff_syst;
        eff_syst_cpl *= tmp_eff_syst_cpl;
      }

      float SF_syst = eff_syst / eff_nominal_unscaled_val;
      float SF_syst_cpl = eff_syst_cpl / eff_nominal_unscaled_val;
      if (activeSysts.size() == 1) SF_err_val = SF_syst - SF_nominal_val;
      else{
        if (verbosity>=TVar::DEBUG) MELAout
          << "\t\t- Adding SF for syst " << asyst << ". SF_syst = " << SF_syst << ", SF_syst_cpl = " << SF_syst_cpl
          << ", SF_nominal=" << SF_nominal_val
          << ", current SF error value = " << SF_err_val
          << endl;
        float SF_syst_eff = (isDownSystematic(syst) ? std::min(SF_syst, SF_syst_cpl) : std::max(SF_syst, SF_syst_cpl));

        SF_err_val = std::sqrt(std::pow(SF_err_val, 2) + std::pow(SF_syst_eff - SF_nominal_val, 2));
        if (SF_syst_eff < SF_nominal_val) SF_err_val *= -1.;
      }
    }
  }

  val = SF_nominal_val + SF_err_val;
  if (effval) *effval = std::max(0.f, std::min(1.f, eff_nominal_unscaled_val*val));

  if (verbosity>=TVar::DEBUG){
    MELAout << "\t- Final SF = " << SF_nominal_val << " + " << SF_err_val << " = " << val << endl;
    if (effval) MELAout << "\t- Final eff = " << eff_nominal_unscaled_val << " x " << val << " = " << *effval << endl;
  }
}
void ElectronScaleFactorHandler::getIdIsoSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, ElectronObject const* obj, float& val, float* effval) const{
  val = 1;
  if (effval) *effval = 1;

  if (!obj) return;
  if (!obj->extras.is_genMatched_prompt) return;
  if (verbosity>=TVar::DEBUG) MELAout << "ElectronScaleFactorHandler::getIdIsoSFAndEff: Electron gen matching flags: " << obj->extras.is_genMatched << ", " << obj->extras.is_genMatched_prompt << endl;

  bool passId = obj->testSelectionBit(ElectronSelectionHelpers::bit_preselectionTight_id);
  bool passLooseIso = passId && obj->testSelectionBit(ElectronSelectionHelpers::kFakeableBaseIso);
  bool passTightIso = passId && obj->testSelectionBit(ElectronSelectionHelpers::bit_preselectionTight_iso);
  if (passTightIso) assert(passLooseIso);

  getIdIsoSFAndEff(syst, obj->pt(), obj->etaSC(), static_cast<unsigned short>(obj->isAnyGap()), passId, passLooseIso, passTightIso, val, effval);
}
