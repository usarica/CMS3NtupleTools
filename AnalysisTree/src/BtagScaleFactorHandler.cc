#include <cassert>
#include "SampleHelpersCore.h"
#include "SamplesCore.h"
#include "BtagScaleFactorHandler.h"
#include "AK4JetSelectionHelpers.h"
#include "ParticleSelectionHelpers.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace SampleHelpers;
using namespace BtagHelpers;
using namespace MELAStreamHelpers;


BtagScaleFactorHandler::BtagScaleFactorHandler() : ScaleFactorHandlerBase()
{
  setup();
}

BtagScaleFactorHandler::~BtagScaleFactorHandler(){ this->reset(); }

void BtagScaleFactorHandler::evalEfficiencyFromHistogram(float& theSF, float const& pt, float const& eta, ExtendedHistogram_2D const& hist, bool etaOnY, bool useAbsEta) const{
  TH2F const* hh = hist.getHistogram();
  if (!hh){
    MELAerr << "BtagScaleFactorHandler::evalScaleFactorFromHistogram: Histogram is null." << endl;
    return;
  }

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
  else if (ix>nbinsx+1) ix=nbinsx+1; // Overflows exist
  if (iy==0) iy=1;
  else if (iy==nbinsy+1) iy=nbinsy;

  theSF = hh->GetBinContent(ix, iy);
}

bool BtagScaleFactorHandler::setup(){
  using namespace SystematicsHelpers;

  bool res = true;
  this->reset();

  if (verbosity>=TVar::INFO) MELAout << "BtagScaleFactorHandler::setup: Setting up efficiency and SF histograms for year " << SampleHelpers::getDataYear() << endl;

  TDirectory* curdir = gDirectory;
  TDirectory* uppermostdir = SampleHelpers::rootTDirectory;

  TString sf_fname_deepCSV = BtagHelpers::getBtagSFFileName(kDeepCSV_Loose);
  WP_calib_map[kDeepCSV_Loose] = new BTagCalibration("DeepCSV", sf_fname_deepCSV.Data());

  TString sf_fname_deepFlavor = BtagHelpers::getBtagSFFileName(kDeepFlav_Loose);
  WP_calib_map[kDeepFlav_Loose] = new BTagCalibration("DeepFlavor", sf_fname_deepFlavor.Data());

  for (int iwp=0; iwp<(int) nBtagWPTypes; iwp++){
    BtagWPType type = (BtagWPType) iwp;
    BTagEntry::OperatingPoint opPoint;
    BTagCalibration* calibration=nullptr;
    switch (type){
    case kDeepCSV_Loose:
      opPoint = BTagEntry::OP_LOOSE;
      calibration = WP_calib_map[kDeepCSV_Loose];
      break;
    case kDeepFlav_Loose:
      opPoint = BTagEntry::OP_LOOSE;
      calibration = WP_calib_map[kDeepFlav_Loose];
      break;
    case kDeepCSV_Medium:
      opPoint = BTagEntry::OP_MEDIUM;
      calibration = WP_calib_map[kDeepCSV_Loose];
      break;
    case kDeepFlav_Medium:
      opPoint = BTagEntry::OP_MEDIUM;
      calibration = WP_calib_map[kDeepFlav_Loose];
      break;
    case kDeepCSV_Tight:
      opPoint = BTagEntry::OP_TIGHT;
      calibration = WP_calib_map[kDeepCSV_Loose];
      break;
    case kDeepFlav_Tight:
      opPoint = BTagEntry::OP_TIGHT;
      calibration = WP_calib_map[kDeepFlav_Loose];
      break;
    default:
      MELAerr << "BtagScaleFactorHandler::setup: No operating point implementation for b tag WP " << type << ". Aborting..." << endl;
      assert(0);
      break;
    }

    WP_calibreader_map_nominal[type] = new BTagCalibrationReader(opPoint, "central");
    loadBTagCalibrations(WP_calibreader_map_nominal[type], calibration);

    WP_calibreader_map_dn[type] = new BTagCalibrationReader(opPoint, "down");
    loadBTagCalibrations(WP_calibreader_map_dn[type], calibration);

    WP_calibreader_map_up[type] = new BTagCalibrationReader(opPoint, "up");
    loadBTagCalibrations(WP_calibreader_map_up[type], calibration);
  }

  std::vector<SystematicsHelpers::SystematicVariationTypes> const allowedSysts{
    sNominal,
    eJECDn, eJECUp,
    eJERDn, eJERUp,
    ePUDn, ePUUp,
    ePUJetIdEffDn, ePUJetIdEffUp
  };
  std::vector< std::pair<BTagEntry::JetFlavor, TString> > const flavpairs{ { BTagEntry::FLAV_B, "b" }, { BTagEntry::FLAV_C, "c" }, { BTagEntry::FLAV_UDSG, "udsg" } };
  std::vector<TString> const strpujetidcats{ "T", "MnT", "LnM", "F" };
  TFile* finput_eff = TFile::Open(BtagHelpers::getBtagEffFileName(), "read"); uppermostdir->cd();
  if (verbosity>=TVar::INFO) MELAout << "BtagScaleFactorHandler::setup: Reading " << finput_eff->GetName() << " to acquire efficiency histograms..." << endl;
  {
    ExtendedHistogram_2D empty_hist; empty_hist.reset();
    TString hname;
    for (auto const& syst:allowedSysts){
      TString systname = SystematicsHelpers::getSystName(syst).data();
      syst_flav_pujetid_WP_mceffhist_map[syst] = std::unordered_map< BTagEntry::JetFlavor, std::vector<std::vector<ExtendedHistogram_2D>> >();
      for (auto const& flavpair:flavpairs){
        BTagEntry::JetFlavor jflav = flavpair.first;
        TString const& strflav = flavpair.second;
        syst_flav_pujetid_WP_mceffhist_map[syst][jflav] = std::vector<std::vector<ExtendedHistogram_2D>>(strpujetidcats.size(), std::vector<ExtendedHistogram_2D>(nBtagWPTypes, empty_hist));
        for (unsigned short ipujetidwp=0; ipujetidwp<strpujetidcats.size(); ipujetidwp++){
          TString const& strpujetidcat = strpujetidcats.at(ipujetidwp);
          for (int iwp=0; iwp<(int) nBtagWPTypes; iwp++){
            BtagWPType wptype = (BtagWPType) iwp;
            hname = BtagHelpers::getBtagEffHistName(wptype, strflav.Data()); hname = hname + "_PUJetId_" + strpujetidcat + "_" + systname;
            if (verbosity>=TVar::DEBUG) MELAout << "\t- Extracting MC efficiency histogram " << hname << "..." << endl;
            bool tmpres = getHistogram<TH2F, ExtendedHistogram_2D>(syst_flav_pujetid_WP_mceffhist_map[syst][jflav].at(ipujetidwp).at(iwp), finput_eff, hname);
            if (!tmpres && verbosity>=TVar::DEBUG) MELAerr << "\t\t- FAILED!" << endl;
            res &= tmpres;
          }
        }
      }
    }
  }
  ScaleFactorHandlerBase::closeFile(finput_eff); curdir->cd();

  return res;
}
void BtagScaleFactorHandler::loadBTagCalibrations(BTagCalibrationReader* const& reader, BTagCalibration* const& calibration){
  reader->load(*calibration, BTagEntry::FLAV_B, "comb");
  reader->load(*calibration, BTagEntry::FLAV_C, "comb");
  reader->load(*calibration, BTagEntry::FLAV_UDSG, "incl");
}

void BtagScaleFactorHandler::reset(){
  for (auto it:WP_calibreader_map_up) delete it.second;
  for (auto it:WP_calibreader_map_dn) delete it.second;
  for (auto it:WP_calibreader_map_nominal) delete it.second;
  for (auto it:WP_calib_map) delete it.second;

  WP_calibreader_map_up.clear();
  WP_calibreader_map_dn.clear();
  WP_calibreader_map_nominal.clear();
  WP_calib_map.clear();
  syst_flav_pujetid_WP_mceffhist_map.clear();
}

float BtagScaleFactorHandler::getSFFromBTagCalibrationReader(
  BTagCalibrationReader const* calibReader, BTagCalibrationReader const* calibReader_Nominal,
  SystematicsHelpers::SystematicVariationTypes const& syst, BTagEntry::JetFlavor const& flav, float const& pt, float const& eta
) const{
  float myPt = pt;
  float const MaxJetEta = (SampleHelpers::getDataYear()<=2016 ? 2.4 : 2.5);
  if (std::abs(eta) > MaxJetEta) return 1; // Do not apply SF for jets with eta higher than the threshold

  std::pair<float, float> pt_min_max = calibReader->min_max_pt(flav, eta);
  if (pt_min_max.second<0.) return 1;
  bool DoubleUncertainty = false;
  if (myPt<pt_min_max.first){
    myPt = pt_min_max.first+1e-5;
    DoubleUncertainty = true;
  }
  else if (myPt>pt_min_max.second){
    myPt = pt_min_max.second-1e-5;
    DoubleUncertainty = true;
  }

  float SF = calibReader->eval(flav, eta, myPt);
  if (DoubleUncertainty && (syst==SystematicsHelpers::eBTagSFDn || syst==SystematicsHelpers::eBTagSFUp)){
    float SFcentral = calibReader_Nominal->eval(flav, eta, myPt);
    SF = 2.f*(SF - SFcentral) + SFcentral;
  }

  return SF;
}

void BtagScaleFactorHandler::getSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, float const& pt, float const& eta, unsigned short const& pujetidcat, BTagEntry::JetFlavor const& flav, float const& btagval, float& val, float* effval) const{
  using namespace SystematicsHelpers;

  val = 1;

  if (this->verbosity>=TVar::DEBUG) MELAout
    << "BtagScaleFactorHandler::getSFAndEff: Calling for jet (pt, eta, PU jet id cat, flav, btagval) = ("
    << pt << ", " << eta << ", " << pujetidcat << ", " << flav << ", " << btagval << "):"
    << endl;

  std::unordered_map<BtagHelpers::BtagWPType, BTagCalibrationReader*> const* WP_calibreader_map = nullptr;
  switch (syst){
  case SystematicsHelpers::eBTagSFDn:
    WP_calibreader_map = &WP_calibreader_map_dn;
    break;
  case SystematicsHelpers::eBTagSFUp:
    WP_calibreader_map = &WP_calibreader_map_up;
    break;
  default:
    WP_calibreader_map = &WP_calibreader_map_nominal;
    break;
  }

  std::vector<SystematicsHelpers::SystematicVariationTypes> const allowedJetSysts{
    sNominal,
    eJECDn, eJECUp,
    eJERDn, eJERUp,
    ePUDn, ePUUp,
    ePUJetIdEffDn, ePUJetIdEffUp
  };
  SystematicsHelpers::SystematicVariationTypes const jetsyst = (HelperFunctions::checkListVariable(allowedJetSysts, syst) ? syst : SystematicsHelpers::sNominal);

  std::vector<BTagCalibrationReader const*> calibReaders;
  std::vector<BTagCalibrationReader const*> calibReaders_Nominal;
  std::vector<float> const btagwps = BtagHelpers::getBtagWPs(false);
  std::vector<ExtendedHistogram_2D> const& effhists = syst_flav_pujetid_WP_mceffhist_map.find(jetsyst)->second.find(flav)->second.at(pujetidcat);
  unsigned short idx_offset_effmc = 0;
  switch (BtagHelpers::btagWPType){
  case kDeepCSV_Loose:
  case kDeepCSV_Medium:
  case kDeepCSV_Tight:
  {
    calibReaders.push_back(WP_calibreader_map->find(BtagHelpers::kDeepCSV_Loose)->second);
    calibReaders.push_back(WP_calibreader_map->find(BtagHelpers::kDeepCSV_Medium)->second);
    calibReaders.push_back(WP_calibreader_map->find(BtagHelpers::kDeepCSV_Tight)->second);

    calibReaders_Nominal.push_back(WP_calibreader_map_nominal.find(BtagHelpers::kDeepCSV_Loose)->second);
    calibReaders_Nominal.push_back(WP_calibreader_map_nominal.find(BtagHelpers::kDeepCSV_Medium)->second);
    calibReaders_Nominal.push_back(WP_calibreader_map_nominal.find(BtagHelpers::kDeepCSV_Tight)->second);

    idx_offset_effmc = (unsigned short) kDeepCSV_Loose;
    break;
  }
  case kDeepFlav_Loose:
  case kDeepFlav_Medium:
  case kDeepFlav_Tight:
  {
    calibReaders.push_back(WP_calibreader_map->find(BtagHelpers::kDeepFlav_Loose)->second);
    calibReaders.push_back(WP_calibreader_map->find(BtagHelpers::kDeepFlav_Medium)->second);
    calibReaders.push_back(WP_calibreader_map->find(BtagHelpers::kDeepFlav_Tight)->second);

    calibReaders_Nominal.push_back(WP_calibreader_map_nominal.find(BtagHelpers::kDeepFlav_Loose)->second);
    calibReaders_Nominal.push_back(WP_calibreader_map_nominal.find(BtagHelpers::kDeepFlav_Medium)->second);
    calibReaders_Nominal.push_back(WP_calibreader_map_nominal.find(BtagHelpers::kDeepFlav_Tight)->second);

    idx_offset_effmc = (unsigned short) kDeepFlav_Loose;
    break;
  }
  default:
    MELAerr << "BtagScaleFactorHandler::getSFAndEff: b tag calibration readers are not implemented." << endl;
    assert(0);
    break;
  }

  std::vector<float> effs_unscaled(calibReaders.size(), 1);
  std::vector<float> effs_scaled(calibReaders.size(), 1);
  std::vector<float> SFs(calibReaders.size(), 1);
  for (unsigned int i=0; i<calibReaders.size(); i++){
    SFs.at(i) = getSFFromBTagCalibrationReader(
      calibReaders.at(i), calibReaders_Nominal.at(i),
      syst, flav, pt, eta
    );
    if (i>0) SFs.at(i) /= SFs.at(i-1);
    evalEfficiencyFromHistogram(effs_unscaled.at(i), pt, eta, effhists.at(i+idx_offset_effmc), true, false);
    effs_scaled.at(i) = std::max(0.f, std::min(1.f, effs_unscaled.at(i) * SFs.at(i)));
  }

  float eff_unscaled = 1;
  float eff_scaled = 1;
  for (unsigned short iwp=0; iwp<calibReaders.size(); iwp++){
    float tmp_eff_unscaled = 1;
    float tmp_eff_scaled = 1;
    bool doContinue = true;
    if (btagval>=btagwps.at(iwp)){
      tmp_eff_unscaled *= effs_unscaled.at(iwp);
      tmp_eff_scaled *= effs_scaled.at(iwp);
      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Passed";
    }
    else{
      tmp_eff_unscaled *= 1.f - effs_unscaled.at(iwp);
      tmp_eff_scaled *= 1.f - effs_scaled.at(iwp);
      doContinue = false;
      if (this->verbosity>=TVar::DEBUG) MELAout << "\t- Failed";
    }
    if (this->verbosity>=TVar::DEBUG) MELAout << " WP " << iwp << " with unscaled, scaled effs = " << tmp_eff_unscaled << ", " << tmp_eff_scaled << endl;
    eff_unscaled *= tmp_eff_unscaled;
    eff_scaled *= tmp_eff_scaled;
    if (!doContinue) break;
  }

  val = (eff_unscaled>0.f ? eff_scaled / eff_unscaled : 0.f);
  if (effval) *effval = (val>0.f ? eff_scaled : 0.f);

  if (this->verbosity>=TVar::DEBUG){
    MELAout << "\t- WPs: " << btagwps << endl;
    MELAout << "\t- Unscaled effs: " << effs_unscaled << endl;
    MELAout << "\t- Scaled effs: " << effs_scaled << endl;
    MELAout << "\t- SFs: " << SFs << endl;
    MELAout << "\t- Final SF: " << val << endl;
    MELAout << "\t- Final eff: " << eff_scaled << endl;
  }
}
void BtagScaleFactorHandler::getSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, AK4JetObject const* obj, float& val, float* effval) const{
  val = 1;
  if (effval) *effval = 1;

  if (!obj) return;
  if (!obj->testSelectionBit(AK4JetSelectionHelpers::kBtaggable_NoPUJetId)) return;

  unsigned short pujetidcat=0;
  if (!obj->testSelectionBit(AK4JetSelectionHelpers::kTightPUJetId)) pujetidcat++;
  if (!obj->testSelectionBit(AK4JetSelectionHelpers::kMediumPUJetId)) pujetidcat++;
  if (!obj->testSelectionBit(AK4JetSelectionHelpers::kLoosePUJetId)) pujetidcat++;

  getSFAndEff(syst, obj->pt(), obj->eta(), pujetidcat, obj->getBTagJetFlavor(), obj->getBtagValue(), val, effval);
}
