#include <cassert>
#include "SampleHelpersCore.h"
#include "SamplesCore.h"
#include "BtagScaleFactorHandler.h"
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
  bool res = true;
  this->reset();

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

  TFile* finput_eff = TFile::Open(BtagHelpers::getBtagEffFileName(), "read"); uppermostdir->cd();
  for (int iwp=0; iwp<(int) nBtagWPTypes; iwp++){
    BtagWPType type = (BtagWPType) iwp;
    TString hname;
    ExtendedHistogram_2D empty_hist; empty_hist.reset();

    WP_flav_mceffhist_map[type] = std::unordered_map<BTagEntry::JetFlavor, ExtendedHistogram_2D>();

    WP_flav_mceffhist_map[type][BTagEntry::FLAV_B] = empty_hist; WP_flav_mceffhist_map[type][BTagEntry::FLAV_B].reset();
    hname = BtagHelpers::getBtagEffHistName(type, "b");
    //MELAout << "BtagScaleFactorHandler::setup: Acquiring histogram " << hname << endl;
    res &= getHistogram<TH2F, ExtendedHistogram_2D>(WP_flav_mceffhist_map[type][BTagEntry::FLAV_B], finput_eff, hname);

    WP_flav_mceffhist_map[type][BTagEntry::FLAV_C] = empty_hist; WP_flav_mceffhist_map[type][BTagEntry::FLAV_C].reset();
    hname = BtagHelpers::getBtagEffHistName(type, "c");
    //MELAout << "BtagScaleFactorHandler::setup: Acquiring histogram " << hname << endl;
    res &= getHistogram<TH2F, ExtendedHistogram_2D>(WP_flav_mceffhist_map[type][BTagEntry::FLAV_C], finput_eff, hname);

    WP_flav_mceffhist_map[type][BTagEntry::FLAV_UDSG] = empty_hist; WP_flav_mceffhist_map[type][BTagEntry::FLAV_UDSG].reset();
    hname = BtagHelpers::getBtagEffHistName(type, "udsg");
    //MELAout << "BtagScaleFactorHandler::setup: Acquiring histogram " << hname << endl;
    res &= getHistogram<TH2F, ExtendedHistogram_2D>(WP_flav_mceffhist_map[type][BTagEntry::FLAV_UDSG], finput_eff, hname);
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
  for (auto it:WP_flav_mceffhist_map){ for (auto it2:it.second) it2.second.reset(); }
  for (auto it:WP_calibreader_map_up) delete it.second;
  for (auto it:WP_calibreader_map_dn) delete it.second;
  for (auto it:WP_calibreader_map_nominal) delete it.second;
  for (auto it:WP_calib_map) delete it.second;

  WP_calibreader_map_up.clear();
  WP_calibreader_map_dn.clear();
  WP_calibreader_map_nominal.clear();
  WP_calib_map.clear();
}

float BtagScaleFactorHandler::getSFFromBTagCalibrationReader(
  BTagCalibrationReader const* calibReader, BTagCalibrationReader const* calibReader_Nominal,
  SystematicsHelpers::SystematicVariationTypes const& syst, BTagEntry::JetFlavor const& flav, float const& pt, float const& eta
) const{
  float myPt = pt;
  float const MaxJetEta = (SampleHelpers::theDataYear<=2016 ? 2.4 : 2.5);
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

void BtagScaleFactorHandler::getSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, BTagEntry::JetFlavor const& flav, float const& btagval, float const& pt, float const& eta, float& val, float* effval) const{
  val = 1;

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

  std::vector<BTagCalibrationReader const*> calibReaders;
  std::vector<BTagCalibrationReader const*> calibReaders_Nominal;
  std::vector<ExtendedHistogram_2D const*> effhists;
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

    effhists.push_back(&(WP_flav_mceffhist_map.find(BtagHelpers::kDeepCSV_Loose)->second.find(flav)->second));
    effhists.push_back(&(WP_flav_mceffhist_map.find(BtagHelpers::kDeepCSV_Medium)->second.find(flav)->second));
    effhists.push_back(&(WP_flav_mceffhist_map.find(BtagHelpers::kDeepCSV_Tight)->second.find(flav)->second));
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

    effhists.push_back(&(WP_flav_mceffhist_map.find(BtagHelpers::kDeepFlav_Loose)->second.find(flav)->second));
    effhists.push_back(&(WP_flav_mceffhist_map.find(BtagHelpers::kDeepFlav_Medium)->second.find(flav)->second));
    effhists.push_back(&(WP_flav_mceffhist_map.find(BtagHelpers::kDeepFlav_Tight)->second.find(flav)->second));
    break;
  }
  default:
    MELAerr << "BtagScaleFactorHandler::getSFAndEff: b tag calibration readers are not implemented." << endl;
    assert(0);
    break;
  }

  std::vector<float> effraw(calibReaders.size(), 1);
  std::vector<float> SFraw(calibReaders.size(), 1);
  for (unsigned int i=0; i<calibReaders.size(); i++){
    SFraw.at(i) = getSFFromBTagCalibrationReader(
      calibReaders.at(i), calibReaders_Nominal.at(i),
      syst, flav, pt, eta
    );
    evalEfficiencyFromHistogram(effraw.at(i), pt, eta, *(effhists.at(i)), true, false);
  }

  std::vector<float> eff(calibReaders.size()+1, 1);
  std::vector<float> effScaled(calibReaders.size()+1, 1);
  for (unsigned int i=0; i<calibReaders.size()+1; i++){
    if (i==0){
      eff.at(i) = std::max(0.f, 1.f - effraw.at(i));
      effScaled.at(i) = std::max(0.f, 1.f - effraw.at(i)*SFraw.at(i));
    }
    else if (i==calibReaders.size()){
      eff.at(i) = std::max(0.f, effraw.at(i-1));
      effScaled.at(i) = std::max(0.f, effraw.at(i-1)*SFraw.at(i-1));
    }
    else{
      eff.at(i) = std::max(0.f, effraw.at(i-1) - effraw.at(i));
      effScaled.at(i) = std::max(0.f, effraw.at(i-1)*SFraw.at(i-1) - effraw.at(i)*SFraw.at(i));
    }
  }

  std::vector<float> btagwps = BtagHelpers::getBtagWPs(false);
  assert(btagwps.size()+1 == eff.size());
  unsigned int idx_tag = btagwps.size();
  for (unsigned int i=0; i<btagwps.size(); i++){
    if (btagval<btagwps.at(i)){
      idx_tag = i;
      break;
    }
  }

  float const& eff_region = eff.at(idx_tag);
  float const& eff_region_corr = effScaled.at(idx_tag);
  if (eff_region == 0.f){
    val = 1;
    if (effval) *effval = 0;
  }
  else{
    val = eff_region_corr / eff_region;
    if (effval) *effval = eff_region_corr;
  }

  if (this->verbosity>=TVar::DEBUG){
    MELAout << "Jet (pt, eta, flav, btagval) = (" << pt << ", " << eta << ", " << flav << ", " << btagval << "):" << endl;
    MELAout << "\t- WPs: " << btagwps << endl;
    MELAout << "\t- Raw effs: " << effraw << endl;
    MELAout << "\t- Raw SFs: " << SFraw << endl;
    MELAout << "\t- Region effs: " << eff << endl;
    MELAout << "\t- Region scaled effs: " << effScaled << endl;
    MELAout << "\t- Final SF: " << val << endl;
    MELAout << "\t- Final eff: " << eff_region_corr << endl;
  }
}
void BtagScaleFactorHandler::getSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, AK4JetObject const* obj, float& val, float* effval) const{
  val = 1;
  if (effval) *effval = 0;

  if (!obj) return;
  if (!ParticleSelectionHelpers::isJetForBtagSF(obj)) return;

  float pt = obj->pt() / obj->currentSystScale * obj->extras.JECNominal;
  getSFAndEff(syst, obj->getBTagJetFlavor(), obj->getBtagValue(), pt, obj->eta(), val, effval);
}
