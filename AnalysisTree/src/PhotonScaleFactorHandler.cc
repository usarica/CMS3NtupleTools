#include <CMSDataTools/AnalysisTree/interface/HostHelpersCore.h>
#include "PhotonScaleFactorHandler.h"
#include "ParticleSelectionHelpers.h"
#include "SampleHelpersCore.h"
#include "SamplesCore.h"
#include "MELAStreamHelpers.hh"


using namespace std;
using namespace SampleHelpers;
using namespace PhotonSelectionHelpers;
using namespace ParticleSelectionHelpers;
using namespace MELAStreamHelpers;


PhotonScaleFactorHandler::PhotonScaleFactorHandler() : ScaleFactorHandlerBase()
{
  this->setup();
}

PhotonScaleFactorHandler::~PhotonScaleFactorHandler(){ this->reset(); }

TString PhotonScaleFactorHandler::getScaleFactorFileName(PhotonSelectionHelpers::SelectionBits const& preselectionBit, int const& year){
  TString res = ANALYSISTREEPKGDATAPATH + Form("ScaleFactors/Photons/%i/", year);
  HostHelpers::ExpandEnvironmentVariables(res);
  if (idType_preselection == kCutBasedId_Fall17V2){
    switch (preselectionBit){
    case kLooseId:
      res += "Fall17V2_PhotonsLoose";
      break;
    case kMediumId:
      res += "Fall17V2_PhotonsMedium";
      break;
    case kTightId:
      res += "Fall17V2_PhotonsTight";
      break;
    default:
      MELAerr << "PhotonScaleFactorHandler::getScaleFactorFileName: Selection bit " << preselectionBit << " is not implemented for cut-based ids." << endl;
      assert(0);
    }
  }
  else if (idType_preselection == kMVAId_Fall17V2){
    switch (preselectionBit){
    case kLooseId:
      res += "Fall17V2_PhotonsMVAwp90";
      break;
    case kMediumId:
    case kTightId:
      res += "Fall17V2_PhotonsMVAwp80";
      break;
    default:
      MELAerr << "PhotonScaleFactorHandler::getScaleFactorFileName: Selection bit " << preselectionBit << " is not implemented for MVA ids." << endl;
      assert(0);
    }
  }
  else{
    MELAerr << "PhotonScaleFactorHandler::getScaleFactorFileName: Id type " << idType_preselection << " is not implemented." << endl;
    assert(0);
  }
  res += ".root";

  if (!HostHelpers::FileReadable(res.Data())){
    MELAerr << "PhotonScaleFactorHandler::getScaleFactorFileName: File " << res << " is not readable." << endl;
    assert(0);
  }

  return res;
}


bool PhotonScaleFactorHandler::setup(){
  bool res = true;
  this->reset();

  TDirectory* curdir = gDirectory;
  TDirectory* uppermostdir = SampleHelpers::rootTDirectory;

  // Get tampon SFs
  {
    TFile* finput = TFile::Open(getScaleFactorFileName(bit_SFTampon_id, SampleHelpers::theDataYear), "read"); uppermostdir->cd();
    res &= getHistogram<TH2F, ExtendedHistogram_2D>(h_eff_mc_tampon, finput, "EGamma_EffMC2D");
    res &= getHistogram<TH2F, ExtendedHistogram_2D>(h_SF_tampon, finput, "EGamma_SF2D");
    ScaleFactorHandlerBase::closeFile(finput); curdir->cd();
  }
  // Get tight SFs
  {
    TFile* finput = TFile::Open(getScaleFactorFileName(bit_preselectionTight_id, SampleHelpers::theDataYear), "read"); uppermostdir->cd();
    res &= getHistogram<TH2F, ExtendedHistogram_2D>(h_eff_mc_tight, finput, "EGamma_EffMC2D");
    res &= getHistogram<TH2F, ExtendedHistogram_2D>(h_SF_tight, finput, "EGamma_SF2D");
    ScaleFactorHandlerBase::closeFile(finput); curdir->cd();
  }

  return res;
}
void PhotonScaleFactorHandler::reset(){
  h_eff_mc_tampon.reset();
  h_SF_tampon.reset();

  h_eff_mc_tight.reset();
  h_SF_tight.reset();
}

void PhotonScaleFactorHandler::evalScaleFactorFromHistogram(float& theSF, float& theSFRelErr, float const& pt, float const& etaSC, ExtendedHistogram_2D const& hist, bool etaOnY, bool useAbsEta) const{
  TH2F const* hh = hist.getHistogram();
  if (!hh){
    MELAerr << "PhotonScaleFactorHandler::evalScaleFactorFromHistogram: Histogram is null." << endl;
    return;
  }

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

void PhotonScaleFactorHandler::getIdIsoSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, float const& pt, float const& etaSC, bool const& isTight, bool const& isTampon, float& val, float* effval) const{
  val = 1;

  float eff_tight=1, eff_relerr_tight=0;
  evalScaleFactorFromHistogram(eff_tight, eff_relerr_tight, pt, etaSC, h_eff_mc_tight, false, false);
  float eff_tampon=1, eff_relerr_tampon=0;
  evalScaleFactorFromHistogram(eff_tampon, eff_relerr_tampon, pt, etaSC, h_eff_mc_tampon, false, false);
  float eff_nonid = std::max(0.f, 1.f - eff_tampon);

  float SF_tight=1, SF_relerr_tight=0;
  evalScaleFactorFromHistogram(SF_tight, SF_relerr_tight, pt, etaSC, h_SF_tight, false, false);
  float SF_tampon=1, SF_relerr_tampon=0;
  evalScaleFactorFromHistogram(SF_tampon, SF_relerr_tampon, pt, etaSC, h_SF_tampon, false, false);

  if (syst == SystematicsHelpers::ePhoEffDn){
    SF_tight *= (1.f - SF_relerr_tight);
    SF_tampon *= (1.f - SF_relerr_tampon);
  }
  else if (syst == SystematicsHelpers::ePhoEffUp){
    SF_tight *= (1.f + SF_relerr_tight);
    SF_tampon *= (1.f + SF_relerr_tampon);
  }

  float eff_region_corr, eff_region;
  if (isTight){
    eff_region_corr = eff_tight * SF_tight;
    eff_region = eff_tight;
  }
  else if (isTampon){
    eff_region_corr = std::max(0.f, eff_tampon*SF_tampon - eff_tight*SF_tight);
    eff_region = std::max(0.f, eff_tampon - eff_tight);
  }
  else{
    eff_region_corr = std::max(0.f, 1.f - eff_tampon*SF_tampon);
    eff_region = eff_nonid;
  }

  /*
  if (eff_region == 0.f || eff_region_corr == 0.f){
    if (this->verbosity >= TVar::ERROR) MELAerr
      << "PhotonScaleFactorHandler::getIdIsoSFAndEff: Photon (pT, etaSC) = (" << pt << ", " << etaSC
      << ") efficiency (uncorrected, corrected) = (" << eff_region << ", " << eff_region_corr
      << ") with id (tight, tampon) = (" << isTight << ", " << isTampon << ")"
      << endl;
  }
  */

  if (eff_region == 0.f){
    val = 1;
    if (effval) *effval = 0;
  }
  else{
    val = eff_region_corr / eff_region;
    if (effval) *effval = eff_region_corr;
  }
}

void PhotonScaleFactorHandler::getIdIsoSFAndEff(SystematicsHelpers::SystematicVariationTypes const& syst, PhotonObject const* obj, float& val, float* effval) const{
  val = 1;
  if (effval) *effval = 0;

  if (!obj) return;

  bool const isTight = ParticleSelectionHelpers::isTightParticle(obj);
  bool const isTampon = obj->testSelectionBit(kSFTampon);

  getIdIsoSFAndEff(syst, obj->pt(), obj->etaSC(), isTight, isTampon, val, effval);
}
